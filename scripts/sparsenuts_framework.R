resolve_project_root <- function() {
  cwd <- getwd()

  candidates <- c(
    cwd,
    normalizePath(file.path(cwd, ".."), mustWork = FALSE)
  )

  for (path in candidates) {
    if (file.exists(file.path(path, "fimsfit.R")) &&
        file.exists(file.path(path, "scripts", "01_build_data.R"))) {
      return(normalizePath(path))
    }
  }

  stop("Could not resolve project root from working directory: ", cwd)
}

default_snuts_control <- function() {
  list(
    adapt_delta = 0.9,
    max_treedepth = 12,
    adapt_window = 50,
    adapt_init_buffer = 75,
    adapt_term_buffer = 50
  )
}

build_fims_inputs_from_payload <- function(payload) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(FIMS)
  })

  data_ebs <- payload$data_ebs
  pm_years <- payload$years
  base <- payload$base

  clear()

  data_4_model <- FIMSFrame(data_ebs)

  cfg <- create_default_configurations(data_4_model) |>
    tidyr::unnest(cols = data) |>
    dplyr::rows_update(
      tibble::tibble(
        module_name = "Selectivity",
        fleet_name = c("fishery", "bts", "ats", "avo"),
        module_type = c("Logistic", "Logistic", "DoubleLogistic", "DoubleLogistic")
      ),
      by = c("module_name", "fleet_name")
    ) |>
    tidyr::nest(.by = c(model_family, module_name, fleet_name))

  pars <- create_default_parameters(cfg, data_4_model) |>
    tidyr::unnest(cols = data) |>
    dplyr::rows_update(
      tibble::tibble(
        module_name = "Recruitment",
        label = c(
          rep("log_devs", length((get_start_year(data_4_model) + 1):get_end_year(data_4_model))),
          "log_sd"
        ),
        time = c(
          (get_start_year(data_4_model) + 1):get_end_year(data_4_model),
          NA_real_
        ),
        value = c(
          rep(0, length((get_start_year(data_4_model) + 1):get_end_year(data_4_model))),
          0.1
        ),
        estimation_type = c(
          rep("fixed_effects", length((get_start_year(data_4_model) + 1):get_end_year(data_4_model))),
          "constant"
        )
      ),
      by = c("module_name", "label", "time")
    ) |>
    dplyr::rows_update(
      tibble::tibble(
        module_name = "Selectivity",
        fleet_name = "avo",
        label = c(
          "inflection_point_asc",
          "slope_asc",
          "inflection_point_desc",
          "slope_desc"
        ),
        value = c(1.5, 2.0, 8.0, 0.1),
        estimation_type = "constant"
      ),
      by = c("module_name", "fleet_name", "label")
    ) |>
    dplyr::mutate(
      selectivity_shared_with = dplyr::case_when(
        module_name == "Selectivity" & fleet_name == "cpue" ~ "fishery",
        TRUE ~ NA_character_
      )
    )

  if (!is.null(base) && !is.null(base$report$F)) {
    pars <- pars |>
      dplyr::rows_update(
        tibble::tibble(
          fleet_name = "fishery",
          label = "log_Fmort",
          time = pm_years,
          value = log(base$report$F)
        ),
        by = c("fleet_name", "label", "time")
      )
  }

  input <- initialize_fims(pars, data_4_model)

  list(
    payload = payload,
    data_4_model = data_4_model,
    configurations = cfg,
    parameters = pars,
    input = input
  )
}

build_fims_inputs <- function(project_root = resolve_project_root()) {
  data_path <- file.path(project_root, "data", "ebs_fims_data.rds")
  if (!file.exists(data_path)) {
    stop("Missing input data. Run scripts/01_build_data.R first: ", data_path)
  }

  payload <- readRDS(data_path)
  build_fims_inputs_from_payload(payload)
}

fit_fims_mode_for_snuts <- function(project_root = resolve_project_root(),
                                    get_sd = FALSE,
                                    number_of_loops = 3,
                                    control = list(eval.max = 10000, iter.max = 10000, trace = 0)) {
  suppressPackageStartupMessages({
    library(FIMS)
  })

  source(file.path(project_root, "fimsfit.R"), local = .GlobalEnv)
  if (!exists("is_fims_verbose", mode = "function")) {
    is_fims_verbose <- function() FALSE
  }

  built <- build_fims_inputs(project_root = project_root)

  fit <- fit_fims(
    built$input,
    get_sd = get_sd,
    number_of_loops = number_of_loops,
    control = control
  )

  list(
    fit = fit,
    built = built
  )
}

make_retro_payload <- function(payload, peel) {
  stopifnot(peel >= 0)

  terminal_year <- max(payload$years) - peel
  peeled_years <- payload$years[payload$years <= terminal_year]

  peeled_data <- payload$data_ebs |>
    dplyr::filter(
      dplyr::case_when(
        type == "weight_at_age" ~ timing <= (terminal_year + 1),
        TRUE ~ timing <= terminal_year
      )
    )

  list(
    data_ebs = peeled_data,
    years = peeled_years,
    ages = payload$ages,
    base = payload$base
  )
}

run_sparse_nuts_fims <- function(fit,
                                 num_samples = 250,
                                 num_warmup = 250,
                                 chains = 4,
                                 cores = 1,
                                 thin = 1,
                                 seed = 123,
                                 metric = "diag",
                                 control = default_snuts_control(),
                                 init = "last.par.best",
                                 model_name = "ebswp_fims") {
  suppressPackageStartupMessages({
    library(SparseNUTS)
  })

  obj <- get_obj(fit)

  SparseNUTS::sample_snuts(
    obj = obj,
    num_samples = num_samples,
    num_warmup = num_warmup,
    chains = chains,
    cores = cores,
    thin = thin,
    seed = seed,
    metric = metric,
    init = init,
    control = control,
    model_name = model_name
  )
}

extract_fims_estimates <- function(fit) {
  obj <- get_obj(fit)
  sdreport <- get_sdreport(fit)
  opt <- get_opt(fit)
  parameter_names <- names(obj[["par"]])

  tmb_output <- FIMS:::reshape_tmb_estimates(
    obj = obj,
    sdreport = sdreport,
    opt = opt,
    parameter_names = parameter_names
  )

  json_output <- FIMS:::reshape_json_estimates(get_model_output(fit))

  estimates <- dplyr::left_join(
    json_output,
    tmb_output |>
      dplyr::filter(!is.na(parameter_id)) |>
      dplyr::select(-initial, -module_name, -module_id, -estimate, -label),
    by = c("parameter_id")
  )

  if (all(c("uncertainty.x", "uncertainty.y") %in% names(estimates))) {
    estimates <- estimates |>
      dplyr::mutate(
        uncertainty = dplyr::coalesce(uncertainty.x, uncertainty.y),
        .after = "estimation_type"
      ) |>
      dplyr::select(-uncertainty.x, -uncertainty.y)
  } else if ("uncertainty.x" %in% names(estimates)) {
    estimates <- estimates |>
      dplyr::rename(uncertainty = uncertainty.x)
  } else if ("uncertainty.y" %in% names(estimates)) {
    estimates <- estimates |>
      dplyr::rename(uncertainty = uncertainty.y)
  }

  estimates
}

save_sparse_nuts_results <- function(snuts_fit,
                                     mle_fit,
                                     built,
                                     out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  out <- list(
    snuts_fit = snuts_fit,
    mle = list(
      max_gradient = get_max_gradient(mle_fit),
      version = get_version(mle_fit),
      timing = get_timing(mle_fit),
      number_of_parameters = get_number_of_parameters(mle_fit)
    ),
    model_input = list(
      years = built$payload$years,
      ages = built$payload$ages
    )
  )

  saveRDS(out, out_path)
  invisible(out)
}
