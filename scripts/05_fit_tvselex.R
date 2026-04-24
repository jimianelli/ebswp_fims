#!/usr/bin/env Rscript

source(file.path("scripts", "sparsenuts_framework.R"))
source("fimsfit.R", local = .GlobalEnv)

if (!exists("is_fims_verbose", mode = "function")) {
  is_fims_verbose <- function() FALSE
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(FIMS)
})

project_root <- resolve_project_root()
data_path <- file.path(project_root, "data", "ebs_fims_data.rds")
if (!file.exists(data_path)) {
  stop("Missing input data. Run scripts/01_build_data.R first: ", data_path)
}

payload <- readRDS(data_path)
data_4_model <- FIMSFrame(payload$data_ebs)

cfg_tvselex <- create_default_configurations(data_4_model) |>
  tidyr::unnest(cols = data) |>
  dplyr::rows_update(
    tibble::tibble(
      module_name = "Selectivity",
      fleet_name = c("fishery", "bts", "ats", "avo"),
      module_type = c("DoubleLogistic", "Logistic", "DoubleLogistic", "DoubleLogistic")
    ),
    by = c("module_name", "fleet_name")
  ) |>
  tidyr::nest(.by = c(model_family, module_name, fleet_name))

pars_tvselex_base <- create_default_parameters(cfg_tvselex, data_4_model) |>
  tidyr::unnest(cols = data) |>
  dplyr::rows_update(
    tibble::tibble(
      module_name = "Recruitment",
      label = "log_devs",
      time = (get_start_year(data_4_model) + 1):get_end_year(data_4_model),
      estimation_type = "fixed_effects"
    ),
    by = c("module_name", "label", "time")
  )

fishery_tvselex <- pars_tvselex_base |>
  dplyr::filter(
    module_name == "Selectivity",
    fleet_name == "fishery",
    label == "inflection_point_asc"
  ) |>
  dplyr::select(-time) |>
  tidyr::crossing(time = payload$years) |>
  dplyr::mutate(estimation_type = "fixed_effects")

pars_tvselex <- pars_tvselex_base |>
  dplyr::filter(
    !(
      module_name == "Selectivity" &
        fleet_name == "fishery" &
        label == "inflection_point_asc"
    )
  ) |>
  dplyr::bind_rows(fishery_tvselex) |>
  dplyr::mutate(
    selectivity_shared_with = dplyr::case_when(
      module_name == "Selectivity" & fleet_name == "cpue" ~ "fishery",
      module_name == "Selectivity" & fleet_name == "avo" ~ "ats",
      TRUE ~ NA_character_
    )
  )

out_path <- file.path(project_root, "outputs", "tvselex_fit_summary.rds")
if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)

input_tvselex <- initialize_fims(pars_tvselex, data_4_model)

fit_summary <- tryCatch(
  {
    fit <- fit_fims(
      input_tvselex,
      get_sd = FALSE,
      number_of_loops = 1,
      control = list(eval.max = 10000, iter.max = 10000, trace = 0)
    )

    list(
      run = "TVselex",
      status = "fit completed",
      estimates = get_estimates(fit),
      report = get_report(fit),
      max_gradient = get_max_gradient(fit),
      version = get_version(fit),
      timing = get_timing(fit),
      number_of_parameters = get_number_of_parameters(fit),
      fishery_selectivity_rows = sum(
        pars_tvselex$module_name == "Selectivity" &
          pars_tvselex$fleet_name == "fishery" &
          pars_tvselex$label == "inflection_point_asc" &
          !is.na(pars_tvselex$time)
      )
    )
  },
  error = function(e) {
    list(
      run = "TVselex",
      status = "fit failed",
      error = conditionMessage(e),
      version = utils::packageVersion("FIMS"),
      number_of_parameters = c(
        fixed_effects = length(input_tvselex$parameters$p),
        random_effects = length(input_tvselex$parameters$re)
      ),
      fishery_selectivity_rows = sum(
        pars_tvselex$module_name == "Selectivity" &
          pars_tvselex$fleet_name == "fishery" &
          pars_tvselex$label == "inflection_point_asc" &
          !is.na(pars_tvselex$time)
      )
    )
  }
)

saveRDS(fit_summary, out_path)
cat("Wrote:", out_path, "\n")
if (!identical(fit_summary$status, "fit completed")) {
  stop(fit_summary$error, call. = FALSE)
}
