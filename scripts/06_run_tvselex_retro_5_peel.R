#!/usr/bin/env Rscript

source(file.path("scripts", "sparsenuts_framework.R"))
source("fimsfit.R", local = .GlobalEnv)

if (!exists("is_fims_verbose", mode = "function")) {
  is_fims_verbose <- function() FALSE
}

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(FIMS)
})

project_root <- resolve_project_root()
data_path <- file.path(project_root, "data", "ebs_fims_data.rds")
out_path <- file.path(project_root, "outputs", "tvselex_retro_5_peel_summary.rds")

if (!file.exists(data_path)) {
  stop("Missing input data. Run scripts/01_build_data.R first: ", data_path)
}

payload <- readRDS(data_path)

build_tvselex_inputs_from_payload <- function(payload) {
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
  apply_age_specific_natural_mortality() |>
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
        TRUE ~ NA_character_
      )
    )

  list(
    payload = payload,
    input = initialize_fims(pars_tvselex, data_4_model),
    parameters = pars_tvselex
  )
}

fit_one_peel <- function(peel) {
  message("Running TVselex retrospective peel ", peel, "...")

  peel_payload <- make_retro_payload(payload, peel = peel)
  built <- build_tvselex_inputs_from_payload(peel_payload)

  tryCatch(
    {
      fit <- fit_fims(
        built$input,
        get_sd = FALSE,
        number_of_loops = 1,
        control = list(eval.max = 10000, iter.max = 10000, trace = 0)
      )

      est <- get_estimates(fit)

      bind_rows(
        est |>
          filter(label == "spawning_biomass") |>
          transmute(
            run = "TVselex",
            peel = peel,
            quantity = "SSB",
            year = built$payload$years[year_i],
            value = estimated,
            max_gradient = get_max_gradient(fit),
            status = "fit completed",
            error = NA_character_
          ),
        est |>
          filter(label == "expected_recruitment") |>
          transmute(
            run = "TVselex",
            peel = peel,
            quantity = "Recruitment",
            year = built$payload$years[year_i],
            value = estimated,
            max_gradient = get_max_gradient(fit),
            status = "fit completed",
            error = NA_character_
          )
      )
    },
    error = function(e) {
      tibble::tibble(
        run = "TVselex",
        peel = peel,
        quantity = c("SSB", "Recruitment"),
        year = NA_integer_,
        value = NA_real_,
        max_gradient = NA_real_,
        status = "fit failed",
        error = conditionMessage(e)
      )
    }
  )
}

retro <- purrr::map_dfr(0:5, fit_one_peel)

if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)
saveRDS(retro, out_path)
message("Wrote: ", out_path)
