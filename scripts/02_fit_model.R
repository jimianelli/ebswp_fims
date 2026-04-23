#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(FIMS)
})

in_path <- "/Users/jim/_mymods/pollock/fims/data/ebs_fims_data.rds"
if (!file.exists(in_path)) {
  stop("Missing input data. Run scripts/01_build_data.R first: ", in_path)
}

payload <- readRDS(in_path)

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
      module_type = "DoubleLogistic"
    ),
    by = c("module_name", "fleet_name")
  ) |>
  tidyr::nest(.by = c(model_family, module_name, fleet_name))

pars <- create_default_parameters(cfg, data_4_model) |>
  tidyr::unnest(cols = data) |>
  dplyr::rows_update(
    tibble::tibble(
      module_name = "Recruitment",
      label = "log_devs",
      time = (get_start_year(data_4_model) + 1):get_end_year(data_4_model),
      estimation_type = "fixed_effects"
    ),
    by = c("module_name", "label", "time")
  ) |>
  # NOTE: selectivity_shared_with is not recognized by FIMS >= 0.9.0;
  # AVO selectivity is estimated independently. Remove if confirmed unsupported.
  dplyr::mutate(
    selectivity_shared_with = dplyr::case_when(
      module_name == "Selectivity" & fleet_name == "avo" ~ "ats",
      TRUE ~ NA_character_
    )
  )

# Optional: initialize F from base.rds time series if available
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

fit <- fit_fims(
  input,
  get_sd = TRUE,
  number_of_loops = 3,
  control = list(eval.max = 10000, iter.max = 10000, trace = 0)
)

fit_summary <- list(
  estimates = get_estimates(fit),
  report = get_report(fit),
  max_gradient = get_max_gradient(fit),
  version = get_version(fit),
  timing = get_timing(fit),
  number_of_parameters = get_number_of_parameters(fit)
)

out_path <- "/Users/jim/_mymods/pollock/fims/outputs/fims_fit_summary.rds"
if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)

saveRDS(fit_summary, out_path)
cat("Wrote:", out_path, "\n")
