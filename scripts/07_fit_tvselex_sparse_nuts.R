#!/usr/bin/env Rscript

source(file.path("scripts", "sparsenuts_framework.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(FIMS)
  library(SparseNUTS)
})

project_root <- resolve_project_root()
out_path <- file.path(project_root, "outputs", "tvselex_sparsenuts_fit.rds")
fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message("Building TVselex mode estimate for SparseNUTS...")
mle_obj <- fit_tvselex_mode_for_snuts(
  project_root = project_root,
  get_sd = FALSE
)

message("Running SparseNUTS on TVselex with package defaults...")
snuts_fit <- tryCatch(
  SparseNUTS::sample_snuts(
    obj = get_obj(mle_obj$fit),
    model_name = "ebswp_fims_tvselex"
  ),
  error = function(e) e
)

if (inherits(snuts_fit, "error")) {
  out <- list(
    status = "SparseNUTS failed",
    error = conditionMessage(snuts_fit),
    sampler_settings = "SparseNUTS::sample_snuts() package defaults; only model_name supplied.",
    mle = list(
      max_gradient = get_max_gradient(mle_obj$fit),
      version = get_version(mle_obj$fit),
      timing = get_timing(mle_obj$fit),
      number_of_parameters = get_number_of_parameters(mle_obj$fit)
    )
  )
  saveRDS(out, out_path)
  message("Wrote failure summary: ", out_path)
  stop(conditionMessage(snuts_fit), call. = FALSE)
}

parameter_lookup <- mle_obj$built$parameters |>
  dplyr::filter(estimation_type %in% c("fixed_effects", "random_effects")) |>
  dplyr::mutate(parameter_index = dplyr::row_number())

key_parameter_lookup <- dplyr::bind_rows(
  parameter_lookup |>
    dplyr::filter(module_name == "Recruitment", label %in% c("log_rzero", "logit_steep", "log_sd")),
  parameter_lookup |>
    dplyr::filter(module_name == "Fleet", label == "log_q"),
  parameter_lookup |>
    dplyr::filter(
      module_name == "Selectivity",
      fleet_name == "fishery",
      label == "inflection_point_asc",
      time %in% stats::quantile(time, probs = c(0, 0.5, 1), na.rm = TRUE, names = FALSE)
    )
) |>
  dplyr::distinct(parameter_index, .keep_all = TRUE) |>
  dplyr::arrange(parameter_index)

pairs_path <- file.path(fig_dir, "tvselex_sparsenuts_pairs_slow.png")
png(pairs_path, width = 1600, height = 1600, res = 180)
pairs(snuts_fit, order = "slow")
dev.off()

marginals_path <- file.path(fig_dir, "tvselex_sparsenuts_key_marginals.png")
png(marginals_path, width = 1800, height = 1400, res = 180)
SparseNUTS::plot_marginals(
  snuts_fit,
  pars = key_parameter_lookup$parameter_index,
  mfrow = c(3, 3)
)
dev.off()

out <- list(
  snuts_fit = snuts_fit,
  mle = list(
    max_gradient = get_max_gradient(mle_obj$fit),
    version = get_version(mle_obj$fit),
    timing = get_timing(mle_obj$fit),
    number_of_parameters = get_number_of_parameters(mle_obj$fit)
  ),
  model_input = list(
    years = mle_obj$built$payload$years,
    ages = mle_obj$built$payload$ages
  ),
  parameter_lookup = parameter_lookup,
  key_parameter_lookup = key_parameter_lookup,
  figures = list(
    pairs_slow = pairs_path,
    key_marginals = marginals_path
  ),
  diagnostics = SparseNUTS::check_snuts_diagnostics(snuts_fit, print = FALSE)
)

saveRDS(out, out_path)
message("Wrote: ", out_path)
message("Wrote: ", pairs_path)
message("Wrote: ", marginals_path)
