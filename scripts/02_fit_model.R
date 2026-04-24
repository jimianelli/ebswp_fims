#!/usr/bin/env Rscript

source(file.path("scripts", "sparsenuts_framework.R"))
source("fimsfit.R", local = .GlobalEnv)

if (!exists("is_fims_verbose", mode = "function")) {
  is_fims_verbose <- function() FALSE
}

suppressPackageStartupMessages({
  library(FIMS)
})

project_root <- resolve_project_root()

built <- build_fims_inputs(project_root = project_root)

fit <- fit_fims(
  built$input,
  get_sd = TRUE,
  number_of_loops = 3,
  control = list(eval.max = 10000, iter.max = 10000, trace = 0)
)

fit_summary <- list(
  estimates = extract_fims_estimates(fit),
  report = get_report(fit),
  max_gradient = get_max_gradient(fit),
  version = get_version(fit),
  timing = get_timing(fit),
  number_of_parameters = get_number_of_parameters(fit)
)

out_path <- file.path(project_root, "outputs", "fims_fit_summary.rds")
if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)

saveRDS(fit_summary, out_path)
cat("Wrote:", out_path, "\n")
