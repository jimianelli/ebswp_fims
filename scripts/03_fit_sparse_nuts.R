#!/usr/bin/env Rscript

parse_arg <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[[1]])
}

args <- commandArgs(trailingOnly = TRUE)

source(file.path("scripts", "sparsenuts_framework.R"))

project_root <- resolve_project_root()

num_samples <- as.integer(parse_arg(args, "num-samples", 250))
num_warmup <- as.integer(parse_arg(args, "num-warmup", 250))
chains <- as.integer(parse_arg(args, "chains", 4))
cores <- as.integer(parse_arg(args, "cores", 1))
thin <- as.integer(parse_arg(args, "thin", 1))
seed <- as.integer(parse_arg(args, "seed", 123))
metric <- parse_arg(args, "metric", "diag")
out_path <- parse_arg(args, "out", file.path(project_root, "outputs", "sparsenuts_fit.rds"))

message("Building FIMS mode estimate for SparseNUTS...")
mle_obj <- fit_fims_mode_for_snuts(
  project_root = project_root,
  get_sd = FALSE
)

message("Running SparseNUTS...")
snuts_fit <- run_sparse_nuts_fims(
  fit = mle_obj$fit,
  num_samples = num_samples,
  num_warmup = num_warmup,
  chains = chains,
  cores = cores,
  thin = thin,
  seed = seed,
  metric = metric
)

save_sparse_nuts_results(
  snuts_fit = snuts_fit,
  mle_fit = mle_obj$fit,
  built = mle_obj$built,
  out_path = out_path
)

message("Wrote: ", out_path)
