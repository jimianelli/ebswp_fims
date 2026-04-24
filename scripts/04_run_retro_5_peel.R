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
out_path <- file.path(project_root, "outputs", "retro_5_peel_summary.rds")

if (!file.exists(data_path)) {
  stop("Missing input data. Run scripts/01_build_data.R first: ", data_path)
}

payload <- readRDS(data_path)

fit_one_peel <- function(peel) {
  message("Running retrospective peel ", peel, "...")

  peel_payload <- make_retro_payload(payload, peel = peel)
  built <- build_fims_inputs_from_payload(peel_payload)

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
        peel = peel,
        quantity = "SSB",
        year = built$payload$years[year_i],
        value = estimated,
        max_gradient = get_max_gradient(fit)
      ),
    est |>
      filter(label == "expected_recruitment") |>
      transmute(
        peel = peel,
        quantity = "Recruitment",
        year = built$payload$years[year_i],
        value = estimated,
        max_gradient = get_max_gradient(fit)
      )
  )
}

retro <- purrr::map_dfr(0:5, fit_one_peel)

saveRDS(retro, out_path)
message("Wrote: ", out_path)
