#!/usr/bin/env Rscript
# Interactive exploration script — loads pre-fitted model from outputs/fims_fit.rds
# To (re)fit: run scripts/01_build_data.R then scripts/02_fit_model.R

suppressPackageStartupMessages({
  library(FIMS)
  library(dplyr)
  library(ggplot2)
})

fit_path <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "outputs/fims_fit.rds")
if (!file.exists(fit_path)) stop("No saved fit found. Run scripts/02_fit_model.R first: ", fit_path)
fit <- readRDS(fit_path)

# ---- Quick diagnostics ----
print(fit)
get_version(fit)
get_timing(fit)
get_max_gradient(fit)

# ---- Parameter estimates ----
est <- get_estimates(fit)
dplyr::count(est, module_name, label)

# ---- Derived quantities from TMB report ----
rep <- get_report(fit)
names(rep)

# ---- Spawning biomass ----
years <- get_input(fit)$data_4_model |> (\(x) get_start_year(x):get_end_year(x))()
ssb <- data.frame(year = years, SSB = unlist(rep$spawning_biomass))

ggplot(ssb, aes(x = year, y = SSB)) +
  geom_line() +
  geom_point() +
  labs(title = "EBS Pollock — Spawning Biomass", y = "SSB (mt)", x = NULL)

# ---- JSON estimates (tidy) ----
json <- get_model_output(fit)
json_est <- reshape_json_estimates(json)

ssb_est <- json_est |>
  dplyr::filter(label %in% c("spawning_biomass", "SSB"))

# ---- Raw TMB objects (for diagnostics) ----
# get_obj(fit)      # TMB object
# get_opt(fit)      # nlminb result
# get_sdreport(fit) # standard errors
