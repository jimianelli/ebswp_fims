# Build a simplified FIMS baseline for EBS pollock from openMSE pm_24.rds
# - No seasonal timing
# - No length comps
# - Recruitment deviations as fixed effects
# - Flexible selectivity (DoubleLogistic) for fishery and surveys

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(FIMS)
})

pm_path <- "pm_24.rds"
base_path <- "base.rds"

pm <- readRDS(pm_path)
base <- readRDS(base_path)

cv_2_sd <- function(x) sqrt(log(x^2 + 1))

# Core dimensions
years <- pm$styr:pm$endyr
ages <- pm$recage:(pm$recage + pm$nages - 1)

# ---- Landings (fishery) ----
# NOTE: No catch uncertainty provided; set a small constant log-SD for now.
catch_log_sd <- 0.05
landings <- tibble::tibble(
  type = "landings",
  name = "fishery",
  age = NA_real_,
  length = NA_real_,
  timing = years,
  value = pm$obs_catch,
  unit = "mt",
  uncertainty = rep(catch_log_sd, length(years))
)

# ---- Indices (3 surveys) ----
make_index <- function(name, years_obs, obs, obs_sd) {
  # Convert sd to log-SD if possible; fallback to mean of available log-SD
  cv <- obs_sd / obs
  cv[!is.finite(cv) | cv <= 0] <- NA_real_
  log_sd <- cv_2_sd(cv)
  fill_sd <- mean(log_sd, na.rm = TRUE)
  if (!is.finite(fill_sd)) {
    fill_sd <- 0.2
  }
  log_sd[!is.finite(log_sd) | log_sd <= 0] <- fill_sd

  full_value <- rep(-999, length(years))
  full_sd <- rep(fill_sd, length(years))
  idx <- match(years_obs, years)
  full_value[idx] <- obs
  full_sd[idx] <- log_sd

  tibble::tibble(
    type = "index",
    name = name,
    age = NA_real_,
    length = NA_real_,
    timing = years,
    value = full_value,
    unit = "mt",
    uncertainty = full_sd
  )
}

index_bts <- make_index("bts", pm$yrs_bts_data, pm$ob_bts, pm$ob_bts_std)
index_ats <- make_index("ats", pm$yrs_ats_data, pm$ob_ats, pm$ob_ats_std)
index_avo <- make_index("avo", pm$yrs_avo, pm$ob_avo, pm$ob_avo_std)

# ---- Age compositions (fishery, BTS, ATS) ----
make_agecomp <- function(name, years_obs, mat, sample_size, already_prop = FALSE) {
  if (!already_prop) {
    rs <- rowSums(mat, na.rm = TRUE)
    rs[rs <= 0 | !is.finite(rs)] <- NA_real_
    mat <- sweep(mat, 1, rs, "/")
  }

  colnames(mat) <- ages
  obs <- as_tibble(mat) |>
    mutate(timing = years_obs) |>
    pivot_longer(-timing, names_to = "age", values_to = "value") |>
    mutate(
      age = as.numeric(age),
      uncertainty = rep(sample_size, each = length(ages))
    )

  full_grid <- expand_grid(timing = years, age = ages)
  out <- full_grid |>
    left_join(obs, by = c("timing", "age")) |>
    mutate(
      value = ifelse(is.na(value), -999, value),
      uncertainty = ifelse(is.na(uncertainty), 0, uncertainty)
    ) |>
    mutate(
      type = "age_comp",
      name = name,
      length = NA_real_,
      unit = "proportion"
    ) |>
    select(type, name, age, length, timing, value, unit, uncertainty)

  out
}

age_fsh <- make_agecomp(
  "fishery",
  pm$yrs_fsh_data,
  pm$oac_fsh_data,
  pm$sam_fsh,
  already_prop = FALSE
)

age_bts <- make_agecomp(
  "bts",
  pm$yrs_bts_data,
  pm$oac_bts,
  pm$sam_bts,
  already_prop = FALSE
)

age_ats <- make_agecomp(
  "ats",
  pm$yrs_ats_data,
  pm$oac_ats,
  pm$sam_ats,
  already_prop = FALSE
)

# ---- Weight-at-age ----
# FIMS >= 0.9 requires weight_at_age for every year × age combination.
# Use ATS WAA where observed; fill remaining years with the age-specific mean.
colnames(pm$wt_ats) <- ages
waa_obs <- as_tibble(pm$wt_ats) |>
  mutate(timing = pm$yrs_ats_data) |>
  pivot_longer(-timing, names_to = "age", values_to = "value") |>
  mutate(age = as.numeric(age))

waa_mean <- waa_obs |>
  group_by(age) |>
  summarize(value_mean = mean(value, na.rm = TRUE), .groups = "drop")

# FIMS initialize_fims uses one extra initialization year before the data start,
# so supply WAA for (styr - 1):endyr rather than styr:endyr.
waa_years <- (min(years) - 1):max(years)
waa <- expand_grid(timing = waa_years, age = ages) |>
  left_join(waa_obs, by = c("timing", "age")) |>
  left_join(waa_mean, by = "age") |>
  mutate(
    value = dplyr::coalesce(value, value_mean),
    type = "weight_at_age",
    name = "fishery",
    length = NA_real_,
    unit = "mt",
    uncertainty = NA_real_
  ) |>
  select(type, name, age, length, timing, value, unit, uncertainty)

# ---- Combine ----
# NOTE: This data frame is ready for FIMSFrame()

data_ebs <- bind_rows(
  landings,
  index_bts, index_ats, index_avo,
  age_fsh, age_bts, age_ats,
  waa
)

# ---- Build FIMS objects ----
#
# data_4_model <- FIMSFrame(data_ebs)
#
# cfg <- create_default_configurations(data_4_model) |>
#   tidyr::unnest(cols = data) |>
#   dplyr::rows_update(
#     tibble::tibble(
#       module_name = "Selectivity",
#       fleet_name = c("fishery", "bts", "ats", "avo"),
#       module_type = "DoubleLogistic"
#     ),
#     by = c("module_name", "fleet_name")
#   ) |>
#   tidyr::nest(.by = c(model_family, module_name, fleet_name))
#
# pars <- create_default_parameters(cfg, data_4_model) |>
#   tidyr::unnest(cols = data) |>
#   dplyr::rows_update(
#     tibble::tibble(
#       module_name = "Recruitment",
#       label = "log_devs",
#       time = (get_start_year(data_4_model) + 1):get_end_year(data_4_model),
#       estimation_type = "fixed_effects"
#     ),
#     by = c("module_name", "label", "time")
#   )
#
# # Share ATS selectivity with AVO
# pars <- pars |>
#   dplyr::mutate(
#     selectivity_shared_with = dplyr::case_when(
#       module_name == "Selectivity" & fleet_name == "avo" ~ "ats",
#       TRUE ~ NA_character_
#     )
#   )
#
# # Optional: initialize log_Fmort with base.rds time series if available
# if (!is.null(base$report$F)) {
#   pars <- pars |>
#     dplyr::rows_update(
#       tibble::tibble(
#         fleet_name = "fishery",
#         label = "log_Fmort",
#         time = years,
#         value = log(base$report$F)
#       ),
#       by = c("fleet_name", "label", "time")
#     )
# }
#
# # Fit
# fit <- pars |>
#   initialize_fims(data_4_model) |>
#   fit_fims(optimize = TRUE)
