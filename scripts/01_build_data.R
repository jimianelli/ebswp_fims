#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

pm_path <- "/Users/jim/_mymods/pollock/openMSE/data/pm_24.rds"
base_path <- "/Users/jim/_mymods/pollock/openMSE/analysis/outputs/base.rds"

if (!file.exists(pm_path)) {
  stop("Missing input: ", pm_path)
}

pm <- readRDS(pm_path)
base <- if (file.exists(base_path)) readRDS(base_path) else NULL

cv_2_sd <- function(x) sqrt(log(x^2 + 1))

years <- pm$styr:pm$endyr
ages <- pm$recage:(pm$recage + pm$nages - 1)

# ---- Landings (fishery) ----
# NOTE: No catch uncertainty provided; set a small constant log-SD.
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
index_cpue <- make_index("cpue", pm$yrs_cpue, pm$obs_cpue, pm$obs_cpue_std)

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
# FIMS uses a single growth weight-at-age vector internally. For SSB, use the
# stock-level spawning weight-at-age series rather than averaging fleet-specific
# weight-at-age inputs.
# FIMS >= 0.9 requires weight_at_age for every age × (n_years + 1).
# FIMSFrame drops years before the modeled data start, so append one terminal
# year rather than prepending one initialization year.
make_waa <- function(name, years_obs, mat) {
  colnames(mat) <- ages

  waa_obs <- as_tibble(mat) |>
    mutate(timing = years_obs) |>
    pivot_longer(-timing, names_to = "age", values_to = "value") |>
    mutate(age = as.numeric(age))

  waa_mean <- waa_obs |>
    group_by(age) |>
    summarize(value_mean = mean(value, na.rm = TRUE), .groups = "drop")

  waa_years <- min(years):(max(years) + 1)

  expand_grid(timing = waa_years, age = ages) |>
    left_join(waa_obs, by = c("timing", "age")) |>
    left_join(waa_mean, by = "age") |>
    mutate(
      value = dplyr::coalesce(value, value_mean),
      type = "weight_at_age",
      name = name,
      length = NA_real_,
      unit = "mt",
      uncertainty = NA_real_
    ) |>
    select(type, name, age, length, timing, value, unit, uncertainty)
}

waa <- make_waa("ssb", years, pm$wt_ssb)

# ---- Combine ----
data_ebs <- bind_rows(
  landings,
  index_bts, index_ats, index_avo, index_cpue,
  age_fsh, age_bts, age_ats,
  waa
)

out_path <- "/Users/jim/_mymods/pollock/fims/data/ebs_fims_data.rds"
if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)

saveRDS(
  list(
    data_ebs = data_ebs,
    years = years,
    ages = ages,
    base = base
  ),
  out_path
)

cat("Wrote:", out_path, "\n")
