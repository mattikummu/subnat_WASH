# 49b_prepare_adm0_historical_sanit.R
# ---------------------------------------------------------------------------
# Prepare ADM0 (country-level) historical sanitation data (1990-2015) from the
# older JMP estimates (1980-2015).
#
# This file contains only s_bas (= san_improved_t) and
# s_od (= san_open_defecation_t).  For s_lim and s_unimp we
# rely entirely on XGBoost predictions, scaled so that
#   s_lim + s_unimp = 100 − (s_bas + s_od)
#
# Method:
#   1. Read data_in/estimates_1980-2015_v1.xlsx, map to ISO3
#   2. Extract ADM0 predictors (population-weighted) for 1990-2015
#   3. Predict all four indicators with saved XGBoost ensembles
#   4. Linearly interpolate observed s_bas & s_od
#   5. Extrapolate observed s_bas & s_od using XGBoost growth ratios
#   6. For s_lim & s_unimp: use XGBoost predictions, then rescale
#      so their sum = 100 − (s_bas + s_od)
#
# Inputs:
#   data_in/estimates_1980-2015_v1.xlsx
#   /Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg  (ADM_0 layer)
#   Predictor rasters (same as script 49, minus aridity & gw_depth)
#   results_sanit/xgboost_sanit_ensemble_*.rds
#
# Outputs:
#   results_sanit/xgb_raster_adm0_historical_sanit.csv
#   results_sanit/sanit_adm0_xgboost_predictions_hist.csv
#   results_sanit/JMP_adm0_sanitation_historical.csv
#   results_sanit/JMP_adm0_sanitation_interpPredExtrap_hist.csv
#   results_sanit/{indicator}_adm0_interpPredExtrap_hist.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(tidymodels)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(purrr)
  library(stringr)
  library(zoo)
  library(countrycode)
})

source("cache_prompt.R")

dir.create("results_sanit", showWarnings = FALSE)

sanit_indicators <- c("s_bas", "s_lim", "s_unimp", "s_od")
predictor_cols   <- c(
  "travel_time", "poverty", "msch", "gnic", "gini",
  "urban", "gdi", "corruption"
)

years_hist <- 1990:2015
gadm_path  <- "/Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg"

clamp_pct <- function(x) pmax(0, pmin(100, x))

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  w_sum <- sum(w[ok], na.rm = TRUE)
  if (is.na(w_sum) || w_sum <= 0) return(NA_real_)
  weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

# ---- 1. Read & clean historical JMP sanitation estimates ----
message("== Reading historical JMP sanitation estimates ==")
xlsx_path <- "data_in/estimates_1980-2015_v1.xlsx"
if (!file.exists(xlsx_path)) stop("Input file not found: ", xlsx_path)

hist_raw <- read_excel(xlsx_path, sheet = "Sheet1")

hist_sanit <- hist_raw |>
  transmute(
    country_file,
    country_clean = country_file |>
      str_replace("_15$", "") |>
      str_replace_all("_", " ") |>
      str_to_title(),
    year  = as.numeric(year),
    s_bas = as.numeric(san_improved_t),
    s_od  = as.numeric(san_open_defecation_t)
  ) |>
  mutate(iso3 = countrycode(country_clean, "country.name", "iso3c", warn = FALSE)) |>
  filter(!is.na(iso3), year >= 1990, year <= 2015) |>
  select(iso3, country_name = country_clean, year, s_bas, s_od)

n_unmapped <- hist_raw |>
  mutate(
    cc = str_replace(country_file, "_15$", "") |>
      str_replace_all("_", " ") |>
      str_to_title()
  ) |>
  filter(is.na(countrycode(cc, "country.name", "iso3c", warn = FALSE))) |>
  distinct(country_file) |>
  nrow()

message(sprintf("  Parsed %d rows, %d countries (%d territories unmapped)",
                nrow(hist_sanit), n_distinct(hist_sanit$iso3), n_unmapped))

write_csv(hist_sanit, "results_sanit/JMP_adm0_sanitation_historical.csv")
message("  Saved: results_sanit/JMP_adm0_sanitation_historical.csv")

message("== Loading global ADM0 predictors and filtering to 1990-2015 ==")
adm0_predictors <- read_csv("results/global_adm0_predictors_1990_2023.csv",
                            show_col_types = FALSE) |>
  filter(year %in% years_hist) |>
  select(iso3, country_name_gadm, year, pop, all_of(predictor_cols)) |>
  arrange(iso3, year)

# ---- 4. XGBoost predictions for all four indicators ----
predict_ensemble <- function(ensemble_path, new_data) {
  ens <- readRDS(ensemble_path)

  out  <- rep(NA_real_, nrow(new_data))
  keep <- complete.cases(new_data[, predictor_cols])

  if (any(keep)) {
    pred_matrix <- map(ens$models, function(mod) {
      predict(mod, new_data = new_data[keep, predictor_cols, drop = FALSE])$.pred
    }) |> do.call(cbind, args = _)

    out[keep] <- apply(pred_matrix, 1, median)
  }

  clamp_pct(out)
}

message("== Predicting ADM0 values with XGBoost ensembles (1990-2015) ==")
adm0_predictions <- adm0_predictors |>
  select(iso3, year, pop, all_of(predictor_cols))

for (indicator in sanit_indicators) {
  ens_path <- sprintf("results_sanit/xgboost_sanit_ensemble_%s.rds", indicator)
  if (!file.exists(ens_path)) {
    stop("Required ensemble model not found: ", ens_path,
         "\nRun 47_train_xgboost_sanit.R first.")
  }
  message("  Predicting ", indicator, " ...")
  adm0_predictions[[paste0("pred_", indicator)]] <-
    predict_ensemble(ens_path, adm0_predictions)
}

write_csv(adm0_predictions, "results_sanit/sanit_adm0_xgboost_predictions_hist.csv")
message(sprintf("  Saved: results_sanit/sanit_adm0_xgboost_predictions_hist.csv (%d rows)",
                nrow(adm0_predictions)))

# ---- 5. Interpolate & extrapolate s_bas and s_od ----
interp_extrap_indicator <- function(indicator, obs_data) {
  message("\n--- Interpolate/extrapolate: ", indicator, " ---")

  pred_col <- paste0("pred_", indicator)

  obs <- obs_data |>
    select(iso3, year, reported = all_of(indicator))

  pred <- adm0_predictions |>
    select(iso3, year, predicted = all_of(pred_col))

  all_iso3 <- union(obs$iso3, pred$iso3)

  full_grid <- expand_grid(iso3 = all_iso3, year = years_hist) |>
    left_join(obs,  by = c("iso3", "year")) |>
    left_join(pred, by = c("iso3", "year")) |>
    arrange(iso3, year)

  interp <- full_grid |>
    group_by(iso3) |>
    mutate(reported_interp = zoo::na.approx(reported, rule = 1, na.rm = FALSE)) |>
    ungroup()

  result <- interp |>
    group_by(iso3) |>
    mutate(
      first_year = suppressWarnings(min(year[!is.na(reported)], na.rm = TRUE)),
      last_year  = suppressWarnings(max(year[!is.na(reported)], na.rm = TRUE)),
      first_year = if_else(is.infinite(first_year), NA_real_, first_year),
      last_year  = if_else(is.infinite(last_year),  NA_real_, last_year),

      anchor_rep_first  = if_else(!is.na(first_year), reported[year == first_year][1], NA_real_),
      anchor_pred_first = if_else(!is.na(first_year), predicted[year == first_year][1], NA_real_),
      anchor_rep_last   = if_else(!is.na(last_year),  reported[year == last_year][1], NA_real_),
      anchor_pred_last  = if_else(!is.na(last_year),  predicted[year == last_year][1], NA_real_),

      filled = case_when(
        !is.na(reported_interp) ~ reported_interp,
        !is.na(anchor_rep_first) & !is.na(anchor_pred_first) &
          anchor_pred_first > 0.01 & !is.na(predicted) & year < first_year ~
          anchor_rep_first * (predicted / anchor_pred_first),
        !is.na(anchor_rep_last) & !is.na(anchor_pred_last) &
          anchor_pred_last > 0.01 & !is.na(predicted) & year > last_year ~
          anchor_rep_last * (predicted / anchor_pred_last),
        (year < first_year | year > last_year) & !is.na(predicted) ~ predicted,
        TRUE ~ reported_interp
      ),

      fill_type = case_when(
        !is.na(reported) ~ "reported",
        is.na(reported) & !is.na(reported_interp) ~ "interpolated",
        is.na(reported) & !is.na(predicted) & !is.na(anchor_rep_first) &
          !is.na(anchor_pred_first) & anchor_pred_first > 0.01 &
          year < first_year ~ "xgb_extrap_before",
        is.na(reported) & !is.na(predicted) & !is.na(anchor_rep_last) &
          !is.na(anchor_pred_last) & anchor_pred_last > 0.01 &
          year > last_year ~ "xgb_extrap_after",
        is.na(reported) & !is.na(predicted) ~ "xgb_direct",
        TRUE ~ NA_character_
      ),

      filled = clamp_pct(filled)
    ) |>
    ungroup()

  result <- result |>
    group_by(iso3) |>
    mutate(
      has_any_reported = any(!is.na(reported)),
      filled = case_when(
        has_any_reported  ~ filled,
        !has_any_reported & !is.na(predicted) ~ clamp_pct(predicted),
        TRUE ~ filled
      ),
      fill_type = case_when(
        has_any_reported  ~ fill_type,
        !has_any_reported & !is.na(predicted) ~ "xgb_direct",
        TRUE ~ fill_type
      )
    ) |>
    ungroup() |>
    select(iso3, year, reported, predicted, filled, fill_type)

  n_filled <- sum(!is.na(result$filled))
  message(sprintf("  Filled %d / %d country-years", n_filled, nrow(result)))
  result
}

interp_list <- list()
for (indicator in c("s_bas", "s_od")) {
  interp_list[[indicator]] <- interp_extrap_indicator(indicator, hist_sanit)
}

# ---- 6. s_lim & s_unimp from XGBoost, scaled ----
message("\n== Scaling s_lim & s_unimp from XGBoost predictions ==")

scaling_df <- interp_list[["s_bas"]] |>
  select(iso3, year, s_bas_filled = filled) |>
  left_join(
    interp_list[["s_od"]] |>
      select(iso3, year, s_od_filled = filled),
    by = c("iso3", "year")
  ) |>
  left_join(
    adm0_predictions |>
      select(iso3, year, pred_s_lim, pred_s_unimp),
    by = c("iso3", "year")
  ) |>
  mutate(
    residual = 100 - clamp_pct(s_bas_filled) - clamp_pct(s_od_filled),
    residual = pmax(0, residual),

    xgb_sum = pred_s_lim + pred_s_unimp,

    s_lim_filled = if_else(
      xgb_sum > 0,
      residual * (pred_s_lim / xgb_sum),
      residual * 0.5
    ),
    s_unimp_filled = if_else(
      xgb_sum > 0,
      residual * (pred_s_unimp / xgb_sum),
      residual * 0.5
    ),

    s_lim_filled   = clamp_pct(s_lim_filled),
    s_unimp_filled = clamp_pct(s_unimp_filled)
  )

for (indicator in c("s_lim", "s_unimp")) {
  pred_col   <- paste0("pred_", indicator)
  filled_col <- paste0(indicator, "_filled")

  interp_list[[indicator]] <- scaling_df |>
    left_join(
      adm0_predictions |> select(iso3, year, predicted = all_of(pred_col)),
      by = c("iso3", "year")
    ) |>
    transmute(
      iso3,
      year,
      reported  = NA_real_,
      predicted,
      filled    = .data[[filled_col]],
      fill_type = if_else(!is.na(filled), "xgb_scaled", NA_character_)
    )
}

# ---- 7. Save outputs ----
message("\n== Saving individual indicator files ==")
for (indicator in sanit_indicators) {
  out_path <- sprintf("results_sanit/%s_adm0_interpPredExtrap_hist.csv", indicator)
  write_csv(interp_list[[indicator]], out_path)
  message("  Saved: ", out_path)
}

message("\n== Building combined ADM0 historical file ==")
country_lookup <- adm0_predictors |>
  distinct(iso3, country_name = country_name_gadm) |>
  bind_rows(hist_sanit |> distinct(iso3, country_name)) |>
  filter(!is.na(iso3), !is.na(country_name)) |>
  distinct(iso3, .keep_all = TRUE)

combined <- expand_grid(
  iso3 = sort(unique(union(hist_sanit$iso3, adm0_predictions$iso3))),
  year = years_hist
) |>
  left_join(country_lookup, by = "iso3")

for (indicator in sanit_indicators) {
  combined <- combined |>
    left_join(
      interp_list[[indicator]] |>
        rename(
          !!paste0(indicator, "_reported")  := reported,
          !!paste0(indicator, "_predicted") := predicted,
          !!paste0(indicator, "_filled")    := filled,
          !!paste0(indicator, "_fill_type") := fill_type
        ),
      by = c("iso3", "year")
    )
}

combined <- combined |>
  mutate(
    total_filled = s_bas_filled + s_lim_filled +
      s_unimp_filled + s_od_filled
  )

n_bad <- sum(!near(combined$total_filled, 100, tol = 0.5) &
               !is.na(combined$total_filled))
if (n_bad > 0) {
  message(sprintf("WARNING: %d rows where sum of filled indicators != 100 (tol=0.5)", n_bad))
} else {
  message("  Sum constraint OK: all filled rows sum to ~100")
}

combined <- combined |> select(-total_filled)

write_csv(combined, "results_sanit/JMP_adm0_sanitation_interpPredExtrap_hist.csv")
message(sprintf("  Saved: results_sanit/JMP_adm0_sanitation_interpPredExtrap_hist.csv (%d rows)",
                nrow(combined)))

message("\n== Coverage summary ==")
for (indicator in sanit_indicators) {
  filled_col <- paste0(indicator, "_filled")
  n_filled <- sum(!is.na(combined[[filled_col]]))
  n_cntry  <- n_distinct(combined$iso3[!is.na(combined[[filled_col]])])
  message(sprintf("  %s: %d country-years | %d countries",
                  indicator, n_filled, n_cntry))
}

message("\n== Fill type summary ==")
for (indicator in sanit_indicators) {
  fill_type_col <- paste0(indicator, "_fill_type")
  message("  ", indicator, ":")
  combined |>
    filter(!is.na(.data[[fill_type_col]])) |>
    count(.data[[fill_type_col]], name = "n") |>
    arrange(desc(n)) |>
    mutate(msg = paste0("    ", .data[[fill_type_col]], ": ", n)) |>
    pull(msg) |>
    walk(message)
}

message("\nDone.")
