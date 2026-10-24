# 09b_prepare_adm0_historical.R
# ---------------------------------------------------------------------------
# Prepare ADM0 (country-level) historical water data (1990-2015) from the
# older JMP estimates (1980-2015).
#
# This file contains only Basic (= wat_total_improved_t) and
# Surface.water (= wat_surface_water_t).  For Limited and Unimproved we
# rely entirely on XGBoost predictions, scaled so that
#   Limited + Unimproved = 100 − (Basic + Surface.water)
#
# Method:
#   1. Read data_in/estimates_1980-2015_v1.xlsx, map to ISO3
#   2. Extract ADM0 predictors (population-weighted) for 1990-2015
#   3. Predict all four indicators with saved XGBoost ensembles
#   4. Linearly interpolate observed Basic & Surface.water
#   5. Extrapolate observed Basic & Surface.water using XGBoost growth ratios
#   6. For Limited & Unimproved: use XGBoost predictions, then rescale
#      so their sum = 100 − (Basic + Surface.water)
#
# Inputs:
#   data_in/estimates_1980-2015_v1.xlsx
#   /Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg  (ADM_0 layer)
#   Predictor rasters (same as script 09)
#   results/xgboost_wash_ensemble_*.rds
#
# Outputs:
#   results/xgb_raster_adm0_historical.csv          (cached predictor extraction)
#   results/wash_adm0_xgboost_predictions_hist.csv   (raw xgboost predictions)
#   results/JMP_adm0_water_historical.csv            (cleaned observed data)
#   results/JMP_adm0_water_interpPredExtrap_hist.csv (combined final output)
#   results/{indicator}_adm0_interpPredExtrap_hist.csv
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

dir.create("results", showWarnings = FALSE)

wash_indicators <- c("Basic", "Limited", "Unimproved", "Surface.water")
predictor_cols  <- c(
  "travel_time", "poverty", "msch", "gnic", "gini",
  "urban", "aridity", "gdi", "corruption", "gw_depth"
)

years_hist <- 1990:2015

clamp_pct <- function(x) pmax(0, pmin(100, x))

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  w_sum <- sum(w[ok], na.rm = TRUE)
  if (is.na(w_sum) || w_sum <= 0) return(NA_real_)
  weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

# ---- 1. Read & clean historical JMP estimates ----
message("== Reading historical JMP estimates ==")
xlsx_path <- "data_in/estimates_1980-2015_v1.xlsx"
if (!file.exists(xlsx_path)) stop("Input file not found: ", xlsx_path)

hist_raw <- read_excel(xlsx_path, sheet = "Sheet1")

hist_water <- hist_raw |>
  transmute(
    country_file,
    country_clean = country_file |>
      str_replace("_15$", "") |>
      str_replace_all("_", " ") |>
      str_to_title(),
    year = as.numeric(year),
    Basic         = as.numeric(wat_total_improved_t),
    Surface.water = as.numeric(wat_surface_water_t)
  ) |>
  mutate(iso3 = countrycode(country_clean, "country.name", "iso3c", warn = FALSE)) |>
  filter(!is.na(iso3), year >= 1990, year <= 2015) |>
  select(iso3, country_name = country_clean, year, Basic, Surface.water)

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
                nrow(hist_water), n_distinct(hist_water$iso3), n_unmapped))

write_csv(hist_water, "results/JMP_adm0_water_historical.csv")
message("  Saved: results/JMP_adm0_water_historical.csv")

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

for (indicator in wash_indicators) {
  ens_path <- sprintf("results/xgboost_wash_ensemble_%s.rds", indicator)
  if (!file.exists(ens_path)) {
    stop("Required ensemble model not found: ", ens_path,
         "\nRun 07_train_xgboost_wash.R first.")
  }
  message("  Predicting ", indicator, " ...")
  adm0_predictions[[paste0("pred_", indicator)]] <-
    predict_ensemble(ens_path, adm0_predictions)
}

write_csv(adm0_predictions, "results/wash_adm0_xgboost_predictions_hist.csv")
message(sprintf("  Saved: results/wash_adm0_xgboost_predictions_hist.csv (%d rows)",
                nrow(adm0_predictions)))

# ---- 5. Interpolate & extrapolate Basic and Surface.water ----
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

  # Linear interpolation between observed points
  interp <- full_grid |>
    group_by(iso3) |>
    mutate(reported_interp = zoo::na.approx(reported, rule = 1, na.rm = FALSE)) |>
    ungroup()

  # Extrapolation using XGBoost growth ratios (anchored to nearest observed)
  result <- interp |>
    group_by(iso3) |>
    mutate(
      first_year = suppressWarnings(min(year[!is.na(reported)], na.rm = TRUE)),
      last_year  = suppressWarnings(max(year[!is.na(reported)], na.rm = TRUE)),

      first_year = if_else(is.infinite(first_year), NA_real_, first_year),
      last_year  = if_else(is.infinite(last_year),  NA_real_, last_year),

      anchor_rep_first  = if_else(!is.na(first_year),
                                  reported[year == first_year][1], NA_real_),
      anchor_pred_first = if_else(!is.na(first_year),
                                  predicted[year == first_year][1], NA_real_),
      anchor_rep_last   = if_else(!is.na(last_year),
                                  reported[year == last_year][1], NA_real_),
      anchor_pred_last  = if_else(!is.na(last_year),
                                  predicted[year == last_year][1], NA_real_),

      filled = case_when(
        !is.na(reported_interp) ~ reported_interp,

        # Extrapolate backwards using ratio to anchor
        !is.na(anchor_rep_first) & !is.na(anchor_pred_first) &
          anchor_pred_first > 0.01 & !is.na(predicted) & year < first_year ~
          anchor_rep_first * (predicted / anchor_pred_first),

        # Extrapolate forwards using ratio to anchor
        !is.na(anchor_rep_last) & !is.na(anchor_pred_last) &
          anchor_pred_last > 0.01 & !is.na(predicted) & year > last_year ~
          anchor_rep_last * (predicted / anchor_pred_last),

        # Fallback: direct XGBoost prediction
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

  # Countries with no observations at all: fall back to direct XGBoost
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

# Interpolate/extrapolate only the two observed indicators
interp_list <- list()
for (indicator in c("Basic", "Surface.water")) {
  interp_list[[indicator]] <- interp_extrap_indicator(indicator, hist_water)
}

# ---- 6. Limited & Unimproved from XGBoost, scaled ----
message("\n== Scaling Limited & Unimproved from XGBoost predictions ==")

# Build a table with filled Basic and Surface.water plus raw XGBoost
# predictions for Limited and Unimproved
scaling_df <- interp_list[["Basic"]] |>
  select(iso3, year, Basic_filled = filled) |>
  left_join(
    interp_list[["Surface.water"]] |>
      select(iso3, year, Surfacewater_filled = filled),
    by = c("iso3", "year")
  ) |>
  left_join(
    adm0_predictions |>
      select(iso3, year, pred_Limited, pred_Unimproved),
    by = c("iso3", "year")
  ) |>
  mutate(
    # Residual that Limited + Unimproved must sum to
    residual = 100 - clamp_pct(Basic_filled) - clamp_pct(Surfacewater_filled),
    residual = pmax(0, residual),

    # Raw xgboost sum for Limited + Unimproved
    xgb_sum = pred_Limited + pred_Unimproved,

    # Scale proportionally
    Limited_filled = if_else(
      xgb_sum > 0,
      residual * (pred_Limited / xgb_sum),
      residual * 0.5    # equal split if both xgboost predictions are zero
    ),
    Unimproved_filled = if_else(
      xgb_sum > 0,
      residual * (pred_Unimproved / xgb_sum),
      residual * 0.5
    ),

    # Clamp
    Limited_filled    = clamp_pct(Limited_filled),
    Unimproved_filled = clamp_pct(Unimproved_filled)
  )

# Build interp_list entries for Limited and Unimproved
# (no observed data — all values are xgb-based)
for (indicator in c("Limited", "Unimproved")) {
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
for (indicator in wash_indicators) {
  out_path <- sprintf("results/%s_adm0_interpPredExtrap_hist.csv", indicator)
  write_csv(interp_list[[indicator]], out_path)
  message("  Saved: ", out_path)
}

message("\n== Building combined ADM0 historical file ==")
country_lookup <- adm0_predictors |>
  distinct(iso3, country_name = country_name_gadm) |>
  bind_rows(hist_water |> distinct(iso3, country_name)) |>
  filter(!is.na(iso3), !is.na(country_name)) |>
  distinct(iso3, .keep_all = TRUE)

combined <- expand_grid(
  iso3 = sort(unique(union(hist_water$iso3, adm0_predictions$iso3))),
  year = years_hist
) |>
  left_join(country_lookup, by = "iso3")

for (indicator in wash_indicators) {
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

# Verify constraint: sum of filled values ≈ 100
combined <- combined |>
  mutate(
    total_filled = Basic_filled + Limited_filled +
      Unimproved_filled + Surface.water_filled
  )

n_bad <- sum(!near(combined$total_filled, 100, tol = 0.5) &
               !is.na(combined$total_filled))
if (n_bad > 0) {
  message(sprintf("WARNING: %d rows where sum of filled indicators != 100 (tol=0.5)", n_bad))
} else {
  message("  Sum constraint OK: all filled rows sum to ~100")
}

combined <- combined |> select(-total_filled)

write_csv(combined, "results/JMP_adm0_water_interpPredExtrap_hist.csv")
message(sprintf("  Saved: results/JMP_adm0_water_interpPredExtrap_hist.csv (%d rows)",
                nrow(combined)))

# ---- Summary ----
message("\n== Coverage summary ==")
for (indicator in wash_indicators) {
  filled_col <- paste0(indicator, "_filled")
  n_filled <- sum(!is.na(combined[[filled_col]]))
  n_cntry  <- n_distinct(combined$iso3[!is.na(combined[[filled_col]])])
  message(sprintf("  %s: %d country-years | %d countries",
                  indicator, n_filled, n_cntry))
}

message("\n== Fill type summary ==")
for (indicator in wash_indicators) {
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
