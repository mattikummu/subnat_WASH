# 50_fill_adm0_sanit_with_xgboost.R
# ---------------------------------------------------------------------------
# Fill missing ADM0 (country-level) JMP sanitation values using the
# XGBoost ensemble models trained in script 46.
#
# Method:
#   1. Extract ADM0 predictors directly from country polygons using
#      population-weighted raster extraction (8 predictors, no aridity/gw_depth)
#   2. Predict ADM0 indicator values with the saved XGBoost ensembles
#   3. Linearly interpolate between observed ADM0 years
#   4. Extrapolate outside observed windows using XGBoost growth ratios
#   5. Fall back to direct XGBoost predictions where no observations exist
#
# Inputs:
#   results_sanit/JMP_adm0_sanitation.csv
#   /Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg (ADM_0 layer)
#   Raw predictor rasters used in script 45
#   results_sanit/xgboost_sanit_ensemble_*.rds
#
# Outputs:
#   results_sanit/xgb_raster_adm0_1990_2023.csv
#   results_sanit/sanit_adm0_xgboost_predictions.csv
#   results_sanit/JMP_adm0_sanitation_interpPredExtrap.csv
#   results_sanit/{indicator}_adm0_interpPredExtrap.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(tidymodels)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(zoo)
})

# Source cache prompt helper
source("cache_prompt.R")

# Source bias correction helper
source("apply_bias_correction.R")

dir.create("results_sanit", showWarnings = FALSE)

sanit_indicators <- c("s_bas", "s_lim", "s_unimp", "s_od")
predictor_cols <- c(
  "travel_time", "poverty", "msch", "gnic", "gini",
  "urban", "gdi", "corruption"
)
years_in <- 1990:2023
years_pred <- 1990:2023

clamp_pct <- function(x) {
  pmax(0, pmin(100, x))
}

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  w_sum <- sum(w[ok], na.rm = TRUE)
  if (is.na(w_sum) || w_sum <= 0) return(NA_real_)
  weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

message("== Loading ADM0 observed data ==")
adm0_obs <- read_csv("results_sanit/JMP_adm0_sanitation.csv", show_col_types = FALSE) |>
  filter(!is.na(iso3), year %in% years_in)

message("== Loading global ADM0 predictors (from script 07) ==")
adm0_predictors <- read_csv("results/global_adm0_predictors_1990_2023.csv",
                            show_col_types = FALSE) |>
  filter(year %in% years_pred) |>
  select(iso3, country_name_gadm, year, pop, all_of(predictor_cols)) |>
  arrange(iso3, year)

message(sprintf("  Loaded %d ADM0 predictor rows (%d countries × %d years)",
                nrow(adm0_predictors),
                n_distinct(adm0_predictors$iso3),
                n_distinct(adm0_predictors$year)))

predict_ensemble <- function(ensemble_path, new_data) {
  ens <- readRDS(ensemble_path)

  out <- rep(NA_real_, nrow(new_data))
  keep <- complete.cases(new_data[, predictor_cols])

  if (any(keep)) {
    pred_matrix <- map(ens$models, function(mod) {
      predict(mod, new_data = new_data[keep, predictor_cols, drop = FALSE])$.pred
    }) |> do.call(cbind, args = _)

    ensemble_pred_raw <- apply(pred_matrix, 1, median)
    out[keep] <- apply_bias_correction(ensemble_pred_raw, ens$bias_correction)
  }

  clamp_pct(out)
}

message("== Predicting ADM0 values with existing XGBoost ensembles ==")
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

write_csv(adm0_predictions, "results_sanit/sanit_adm0_xgboost_predictions.csv")
message(sprintf("Saved: results_sanit/sanit_adm0_xgboost_predictions.csv (%d rows)",
                nrow(adm0_predictions)))

interp_extrap_indicator <- function(indicator) {
  message("\n--- Processing: ", indicator, " ---")

  pred_col <- paste0("pred_", indicator)

  obs <- adm0_obs |>
    select(iso3, year, reported = all_of(indicator))

  pred <- adm0_predictions |>
    select(iso3, year, predicted = all_of(pred_col))

  all_iso3 <- union(obs$iso3, pred$iso3)

  full_grid <- expand_grid(iso3 = all_iso3, year = years_in) |>
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
      last_year  = if_else(is.infinite(last_year), NA_real_, last_year),

      anchor_rep_first  = if_else(!is.na(first_year), reported[year == first_year][1], NA_real_),
      anchor_pred_first = if_else(!is.na(first_year), predicted[year == first_year][1], NA_real_),
      anchor_rep_last   = if_else(!is.na(last_year), reported[year == last_year][1], NA_real_),
      anchor_pred_last  = if_else(!is.na(last_year), predicted[year == last_year][1], NA_real_),

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
          !is.na(anchor_pred_first) & anchor_pred_first > 0.01 & year < first_year ~ "xgb_extrap_before",
        is.na(reported) & !is.na(predicted) & !is.na(anchor_rep_last) &
          !is.na(anchor_pred_last) & anchor_pred_last > 0.01 & year > last_year ~ "xgb_extrap_after",
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
        has_any_reported ~ filled,
        !has_any_reported & !is.na(predicted) ~ clamp_pct(predicted),
        TRUE ~ filled
      ),
      fill_type = case_when(
        has_any_reported ~ fill_type,
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
for (indicator in sanit_indicators) {
  interp_list[[indicator]] <- interp_extrap_indicator(indicator)
}

message("\n== Saving individual indicator files ==")
for (indicator in sanit_indicators) {
  out_path <- sprintf("results_sanit/%s_adm0_interpPredExtrap.csv", indicator)
  write_csv(interp_list[[indicator]], out_path)
  message("  Saved: ", out_path)
}

message("\n== Building combined ADM0 file ==")
country_lookup <- adm0_predictors |>
  distinct(iso3, country_name = country_name_gadm) |>
  bind_rows(adm0_obs |> distinct(iso3, country_name)) |>
  filter(!is.na(iso3), !is.na(country_name)) |>
  distinct(iso3, .keep_all = TRUE)

combined <- expand_grid(
  iso3 = sort(unique(union(adm0_obs$iso3, adm0_predictions$iso3))),
  year = years_in
) |>
  left_join(country_lookup, by = "iso3")

for (indicator in sanit_indicators) {
  combined <- combined |>
    left_join(
      interp_list[[indicator]] |>
        rename(
          !!paste0(indicator, "_reported") := reported,
          !!paste0(indicator, "_predicted") := predicted,
          !!paste0(indicator, "_filled") := filled,
          !!paste0(indicator, "_fill_type") := fill_type
        ),
      by = c("iso3", "year")
    )
}

write_csv(combined, "results_sanit/JMP_adm0_sanitation_interpPredExtrap.csv")
message(sprintf("Saved: results_sanit/JMP_adm0_sanitation_interpPredExtrap.csv (%d rows)",
                nrow(combined)))

message("\n== Coverage summary ==")
for (indicator in sanit_indicators) {
  filled_col <- paste0(indicator, "_filled")
  n_filled <- sum(!is.na(combined[[filled_col]]))
  n_cntry  <- n_distinct(combined$iso3[!is.na(combined[[filled_col]])])
  message(sprintf("  %s: %d country-years | %d countries", indicator, n_filled, n_cntry))
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
