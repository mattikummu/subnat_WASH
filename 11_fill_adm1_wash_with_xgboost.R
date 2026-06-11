# 11_fill_adm1_wash_with_xgboost.R
# ---------------------------------------------------------------------------
# Fill missing ADM1 (subnational) JMP water values using existing XGBoost
# ensemble models (script 06) with one global ADM1 geometry base.
#
# Workflow:
#   1) Build global ADM1 geometry = best division + GADM fallback countries
#   2) Ensure all polygons have stable adm_id
#   3) Extract predictors for all ADM1 polygons at once (population-weighted)
#   4) Predict all indicators with saved ensembles
#   5) Interpolate/extrapolate per adm_id (same logic as ADM0 step)
#
# Inputs:
#   results/JMP_adm0_water_interpPredExtrap.csv
#   results/wash_best_division_data.csv
#   results/wash_best_division.gpkg
#   /Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg (ADM_1)
#   results/xgboost_wash_ensemble_*.rds
#
# Outputs:
#   results/wash_adm1_global.gpkg
#   results/adm1_global_lookup.csv
#   results/xgb_raster_adm1_1990_2023.csv
#   results/wash_adm1_xgboost_predictions.csv
#   results/JMP_adm1_water_interpPredExtrap.csv
#   results/{indicator}_adm1_interpPredExtrap.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(tidymodels)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(zoo)
  library(countrycode)
})

# Source cache prompt helper
source("cache_prompt.R")

# Source bias correction helper
source("apply_bias_correction.R")

dir.create("results", showWarnings = FALSE)

wash_indicators <- c("Basic", "Limited", "Unimproved", "Surface.water")
predictor_cols <- c(
  "travel_time", "poverty", "msch", "gnic", "gini",
  "urban", "aridity", "gdi", "corruption", "gw_depth"
)
years_pred <- 1990:2023
gadm_path <- "/Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg"

clamp_pct <- function(x) {
  pmax(0, pmin(100, x))
}

mean_or_na <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    mean(x, na.rm = TRUE)
  }
}

message("== Loading ADM0/ADM1 inputs ==")
adm0_filled <- read_csv("results/JMP_adm0_water_interpPredExtrap.csv", show_col_types = FALSE)
adm1_best_obs <- read_csv("results/wash_best_division_data.csv", show_col_types = FALSE) |>
  filter(year %in% years_pred) |>
  group_by(adm_id, year) |>
  summarise(
    iso3 = first(iso3),
    region_name = first(region_name),
    Basic = mean_or_na(Basic),
    Limited = mean_or_na(Limited),
    Unimproved = mean_or_na(Unimproved),
    Surface.water = mean_or_na(Surface.water),
    .groups = "drop"
  )
best_poly <- vect("results/wash_best_division.gpkg")

message(sprintf("  best polygons: %d | best obs rows: %d", nrow(best_poly), nrow(adm1_best_obs)))

# ---------------------------------------------------------------------------
# 1) Load global ADM1 geometry and predictors
# ---------------------------------------------------------------------------
message("== Loading global ADM1 geometry and predictors (from script 07) ==")
adm1_global <- vect("results/wash_adm1_global.gpkg")
adm1_ids <- as.data.frame(adm1_global) |>
  distinct(adm_id) |>
  pull(adm_id)

adm1_predictors <- read_csv("results/global_adm1_predictors_1990_2023.csv",
                            show_col_types = FALSE) |>
  filter(adm_id %in% adm1_ids) |>
  filter(year %in% years_pred) |>
  arrange(adm_id, year)

message(sprintf("  Loaded %d ADM1 predictor rows (%d adm_ids × %d years)",
                nrow(adm1_predictors),
                n_distinct(adm1_predictors$adm_id),
                n_distinct(adm1_predictors$year)))

# ---------------------------------------------------------------------------
# 3) Predict all indicators
# ---------------------------------------------------------------------------
message("== Predicting ADM1 values with existing XGBoost ensembles ==")

attrs <- as.data.frame(adm1_global) |>
  transmute(poly_idx = row_number(), iso3, region_name, adm_id, source)

pred_input <- adm1_predictors |>
  select(adm_id, iso3, region_name, source, year, pop, all_of(predictor_cols))

pred_cache_path <- "results/wash_adm1_xgboost_predictions.csv"
pred_meta_path <- "results/wash_adm1_xgboost_predictions_meta.rds"

pred_key <- pred_input |>
  arrange(adm_id, year) |>
  transmute(adm_id = as.numeric(adm_id), year = as.numeric(year))

pred_sig <- list(
  n_rows = nrow(pred_key),
  sum_adm_id = sum(pred_key$adm_id, na.rm = TRUE),
  sum_year = sum(pred_key$year, na.rm = TRUE),
  sum_cross = sum(pred_key$adm_id * pred_key$year, na.rm = TRUE)
)

predict_ensemble <- function(ensemble_path, new_data, pred_cols) {
  if (!file.exists(ensemble_path)) return(rep(NA_real_, nrow(new_data)))

  ens <- readRDS(ensemble_path)
  out <- rep(NA_real_, nrow(new_data))
  # Only skip rows where ALL predictors are NA; XGBoost handles partial NAs
  ok <- rowSums(!is.na(new_data[, pred_cols, drop = FALSE])) > 0

  if (any(ok)) {
    pred_matrix <- map(ens$models, function(mod) {
      predict(mod, new_data = new_data[ok, pred_cols, drop = FALSE])$.pred
    }) |>
      do.call(cbind, args = _)

    ensemble_pred_raw <- apply(pred_matrix, 1, median)
    out[ok] <- apply_bias_correction(ensemble_pred_raw, ens$bias_correction)
  }

  clamp_pct(out)
}

use_pred_cache <- use_existing_cache(pred_cache_path, "ADM1 prediction cache")

if (use_pred_cache && file.exists(pred_meta_path)) {
  cached_sig <- readRDS(pred_meta_path)
  if (isTRUE(all.equal(cached_sig, pred_sig))) {
    pred_tbl <- read_csv(pred_cache_path, show_col_types = FALSE)
    message(sprintf("Loaded: %s (%d rows)", pred_cache_path, nrow(pred_tbl)))
  } else {
    message("Prediction cache signature mismatch; re-predicting.")
    pred_tbl <- pred_input |>
      select(iso3, region_name, adm_id, source, year)

    for (ind in wash_indicators) {
      message(sprintf("  Predicting %s ...", ind))
      ens_path <- sprintf("results/xgboost_wash_ensemble_%s.rds", ind)
      pred_tbl[[paste0("pred_", ind)]] <- predict_ensemble(ens_path, pred_input, predictor_cols)
    }

    write_csv(pred_tbl, pred_cache_path)
    saveRDS(pred_sig, pred_meta_path)
    message(sprintf("Saved: %s (%d rows)", pred_cache_path, nrow(pred_tbl)))
  }
} else if (use_pred_cache && !file.exists(pred_meta_path)) {
  message("Prediction cache metadata missing; re-predicting.")
  pred_tbl <- pred_input |>
    select(iso3, region_name, adm_id, source, year)

  for (ind in wash_indicators) {
    message(sprintf("  Predicting %s ...", ind))
    ens_path <- sprintf("results/xgboost_wash_ensemble_%s.rds", ind)
    pred_tbl[[paste0("pred_", ind)]] <- predict_ensemble(ens_path, pred_input, predictor_cols)
  }

  write_csv(pred_tbl, pred_cache_path)
  saveRDS(pred_sig, pred_meta_path)
  message(sprintf("Saved: %s (%d rows)", pred_cache_path, nrow(pred_tbl)))
} else {
  pred_tbl <- pred_input |>
    select(iso3, region_name, adm_id, source, year)

  for (ind in wash_indicators) {
    message(sprintf("  Predicting %s ...", ind))
    ens_path <- sprintf("results/xgboost_wash_ensemble_%s.rds", ind)
    pred_tbl[[paste0("pred_", ind)]] <- predict_ensemble(ens_path, pred_input, predictor_cols)
  }

  write_csv(pred_tbl, pred_cache_path)
  saveRDS(pred_sig, pred_meta_path)
  message(sprintf("Saved: %s (%d rows)", pred_cache_path, nrow(pred_tbl)))
}

# ---------------------------------------------------------------------------
# 4) Interpolate + extrapolate per ADM1
# ---------------------------------------------------------------------------
message("== Interpolating and extrapolating ADM1 series ==")

base_grid <- expand_grid(
  adm_id = unique(attrs$adm_id),
  year = years_pred
) |>
  left_join(attrs |> select(adm_id, iso3, region_name, source) |> distinct(), by = "adm_id")

fill_one_indicator <- function(indicator) {
  message(sprintf("  Processing: %s", indicator))

  reported_col <- indicator
  pred_col <- paste0("pred_", indicator)

  tmp <- base_grid |>
    left_join(
      adm1_best_obs |>
        select(adm_id, year, reported = all_of(reported_col)),
      by = c("adm_id", "year")
    ) |>
    left_join(
      pred_tbl |>
        select(adm_id, year, predicted = all_of(pred_col)),
      by = c("adm_id", "year")
    ) |>
    arrange(adm_id, year)

  filled <- tmp |>
    group_by(adm_id) |>
    group_modify(~{
      d <- .x |>
        group_by(year) |>
        summarise(
          reported = mean_or_na(reported),
          predicted = mean_or_na(predicted),
          .groups = "drop"
        ) |>
        arrange(year)
      repv <- d$reported
      predv <- d$predicted

      # Only interpolate if there are at least 2 non-NA observed values
      n_obs <- sum(!is.na(repv))
      if (n_obs >= 2) {
        interpv <- na.approx(repv, x = d$year, na.rm = FALSE, rule = 1, f = 0.5)
      } else {
        interpv <- repv  # No interpolation if < 2 observations
      }


      # Initialize filled values with reported; keep track of fill types
      filledv <- repv
      ftype <- ifelse(!is.na(repv), "reported", NA_character_)

      # Only apply interpolation if we have 2+ observations AND computed interpv differs from repv
      if (n_obs >= 2) {
        idx_interp <- is.na(filledv) & !is.na(interpv)
        if (any(idx_interp)) {
          filledv[idx_interp] <- interpv[idx_interp]
          ftype[idx_interp] <- "interpolated"
        }
      }
      obs_idx <- which(!is.na(repv))
      if (length(obs_idx) > 0) {
        first_i <- obs_idx[1]
        last_i <- obs_idx[length(obs_idx)]

        first_pred <- predv[first_i]
        last_pred <- predv[last_i]

        idx_before <- which(
          is.na(filledv) & seq_along(filledv) < first_i &
            !is.na(predv) & !is.na(first_pred) & first_pred > 0.01
        )
        if (length(idx_before) > 0) {
          filledv[idx_before] <- repv[first_i] * (predv[idx_before] / first_pred)
          ftype[idx_before] <- "xgb_extrap_before"
        }

        idx_after <- which(
          is.na(filledv) & seq_along(filledv) > last_i &
            !is.na(predv) & !is.na(last_pred) & last_pred > 0.01
        )
        if (length(idx_after) > 0) {
          filledv[idx_after] <- repv[last_i] * (predv[idx_after] / last_pred)
          ftype[idx_after] <- "xgb_extrap_after"
        }
      }

      idx_direct <- which(is.na(filledv) & !is.na(predv))
      if (length(idx_direct) > 0) {
        filledv[idx_direct] <- predv[idx_direct]
        ftype[idx_direct] <- "xgb_direct"
      }

      tibble(
        year = d$year,
        reported = repv,
        predicted = predv,
        filled = clamp_pct(filledv),
        fill_type = ftype
      )
    }) |>
    ungroup()

  tmp |>
    select(iso3, adm_id, region_name, source, year) |>
    bind_cols(
      filled |>
        select(reported, predicted, filled, fill_type) |>
        rename(
          !!paste0(indicator, "_reported") := reported,
          !!paste0(indicator, "_predicted") := predicted,
          !!paste0(indicator, "_filled") := filled,
          !!paste0(indicator, "_fill_type") := fill_type
        )
    )
}

interp_list <- map(wash_indicators, fill_one_indicator)
names(interp_list) <- wash_indicators

combined <- reduce(
  interp_list,
  ~left_join(.x, .y, by = c("iso3", "adm_id", "region_name", "source", "year"))
)

write_csv(combined, "results/JMP_adm1_water_interpPredExtrap.csv")
message(sprintf("Saved: results/JMP_adm1_water_interpPredExtrap.csv (%d rows)", nrow(combined)))

walk(wash_indicators, function(ind) {
  out_path <- sprintf("results/%s_adm1_interpPredExtrap.csv", ind)
  one <- combined |>
    select(iso3, adm_id, region_name, source, year,
           all_of(paste0(ind, "_reported")),
           all_of(paste0(ind, "_predicted")),
           all_of(paste0(ind, "_filled")),
           all_of(paste0(ind, "_fill_type")))
  write_csv(one, out_path)
  message(sprintf("  Saved: %s", out_path))
})

message("\n== Coverage summary ==")
walk(wash_indicators, function(ind) {
  filled_col <- paste0(ind, "_filled")
  n_filled <- sum(!is.na(combined[[filled_col]]))
  n_cntry <- n_distinct(combined$iso3[!is.na(combined[[filled_col]])])
  message(sprintf("  %s: %d adm1-years | %d countries", ind, n_filled, n_cntry))
})

message("\nDone.")
