# 46_prepare_xgboost_data_sanit.R
# ---------------------------------------------------------------------------
# Pop-weighted extraction of independent variables for XGBoost modelling
# of subnational SANITATION indicators.
#
# Mirrors script 05 (water) but:
#   - Extracts sanitation indicators (s_bas, s_lim, s_unimp, s_od) instead
#     of water indicators
#   - Excludes aridity and groundwater depth from predictors
#   - Outputs to results_sanit/
#
# Key optimisation vs old version:
#   The polygons (best + second + ADM0) and predictor rasters are IDENTICAL
#   to the water pipeline — only the response (WASH) data differs. So
#   instead of extracting predictors from scratch we reuse the water
#   pipeline's `results/xgb_raster_all_years.csv` for any polygon that
#   appears in both. Only new ADM0 countries unique to sanitation need
#   fresh extraction.
#
# Other changes mirror the water pipeline:
#   - get_adm0_eligible_countries() instead of size-based eligibility
#   - HKG / MAC NA adm_id guard
#   - Drop all-NA response rows in subnational AND ADM0
#
# Outputs:
#   results_sanit/JMP_regions_sanitation.csv
#   results_sanit/sanit_best_division_data.csv
#   results_sanit/sanit_second_division_data.csv
#   results_sanit/xgb_raster_all_years_sanit.csv
#   results_sanit/xgb_training_data_sanit.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(purrr)
  library(tibble)
})

source("cache_prompt.R")

dir.create("results_sanit", showWarnings = FALSE)
source("functions/f_small_country_training.R")

years   <- 2000:2023
n_years <- length(years)

# == 0. Extract sanitation data from JMP subnational Excel files ==============

message("== Extracting sanitation data from JMP subnational files ==")

target_vars_sanit <- tibble(
  var1 = c("Sanitation",  "Sanitation",  "Sanitation ladder", "Sanitation ladder",
           "Sanitation ladder", "Sanitation ladder"),
  var2 = c("Improved",    "Not Improved", "Basic",             "Limited",
           "Unimproved",  "Open defecation")
)

in_dir <- "data_in/JMP_subnational"
files  <- list.files(in_dir, pattern = "\\.xlsx$", full.names = TRUE)

read_country_sanit <- function(path) {
  ds <- read_excel(path, sheet = "data_summary")

  df <- ds |>
    semi_join(target_vars_sanit, by = c("var1", "var2")) |>
    select(iso3, year, data_source, var2,
           matches("^region_\\d+_name$"),
           matches("^region_\\d+$"))

  names_long <- df |>
    select(iso3, year, data_source, var2, matches("^region_\\d+_name$")) |>
    pivot_longer(
      cols      = matches("^region_\\d+_name$"),
      names_to  = "region_num",
      values_to = "region_name",
      names_pattern = "region_(\\d+)_name"
    )

  values_long <- df |>
    select(iso3, year, data_source, var2, matches("^region_\\d+$")) |>
    pivot_longer(
      cols      = matches("^region_\\d+$"),
      names_to  = "region_num",
      values_to = "value",
      names_pattern = "region_(\\d+)"
    )

  left_join(names_long, values_long,
            by = c("iso3", "year", "data_source", "var2", "region_num")) |>
    filter(!is.na(region_name)) |>
    select(-region_num) |>
    mutate(value = as.numeric(value)) |>
    summarise(value = mean(value, na.rm = TRUE),
              .by = c(iso3, year, data_source, region_name, var2)) |>
    mutate(value = if_else(is.nan(value), NA_real_, value)) |>
    pivot_wider(names_from = var2, values_from = value) |>
    select(iso3, region_name, year, data_source,
           any_of(c("Improved", "Not Improved", "Basic",
                    "Limited", "Unimproved", "Open defecation")))
}

sanit_cache <- "results_sanit/JMP_regions_sanitation.csv"
if (file.exists(sanit_cache)) {
  message("  Loading cached sanitation extraction ...")
  all_sanit <- read_csv(sanit_cache, show_col_types = FALSE)
} else {
  all_sanit <- map(files, \(f) {
    tryCatch(
      read_country_sanit(f),
      error = function(e) {
        warning("Failed: ", basename(f), " - ", conditionMessage(e))
        NULL
      }
    )
  }, .progress = TRUE) |>
    list_rbind()

  all_sanit <- all_sanit |>
    rename(
      s_bas   = Basic,
      s_lim   = Limited,
      s_unimp = Unimproved,
      s_od    = `Open defecation`
    )

  n_before <- nrow(all_sanit)
  all_sanit <- all_sanit |>
    filter(
      if_all(c(s_bas, s_lim, s_unimp, s_od),
             \(x) is.na(x) | (x >= 0 & x <= 100))
    )
  n_dropped <- n_before - nrow(all_sanit)
  if (n_dropped > 0) {
    message(sprintf("  Dropped %d rows with indicator values outside [0, 100]", n_dropped))
  }

  # Floor any negative numerical noise to 0 (mirrors water pipeline)
  n_neg <- all_sanit |>
    summarise(across(c(s_bas, s_lim, s_unimp, s_od), ~sum(.x < 0, na.rm = TRUE))) |>
    rowSums()
  if (n_neg > 0) {
    message(sprintf("  Flooring %d negative sanitation values to 0", n_neg))
    all_sanit <- all_sanit |>
      mutate(across(c(s_bas, s_lim, s_unimp, s_od), ~pmax(.x, 0)))
  }

  write_csv(all_sanit, sanit_cache)
  message(sprintf("  Saved: %s (%d rows, %d countries)",
                  sanit_cache, nrow(all_sanit), n_distinct(all_sanit$iso3)))
}

message(sprintf("  Sanitation data: %d rows, %d countries",
                nrow(all_sanit), n_distinct(all_sanit$iso3)))

# == 0b. Map sanitation data to matched boundaries ============================
#
# Uses the same poly_id / name_poly_lookup approach as script 05, so that
# region-name variants (capitalisation, splits, etc.) are resolved via the
# centroid fingerprint rather than a fragile direct region_name join.

message("== Mapping sanitation data to matched boundaries ==")

wash_sv_full <- vect("results/JMP_regions_matched.gpkg")
full_df      <- as.data.frame(wash_sv_full)
full_cents   <- crds(centroids(wash_sv_full))
full_df$poly_id <- paste(round(full_cents[, 1], 3), round(full_cents[, 2], 3))

# Canonicalise boundary_source using the same fingerprint logic as script 05
source_fp <- full_df |>
  group_by(iso3, boundary_source) |>
  summarise(fp = paste(sort(unique(poly_id)), collapse = "|"), .groups = "drop")
src_rows <- full_df |> count(iso3, boundary_source, name = "n_rows")
canon_map <- source_fp |>
  left_join(src_rows, by = c("iso3", "boundary_source")) |>
  mutate(n_rows = coalesce(n_rows, 0L)) |>
  group_by(iso3, fp) |>
  mutate(canon = boundary_source[which.max(n_rows)]) |>
  ungroup() |>
  select(iso3, boundary_source, canon)

full_df <- full_df |>
  left_join(canon_map, by = c("iso3", "boundary_source")) |>
  mutate(boundary_source = coalesce(canon, boundary_source)) |>
  select(-canon)

boundary_meta <- full_df |>
  select(iso3, region_name, year, data_source, survey_type, survey_family,
         boundary_source, match_type) |>
  distinct()

sanit_matched <- all_sanit |>
  inner_join(boundary_meta,
             by = c("iso3", "region_name", "year", "data_source"))

message(sprintf("  Matched sanitation rows: %d / %d (%.1f%%)",
                nrow(sanit_matched), nrow(all_sanit),
                100 * nrow(sanit_matched) / nrow(all_sanit)))

# Load division polygons and compute their poly_id from centroids
best_sv_tmp    <- vect("results/wash_best_division.gpkg")
second_sv_tmp  <- vect("results/wash_second_division.gpkg")
best_attrs_tmp   <- as.data.frame(best_sv_tmp)
second_attrs_tmp <- as.data.frame(second_sv_tmp)
best_cents_tmp   <- crds(centroids(best_sv_tmp))
second_cents_tmp <- crds(centroids(second_sv_tmp))
best_attrs_tmp$poly_id   <- paste(round(best_cents_tmp[, 1], 3),
                                  round(best_cents_tmp[, 2], 3))
second_attrs_tmp$poly_id <- paste(round(second_cents_tmp[, 1], 3),
                                  round(second_cents_tmp[, 2], 3))
rm(best_sv_tmp, second_sv_tmp, best_cents_tmp, second_cents_tmp)

# Build name_poly_lookup: (iso3, boundary_source, region_name) -> poly_id -> adm_id
# Prefers the polygon matched most often (handles split-polygon ambiguity).
name_poly_lookup <- full_df |>
  count(iso3, boundary_source, region_name, poly_id, name = "n_match") |>
  group_by(iso3, boundary_source, region_name) |>
  slice_max(n_match, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(iso3, boundary_source, region_name, poly_id) |>
  left_join(
    bind_rows(
      best_attrs_tmp   |> select(iso3, boundary_source, poly_id, adm_id),
      second_attrs_tmp |> select(iso3, boundary_source, poly_id, adm_id)
    ) |> distinct(iso3, boundary_source, poly_id, .keep_all = TRUE),
    by = c("iso3", "boundary_source", "poly_id")
  )

message(sprintf("  %d name-to-polygon mappings", nrow(name_poly_lookup)))

# Assign adm_id via lookup (robust to name variants)
sanit_matched <- sanit_matched |>
  left_join(
    name_poly_lookup |> select(iso3, boundary_source, region_name, adm_id),
    by = c("iso3", "boundary_source", "region_name")
  )

n_no_adm <- sum(is.na(sanit_matched$adm_id))
if (n_no_adm > 0) {
  message(sprintf("  %d sanit rows have no adm_id after lookup (mostly outside selected best/second boundaries)",
                  n_no_adm))
}

best_bs   <- best_attrs_tmp   |> distinct(iso3, boundary_source)
second_bs <- second_attrs_tmp |> distinct(iso3, boundary_source)

sanit_best   <- sanit_matched |> semi_join(best_bs,   by = c("iso3", "boundary_source"))
sanit_second <- sanit_matched |> semi_join(second_bs, by = c("iso3", "boundary_source"))

n_no_adm_best <- sum(is.na(sanit_best$adm_id))
n_no_adm_second <- sum(is.na(sanit_second$adm_id))
if ((n_no_adm_best + n_no_adm_second) > 0) {
  message(sprintf("  WARNING: %d selected-division rows still have no adm_id (best: %d, second: %d)",
                  n_no_adm_best + n_no_adm_second, n_no_adm_best, n_no_adm_second))
}

write_csv(sanit_best,   "results_sanit/sanit_best_division_data.csv")
write_csv(sanit_second, "results_sanit/sanit_second_division_data.csv")
message(sprintf("  Best division: %d rows | Second division: %d rows",
                nrow(sanit_best), nrow(sanit_second)))
message(sprintf("  NA adm_id — best: %d | second: %d",
                sum(is.na(sanit_best$adm_id)), sum(is.na(sanit_second$adm_id))))

rm(wash_sv_full, full_df, full_cents, boundary_meta,
   best_attrs_tmp, second_attrs_tmp, source_fp, src_rows, canon_map, name_poly_lookup)
gc(verbose = FALSE)

# == 1. Reuse water pipeline's extracted predictors ===========================
#
# The polygons and predictor rasters are identical between water and
# sanitation pipelines. Reuse the water pipeline's extraction directly
# instead of re-running the expensive raster × polygon operation.

message("== Reusing predictor extraction from water pipeline ==")

water_raster_path <- "results/xgb_raster_all_years.csv"
if (!file.exists(water_raster_path)) {
  stop(sprintf(
    "Cannot reuse water predictor extraction: %s does not exist.\n  Run script 05 first.",
    water_raster_path
  ))
}

# The water file contains 10 predictors: travel_time, poverty, msch, gnic,
# gini, urban, aridity, gdi, corruption, gw_depth.
# Sanitation uses 8 (no aridity, no gw_depth) — drop them on read.
all_extracted <- read_csv(water_raster_path, show_col_types = FALSE) |>
  select(division_rank, adm_id, iso3, region_name, boundary_source, poly_id, year,
         pop, poverty, msch, gnic, gini, urban, travel_time, gdi, corruption)

message(sprintf("  Loaded %d rows from water pipeline (%d best + %d second + %d ADM0)",
                nrow(all_extracted),
                sum(all_extracted$division_rank == 1L),
                sum(all_extracted$division_rank == 2L),
                sum(all_extracted$division_rank == 0L)))

write_csv(all_extracted, "results_sanit/xgb_raster_all_years_sanit.csv")
message(sprintf("Saved: results_sanit/xgb_raster_all_years_sanit.csv (%d rows)",
                nrow(all_extracted)))

# == 2. Determine ADM0 eligibility for sanitation =============================
#
# ADM0 eligibility for sanitation may differ from water — a country can
# have JMP ADM0 sanitation data without ADM0 water data, or vice versa.
# We compute eligibility independently here.

predictor_cols_adm0 <- c(
  "poverty", "msch", "gnic", "gini", "urban",
  "travel_time", "gdi", "corruption"
)

countries_with_subnat <- sort(unique(c(sanit_best$iso3, sanit_second$iso3)))

adm0_obs_full <- read_csv("results_sanit/JMP_adm0_sanitation.csv",
                          show_col_types = FALSE)

adm0_eligible <- get_adm0_eligible_countries(
  countries_with_subnat = countries_with_subnat,
  adm0_obs_iso3         = adm0_obs_full$iso3
)

message("\n== Determining ADM0 eligibility (sanitation) ==")
message(sprintf("  Countries with subnational sanitation data: %d",
                length(countries_with_subnat)))
message(sprintf("  Countries with JMP ADM0 sanit obs:           %d",
                n_distinct(adm0_obs_full$iso3)))
message(sprintf("  Excluded as too large:                       %d (%s)",
                length(haqi_large_countries()),
                paste(haqi_large_countries(), collapse = ", ")))
message(sprintf("  ADM0-eligible (sanitation):                  %d",
                length(adm0_eligible)))

# == 2c. Read sanitation-only ADM0 predictors from global lookup =============
message("\n== Reading sanitation ADM0 predictors from global lookup ==")

sanit_pred_cols <- c("poverty", "msch", "gnic", "gini", "urban",
                     "travel_time", "gdi", "corruption")

adm0_pred_global <- read_csv("results/global_adm0_predictors_1990_2023.csv",
                             show_col_types = FALSE) |>
  filter(year %in% years, iso3 %in% adm0_eligible)

adm0_long_extra <- adm0_pred_global |>
  mutate(division_rank   = 0L,
         region_name     = NA_character_,
         boundary_source = "ADM0_national",
         poly_id         = NA_character_) |>
  select(division_rank, adm_id, iso3, region_name, boundary_source, poly_id,
         year, pop, all_of(sanit_pred_cols))

# == 3. Build full predictor pool: water file + global ADM0 rows =============
all_extracted <- all_extracted |>
  filter(division_rank != 0L)

if (nrow(adm0_long_extra) > 0) {
  all_extracted <- bind_rows(all_extracted, adm0_long_extra)
  message(sprintf("\n  Added %d ADM0 rows from global lookup", nrow(adm0_long_extra)))
}

write_csv(all_extracted, "results_sanit/xgb_raster_all_years_sanit.csv")
message(sprintf("Saved (final): results_sanit/xgb_raster_all_years_sanit.csv (%d rows)",
                nrow(all_extracted)))

# == 4. Join predictors with sanitation observations ==========================

message("\n== Joining predictors with sanitation observations ==")

sanit_cols <- c("s_bas", "s_lim", "s_unimp", "s_od")

sanit_all <- bind_rows(
  sanit_best   |> mutate(division_rank = 1L),
  sanit_second |> mutate(division_rank = 2L)
)

# Join raster values on adm_id + year (subnational only here)
training_data <- sanit_all |>
  select(division_rank, adm_id, iso3, region_name, boundary_source,
         year, data_source, survey_type, survey_family,
         all_of(sanit_cols)) |>
  left_join(
    all_extracted |>
      filter(division_rank > 0L) |>
      select(adm_id, year, pop, poverty, msch, gnic, gini,
             urban, travel_time, gdi, corruption),
    by = c("adm_id", "year")
  ) |>
  mutate(poly_id = NA_character_, .after = boundary_source)

n_joined <- sum(!is.na(training_data$pop))
message(sprintf("  %d / %d sanitation rows matched with raster values (%.1f%%)",
                n_joined, nrow(training_data),
                100 * n_joined / nrow(training_data)))

# Drop subnational rows where ALL FOUR sanitation indicators are NA.
# Mirrors the same filter in the water pipeline.
sanit_any_nonNA_subnat <- !is.na(training_data$s_bas) |
                          !is.na(training_data$s_lim) |
                          !is.na(training_data$s_unimp) |
                          !is.na(training_data$s_od)
n_allna_subnat <- sum(!sanit_any_nonNA_subnat)
if (n_allna_subnat > 0) {
  iso_dropped <- sort(unique(training_data$iso3[!sanit_any_nonNA_subnat]))
  message(sprintf("  Dropping %d subnational all-NA sanit rows (countries: %s)",
                  n_allna_subnat, paste(iso_dropped, collapse = ", ")))
  training_data <- training_data[sanit_any_nonNA_subnat, ]
}

# == 5. Append ADM0 sanitation rows ===========================================

message("\n== Appending ADM0 sanitation rows ==")

adm0_obs <- adm0_obs_full |>
  filter(year %in% years, iso3 %in% adm0_eligible)

if (length(adm0_eligible) > 0 && nrow(adm0_obs) > 0) {

  adm0_pred <- all_extracted |>
    filter(division_rank == 0L) |>
    select(adm_id, iso3, year, pop, all_of(predictor_cols_adm0))

  adm0_training <- adm0_obs |>
    left_join(adm0_pred, by = c("iso3", "year")) |>
    transmute(
      division_rank   = 0L,
      adm_id,
      iso3,
      region_name     = country_name,
      boundary_source = "ADM0_national",
      poly_id         = NA_character_,
      year,
      data_source     = "JMP_ADM0",
      survey_type     = NA_character_,
      survey_family   = NA_character_,
      s_bas, s_lim, s_unimp, s_od,
      pop, poverty, msch, gnic, gini, urban, travel_time,
      gdi, corruption
    )

  # Drop ADM0 rows where ALL FOUR sanit indicators are NA
  sanit_any_nonNA_adm0 <- !is.na(adm0_training$s_bas) |
                          !is.na(adm0_training$s_lim) |
                          !is.na(adm0_training$s_unimp) |
                          !is.na(adm0_training$s_od)
  n_allna_adm0 <- sum(!sanit_any_nonNA_adm0)
  if (n_allna_adm0 > 0) {
    iso_adm0_dropped <- sort(unique(adm0_training$iso3[!sanit_any_nonNA_adm0]))
    message(sprintf("  Dropping %d all-NA sanit ADM0 rows (%d countries: %s)",
                    n_allna_adm0,
                    length(iso_adm0_dropped),
                    paste(iso_adm0_dropped, collapse = ", ")))
    adm0_training <- adm0_training[sanit_any_nonNA_adm0, ]
  }

  # Defensive overlap check
  overlap_iso3 <- intersect(unique(training_data$iso3),
                            unique(adm0_training$iso3))
  if (length(overlap_iso3) > 0) {
    stop(sprintf("ADM0/ADM1 overlap detected for %d country(ies): %s",
                 length(overlap_iso3), paste(overlap_iso3, collapse = ", ")))
  }

  n_adm0_joined <- sum(!is.na(adm0_training$pop))
  message(sprintf("  Added ADM0 rows: %d (%d with predictor match, across %d countries)",
                  nrow(adm0_training), n_adm0_joined,
                  n_distinct(adm0_training$iso3)))

  training_data <- bind_rows(training_data, adm0_training)
} else {
  message("  No ADM0 rows to add.")
}

message(sprintf("  Final training rows (ADM1 + ADM0): %d", nrow(training_data)))

write_csv(training_data, "results_sanit/xgb_training_data_sanit.csv")
message(sprintf("Saved: results_sanit/xgb_training_data_sanit.csv (%d rows)",
                nrow(training_data)))

message("\nDone!")
