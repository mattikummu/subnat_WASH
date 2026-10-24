# prepare_xgboost_data.R
# ---------------------------------------------------------------------------
# Pop-weighted extraction of independent variables for XGBoost modelling
# of subnational WASH indicators.
#
#  1. Build 24-layer raster stacks (2000-2023) for each variable
#  2. Pop-weighted extraction for all admin polygons (best + second division)
#  3. Join with WASH survey data
#  4. Save ready-to-use training dataset
#
# Outputs:
#   results/xgb_raster_all_years.csv  -- raster values for ALL polygons x years
#   results/xgb_training_data.csv     -- WASH + raster values (training-ready)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
})

dir.create("results", showWarnings = FALSE)
source("functions/f_small_country_training.R")

years   <- 2000:2023
n_years <- length(years)

# == 1. Load admin polygons ===================================================

message("== Loading admin polygons ==")
best_sv   <- vect("results/wash_best_division.gpkg")
second_sv <- vect("results/wash_second_division.gpkg")

# adm_id is stored as an attribute in the gpkg (written by select_best_divisions.R)
add_poly_id <- function(sv) {
  cents <- crds(centroids(sv))
  df <- as.data.frame(sv)           # includes adm_id from gpkg attributes
  df$poly_idx <- seq_len(nrow(df))
  df$poly_id  <- paste(round(cents[, 1], 3), round(cents[, 2], 3))
  df
}

best_attrs   <- add_poly_id(best_sv)
second_attrs <- add_poly_id(second_sv)
message(sprintf("  Best: %d polygons | Second: %d polygons",
                nrow(best_sv), nrow(second_sv)))

message(sprintf("  adm_id range: %d – %d",
                min(best_attrs$adm_id), max(best_attrs$adm_id)))

# == 1b. Determine ADM0 eligibility and load ADM0 polygons EARLY ==============
# We load ADM0 polygons here (before extraction) so that all three polygon
# sets (best + second + adm0) can be extracted in a single pass — much
# faster than three separate extraction passes.

message("\n== Determining ADM0 eligibility ==")

source("functions/f_small_country_training.R")

predictor_cols_adm0 <- c(
  "poverty", "msch", "gnic", "gini", "urban",
  "travel_time", "aridity", "gdi", "corruption", "gw_depth"
)

# Countries already in subnational training set come from polygon attributes
countries_with_subnat <- sort(unique(c(best_attrs$iso3, second_attrs$iso3)))

adm0_obs_full <- read_csv("results/JMP_adm0_water.csv", show_col_types = FALSE)

adm0_eligible <- get_adm0_eligible_countries(
  countries_with_subnat = countries_with_subnat,
  adm0_obs_iso3         = adm0_obs_full$iso3
)

message(sprintf("  Countries with subnational data: %d", length(countries_with_subnat)))
message(sprintf("  Countries with JMP ADM0 obs:     %d", n_distinct(adm0_obs_full$iso3)))
message(sprintf("  Excluded as too large:           %d (%s)",
                length(haqi_large_countries()),
                paste(haqi_large_countries(), collapse = ", ")))
message(sprintf("  ADM0-eligible:                   %d", length(adm0_eligible)))

# Load ADM0 polygons for eligible countries (or empty SpatVector if none)
if (length(adm0_eligible) > 0) {
  adm0_attrs <- read_csv("results/global_adm0_predictors_1990_2023.csv",
                         show_col_types = FALSE) |>
    filter(iso3 %in% adm0_eligible) |>
    distinct(iso3, adm_id) |>
    mutate(boundary_source = "ADM0_national",
           region_name     = NA_character_,
           poly_id         = NA_character_,
           poly_idx        = row_number())
  message(sprintf("  ADM0 attrs: %d rows", nrow(adm0_attrs)))
} else {
  adm0_attrs <- tibble(
    poly_idx = integer(), iso3 = character(),
    region_name = character(), boundary_source = character(),
    poly_id = character(), adm_id = integer()
  )
}

# == 2. Build poly_id lookup for WASH data joining ============================
# The matched gpkg has all (region_name, year, survey) variants pointing to
# the same physical polygon.  Map each variant to its poly_id so we can
# join WASH data (variant names) to extracted rasters (canonical name from
# dedup'd gpkg).  Must also apply the same boundary_source merging that
# select_best_divisions.R uses (fingerprint approach).

message("Building poly_id lookup from matched gpkg ...")
wash_sv_full <- vect("results/JMP_regions_matched.gpkg")
full_df      <- as.data.frame(wash_sv_full)
full_cents   <- crds(centroids(wash_sv_full))
full_df$poly_id <- paste(round(full_cents[, 1], 3), round(full_cents[, 2], 3))

# Recompute canonical boundary_source (same fingerprint logic)
source_fp <- full_df |>
  group_by(iso3, boundary_source) |>
  summarise(fp = paste(sort(unique(poly_id)), collapse = "|"), .groups = "drop")
src_rows  <- full_df |> count(iso3, boundary_source, name = "n_rows")
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

# Some boundary files have multiple polygons with the same region name (e.g.
# COD "Equateur" split into two).  Keep the polygon with the most rows in the
# matched gpkg (= the one most surveys matched to).
name_poly_lookup <- full_df |>
  count(iso3, boundary_source, region_name, poly_id, name = "n_match") |>
  group_by(iso3, boundary_source, region_name) |>
  slice_max(n_match, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(iso3, boundary_source, region_name, poly_id) |>
  left_join(
    bind_rows(
      best_attrs   |> select(iso3, boundary_source, poly_id, adm_id),
      second_attrs |> select(iso3, boundary_source, poly_id, adm_id)
    ) |> distinct(iso3, boundary_source, poly_id, .keep_all = TRUE),
    by = c("iso3", "boundary_source", "poly_id")
  )
message(sprintf("  %d name-to-polygon mappings", nrow(name_poly_lookup)))

rm(wash_sv_full, full_df, full_cents, source_fp, src_rows, canon_map)
gc(verbose = FALSE)

# == 4. Read predictors from global lookup (written by script 07) =============
message("\n== Reading global predictors (from script 07) ==")

adm1_pred <- read_csv("results/global_adm1_predictors_1990_2023.csv",
                      show_col_types = FALSE) |>
  filter(year %in% years)

adm0_pred <- read_csv("results/global_adm0_predictors_1990_2023.csv",
                      show_col_types = FALSE) |>
  filter(year %in% years)

training_ids <- c(best_attrs$adm_id, second_attrs$adm_id)
adm1_pred_train <- adm1_pred |> filter(adm_id %in% training_ids)

stopifnot(setequal(unique(adm1_pred_train$adm_id), unique(training_ids)))

adm0_pred_train <- adm0_pred |> filter(iso3 %in% adm0_eligible)

best_long <- adm1_pred_train |>
  inner_join(best_attrs |> select(adm_id, boundary_source, poly_id),
             by = "adm_id") |>
  mutate(division_rank = 1L) |>
  select(division_rank, adm_id, iso3, region_name, boundary_source, poly_id,
         year, pop, all_of(c("poverty", "msch", "gnic", "gini", "urban",
                             "travel_time", "aridity", "gdi", "corruption", "gw_depth")))

second_long <- adm1_pred_train |>
  inner_join(second_attrs |> select(adm_id, boundary_source, poly_id),
             by = "adm_id") |>
  mutate(division_rank = 2L) |>
  select(division_rank, adm_id, iso3, region_name, boundary_source, poly_id,
         year, pop, all_of(c("poverty", "msch", "gnic", "gini", "urban",
                             "travel_time", "aridity", "gdi", "corruption", "gw_depth")))

adm0_long <- if (length(adm0_eligible) > 0) {
  adm0_pred_train |>
    mutate(division_rank   = 0L,
           region_name     = NA_character_,
           boundary_source = "ADM0_national",
           poly_id         = NA_character_) |>
    select(division_rank, adm_id, iso3, region_name, boundary_source, poly_id,
           year, pop, all_of(c("poverty", "msch", "gnic", "gini", "urban",
                               "travel_time", "aridity", "gdi", "corruption", "gw_depth")))
} else {
  tibble()
}

all_extracted <- bind_rows(best_long, second_long, adm0_long)

message(sprintf("  Total: %d rows (%d best + %d second + %d adm0)",
                nrow(all_extracted), nrow(best_long), nrow(second_long),
                nrow(adm0_long)))

write_csv(all_extracted, "results/xgb_raster_all_years.csv")
message(sprintf("Saved: results/xgb_raster_all_years.csv (%d rows)",
                nrow(all_extracted)))

# == 6. Join with WASH survey data ============================================

message("\n== Joining with WASH data ==")

wash_best   <- read_csv("results/wash_best_division_data.csv",   show_col_types = FALSE)
wash_second <- read_csv("results/wash_second_division_data.csv", show_col_types = FALSE)

wash_cols <- intersect(
  c("Improved", "Not.Improved", "Basic", "Limited", "Unimproved", "Surface.water"),
  names(wash_best)
)

wash_all <- bind_rows(
  wash_best   |> mutate(division_rank = 1L),
  wash_second |> mutate(division_rank = 2L)
)

# Add poly_id using the name lookup (handles name variants + canonical BS)
wash_all <- wash_all |>
  left_join(
    name_poly_lookup |> select(iso3, boundary_source, region_name, poly_id),
    by = c("iso3", "boundary_source", "region_name")
  )

n_no_poly <- sum(is.na(wash_all$poly_id))
if (n_no_poly > 0) {
  message(sprintf("  WARNING: %d WASH rows could not be matched to a polygon",
                  n_no_poly))
  unmatched <- wash_all |>
    filter(is.na(poly_id)) |>
    distinct(iso3, boundary_source, region_name)
  message("  Unmatched (iso3 / boundary_source / region_name):")
  for (i in seq_len(min(nrow(unmatched), 20))) {
    message(sprintf("    %s | %s | %s",
                    unmatched$iso3[i], unmatched$boundary_source[i],
                    unmatched$region_name[i]))
  }
}

# Join raster values on adm_id + year (stable numeric key -- no name matching)
training_data <- wash_all |>
  select(division_rank, adm_id, iso3, region_name, boundary_source, poly_id,
         year, data_source, survey_type, survey_family,
         all_of(wash_cols)) |>
  left_join(
    all_extracted |>
      select(adm_id, year, pop, poverty, msch, gnic, gini,
             urban, travel_time, aridity, gdi, corruption, gw_depth),
    by = c("adm_id", "year")
  )

n_joined <- sum(!is.na(training_data$pop))
message(sprintf("  %d / %d WASH rows matched with raster values (%.1f%%)",
                n_joined, nrow(training_data),
                100 * n_joined / nrow(training_data)))

# Drop subnational rows where all four WASH indicators are NA.
# These rows (e.g. EGY, HND, SOM, VNM subnational polygons) entered via
# wash_best/second_division_data.csv but carry no trainable signal.
# They waste year-sampler slots in script 06 and can cause countries to
# silently drop from training when sampling draws only these empty rows.
wash_any_nonNA_subnat <- !is.na(training_data$Basic) |
                         !is.na(training_data$Limited) |
                         !is.na(training_data$Unimproved) |
                         !is.na(training_data$Surface.water)
n_allna_subnat <- sum(!wash_any_nonNA_subnat)
if (n_allna_subnat > 0) {
  iso_dropped <- sort(unique(training_data$iso3[!wash_any_nonNA_subnat]))
  message(sprintf("  Dropping %d subnational all-NA WASH rows (countries: %s)",
                  n_allna_subnat, paste(iso_dropped, collapse = ", ")))
  training_data <- training_data[wash_any_nonNA_subnat, ]
}

# == 7. Append ADM0 rows for countries without subnational training data =====
#
# ADM0 raster predictors were already extracted in section 4b (combined
# extraction). Here we just join those predictors with JMP ADM0 observations
# to construct the ADM0 training rows.

message("\n== Appending ADM0 rows (countries without subnational training data) ==")

# adm0_obs was loaded earlier (section 1b) but not yet filtered to year range
adm0_obs <- adm0_obs_full |>
  filter(year %in% years, iso3 %in% adm0_eligible)

if (length(adm0_eligible) > 0 && nrow(adm0_long) > 0) {

  # adm0_long has division_rank = 0L; pull predictors per (iso3, year).
  # Each ADM0 row already has a unique adm_id (from make_country_adm_id_lookup
  # in section 1b), so we can join straight on adm_id+year via adm0_attrs.
  adm0_pred <- adm0_long |>
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
      Improved        = NA_real_,
      Not.Improved    = NA_real_,
      Basic,
      Limited,
      Unimproved,
      Surface.water,
      pop,
      poverty,
      msch,
      gnic,
      gini,
      urban,
      travel_time,
      aridity,
      gdi,
      corruption,
      gw_depth
    )

  # Drop ADM0 rows where all four WASH indicators are NA — no response to
  # train on. Keeps rows where at least one indicator is non-NA (partial
  # coverage is still useful since the model trains per-indicator).
  wash_any_nonNA_adm0 <- !is.na(adm0_training$Basic) |
                         !is.na(adm0_training$Limited) |
                         !is.na(adm0_training$Unimproved) |
                         !is.na(adm0_training$Surface.water)
  n_allna_adm0 <- sum(!wash_any_nonNA_adm0)
  if (n_allna_adm0 > 0) {
    iso_adm0_dropped <- sort(unique(adm0_training$iso3[!wash_any_nonNA_adm0]))
    message(sprintf("  Dropping %d all-NA WASH ADM0 rows (%d countries: %s)",
                    n_allna_adm0,
                    length(iso_adm0_dropped),
                    paste(iso_adm0_dropped, collapse = ", ")))
    adm0_training <- adm0_training[wash_any_nonNA_adm0, ]
  }

  # Defensive: this should never trigger because adm0_eligible was built via
  # setdiff(adm0_obs_iso3, countries_with_subnat). If it does, something
  # upstream has gone wrong.
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

write_csv(training_data, "results/xgb_training_data.csv")
message(sprintf("Saved: results/xgb_training_data.csv (%d rows)",
                nrow(training_data)))

message("\nDone!")
