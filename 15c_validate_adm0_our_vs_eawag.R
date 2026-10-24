# 15c_validate_adm0_our_vs_eawag.R
# ---------------------------------------------------------------------------
# ADM0 validation for Eawag product (year 2020), reported analogously to
# the IHME workflow in script 15:
#   - Compare Eawag improved vs reported JMP Basic+Limited
#   - Compare OUR Basic+Limited vs reported JMP Basic+Limited
#   - Use the same country set: Eawag 135-country scope only
#
# Inputs:
#   data_in/EawagPHH-mapping-safe-drinking-water-use-LMICs-00a413e/
#     5_Predictions/prediction_results/Fullprediction_improved.csv
#   data_in/EawagPHH-mapping-safe-drinking-water-use-LMICs-00a413e/
#     5_Predictions/prediction_results/SMDWPredictions_GADMlevel1.csv
#   data_in/JMP_2025_WLD.xlsx (sheet 'wat')
#   /Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg
#   /Users/mkummu/R/GIS_data_common/GHS_POP/r_pop_GHS_1985_2025_5arcmin.tif
#   results_gis/rast_wash_Basic_1990_2023.tif
#   results_gis/rast_wash_Limited_1990_2023.tif
#
# Outputs:
#   results/validation_eawag_scope_iso3.csv
#   results/validation_eawag_adm1_join_status.csv
#   results/validation_eawag_adm1_values.csv
#   results/validation_adm0_eawag_extracted_2020.csv
#   results/validation_adm0_our_extracted_eawag_scope_2020.csv
#   results/validation_adm0_reported_basic_limited_eawag_scope_2020.csv
#   results/validation_adm0_eawag_reported_common_2020.csv
#   results/validation_adm0_eawag_vs_reported_2020.csv
#   results/validation_adm0_our_vs_reported_eawag_scope_2020.csv
#   results/validation_adm0_eawag_reported_summary_2020.csv
#   results_gis/polyg_eawag_improved_adm1.gpkg
#   results_gis/rast_eawag_adm1_template_5arcmin.tif
#   results_gis/rast_eawag_improved_adm1_5arcmin.tif
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(stringi)
})

dir.create("results", showWarnings = FALSE)
dir.create("results_gis", showWarnings = FALSE)

### 0. Paths ----
source("config.R")
paths <- subnat_wash_paths()
eawag_csv <- "data_in/EawagPHH-mapping-safe-drinking-water-use-LMICs-00a413e/5_Predictions/prediction_results/Fullprediction_improved.csv"
eawag_scope_csv <- "data_in/EawagPHH-mapping-safe-drinking-water-use-LMICs-00a413e/5_Predictions/prediction_results/SMDWPredictions_GADMlevel1.csv"
xlsx_path <- "data_in/JMP_2025_WLD.xlsx"
gadm_path <- paths$gadm
pop_path <- paths$population
our_basic_tif <- "results_gis/rast_wash_Basic_1990_2023.tif"
our_limited_tif <- "results_gis/rast_wash_Limited_1990_2023.tif"

year_target <- 2020L

for (p in c(eawag_csv, eawag_scope_csv, xlsx_path, gadm_path, pop_path, our_basic_tif, our_limited_tif)) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}

normalise_name <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    tolower() |>
    str_replace_all("[^a-z0-9 ]", " ") |>
    str_squish()
}

### 1. Eawag scope (135 countries) ----
message("== Building Eawag 135-country scope ==")

eawag_scope <- read_csv(eawag_scope_csv, show_col_types = FALSE) |>
  transmute(iso3 = toupper(as.character(GID_0))) |>
  filter(!is.na(iso3), iso3 != "") |>
  distinct() |>
  arrange(iso3)

write_csv(eawag_scope, "results/validation_eawag_scope_iso3.csv")
message(sprintf("  Scope countries: %d", nrow(eawag_scope)))

### 2. Load Eawag predictions and GADM ADM1 ----
message("== Loading Eawag improved predictions and GADM ADM1 ==")

eawag <- read_csv(eawag_csv, show_col_types = FALSE) |>
  transmute(
    GID_1 = as.character(GID_1),
    iso3 = toupper(as.character(GID_0)),
    country = as.character(country),
    region = as.character(region),
    wb_region = as.character(`World Bank Region`),
    income_group = as.character(Income.group),
    value = as.numeric(predict) * 100
  ) |>
  semi_join(eawag_scope, by = "iso3")

gadm1 <- vect(gadm_path, layer = "ADM_1")
gadm1_df <- as.data.frame(gadm1) |>
  transmute(
    GID_1 = as.character(GID_1),
    iso3 = as.character(GID_0),
    adm1_name = as.character(NAME_1)
  )

# Exact GID_1 match
e_join <- eawag |>
  left_join(gadm1_df, by = c("GID_1", "iso3"), relationship = "one-to-one") |>
  mutate(match_type = if_else(!is.na(adm1_name), "gid1_exact", NA_character_))

# Fallback by name for unmatched
unmatched <- e_join |>
  filter(is.na(adm1_name)) |>
  mutate(region_norm = normalise_name(region))

if (nrow(unmatched) > 0) {
  gadm_name_pool <- gadm1_df |>
    mutate(adm1_norm = normalise_name(adm1_name))

  fallback_exact <- unmatched |>
    left_join(
      gadm_name_pool |>
        select(iso3, GID_1_fb = GID_1, adm1_name_fb = adm1_name, adm1_norm),
      by = c("iso3", "region_norm" = "adm1_norm"),
      relationship = "many-to-many"
    ) |>
    group_by(GID_1, iso3) |>
    slice_head(n = 1) |>
    ungroup() |>
    mutate(match_type_fb = if_else(!is.na(GID_1_fb), "name_exact", NA_character_))

  still_unmatched <- fallback_exact |> filter(is.na(GID_1_fb))

  if (nrow(still_unmatched) > 0) {
    fuzzy_best <- function(region_norm, pool_norm, pool_gid, pool_name) {
      if (length(pool_norm) == 0) {
        return(list(gid = NA_character_, nm = NA_character_, dist = NA_real_))
      }
      d <- stringdist::stringdist(region_norm, pool_norm, method = "jw")
      i <- which.min(d)
      list(gid = pool_gid[i], nm = pool_name[i], dist = d[i])
    }

    fuzzy_tbl <- still_unmatched |>
      rowwise() |>
      mutate(
        fuzzy = list({
          pool <- gadm_name_pool |> filter(iso3 == .data$iso3)
          fuzzy_best(
            region_norm = .data$region_norm,
            pool_norm = pool$adm1_norm,
            pool_gid = pool$GID_1,
            pool_name = pool$adm1_name
          )
        }),
        GID_1_fz = fuzzy$gid,
        adm1_name_fz = fuzzy$nm,
        dist_fz = as.numeric(fuzzy$dist),
        match_type_fz = if_else(!is.na(dist_fz) & dist_fz <= 0.12, "name_fuzzy", NA_character_)
      ) |>
      ungroup() |>
      select(GID_1, iso3, GID_1_fz, adm1_name_fz, match_type_fz)

    fallback_exact <- fallback_exact |>
      left_join(fuzzy_tbl, by = c("GID_1", "iso3")) |>
      mutate(
        GID_1_final = coalesce(GID_1_fb, GID_1_fz),
        adm1_name_final = coalesce(adm1_name_fb, adm1_name_fz),
        match_type_final = coalesce(match_type_fb, match_type_fz)
      )
  } else {
    fallback_exact <- fallback_exact |>
      mutate(
        GID_1_final = GID_1_fb,
        adm1_name_final = adm1_name_fb,
        match_type_final = match_type_fb
      )
  }

  e_join <- e_join |>
    left_join(
      fallback_exact |>
        select(GID_1, iso3, GID_1_final, adm1_name_final, match_type_final),
      by = c("GID_1", "iso3")
    ) |>
    mutate(
      GID_1_use = coalesce(GID_1_final, GID_1),
      adm1_name = coalesce(adm1_name, adm1_name_final),
      match_type = coalesce(match_type, match_type_final)
    ) |>
    select(-GID_1_final, -adm1_name_final, -match_type_final)
} else {
  e_join <- e_join |> mutate(GID_1_use = GID_1)
}

join_status <- e_join |>
  mutate(is_matched = !is.na(adm1_name)) |>
  count(match_type, is_matched, name = "n") |>
  arrange(desc(is_matched), match_type)

write_csv(join_status, "results/validation_eawag_adm1_join_status.csv")

e_matched <- e_join |>
  filter(!is.na(adm1_name), !is.na(value)) |>
  distinct(GID_1_use, .keep_all = TRUE)

message(sprintf("  Eawag rows in scope: %d | matched to GADM ADM1: %d", nrow(eawag), nrow(e_matched)))

### 3. Build Eawag ADM1 polygons + rasters ----
message("== Building Eawag ADM1 polygons and rasters ==")

gadm1_use <- gadm1[gadm1$GID_1 %in% e_matched$GID_1_use]
e_matched_join <- e_matched |>
  select(-GID_1) |>
  rename(GID_1 = GID_1_use)

gadm1_use <- merge(gadm1_use, e_matched_join, by = "GID_1", all.x = TRUE)
gadm1_use$eawag_id <- as.integer(seq_len(nrow(gadm1_use)))

writeVector(gadm1_use, "results_gis/polyg_eawag_improved_adm1.gpkg", overwrite = TRUE)

write_csv(
  as.data.frame(gadm1_use) |>
    transmute(
      eawag_id,
      GID_1,
      iso3,
      country,
      region,
      wb_region,
      income_group,
      eawag_improved = value,
      match_type
    ),
  "results/validation_eawag_adm1_values.csv"
)

r_ref_5 <- rast(pop_path, lyrs = 1)
r_ref_1 <- rast(
  nrows = nrow(r_ref_5) * 5,
  ncols = ncol(r_ref_5) * 5,
  ext = ext(r_ref_5),
  crs = crs(r_ref_5)
)

r_eawag_zone_1 <- rasterize(gadm1_use, r_ref_1, field = "eawag_id")
r_eawag_zone_5 <- aggregate(r_eawag_zone_1, fact = 5, fun = "modal", na.rm = TRUE)
names(r_eawag_zone_5) <- "eawag_id"

r_eawag_value_5 <- classify(r_eawag_zone_5, cbind(gadm1_use$eawag_id, gadm1_use$value))
names(r_eawag_value_5) <- "eawag_improved"

writeRaster(
  r_eawag_zone_5,
  "results_gis/rast_eawag_adm1_template_5arcmin.tif",
  datatype = "INT4S",
  gdal = "COMPRESS=LZW",
  overwrite = TRUE
)

writeRaster(
  r_eawag_value_5,
  "results_gis/rast_eawag_improved_adm1_5arcmin.tif",
  gdal = "COMPRESS=LZW",
  overwrite = TRUE
)

### 4. Build ADM0 template + population 2020 ----
message("== Building ADM0 template and population weights ==")

gadm0 <- vect(gadm_path, layer = "ADM_0")
gadm0$adm0_id <- as.integer(seq_len(nrow(gadm0)))

r_adm0_1 <- rasterize(gadm0, r_ref_1, field = "adm0_id")
r_adm0_5 <- aggregate(r_adm0_1, fact = 5, fun = "modal", na.rm = TRUE)
names(r_adm0_5) <- "adm0_id"

adm0_lookup <- tibble(
  adm0_id = gadm0$adm0_id,
  iso3 = as.character(gadm0$GID_0),
  country_name = as.character(gadm0$COUNTRY)
) |>
  semi_join(eawag_scope, by = "iso3")

pop_years <- 1985:2025
pop_idx <- match(year_target, pop_years)
r_pop_2020 <- rast(pop_path)[[pop_idx]]
names(r_pop_2020) <- paste0("pop_", year_target)

if (!compareGeom(r_pop_2020, r_adm0_5, stopOnError = FALSE)) {
  r_pop_2020 <- resample(r_pop_2020, r_adm0_5, method = "bilinear")
}

### 5. Aggregate Eawag improved to ADM0 ----
message("== Aggregating Eawag improved to ADM0 ==")

if (!compareGeom(r_eawag_value_5, r_adm0_5, stopOnError = FALSE)) {
  r_eawag_value_5 <- resample(r_eawag_value_5, r_adm0_5, method = "bilinear")
}

r_e_prod <- r_eawag_value_5 * r_pop_2020
z_e_prod <- zonal(r_e_prod, r_adm0_5, fun = "sum", na.rm = TRUE)
z_e_pop <- zonal(r_pop_2020, r_adm0_5, fun = "sum", na.rm = TRUE)

eawag_adm0 <- tibble(
  adm0_id = as.integer(z_e_prod[[1]]),
  value = z_e_prod[[2]] / z_e_pop[[2]],
  source = "eawag",
  indicator = "Basic+Limited",
  year = year_target
) |>
  mutate(value = if_else(is.nan(value), NA_real_, value)) |>
  left_join(adm0_lookup, by = "adm0_id") |>
  filter(!is.na(iso3)) |>
  relocate(source, indicator, iso3, country_name, adm0_id, year, value)

write_csv(eawag_adm0, "results/validation_adm0_eawag_extracted_2020.csv")

### 6. Aggregate OUR Basic+Limited to ADM0 (same scope) ----
message("== Aggregating OUR Basic+Limited to ADM0 (same countries) ==")

get_year_layer <- function(r, year) {
  lyr_years <- suppressWarnings(as.integer(str_extract(names(r), "\\d{4}$")))
  idx <- which(lyr_years == year)
  if (length(idx) == 0) stop("Year ", year, " not found in raster")
  r[[idx[1]]]
}

r_basic <- rast(our_basic_tif)
r_limited <- rast(our_limited_tif)
r_our_bl <- get_year_layer(r_basic, year_target) + get_year_layer(r_limited, year_target)
names(r_our_bl) <- "our_basic_limited"

if (!compareGeom(r_our_bl, r_adm0_5, stopOnError = FALSE)) {
  r_our_bl <- resample(r_our_bl, r_adm0_5, method = "bilinear")
}

r_o_prod <- r_our_bl * r_pop_2020
z_o_prod <- zonal(r_o_prod, r_adm0_5, fun = "sum", na.rm = TRUE)
z_o_pop <- zonal(r_pop_2020, r_adm0_5, fun = "sum", na.rm = TRUE)

our_adm0 <- tibble(
  adm0_id = as.integer(z_o_prod[[1]]),
  value = z_o_prod[[2]] / z_o_pop[[2]],
  source = "our",
  indicator = "Basic+Limited",
  year = year_target
) |>
  mutate(value = if_else(is.nan(value), NA_real_, value)) |>
  left_join(adm0_lookup, by = "adm0_id") |>
  filter(!is.na(iso3)) |>
  relocate(source, indicator, iso3, country_name, adm0_id, year, value)

write_csv(our_adm0, "results/validation_adm0_our_extracted_eawag_scope_2020.csv")

### 7. Read reported ADM0 JMP Basic+Limited (2020) ----
message("== Reading reported ADM0 JMP Basic+Limited (2020) ==")

wat_raw <- read_excel(xlsx_path, sheet = "wat", skip = 1, col_names = FALSE)
wat_data <- wat_raw[-1, ]

reported_adm0 <- tibble(
  iso3 = toupper(as.character(wat_data[[46]])),
  year = as.integer(as.numeric(wat_data[[2]])),
  Basic = as.numeric(wat_data[[15]]),
  Limited = as.numeric(wat_data[[16]])
) |>
  filter(!is.na(iso3), year == year_target) |>
  mutate(reported_value = Basic + Limited,
         indicator = "Basic+Limited") |>
  select(iso3, year, indicator, reported_value) |>
  semi_join(eawag_scope, by = "iso3")

write_csv(reported_adm0, "results/validation_adm0_reported_basic_limited_eawag_scope_2020.csv")

### 8. Build IHME-analog comparison tables ----
message("== Building common comparison tables (source vs reported) ==")

compare_common <- bind_rows(eawag_adm0, our_adm0) |>
  left_join(reported_adm0, by = c("iso3", "year", "indicator")) |>
  mutate(
    diff = value - reported_value,
    abs_diff = abs(diff)
  ) |>
  filter(!is.na(value), !is.na(reported_value))

write_csv(compare_common, "results/validation_adm0_eawag_reported_common_2020.csv")

write_csv(
  compare_common |> filter(source == "eawag"),
  "results/validation_adm0_eawag_vs_reported_2020.csv"
)

write_csv(
  compare_common |> filter(source == "our"),
  "results/validation_adm0_our_vs_reported_eawag_scope_2020.csv"
)

summary_tbl <- compare_common |>
  group_by(source, indicator) |>
  summarise(
    n = n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mae = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(diff^2, na.rm = TRUE)),
    mape_gt1 = mean(abs(diff) / if_else(abs(reported_value) <= 1, NA_real_, abs(reported_value)), na.rm = TRUE) * 100,
    smape = mean(2 * abs(diff) / (abs(value) + abs(reported_value)), na.rm = TRUE) * 100,
    pearson_r = suppressWarnings(cor(value, reported_value, use = "complete.obs")),
    r_squared = pearson_r^2,
    .groups = "drop"
  ) |>
  arrange(source)

write_csv(summary_tbl, "results/validation_adm0_eawag_reported_summary_2020.csv")

message("  Saved: results/validation_adm0_eawag_reported_common_2020.csv")
message("  Saved: results/validation_adm0_eawag_vs_reported_2020.csv")
message("  Saved: results/validation_adm0_our_vs_reported_eawag_scope_2020.csv")
message("  Saved: results/validation_adm0_eawag_reported_summary_2020.csv")
message("\nScript 15c completed successfully.")
