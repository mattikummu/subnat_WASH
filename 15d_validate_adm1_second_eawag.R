# 15d_validate_adm1_second_eawag.R
# ---------------------------------------------------------------------------
# ADM1 (second division) validation for Eawag product (years 2017-2019),
# analogous to IHME-style common-source reporting:
#   - Compare Eawag improved vs reported Basic+Limited
#   - Compare OUR Basic+Limited vs reported Basic+Limited
#   - Use Eawag 135-country scope only
#
# Inputs:
#   results/wash_second_division.gpkg
#   results/wash_second_division_data.csv
#   results/validation_eawag_scope_iso3.csv
#     (fallback: Eawag SMDWPredictions_GADMlevel1.csv)
#   /Users/mkummu/R/GIS_data_common/GHS_POP/r_pop_GHS_1985_2025_5arcmin.tif
#   results_gis/rast_wash_Basic_1990_2023.tif
#   results_gis/rast_wash_Limited_1990_2023.tif
#   results_gis/rast_eawag_improved_adm1_5arcmin.tif
#
# Outputs:
#   results/validation_adm1_second_eawag_scope_iso3.csv
#   results/validation_adm1_second_template_5arcmin.tif
#   results/validation_adm1_second_lookup.csv
#   results/validation_adm1_second_eawag_extracted_2017_2019.csv
#   results/validation_adm1_second_our_extracted_eawag_scope_2017_2019.csv
#   results/validation_adm1_second_reported_basic_limited_eawag_scope_2017_2019.csv
#   results/validation_adm1_second_eawag_reported_common_2017_2019.csv
#   results/validation_adm1_second_eawag_vs_reported_2017_2019.csv
#   results/validation_adm1_second_our_vs_reported_eawag_scope_2017_2019.csv
#   results/validation_adm1_second_eawag_reported_summary_2017_2019.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

dir.create("results", showWarnings = FALSE)

### 0. Paths ----
second_gpkg_path <- "results/wash_second_division.gpkg"
second_data_path <- "results/wash_second_division_data.csv"
scope_path <- "results/validation_eawag_scope_iso3.csv"
scope_fallback <- "data_in/EawagPHH-mapping-safe-drinking-water-use-LMICs-00a413e/5_Predictions/prediction_results/SMDWPredictions_GADMlevel1.csv"
pop_path <- "/Users/mkummu/R/GIS_data_common/GHS_POP/r_pop_GHS_1985_2025_5arcmin.tif"
our_basic_tif <- "results_gis/rast_wash_Basic_1990_2023.tif"
our_limited_tif <- "results_gis/rast_wash_Limited_1990_2023.tif"
eawag_tif <- "results_gis/rast_eawag_improved_adm1_5arcmin.tif"

years_target <- 2017:2019

for (p in c(second_gpkg_path, second_data_path, pop_path, our_basic_tif, our_limited_tif, eawag_tif)) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}

### 1. Eawag scope (135 countries) ----
message("== Loading Eawag scope ==")

if (file.exists(scope_path)) {
  eawag_scope <- read_csv(scope_path, show_col_types = FALSE) |>
    transmute(iso3 = toupper(as.character(iso3))) |>
    filter(!is.na(iso3), iso3 != "") |>
    distinct() |>
    arrange(iso3)
} else {
  if (!file.exists(scope_fallback)) stop("Missing scope file: ", scope_fallback)
  eawag_scope <- read_csv(scope_fallback, show_col_types = FALSE) |>
    transmute(iso3 = toupper(as.character(GID_0))) |>
    filter(!is.na(iso3), iso3 != "") |>
    distinct() |>
    arrange(iso3)
}

write_csv(eawag_scope, "results/validation_adm1_second_eawag_scope_iso3.csv")
message(sprintf("  Scope countries: %d", nrow(eawag_scope)))

reported_years_available <- read_csv(second_data_path, show_col_types = FALSE) |>
  transmute(year = as.integer(year)) |>
  filter(!is.na(year)) |>
  distinct() |>
  arrange(year)

years_use <- intersect(years_target, reported_years_available$year)
if (length(years_use) == 0) {
  stop("None of requested years found in second-division reported data: ", paste(years_target, collapse = ", "))
}
message(sprintf("  Using years: %s", paste(years_use, collapse = ", ")))

### 2. Build second-ADM1 zone template (1 arcmin -> 5 arcmin) ----
message("== Building ADM1-second zone template ==")

second_sv <- vect(second_gpkg_path)
if (!"adm_id" %in% names(as.data.frame(second_sv))) {
  stop("Second-division geometry missing required column: adm_id")
}

second_sv$adm_id <- as.integer(second_sv$adm_id)
if (anyDuplicated(second_sv$adm_id) > 0) {
  stop("Duplicate adm_id found in second-division geometry")
}

second_lookup <- as.data.frame(second_sv) |>
  transmute(
    adm_id = as.integer(adm_id),
    iso3 = toupper(as.character(iso3)),
    region_name = as.character(region_name),
    boundary_source = as.character(boundary_source)
  ) |>
  semi_join(eawag_scope, by = "iso3")

r_pop_full <- rast(pop_path)
r_ref_5 <- r_pop_full[[1]]

r_ref_1 <- rast(
  nrows = nrow(r_ref_5) * 5,
  ncols = ncol(r_ref_5) * 5,
  ext = ext(r_ref_5),
  crs = crs(r_ref_5)
)

r_adm1_1 <- rasterize(second_sv, r_ref_1, field = "adm_id")
r_adm1_5 <- aggregate(r_adm1_1, fact = 5, fun = "modal", na.rm = TRUE)
names(r_adm1_5) <- "adm_id"

writeRaster(
  r_adm1_5,
  "results/validation_adm1_second_template_5arcmin.tif",
  datatype = "INT4S",
  gdal = "COMPRESS=LZW",
  overwrite = TRUE
)

write_csv(second_lookup, "results/validation_adm1_second_lookup.csv")

### 3. Population helpers ----
message("== Preparing population weights ==")

pop_years <- 1985:2025
missing_pop <- years_use[is.na(match(years_use, pop_years))]
if (length(missing_pop) > 0) {
  stop("Population layers missing for years: ", paste(missing_pop, collapse = ", "))
}

extract_with_pop <- function(value_raster, year) {
  pop_idx <- match(year, pop_years)
  r_pop <- r_pop_full[[pop_idx]]
  names(r_pop) <- paste0("pop_", year)

  if (!compareGeom(value_raster, r_adm1_5, stopOnError = FALSE)) {
    value_raster <- resample(value_raster, r_adm1_5, method = "bilinear")
  }

  if (!compareGeom(r_pop, r_adm1_5, stopOnError = FALSE)) {
    r_pop <- resample(r_pop, r_adm1_5, method = "bilinear")
  }

  r_prod <- value_raster * r_pop
  z_prod <- zonal(r_prod, r_adm1_5, fun = "sum", na.rm = TRUE)
  z_pop <- zonal(r_pop, r_adm1_5, fun = "sum", na.rm = TRUE)

  tibble(
    adm_id = as.integer(z_prod[[1]]),
    year = year,
    value = z_prod[[2]] / z_pop[[2]]
  ) |>
    mutate(value = if_else(is.nan(value), NA_real_, value))
}

### 4. Extract Eawag to ADM1-second ----
message("== Extracting Eawag to ADM1-second ==")

r_eawag <- rast(eawag_tif)

eawag_long <- bind_rows(lapply(years_use, function(yr) {
  extract_with_pop(r_eawag, yr) |>
    mutate(source = "eawag", indicator = "Basic+Limited")
})) |>
  left_join(second_lookup, by = "adm_id") |>
  filter(!is.na(iso3)) |>
  relocate(source, indicator, iso3, region_name, adm_id, year, value, boundary_source)

write_csv(eawag_long, "results/validation_adm1_second_eawag_extracted_2017_2019.csv")

### 5. Extract OUR Basic+Limited to ADM1-second ----
message("== Extracting OUR Basic+Limited to ADM1-second ==")

get_year_layer <- function(r, year) {
  lyr_years <- suppressWarnings(as.integer(str_extract(names(r), "\\d{4}$")))
  idx <- which(lyr_years == year)
  if (length(idx) == 0) stop("Year ", year, " not found in raster")
  r[[idx[1]]]
}

r_basic <- rast(our_basic_tif)
r_limited <- rast(our_limited_tif)

our_long <- bind_rows(lapply(years_use, function(yr) {
  r_our_bl <- get_year_layer(r_basic, yr) + get_year_layer(r_limited, yr)
  names(r_our_bl) <- "our_basic_limited"
  extract_with_pop(r_our_bl, yr) |>
    mutate(source = "our", indicator = "Basic+Limited")
})) |>
  left_join(second_lookup, by = "adm_id") |>
  filter(!is.na(iso3)) |>
  relocate(source, indicator, iso3, region_name, adm_id, year, value, boundary_source)

write_csv(our_long, "results/validation_adm1_second_our_extracted_eawag_scope_2017_2019.csv")

### 6. Reported ADM1-second Basic+Limited (2017-2019) ----
message("== Reading reported ADM1-second Basic+Limited ==")

reported_long <- read_csv(second_data_path, show_col_types = FALSE) |>
  transmute(
    adm_id = as.integer(adm_id),
    iso3 = toupper(as.character(iso3)),
    region_name = as.character(region_name),
    year = as.integer(year),
    Basic = as.numeric(Basic),
    Limited = as.numeric(Limited)
  ) |>
  filter(year %in% years_use) |>
  mutate(
    indicator = "Basic+Limited",
    reported_value = Basic + Limited
  ) |>
  group_by(adm_id, year, indicator) |>
  summarise(
    iso3 = first(iso3),
    region_name = first(region_name),
    reported_value = mean(reported_value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  semi_join(eawag_scope, by = "iso3")

write_csv(reported_long, "results/validation_adm1_second_reported_basic_limited_eawag_scope_2017_2019.csv")

### 7. Build common comparison tables ----
message("== Building common comparison tables ==")

compare_common <- bind_rows(eawag_long, our_long) |>
  left_join(
    reported_long |>
      select(adm_id, year, indicator, reported_value),
    by = c("adm_id", "year", "indicator")
  ) |>
  mutate(
    diff = value - reported_value,
    abs_diff = abs(diff)
  ) |>
  filter(!is.na(value), !is.na(reported_value))

write_csv(compare_common, "results/validation_adm1_second_eawag_reported_common_2017_2019.csv")

write_csv(
  compare_common |> filter(source == "eawag"),
  "results/validation_adm1_second_eawag_vs_reported_2017_2019.csv"
)

write_csv(
  compare_common |> filter(source == "our"),
  "results/validation_adm1_second_our_vs_reported_eawag_scope_2017_2019.csv"
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
    pearson_r = {
      cc <- complete.cases(value, reported_value)
      if (sum(cc) >= 2) suppressWarnings(cor(value[cc], reported_value[cc])) else NA_real_
    },
    r_squared = pearson_r^2,
    .groups = "drop"
  ) |>
  arrange(source)

write_csv(summary_tbl, "results/validation_adm1_second_eawag_reported_summary_2017_2019.csv")

message("  Saved: results/validation_adm1_second_eawag_reported_common_2017_2019.csv")
message("  Saved: results/validation_adm1_second_eawag_vs_reported_2017_2019.csv")
message("  Saved: results/validation_adm1_second_our_vs_reported_eawag_scope_2017_2019.csv")
message("  Saved: results/validation_adm1_second_eawag_reported_summary_2017_2019.csv")
message("\nScript 15d completed successfully.")
