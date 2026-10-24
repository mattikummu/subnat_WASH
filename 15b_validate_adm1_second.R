# 15b_validate_adm1_second.R
# ---------------------------------------------------------------------------
# Validate ADM1 (second-best reported division) water indicators by extracting
# raster values as population-weighted means and comparing against reported
# second-division ADM1 values.
#
# Coverage:
#   1) Our raster products (1990-2023): Basic, Limited, Unimproved, Surface.water
#   2) IHME external rasters (2000-2017): Unimproved, Surface.water
#
# Method:
#   - Use second-best ADM1 polygons from results/wash_second_division.gpkg
#   - Build zone raster from those polygons:
#       1 arcmin rasterization -> 5 arcmin modal aggregation
#   - Extract zonal population-weighted means:
#       sum(value * pop) / sum(pop)
#
# Inputs:
#   results/wash_second_division.gpkg
#   results/wash_second_division_data.csv
#   /Users/mkummu/R/GIS_data_common/GHS_POP/r_pop_GHS_1985_2025_5arcmin.tif
#   results_gis/rast_wash_{indicator}_1990_2023.tif
#   data_in/ihme_wat_access/unimproved/*.TIF
#   data_in/ihme_wat_access/surf_water/*.TIF
#
# Outputs:
#   results/validation_adm1_second_template_5arcmin.tif
#   results/validation_adm1_second_lookup.csv
#   results/validation_adm1_second_our_extracted.csv
#   results/validation_adm1_second_ihme_extracted.csv
#   results/validation_adm1_second_reported.csv
#   results/validation_adm1_second_compare_common_2000_2017.csv
#   results/validation_adm1_second_compare_our_full.csv
#   results/validation_adm1_second_compare_summary.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

dir.create("results", showWarnings = FALSE)

### 0. Paths and setup ----
source("config.R")
paths <- subnat_wash_paths()
second_gpkg_path <- "results/wash_second_division.gpkg"
second_data_path <- "results/wash_second_division_data.csv"
pop_path <- paths$population

our_rasters <- c(
  Basic = "results_gis/rast_wash_Basic_1990_2023.tif",
  Limited = "results_gis/rast_wash_Limited_1990_2023.tif",
  Unimproved = "results_gis/rast_wash_Unimproved_1990_2023.tif",
  Surface.water = "results_gis/rast_wash_Surface.water_1990_2023.tif"
)

ihme_dirs <- c(
  Unimproved = "data_in/ihme_wat_access/unimproved",
  Surface.water = "data_in/ihme_wat_access/surf_water"
)

for (p in c(second_gpkg_path, second_data_path, pop_path, unname(our_rasters))) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}

### 1. Build second-ADM1 template raster (1 arcmin -> 5 arcmin) ----
message("== Building ADM1-second zone raster template ==")

second_sv <- vect(second_gpkg_path)

if (!"adm_id" %in% names(as.data.frame(second_sv))) {
  stop("Second-division geometry missing required column: adm_id")
}

second_sv$adm_id <- as.integer(second_sv$adm_id)
if (anyDuplicated(second_sv$adm_id) > 0) {
  stop("Duplicate adm_id found in second-division geometry")
}

second_lookup <- as.data.frame(second_sv) |>
  select(adm_id, iso3, region_name, boundary_source) |>
  mutate(iso3 = toupper(as.character(iso3)))

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

template_path <- "results/validation_adm1_second_template_5arcmin.tif"
writeRaster(
  r_adm1_5,
  template_path,
  datatype = "INT4S",
  gdal = "COMPRESS=LZW",
  overwrite = TRUE
)

write_csv(second_lookup, "results/validation_adm1_second_lookup.csv")
message(sprintf("  Saved: %s", template_path))
message("  Saved: results/validation_adm1_second_lookup.csv")

### 2. Prepare population stack aligned to template ----
message("== Preparing population weights ==")

pop_years <- 1985:2025
r_pop <- r_pop_full
names(r_pop) <- paste0("pop_", pop_years)

if (!compareGeom(r_pop, r_adm1_5, stopOnError = FALSE)) {
  message("  Resampling population stack to ADM1-second template grid...")
  r_pop <- resample(r_pop, r_adm1_5, method = "bilinear")
}

weighted_zonal_mean <- function(value_stack, years, zone_raster, pop_stack, indicator, source) {
  pop_idx <- match(years, pop_years)
  if (any(is.na(pop_idx))) {
    stop(
      "Population years missing for indicator ", indicator, ": ",
      paste(years[is.na(pop_idx)], collapse = ", ")
    )
  }

  pop_use <- pop_stack[[pop_idx]]

  if (!compareGeom(value_stack, zone_raster, stopOnError = FALSE)) {
    value_stack <- resample(value_stack, zone_raster, method = "bilinear")
  }

  if (!compareGeom(value_stack, pop_use, stopOnError = FALSE)) {
    pop_use <- resample(pop_use, value_stack, method = "bilinear")
  }

  names(pop_use) <- paste0("pop_", years)

  r_prod <- value_stack * pop_use

  z_prod <- zonal(r_prod, zone_raster, fun = "sum", na.rm = TRUE)
  z_pop <- zonal(pop_use, zone_raster, fun = "sum", na.rm = TRUE)

  adm_col <- z_prod[[1]]

  bind_rows(lapply(seq_along(years), function(i) {
    tibble(
      adm_id = as.integer(adm_col),
      year = years[i],
      indicator = indicator,
      source = source,
      value = z_prod[[i + 1]] / z_pop[[i + 1]]
    )
  })) |>
    mutate(value = if_else(is.nan(value), NA_real_, value))
}

### 3. Extract OUR rasters to ADM1-second ----
message("== Extracting OUR water rasters to ADM1-second ==")

extract_our_indicator <- function(indicator, tif_path) {
  message(sprintf("  %s", indicator))

  r <- rast(tif_path)
  lyr_years <- suppressWarnings(as.integer(str_extract(names(r), "\\d{4}$")))

  if (all(is.na(lyr_years))) {
    yr_match <- str_extract(basename(tif_path), "\\d{4}_\\d{4}")
    if (is.na(yr_match)) stop("Could not parse year range from: ", tif_path)
    yr_pair <- as.integer(str_split(yr_match, "_", simplify = TRUE))
    lyr_years <- seq.int(yr_pair[1], yr_pair[2])
    if (length(lyr_years) != nlyr(r)) {
      stop("Layer count does not match parsed year range for: ", tif_path)
    }
  }

  keep <- !is.na(lyr_years) & (lyr_years %in% pop_years)
  r <- r[[which(keep)]]
  years <- lyr_years[keep]

  names(r) <- paste0(indicator, "_", years)

  weighted_zonal_mean(
    value_stack = r,
    years = years,
    zone_raster = r_adm1_5,
    pop_stack = r_pop,
    indicator = indicator,
    source = "our"
  )
}

our_long <- bind_rows(mapply(extract_our_indicator, names(our_rasters), our_rasters, SIMPLIFY = FALSE)) |>
  left_join(second_lookup, by = "adm_id") |>
  relocate(source, indicator, iso3, region_name, adm_id, year, value, boundary_source)

write_csv(our_long, "results/validation_adm1_second_our_extracted.csv")
message(sprintf("  Saved: results/validation_adm1_second_our_extracted.csv (%d rows)", nrow(our_long)))

### 4. Stack IHME rasters (2000-2017) and extract to ADM1-second ----
message("== Extracting IHME rasters to ADM1-second ==")

read_ihme_stack <- function(indicator, dir_path) {
  tif_files <- list.files(dir_path, pattern = "\\.TIF$", full.names = TRUE)
  if (length(tif_files) == 0) stop("No IHME TIF files found in: ", dir_path)

  file_tbl <- tibble(path = tif_files) |>
    mutate(
      year = suppressWarnings(
        as.integer(sub(
          ".*_(20[0-9]{2})_Y[0-9]{4}M[0-9]{2}D[0-9]{2}\\.TIF$",
          "\\1",
          basename(path)
        ))
      )
    ) |>
    filter(!is.na(year), year >= 2000, year <= 2017) |>
    arrange(year)

  if (nrow(file_tbl) == 0) stop("No IHME years 2000-2017 found in: ", dir_path)

  r <- rast(file_tbl$path)
  years <- file_tbl$year
  names(r) <- paste0(indicator, "_", years)

  list(r = r, years = years)
}

extract_ihme_indicator <- function(indicator, dir_path) {
  message(sprintf("  %s", indicator))
  ih <- read_ihme_stack(indicator, dir_path)

  weighted_zonal_mean(
    value_stack = ih$r,
    years = ih$years,
    zone_raster = r_adm1_5,
    pop_stack = r_pop,
    indicator = indicator,
    source = "ihme"
  )
}

ihme_long <- bind_rows(mapply(extract_ihme_indicator, names(ihme_dirs), ihme_dirs, SIMPLIFY = FALSE)) |>
  left_join(second_lookup, by = "adm_id") |>
  relocate(source, indicator, iso3, region_name, adm_id, year, value, boundary_source)

write_csv(ihme_long, "results/validation_adm1_second_ihme_extracted.csv")
message(sprintf("  Saved: results/validation_adm1_second_ihme_extracted.csv (%d rows)", nrow(ihme_long)))

### 5. Read reported second-division ADM1 data ----
message("== Reading reported ADM1 second-division data ==")

reported_long <- read_csv(second_data_path, show_col_types = FALSE) |>
  transmute(
    adm_id = as.integer(adm_id),
    iso3 = toupper(as.character(iso3)),
    region_name = as.character(region_name),
    year = as.integer(year),
    Basic,
    Limited,
    Unimproved,
    Surface.water
  ) |>
  pivot_longer(
    cols = c(Basic, Limited, Unimproved, Surface.water),
    names_to = "indicator",
    values_to = "reported_value"
  ) |>
  group_by(adm_id, year, indicator) |>
  summarise(
    iso3 = first(iso3),
    region_name = first(region_name),
    reported_value = mean(reported_value, na.rm = TRUE),
    n_reported = n(),
    .groups = "drop"
  )

write_csv(reported_long, "results/validation_adm1_second_reported.csv")
message(sprintf("  Saved: results/validation_adm1_second_reported.csv (%d rows)", nrow(reported_long)))

### 6. Build comparison tables ----
message("== Building comparison outputs ==")

all_extracted <- bind_rows(our_long, ihme_long)

compare_all <- all_extracted |>
  left_join(
    reported_long |> select(adm_id, year, indicator, reported_value),
    by = c("adm_id", "year", "indicator")
  ) |>
  mutate(
    diff = value - reported_value,
    abs_diff = abs(diff)
  )

compare_common <- compare_all |>
  filter(
    year >= 2000,
    year <= 2017,
    indicator %in% c("Unimproved", "Surface.water"),
    source %in% c("our", "ihme"),
    !is.na(value),
    !is.na(reported_value)
  )

common_keys <- compare_common |>
  distinct(source, adm_id, year, indicator) |>
  mutate(has_data = TRUE) |>
  tidyr::pivot_wider(
    names_from = source,
    values_from = has_data,
    values_fill = FALSE
  ) |>
  filter(our, ihme) |>
  select(adm_id, year, indicator)

compare_common <- compare_common |>
  semi_join(common_keys, by = c("adm_id", "year", "indicator"))

compare_our_full <- compare_all |>
  filter(source == "our")

write_csv(compare_common, "results/validation_adm1_second_compare_common_2000_2017.csv")
write_csv(compare_our_full, "results/validation_adm1_second_compare_our_full.csv")

summary_tbl <- compare_all |>
  filter(!is.na(reported_value), !is.na(value)) |>
  group_by(source, indicator) |>
  summarise(
    n = n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mae = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(diff^2, na.rm = TRUE)),
    mape_gt1 = mean(abs(diff) / if_else(abs(reported_value) <= 1, NA_real_, abs(reported_value)), na.rm = TRUE) * 100,
    smape = mean(2 * abs(diff) / (abs(value) + abs(reported_value)), na.rm = TRUE) * 100,
    .groups = "drop"
  ) |>
  arrange(source, indicator)

write_csv(summary_tbl, "results/validation_adm1_second_compare_summary.csv")

message("  Saved: results/validation_adm1_second_compare_common_2000_2017.csv")
message("  Saved: results/validation_adm1_second_compare_our_full.csv")
message("  Saved: results/validation_adm1_second_compare_summary.csv")

message("\nScript 15b completed successfully.")
