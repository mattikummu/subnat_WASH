# 15_validate_adm0_wash.R
# ---------------------------------------------------------------------------
# Validate country-level (ADM0) water indicators by extracting raster values
# as population-weighted means and comparing against JMP reported ADM0 values.
#
# Coverage:
#   1) Our raster products (1990-2023): Basic, Limited, Unimproved, Surface.water
#   2) IHME external rasters (2000-2017): Unimproved, Surface.water
#
# Method:
#   - Build ADM0 zone raster from GADM level 0:
#       1 arcmin rasterization -> 5 arcmin modal aggregation
#   - Extract zonal population-weighted means:
#       sum(value * pop) / sum(pop)
#
# Inputs:
#   /Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg
#   /Users/mkummu/R/GIS_data_common/GHS_POP/r_pop_GHS_1985_2025_5arcmin.tif
#   results_gis/rast_wash_{indicator}_1990_2023.tif
#   data_in/ihme_wat_access/unimproved/*.TIF
#   data_in/ihme_wat_access/surf_water/*.TIF
#   data_in/JMP_2025_WLD.xlsx (sheet 'wat')
#
# Outputs:
#   results/validation_adm0_template_5arcmin.tif
#   results/validation_adm0_lookup.csv
#   results/validation_adm0_our_extracted.csv
#   results/validation_adm0_ihme_extracted.csv
#   results/validation_adm0_jmp_water.csv
#   results/validation_adm0_compare_common_2000_2017.csv
#   results/validation_adm0_compare_our_full.csv
#   results/validation_adm0_compare_summary.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
})

dir.create("results", showWarnings = FALSE)

### 0. Paths and setup ----
source("config.R")
paths <- subnat_wash_paths()
gadm_path <- paths$gadm
pop_path <- paths$population
xlsx_path <- "data_in/JMP_2025_WLD.xlsx"

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

for (p in c(gadm_path, pop_path, xlsx_path, unname(our_rasters))) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}

### 1. Build ADM0 template raster (1 arcmin -> 5 arcmin) ----
message("== Building ADM0 zone raster template ==")

lay_info <- st_layers(gadm_path)
adm0_layer <- lay_info$name[grepl("^ADM_0$|adm0|level0|level_0", lay_info$name, ignore.case = TRUE)][1]
if (is.na(adm0_layer)) stop("Could not identify ADM0 layer in: ", gadm_path)

adm0_sf <- st_read(gadm_path, layer = adm0_layer, quiet = TRUE) |>
  select(GID_0, COUNTRY)

if (!all(c("GID_0", "COUNTRY") %in% names(adm0_sf))) {
  stop("ADM0 layer must contain fields 'GID_0' and 'COUNTRY'")
}

adm0_sf <- adm0_sf |>
  mutate(
    adm0_id = as.integer(row_number()),
    iso3 = as.character(GID_0),
    country_name = as.character(COUNTRY)
  ) |>
  select(adm0_id, iso3, country_name)

adm0_sv <- vect(adm0_sf)

r_pop_full <- rast(pop_path)
r_ref_5 <- r_pop_full[[1]]

r_ref_1 <- rast(
  nrows = nrow(r_ref_5) * 5,
  ncols = ncol(r_ref_5) * 5,
  ext = ext(r_ref_5),
  crs = crs(r_ref_5)
)

r_adm0_1 <- rasterize(adm0_sv, r_ref_1, field = "adm0_id")
r_adm0_5 <- aggregate(r_adm0_1, fact = 5, fun = "modal", na.rm = TRUE)
names(r_adm0_5) <- "adm0_id"

adm0_template_path <- "results/validation_adm0_template_5arcmin.tif"
writeRaster(
  r_adm0_5,
  adm0_template_path,
  datatype = "INT4S",
  gdal = "COMPRESS=LZW",
  overwrite = TRUE
)

adm0_lookup <- st_drop_geometry(adm0_sf) |>
  mutate(iso3 = toupper(iso3))

write_csv(adm0_lookup, "results/validation_adm0_lookup.csv")
message(sprintf("  Saved: %s", adm0_template_path))
message("  Saved: results/validation_adm0_lookup.csv")

### 2. Prepare population stack aligned to template ----
message("== Preparing population weights ==")

pop_years <- 1985:2025
r_pop <- r_pop_full
names(r_pop) <- paste0("pop_", pop_years)

if (!compareGeom(r_pop, r_adm0_5, stopOnError = FALSE)) {
  message("  Resampling population stack to ADM0 template grid...")
  r_pop <- resample(r_pop, r_adm0_5, method = "bilinear")
}

weighted_zonal_mean <- function(value_stack, years, zone_raster, pop_stack, indicator, source) {
  pop_idx <- match(years, pop_years)
  if (any(is.na(pop_idx))) {
    stop("Population years missing for indicator ", indicator, ": ",
         paste(years[is.na(pop_idx)], collapse = ", "))
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
  z_pop  <- zonal(pop_use, zone_raster, fun = "sum", na.rm = TRUE)

  adm_col <- z_prod[[1]]

  bind_rows(lapply(seq_along(years), function(i) {
    tibble(
      adm0_id = as.integer(adm_col),
      year = years[i],
      indicator = indicator,
      source = source,
      value = z_prod[[i + 1]] / z_pop[[i + 1]]
    )
  })) |>
    mutate(value = if_else(is.nan(value), NA_real_, value))
}

### 3. Extract OUR rasters to ADM0 ----
message("== Extracting OUR water rasters to ADM0 ==")

extract_our_indicator <- function(indicator, tif_path) {
  message(sprintf("  %s", indicator))

  r <- rast(tif_path)
  lyr_years <- suppressWarnings(as.integer(str_extract(names(r), "\\d{4}$")))

  if (all(is.na(lyr_years))) {
    yr_match <- str_extract(basename(tif_path), "\\d{4}_\\d{4}")
    if (is.na(yr_match)) stop("Could not parse year range from: ", tif_path)
    yr_pair <- as.integer(str_split(yr_match, "_", simplify = TRUE))
    lyr_years <- seq.int(yr_pair[1], yr_pair[2])
    if (length(lyr_years) != nlyr(r)) stop("Layer count does not match parsed year range for: ", tif_path)
  }

  keep <- !is.na(lyr_years) & (lyr_years %in% pop_years)
  r <- r[[which(keep)]]
  years <- lyr_years[keep]

  names(r) <- paste0(indicator, "_", years)

  weighted_zonal_mean(
    value_stack = r,
    years = years,
    zone_raster = r_adm0_5,
    pop_stack = r_pop,
    indicator = indicator,
    source = "our"
  )
}

our_long <- bind_rows(mapply(extract_our_indicator, names(our_rasters), our_rasters, SIMPLIFY = FALSE)) |>
  left_join(adm0_lookup, by = "adm0_id") |>
  relocate(source, indicator, iso3, country_name, adm0_id, year, value)

write_csv(our_long, "results/validation_adm0_our_extracted.csv")
message(sprintf("  Saved: results/validation_adm0_our_extracted.csv (%d rows)", nrow(our_long)))

### 4. Stack IHME rasters (2000-2017) and extract to ADM0 ----
message("== Extracting IHME rasters to ADM0 ==")

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
    zone_raster = r_adm0_5,
    pop_stack = r_pop,
    indicator = indicator,
    source = "ihme"
  )
}

ihme_long <- bind_rows(mapply(extract_ihme_indicator, names(ihme_dirs), ihme_dirs, SIMPLIFY = FALSE)) |>
  left_join(adm0_lookup, by = "adm0_id") |>
  relocate(source, indicator, iso3, country_name, adm0_id, year, value)

write_csv(ihme_long, "results/validation_adm0_ihme_extracted.csv")
message(sprintf("  Saved: results/validation_adm0_ihme_extracted.csv (%d rows)", nrow(ihme_long)))

### 5. Read JMP ADM0 reported data ----
message("== Reading JMP ADM0 reported data ==")

wat_raw <- read_excel(xlsx_path, sheet = "wat", skip = 1, col_names = FALSE)
wat_data <- wat_raw[-1, ]

jmp_adm0 <- tibble(
  country_name = as.character(wat_data[[1]]),
  year = as.integer(as.numeric(wat_data[[2]])),
  iso3 = toupper(as.character(wat_data[[46]])),
  Basic = as.numeric(wat_data[[15]]),
  Limited = as.numeric(wat_data[[16]]),
  Unimproved = as.numeric(wat_data[[17]]),
  Surface.water = as.numeric(wat_data[[18]])
) |>
  filter(!is.na(iso3), !is.na(year), !is.na(country_name)) |>
  pivot_longer(
    cols = c(Basic, Limited, Unimproved, Surface.water),
    names_to = "indicator",
    values_to = "jmp_value"
  )

write_csv(jmp_adm0, "results/validation_adm0_jmp_water.csv")
message(sprintf("  Saved: results/validation_adm0_jmp_water.csv (%d rows)", nrow(jmp_adm0)))

### 6. Build comparison tables ----
message("== Building comparison outputs ==")

all_extracted <- bind_rows(our_long, ihme_long)

compare_all <- all_extracted |>
  left_join(jmp_adm0 |> select(iso3, year, indicator, jmp_value),
            by = c("iso3", "year", "indicator")) |>
  mutate(
    diff = value - jmp_value,
    abs_diff = abs(diff)
  )

# Common overlap: both OUR and IHME available against JMP for 2000-2017, two IHME indicators
common_base <- compare_all |>
  filter(
    year >= 2000,
    year <= 2017,
    indicator %in% c("Unimproved", "Surface.water"),
    source %in% c("our", "ihme"),
    !is.na(value),
    !is.na(jmp_value)
  )

common_keys <- common_base |>
  distinct(source, adm0_id, iso3, year, indicator) |>
  mutate(has_data = TRUE) |>
  tidyr::pivot_wider(
    names_from = source,
    values_from = has_data,
    values_fill = FALSE
  ) |>
  filter(our, ihme) |>
  select(adm0_id, iso3, year, indicator)

compare_common <- common_base |>
  semi_join(common_keys, by = c("adm0_id", "iso3", "year", "indicator"))

# Full coverage for OUR rasters against JMP (all four indicators; years where JMP exists)
compare_our_full <- compare_all |>
  filter(source == "our")

write_csv(compare_common, "results/validation_adm0_compare_common_2000_2017.csv")
write_csv(compare_our_full, "results/validation_adm0_compare_our_full.csv")

summary_tbl <- compare_all |>
  filter(!is.na(jmp_value), !is.na(value)) |>
  group_by(source, indicator) |>
  summarise(
    n = n(),
    mean_diff = mean(diff, na.rm = TRUE),
    mae = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(diff^2, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  arrange(source, indicator)

write_csv(summary_tbl, "results/validation_adm0_compare_summary.csv")

message("  Saved: results/validation_adm0_compare_common_2000_2017.csv")
message("  Saved: results/validation_adm0_compare_our_full.csv")
message("  Saved: results/validation_adm0_compare_summary.csv")

message("\nScript 15 completed successfully.")
