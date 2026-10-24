# 21_extract_water_quality.R
# ---------------------------------------------------------------------------
# Extract DynQual water-quality class rasters (produced by 20_classify_pathogen.R)
# to ADM1 polygons using population-weighted zonal means.
# Scale: 1–5  (1 = best quality, 5 = worst quality).
# Data covers ADM1 units × years 1980–2019.
#
# Inputs:
#   results_gis/pathogen_class_annual_1980_2019.tif  – pathogen classes
#   results_gis/organic_class_annual_1980_2019.tif   – organic classes
#   results_gis/salinity_class_annual_1980_2019.tif  – salinity classes
#   (overall computed here as geometric mean of the three above)
#   results_gis/wash_adm1_template_5arcmin.tif       – ADM1 zone raster
#   results_gis/polyg_wash_adm1_1990_2023_geomonly_simpl.gpkg – ADM1 geometry
#   /Users/mkummu/R/GIS_data_common/GHS_POP/r_pop_GHS_1985_2025_5arcmin.tif
#
# Outputs:
#   results/wq_adm1_extracted.csv     – long table: adm_id, year, indicator, wq_mean
#   results_gis/polyg_wq_adm1.gpkg   – simplified gpkg with WQ columns (year 2019)
# ---------------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(readr)

dir.create("results",     showWarnings = FALSE)
dir.create("results_gis", showWarnings = FALSE)

### 0. Setup ----
source("config.R")
paths <- subnat_wash_paths()
pop_path <- paths$population
gpkg_simpl_path <- "results_gis/polyg_wash_adm1_1990_2023_geomonly_simpl.gpkg"

wq_rasters <- list(
  pathogen = "results_gis/pathogen_class_annual_1980_2019.tif",
  organic  = "results_gis/organic_class_annual_1980_2019.tif",
  salinity = "results_gis/salinity_class_annual_1980_2019.tif"
)

wq_years_all <- 1980:2019
pop_years    <- 1985:2025
years_use    <- intersect(wq_years_all, pop_years)

### 1. Load ADM1 template raster ----
message("== Loading ADM1 template raster ==")
r_template <- rast("results_gis/wash_adm1_template_5arcmin.tif")
message(sprintf("  Template: %d x %d", nrow(r_template), ncol(r_template)))

### 2. Load population raster (for weighting) ----
message("== Loading population raster ==")
r_pop_full <- rast(pop_path)
# GHS-POP layers: 1985 = layer 1; select years overlapping WQ and population
pop_idx <- which(pop_years %in% years_use)
r_pop   <- r_pop_full[[pop_idx]]
names(r_pop) <- paste0("pop_", years_use)

if (!compareGeom(r_pop, r_template, stopOnError = FALSE)) {
  message("  Resampling population to template grid...")
  r_pop <- resample(r_pop, r_template, method = "bilinear")
}

### 3. Compute overall = geometric mean of the three indicator rasters ----
message("== Computing overall WQ raster (geometric mean) ==")

overall_tif <- "results_gis/overall_class_annual_1980_2019.tif"

r_path <- rast(wq_rasters$pathogen)
r_orga <- rast(wq_rasters$organic)
r_sali <- rast(wq_rasters$salinity)

# Resample all to template grid if needed
for (r_obj in c("r_path", "r_orga", "r_sali")) {
  r <- get(r_obj)
  if (!compareGeom(r, r_template, stopOnError = FALSE)) {
    assign(r_obj, resample(r, r_template, method = "bilinear"))
  }
}

# Flip each indicator so that 1 = best, 5 = worst -> 6 - raster inverts scale
r_path_f <- 6 - r_path
r_orga_f <- 6 - r_orga
r_sali_f <- 6 - r_sali

# Geometric mean of flipped rasters: (p * o * s)^(1/3)
r_overall_f <- (r_path_f * r_orga_f * r_sali_f) ^ (1/3)

# Flip back so that 5 = worst (original scale)
r_overall <- 6 - r_overall_f
names(r_overall) <- paste0("overall_class_", wq_years_all)

writeRaster(r_overall, overall_tif,
            gdal = c("COMPRESS=LZW", "PREDICTOR=3"), overwrite = TRUE)
message(sprintf("  Saved: %s", overall_tif))

wq_rasters$overall <- overall_tif

rm(r_path, r_orga, r_sali, r_overall)

### 4. Extract WQ for each indicator ----
message("== Extracting WQ to ADM1 (population-weighted mean via zonal) ==")

extract_wq_indicator <- function(indicator, tif_path) {
  message(sprintf("  %s", indicator))

  r_wq <- rast(tif_path)
  yr_idx <- which(wq_years_all %in% years_use)
  r_wq <- r_wq[[yr_idx]]
  names(r_wq) <- paste0(indicator, "_", years_use)

  if (!compareGeom(r_wq, r_template, stopOnError = FALSE)) {
    message("    Resampling WQ to template grid...")
    r_wq <- resample(r_wq, r_template, method = "bilinear")
  }

  # Population-weighted product
  r_prod <- r_wq * r_pop

  z_prod <- zonal(r_prod, r_template, fun = "sum", na.rm = TRUE)
  z_pop  <- zonal(r_pop,  r_template, fun = "sum", na.rm = TRUE)

  adm_col <- z_prod[[1]]

  results <- lapply(seq_along(years_use), function(i) {
    tibble(
      adm_id    = adm_col,
      year      = years_use[i],
      indicator = indicator,
      wq_mean   = z_prod[[i + 1]] / z_pop[[i + 1]]
    )
  })

  bind_rows(results)
}

wq_long <- bind_rows(mapply(extract_wq_indicator,
                            names(wq_rasters), wq_rasters,
                            SIMPLIFY = FALSE))

wq_long <- wq_long |>
  mutate(wq_mean = if_else(is.nan(wq_mean), NA_real_, wq_mean))

message(sprintf("  Total records:  %d", nrow(wq_long)))
message(sprintf("  ADM1 units:     %d", n_distinct(wq_long$adm_id)))
message(sprintf("  Indicators:     %s", paste(sort(unique(wq_long$indicator)), collapse = ", ")))
message(sprintf("  Year range:     %d – %d", min(wq_long$year), max(wq_long$year)))
message(sprintf("  WQ class range: %.2f – %.2f  (1=best, 5=worst)",
                min(wq_long$wq_mean, na.rm = TRUE),
                max(wq_long$wq_mean, na.rm = TRUE)))

### 5. Save long-format CSV ----
write_csv(wq_long, "results/wq_adm1_extracted.csv")
message("  Saved: results/wq_adm1_extracted.csv")

### 6. Build spatial layer and save GeoPackage ----
message("\n== Building spatial layer (all indicators, 2019) ==")

gpkg_sf <- st_read(gpkg_simpl_path, quiet = TRUE) |>
  select(adm_id)

wq_2019 <- wq_long |>
  filter(year == 2019) |>
  pivot_wider(names_from = indicator, values_from = wq_mean,
              names_prefix = "wq_")

polyg_wq <- gpkg_sf |>
  left_join(wq_2019, by = "adm_id")

st_write(polyg_wq, "results_gis/polyg_wq_adm1.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
message("  Saved: results_gis/polyg_wq_adm1.gpkg")

message("\nScript 21 completed successfully.")
