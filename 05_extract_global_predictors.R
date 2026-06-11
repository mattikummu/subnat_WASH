# 05_extract_global_predictors.R
# ---------------------------------------------------------------------------
# Single-pass extraction of population-weighted predictors for ALL admin
# polygons used in the global pipeline:
#   - ADM0 country polygons
#   - ADM1 best polygons
#   - ADM1 second-best polygons (training only)
#   - ADM1 GADM fallback polygons (for countries without best)
# years 1990-2023.
#
# Outputs (consumed by downstream scripts):
#   results/global_adm0_predictors_1990_2023.csv
#   results/global_adm1_predictors_1990_2023.csv
#   results/wash_adm1_global.gpkg
#   results/adm1_global_lookup.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(countrycode)
  library(exactextractr)
  library(sf)
  library(future)
  library(furrr)
})

source("cache_prompt.R")
source("functions/f_small_country_training.R")

# Increase GDAL cache for faster raster I/O (MB).
# Total cache footprint is roughly n_cores * cache_size_mb.
cache_size_mb <- 4000L
gdalCache(cache_size_mb)
message(sprintf("GDAL cache set to %.1f GB", cache_size_mb / 1024))

n_cores <- 5
total_cache_gb <- (n_cores * cache_size_mb) / 1024
message(sprintf("Total GDAL cache across %d workers: %.1f GB", n_cores, total_cache_gb))
if (total_cache_gb > 12) {
  message("  WARNING: High memory usage. Consider reducing n_cores or cache_size_mb.")
}

plan(multisession, workers = n_cores)
on.exit(plan(sequential), add = TRUE)
message(sprintf("Using %d cores for parallel extraction", n_cores))

dir.create("results", showWarnings = FALSE)

predictor_cols <- c(
  "travel_time", "poverty", "msch", "gnic", "gini",
  "urban", "aridity", "gdi", "corruption", "gw_depth"
)
years_pred <- 1990:2023
n_years <- length(years_pred)
gadm_path <- "/Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg"
pop_path <- "/Users/mkummu/R/GIS_data_common/GHS_POP/r_pop_GHS_1985_2025_5arcmin.tif"

adm0_cache <- "results/global_adm0_predictors_1990_2023.csv"
adm1_cache <- "results/global_adm1_predictors_1990_2023.csv"
combined_cache <- "results/global_all_admin_predictors_1990_2023.csv"
gpkg_path <- "results/wash_adm1_global.gpkg"
lookup_path <- "results/adm1_global_lookup.csv"

message("== Loading ADM0 country polygons ==")
gadm0 <- vect(gadm_path, layer = "ADM_0")
gadm0 <- gadm0[!is.na(gadm0$GID_0), ]
gadm0 <- gadm0[order(gadm0$GID_0), ]

adm0_lookup <- make_country_adm_id_lookup(gadm0$GID_0)

adm0_attrs <- tibble(
  poly_idx = seq_len(nrow(gadm0)),
  iso3 = as.character(gadm0$GID_0),
  country_name_gadm = as.character(gadm0$COUNTRY)
) |>
  left_join(adm0_lookup, by = "iso3")

na_idx <- which(is.na(adm0_attrs$adm_id))
if (length(na_idx) > 0) {
  message(sprintf("  WARNING: %d ADM0 polygons have NA adm_id: %s",
                  length(na_idx),
                  paste(adm0_attrs$iso3[na_idx], collapse = ", ")))
  used_ids <- unique(stats::na.omit(adm0_attrs$adm_id))
  fallback_ids <- setdiff(901L:999L, used_ids)
  if (length(fallback_ids) < length(na_idx)) {
    stop("Unable to assign ADM0 fallback IDs below 1000; fallback range 901-999 exhausted")
  }
  adm0_attrs$adm_id[na_idx] <- fallback_ids[seq_along(na_idx)]
}

message(sprintf("  %d ADM0 country polygons", nrow(gadm0)))

message("== Building global ADM1 geometry (best + GADM fallback) ==")
best_poly <- vect("results/wash_best_division.gpkg")
second_poly <- vect("results/wash_second_division.gpkg")

countries_all <- sort(unique(gadm0$GID_0))
countries_with_best <- sort(unique(best_poly$iso3))
countries_need_gadm <- setdiff(countries_all, countries_with_best)

message(sprintf("  Countries total: %d", length(countries_all)))
message(sprintf("  Countries with best ADM1: %d", length(countries_with_best)))
message(sprintf("  Countries needing GADM fallback: %d", length(countries_need_gadm)))

best_core <- best_poly
best_core$source <- "best"
second_core <- second_poly
second_core$source <- "second"

gadm1_all <- vect(gadm_path, layer = "ADM_1")
gadm1_fb <- gadm1_all[gadm1_all$GID_0 %in% countries_need_gadm]

if (nrow(gadm1_fb) > 0) {
  gadm_df <- tibble(
    row_id = seq_len(nrow(gadm1_fb)),
    iso3 = gadm1_fb$GID_0,
    region_name = gadm1_fb$NAME_1,
    region_code = gadm1_fb$GID_1
  ) |>
    mutate(
      cntry_id = suppressWarnings(countrycode(iso3, "iso3c", "iso3n")),
      cntry_id = as.integer(cntry_id),
      cntry_id = if_else(is.na(cntry_id), 900L + dense_rank(iso3), cntry_id)
    ) |>
    arrange(iso3, region_code) |>
    group_by(iso3) |>
    mutate(adm_id = as.integer(cntry_id) * 1000L + row_number()) |>
    ungroup() |>
    arrange(row_id)

  gadm1_fb$iso3 <- gadm_df$iso3
  gadm1_fb$region_name <- gadm_df$region_name
  gadm1_fb$adm_id <- gadm_df$adm_id
  gadm1_fb$source <- "gadm_fallback"

  gadm_core <- gadm1_fb
} else {
  gadm_core <- NULL
}

best_core <- best_core[, c("iso3", "region_name", "adm_id", "source")]
second_core <- second_core[, c("iso3", "region_name", "adm_id", "source")]

if (!is.null(gadm_core)) {
  gadm_core <- gadm_core[, c("iso3", "region_name", "adm_id", "source")]
  adm1_global <- rbind(best_core, gadm_core)
} else {
  adm1_global <- best_core
}

adm1_extract <- if (!is.null(gadm_core)) {
  rbind(best_core, second_core, gadm_core)
} else {
  rbind(best_core, second_core)
}

extract_df <- as.data.frame(adm1_extract) |>
  mutate(poly_idx = row_number()) |>
  group_by(adm_id) |>
  slice_head(n = 1) |>
  ungroup()
adm1_extract <- adm1_extract[extract_df$poly_idx, ]

global_df <- as.data.frame(adm1_global) |>
  mutate(poly_idx = row_number()) |>
  group_by(adm_id) |>
  slice_head(n = 1) |>
  ungroup()
adm1_global <- adm1_global[global_df$poly_idx, ]

writeVector(adm1_global, gpkg_path, overwrite = TRUE)
lookup <- as.data.frame(adm1_global) |>
  select(iso3, region_name, adm_id, source) |>
  distinct()
write_csv(lookup, lookup_path)

message(sprintf("  Global ADM1 polygons: %d", nrow(adm1_global)))
message(sprintf("  ADM1 polygons for predictor extraction: %d (best + second + GADM)", nrow(adm1_extract)))
message(sprintf("Saved: %s", gpkg_path))
message(sprintf("Saved: %s", lookup_path))

worker_extract_predictor <- function(spec,
                                     pop_path,
                                     years_pred,
                                     n_years,
                                     pop_mat,
                                     polygons_path,
                                     polygons_layer) {
  polygons_sf_worker <- sf::st_read(polygons_path, layer = polygons_layer, quiet = TRUE)
  pop_stack_worker <- terra::rast(pop_path)[[6:39]]
  names(pop_stack_worker) <- paste0("y", years_pred)

  message(sprintf("  Extracting %s (%d layers x %d polygons) ...",
                  spec$name, terra::nlyr(pop_stack_worker), nrow(polygons_sf_worker)))

  r <- terra::rast(spec$path)
  if (is.na(spec$start)) {
    r1 <- terra::resample(r, pop_stack_worker[[1]], method = "bilinear")
    var_stack <- r1[[rep(1, n_years)]]
    rm(r1)
  } else {
    idx_start <- min(years_pred) - spec$start + 1
    idx_end <- max(years_pred) - spec$start + 1
    n_avail <- terra::nlyr(r)
    idx_end_clamped <- min(idx_end, n_avail)
    var_stack <- r[[idx_start:idx_end_clamped]]
    if (terra::nlyr(var_stack) < n_years) {
      n_pad <- n_years - terra::nlyr(var_stack)
      message(sprintf("    %s: padding last layer for %d extra year(s)",
                      spec$name, n_pad))
      var_stack <- c(var_stack, r[[rep(n_avail, n_pad)]])
    }
  }
  names(var_stack) <- paste0("y", years_pred)
  rm(r)

  if (isTRUE(spec$log10)) {
    var_stack <- log10(var_stack)
    var_stack[!is.finite(var_stack)] <- NA
  }
  if (isTRUE(spec$log)) {
    var_stack <- log(var_stack)
    var_stack[!is.finite(var_stack)] <- NA
  }

  numerator <- exactextractr::exact_extract(
    var_stack * pop_stack_worker,
    polygons_sf_worker,
    fun = "sum",
    progress = FALSE
  )
  result <- as.data.frame(as.matrix(numerator) / pop_mat)
  colnames(result) <- paste0(spec$name, "_", years_pred)

  rm(var_stack, pop_stack_worker, polygons_sf_worker)
  gc(verbose = FALSE)
  result
}

var_specs <- list(
  list(path = "/Users/mkummu/R/subnat_poverty/results_gis/rast_poverty_combined_income_weighted_1990_2022.tif",
       name = "poverty", start = 1990),
  list(path = "/Users/mkummu/R/hdi_subnat/results/rast_msch_1990_2023.tif",
       name = "msch", start = 1990),
  list(path = "/Users/mkummu/R/subnatGNI/results/rast_gnic_1990_2023.tif",
       name = "gnic", start = 1990, log = TRUE, log10 = FALSE),
  list(path = "/Users/mkummu/R/subnatGini/results/rast_gini_disp_1990_2023.tif",
       name = "gini", start = 1990),
  list(path = "/Users/mkummu/R/misc/percentileNormalisation/output/urbanisationCntryWise_GHS2023a.tif",
       name = "urban", start = 1985),
  list(path = "/Users/mkummu/R/GIS_data_common/travel_time_cities/travel_time_to_cities_11.tif",
       name = "travel_time", start = NA, log10 = TRUE),
  list(path = "/Users/mkummu/R/migration_data_bee/results/aridity_index.tif",
       name = "aridity", start = NA),
  list(path = "/Users/mkummu/R/hdi_subnat/results/rast_gdi_1990_2023.tif",
       name = "gdi", start = 1990),
  list(path = "/Users/mkummu/R/vdem_democracy_dataPrep/results_gis/rast_vdem_adm1_v2x_corr_1990_2024.tif",
       name = "corruption", start = 1990),
  list(path = "/Users/mkummu/R/GIS_data_common/GLOBGM_groundwater/globgm-wtd-ss.tif",
       name = "gw_depth", start = NA)
)

extract_all_vars_combined <- function(polygons_sf, set_id, pop_stack, cache_path, key_label) {
  cache_meta_path <- sub("\\.csv$", "_meta.rds", cache_path)
  cache_sig <- list(
    n_polygons = nrow(polygons_sf),
    n_adm0 = sum(set_id == "adm0"),
    n_adm1 = sum(set_id == "adm1")
  )

  message(sprintf("== Extracting %s predictors from rasters ==", key_label))
  message("  Extracting population (34 layers x polygons) ...")
  pop_df <- exact_extract(pop_stack, polygons_sf, fun = "sum", progress = FALSE)
  pop_mat <- as.matrix(pop_df)
  colnames(pop_df) <- paste0("pop_", years_pred)
  result <- as.data.frame(pop_df)

  polygons_path <- tempfile(pattern = "combined_polygons_", fileext = ".gpkg")
  polygons_layer <- "combined_polygons"
  sf::st_write(polygons_sf, polygons_path, layer = polygons_layer,
               delete_dsn = TRUE, quiet = TRUE)
  on.exit(unlink(polygons_path), add = TRUE)

  message(sprintf("  Extracting %d predictors in parallel (%d cores) ...",
                  length(var_specs), n_cores))

  predictor_results <- future_map(
    var_specs,
    worker_extract_predictor,
    pop_path = pop_path,
    years_pred = years_pred,
    n_years = n_years,
    pop_mat = pop_mat,
    polygons_path = polygons_path,
    polygons_layer = polygons_layer,
    .options = furrr_options(
      seed = TRUE,
      packages = c("terra", "sf", "exactextractr"),
      globals = list(
        worker_extract_predictor = worker_extract_predictor,
        pop_path = pop_path,
        years_pred = years_pred,
        n_years = n_years,
        pop_mat = pop_mat,
        polygons_path = polygons_path,
        polygons_layer = polygons_layer
      )
    )
  )

  result <- bind_cols(result, predictor_results)

  result$set_id <- set_id
  saveRDS(cache_sig, cache_meta_path)
  write_csv(result, cache_path)

  result
}

message("== Loading population stack (1990-2023) ==")
r_pop_full <- rast(pop_path)
r_pop <- r_pop_full[[6:39]]
names(r_pop) <- paste0("y", years_pred)
rm(r_pop_full)
gc(verbose = FALSE)

message("== Building combined admin polygon object (ADM0 + ADM1 extract set) ==")
adm0_core <- gadm0
adm0_core$iso3 <- as.character(gadm0$GID_0)
adm0_core$region_name <- as.character(gadm0$COUNTRY)
adm0_core$adm_id <- adm0_attrs$adm_id
adm0_core$source <- "adm0"
adm0_core <- adm0_core[, c("iso3", "region_name", "adm_id", "source")]

combined_polygons <- rbind(adm0_core, adm1_extract)
set_id <- c(rep("adm0", nrow(adm0_core)), rep("adm1", nrow(adm1_extract)))
combined_polygons_sf <- st_as_sf(combined_polygons)

message(sprintf("  Combined polygons: %d total (%d ADM0 + %d ADM1)",
                nrow(combined_polygons), nrow(adm0_core), nrow(adm1_extract)))

combined_wide <- extract_all_vars_combined(
  polygons_sf = combined_polygons_sf,
  set_id = set_id,
  pop_stack = r_pop,
  cache_path = combined_cache,
  key_label = "combined admin"
)

adm0_wide <- combined_wide |>
  filter(set_id == "adm0") |>
  select(-set_id)
adm1_wide <- combined_wide |>
  filter(set_id == "adm1") |>
  select(-set_id)

adm0_long <- adm0_wide |>
  mutate(poly_idx = seq_len(n())) |>
  pivot_longer(cols = -poly_idx,
               names_to = c("variable", "year"),
               names_pattern = "(.+)_(\\d{4})",
               values_to = "value") |>
  mutate(year = as.integer(year)) |>
  pivot_wider(names_from = variable, values_from = value) |>
  left_join(adm0_attrs, by = "poly_idx") |>
  select(adm_id, iso3, country_name_gadm, year, pop, all_of(predictor_cols)) |>
  arrange(iso3, year)

write_csv(adm0_long, adm0_cache)
message(sprintf("Saved: %s (%d rows)", adm0_cache, nrow(adm0_long)))

adm1_attrs <- as.data.frame(adm1_extract) |>
  transmute(poly_idx = row_number(),
            adm_id, iso3, region_name, source)

adm1_long <- adm1_wide |>
  mutate(poly_idx = seq_len(n())) |>
  pivot_longer(cols = -poly_idx,
               names_to = c("variable", "year"),
               names_pattern = "(.+)_(\\d{4})",
               values_to = "value") |>
  mutate(year = as.integer(year)) |>
  pivot_wider(names_from = variable, values_from = value) |>
  left_join(adm1_attrs, by = "poly_idx") |>
  select(adm_id, iso3, region_name, source, year, pop, all_of(predictor_cols)) |>
  arrange(adm_id, year)

write_csv(adm1_long, adm1_cache)
message(sprintf("Saved: %s (%d rows)", adm1_cache, nrow(adm1_long)))

stopifnot(nrow(adm0_long) == nrow(gadm0) * n_years)
stopifnot(nrow(adm1_long) == nrow(adm1_extract) * n_years)
stopifnot(all(predictor_cols %in% names(adm0_long)))
stopifnot(all(predictor_cols %in% names(adm1_long)))

stopifnot(all(adm0_long$adm_id < 1000L))
stopifnot(all(adm1_long$adm_id >= 1000L))
stopifnot(!any(is.na(adm0_long$adm_id)))
stopifnot(!any(is.na(adm1_long$adm_id)))
stopifnot(n_distinct(adm0_long$adm_id) == nrow(gadm0))
stopifnot(n_distinct(adm1_long$adm_id) == nrow(adm1_extract))

message("All sanity checks passed.")

plan(sequential)
message("Parallel processing complete.")
