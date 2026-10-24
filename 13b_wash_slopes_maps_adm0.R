# 13b_wash_slopes_maps_adm0.R
# ---------------------------------------------------------------------------
# Theil-Sen slope calculation and tmap maps for WASH water indicators
# at ADM0 resolution (1990–2023).
#
# Part 1 (sections 1–3): Slope calculation
#   - Reads merged ADM0 filled data from script 10b
#   - Computes mblm Theil-Sen slope per iso3 for each indicator
#   - Saves slope rasters on ADM0 5-arcmin template
#
# Part 2 (sections 4–6): tmap maps
#   - First-year (1990) and last-year (2023) level maps
#   - Slope maps (diverging palette, non-significant -> NA/white)
#
# Inputs:
#   results/JMP_adm0_water_merged_1990_2023.csv
#   /Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg
#   /Users/mkummu/R/GIS_data_common/GHS_POP/r_pop_GHS_1985_2025_5arcmin.tif
#
# Outputs:
#   results/wash_slopes_adm0_1990_2023.csv
#   results_gis/wash_adm0_template_5arcmin.tif
#   results_gis/wash_adm0_lookup.csv
#   results_gis/wash_adm0_values_1990_2023_long.csv
#   results_gis/wash_adm0_values_1990_2023_wide.csv
#   results_gis/rast_slope_wash_adm0_{indicator}_1990_2023.tif
#   figures/maps/wash_adm0_{indicator}_{year}.png
#   figures/maps/wash_adm0_{indicator}_slope.png
#   figures/legends/wash_adm0_legend_{indicator}_{year}.png
#   figures/legends/wash_adm0_legend_{indicator}_slope.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(tmap)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(broom)
  library(mblm)
  library(cols4all)
})

### 0. Setup ----
source("config.R")
paths <- subnat_wash_paths()
years_in <- 1990:2023
first_year <- years_in[1]
last_year <- years_in[length(years_in)]
ref_year <- 2019
indicators <- c("Basic", "Limited", "Unimproved", "Surface.water")

gadm_path <- paths$gadm
pop_path <- paths$population
adm0_data_path <- "results/JMP_adm0_water_merged_1990_2023.csv"
adm0_simpl_path <- "results_gis/polyg_wash_adm0_geomonly_simpl.gpkg"

dir.create("results", showWarnings = FALSE)
dir.create("results_gis", showWarnings = FALSE)
dir.create("figures/maps", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/legends", recursive = TRUE, showWarnings = FALSE)

if (!file.exists(adm0_data_path)) {
  stop("Missing input: ", adm0_data_path,
       "\nRun 10b_merge_adm0_historical_modern.R first.")
}
if (!file.exists(gadm_path)) stop("Missing input: ", gadm_path)
if (!file.exists(pop_path)) stop("Missing input: ", pop_path)
if (!file.exists(adm0_simpl_path)) {
  stop("Missing input: ", adm0_simpl_path,
       "\nRun 05b_build_adm1_shared_assets.R first.")
}

### 1. Load ADM0 data and geometry ----
message("== Loading ADM0 data and geometry ==")

wash_data <- read_csv(adm0_data_path, show_col_types = FALSE) |>
  filter(year %in% years_in) |>
  transmute(
    iso3 = toupper(as.character(iso3)),
    country_name = as.character(country_name),
    year = as.integer(year),
    Basic = as.numeric(Basic_filled),
    Limited = as.numeric(Limited_filled),
    Unimproved = as.numeric(Unimproved_filled),
    Surface.water = as.numeric(Surface.water_filled)
  )

lay_info <- st_layers(gadm_path)
adm0_layer <- lay_info$name[grepl("^ADM_0$|adm0|level0|level_0", lay_info$name, ignore.case = TRUE)][1]
if (is.na(adm0_layer)) stop("Could not identify ADM0 layer in: ", gadm_path)

adm0_sf <- st_read(gadm_path, layer = adm0_layer, quiet = TRUE) |>
  select(GID_0, COUNTRY) |>
  mutate(
    iso3 = toupper(as.character(GID_0)),
    country_name = as.character(COUNTRY),
    adm0_id = row_number()
  ) |>
  select(adm0_id, iso3, country_name)

adm0_lookup <- adm0_sf |>
  st_drop_geometry() |>
  select(adm0_id, iso3, country_name)

message(sprintf("  Country polygons: %d", nrow(adm0_sf)))
message(sprintf("  Countries in data: %d", n_distinct(wash_data$iso3)))
message(sprintf("  Years: %d-%d", min(wash_data$year, na.rm = TRUE), max(wash_data$year, na.rm = TRUE)))

### 2. Build ADM0 5-arcmin template raster ----
message("== Building ADM0 5-arcmin template ==")

r_pop_full <- rast(pop_path)
r_ref_5 <- r_pop_full[[1]]

r_ref_1 <- rast(
  nrows = nrow(r_ref_5) * 5,
  ncols = ncol(r_ref_5) * 5,
  ext = ext(r_ref_5),
  crs = crs(r_ref_5)
)

adm0_sv <- vect(adm0_sf)
r_adm0_1 <- rasterize(adm0_sv, r_ref_1, field = "adm0_id")
r_template <- aggregate(r_adm0_1, fact = 5, fun = "modal", na.rm = TRUE)
names(r_template) <- "adm0_id"

template_path <- "results_gis/wash_adm0_template_5arcmin.tif"
writeRaster(
  r_template,
  template_path,
  datatype = "INT4S",
  gdal = "COMPRESS=LZW",
  overwrite = TRUE
)

write_csv(adm0_lookup, "results_gis/wash_adm0_lookup.csv")
message(sprintf("  Saved: %s", template_path))
message("  Saved: results_gis/wash_adm0_lookup.csv")

### 3. Slope calculation function ----
calc_slopes <- function(indicator, data) {
  message(sprintf("\n--- Slopes: %s ---", indicator))

  vals <- data |>
    select(iso3, year, value = all_of(indicator)) |>
    filter(!is.na(value)) |>
    distinct(iso3, year, .keep_all = TRUE)

  slopes <- vals |>
    group_by(iso3) |>
    mutate(time = row_number()) |>
    ungroup() |>
    nest(.by = iso3) |>
    mutate(
      model = map(data, \(d) {
        if (nrow(d) < 3) return(NULL)
        tryCatch(
          mblm(value ~ time, data = d, repeated = FALSE),
          error = \(e) NULL
        )
      }),
      tidy = map(model, \(m) {
        if (is.null(m)) return(NULL)
        tryCatch(tidy(m), error = \(e) NULL)
      })
    ) |>
    filter(!map_lgl(tidy, is.null)) |>
    unnest(tidy) |>
    filter(term == "time") |>
    select(iso3, slope_est = estimate, p.value)

  message(sprintf("  Slopes computed for %d countries", nrow(slopes)))

  slopes |>
    mutate(
      p.value_sig = p.value < 0.05,
      slope = if_else(p.value_sig, slope_est, NA_real_)
    )
}

### 4. Compute slopes for all indicators ----
message("== Computing Theil-Sen slopes ==")

slope_list <- list()
for (ind in indicators) {
  slope_list[[ind]] <- calc_slopes(ind, wash_data)
}

slope_years <- years_in[years_in <= ref_year]
wash_data_ref <- wash_data |>
  filter(year %in% slope_years)

slope_list_ref <- list()
for (ind in indicators) {
  slope_list_ref[[ind]] <- calc_slopes(ind, wash_data_ref)
}

slope_df <- bind_rows(slope_list, .id = "indicator")
write_csv(slope_df, sprintf("results/wash_slopes_adm0_%d_%d.csv", first_year, last_year))
write_csv(
  bind_rows(slope_list_ref, .id = "indicator"),
  sprintf("results/wash_slopes_adm0_%d_%d.csv", first_year, ref_year)
)
message(sprintf("  Saved: results/wash_slopes_adm0_%d_%d.csv", first_year, last_year))

### 5. Save ADM0 value attributes (long + wide) ----
message("== Saving ADM0 attributes ==")

adm0_long <- wash_data |>
  arrange(iso3, year)

adm0_wide <- wash_data |>
  pivot_longer(
    cols = all_of(indicators),
    names_to = "indicator",
    values_to = "value"
  ) |>
  mutate(var = paste0(indicator, "_", year)) |>
  select(iso3, var, value) |>
  pivot_wider(names_from = var, values_from = value)

write_csv(adm0_long, sprintf("results_gis/wash_adm0_values_%d_%d_long.csv", first_year, last_year))
write_csv(adm0_wide, sprintf("results_gis/wash_adm0_values_%d_%d_wide.csv", first_year, last_year))

### 6. Create slope rasters ----
message("== Creating slope rasters ==")

for (ind in indicators) {
  slopes_ind <- slope_list[[ind]] |>
    left_join(adm0_lookup |> select(adm0_id, iso3), by = "iso3")

  tmpl <- r_template
  all_ids <- as.numeric(na.omit(unique(values(tmpl))))
  used_ids <- as.numeric(unique(slopes_ind$adm0_id))
  missing <- all_ids[!all_ids %in% used_ids]
  if (length(missing) > 0) tmpl[tmpl %in% missing] <- NA

  r_slope <- classify(
    tmpl,
    cbind(as.numeric(slopes_ind$adm0_id), as.numeric(slopes_ind$slope_est)),
    others = NA
  )
  names(r_slope) <- paste0(ind, "_slope_", first_year, "_", last_year)

  rast_path <- sprintf("results_gis/rast_slope_wash_adm0_%s_%d_%d.tif", ind, first_year, last_year)
  writeRaster(r_slope, rast_path, gdal = "COMPRESS=LZW", overwrite = TRUE)
  message(sprintf("  Saved: %s", rast_path))
}

### 7. Define tmap plotting functions ----
# Country borders overlay
adm0_border_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
  filter(iso_a3 != "ATA")

# Sequential palette for WASH levels (low coverage = dark)
pal_level <- cols4all::c4a("bilbao", reverse = TRUE)
# Diverging palette for slopes
pal_slope <- cols4all::c4a("bu_br_div", reverse = FALSE)

f_map_level <- function(sf_in, data_col, mybreaks, pal, title_label) {
  tm_shape(sf_in) +
    tm_fill(
      fill = data_col,
      fill.scale = tm_scale(
        values = pal,
        breaks = mybreaks,
        value.na = "white"
      ),
      fill.legend = tm_legend(
        orientation = "landscape",
        position = tm_pos_out("center", "bottom"),
        frame = FALSE,
        title = title_label,
        title.size = 1.5,
        text.size = 1.2
      )
    ) +
    tm_shape(adm0_border_sf) +
    tm_borders(col = "grey35", lwd = 0.025) +
    tm_crs("+proj=robin") +
    tm_layout(main.title.position = "center", frame = FALSE)
}

f_map_slope <- function(sf_in, slope_col, mybreaks, pal, title_label) {
  tm_shape(sf_in) +
    tm_fill(
      fill = slope_col,
      fill.scale = tm_scale(
        values = pal,
        breaks = mybreaks,
        midpoint = 0,
        value.na = "white"
      ),
      fill.legend = tm_legend(
        orientation = "landscape",
        position = tm_pos_out("center", "bottom"),
        frame = FALSE,
        title = title_label
      )
    ) +
    tm_shape(adm0_border_sf) +
    tm_borders(col = "grey35", lwd = 0.025) +
    tm_crs("+proj=robin") +
    tm_layout(main.title.position = "center", frame = FALSE)
}

save_map <- function(p_map, filename, width = 80, units = "mm", dpi = 450) {
  p_no_legend <- p_map + tm_layout(legend.show = FALSE)
  tmap_save(p_no_legend, filename = filename, width = width, units = units, dpi = dpi)
}

save_legend <- function(p_map, filename, width = 120, height = 20,
                        units = "mm", dpi = 600) {
  p_leg <- p_map +
    tm_layout(
      legend.only = TRUE,
      frame = FALSE,
      outer.margins = 0,
      inner.margins = c(0, 0, 0, 0),
      font.size = 1.5
    )
  tmap_save(p_leg, filename = filename, width = width, height = height,
            units = units, dpi = dpi)
}

### 8. Prepare map sf with levels and slopes ----
message("== Preparing map attributes ==")

# Use simplified polygons for fast mapping
gpkg_sf <- st_read(adm0_simpl_path, quiet = TRUE) |>
  select(adm0_id, iso3, country_name)

wide_attrs <- read_csv(
  sprintf("results_gis/wash_adm0_values_%d_%d_wide.csv", first_year, last_year),
  show_col_types = FALSE
)

gpkg_sf <- gpkg_sf |>
  left_join(wide_attrs, by = "iso3")

for (ind in indicators) {
  slope_col <- paste0(ind, "_slope")
  pval_col <- paste0(ind, "_p.value")
  sig_col <- paste0(ind, "_sig")

  join_df <- slope_list[[ind]] |>
    transmute(
      iso3,
      !!slope_col := slope,
      !!pval_col := p.value,
      !!sig_col := p.value_sig
    )

  gpkg_sf <- gpkg_sf |>
    select(-any_of(c(slope_col, pval_col, sig_col))) |>
    left_join(join_df, by = "iso3")
}

### 9. Level maps (2019 and 2023) ----
message("== Creating ADM0 level maps ==")

mybreaks_level <- c(0, 1, 2, 4, 8, 16, 32, 64)

ind_titles <- c(
  Basic = "Basic water service (%)",
  Limited = "Limited water service (%)",
  Unimproved = "Unimproved water source (%)",
  Surface.water = "Surface water use (%)"
)

for (ind in indicators) {
  message(sprintf("  Level maps: %s", ind))

  map_years <- c(2019, 2023)
  map_years <- map_years[map_years %in% years_in]
  map_years <- unique(c(map_years, last_year))

  for (yr in map_years) {
    col_name <- paste0(ind, "_", yr)

    p <- f_map_level(
      sf_in = gpkg_sf,
      data_col = col_name,
      mybreaks = mybreaks_level,
      pal = pal_level,
      title_label = ind_titles[[ind]]
    )

    save_map(p, sprintf("figures/maps/wash_adm0_%s_%d.png", ind, yr))
    save_legend(p, sprintf("figures/legends/wash_adm0_legend_%s_%d.png", ind, yr))

    message(sprintf("    Saved: wash_adm0_%s_%d.png", ind, yr))
  }
}

### 10. Slope maps ----
message("== Creating ADM0 slope maps ==")

ind_slope_titles <- c(
  Basic = "Trend, basic water service (pp per year)",
  Limited = "Trend, limited water service (pp per year)",
  Unimproved = "Trend, unimproved water source (pp per year)",
  Surface.water = "Trend, surface water use (pp per year)"
)

all_slopes <- numeric(0)
for (ind in indicators) {
  all_slopes <- c(all_slopes, gpkg_sf[[paste0(ind, "_slope")]])
}

min_s <- quantile(all_slopes, 0.05, na.rm = TRUE)
max_s <- quantile(all_slopes, 0.95, na.rm = TRUE)
abs_max <- max(abs(c(min_s, max_s)))
accur <- round((abs_max * 2) / 20, digits = 5)

slope_breaks <- seq(
  plyr::round_any(-abs_max, accuracy = accur, f = floor),
  plyr::round_any(abs_max, accuracy = accur, f = ceiling),
  by = plyr::round_any(abs_max * 2 / 20, accuracy = 0.00001, f = floor)
)

message(sprintf("  Slope range (5th-95th pct): %.4f to %.4f", -abs_max, abs_max))

for (ind in indicators) {
  message(sprintf("  Slope map: %s", ind))

  slope_col <- paste0(ind, "_slope")

  p <- f_map_slope(
    sf_in = gpkg_sf,
    slope_col = slope_col,
    mybreaks = slope_breaks,
    pal = pal_slope,
    title_label = ind_slope_titles[[ind]]
  )

  save_map(p, sprintf("figures/maps/wash_adm0_%s_slope.png", ind))
  save_legend(
    p,
    sprintf("figures/legends/wash_adm0_legend_%s_slope.png", ind),
    width = 150,
    height = 30
  )

  message(sprintf("    Saved: wash_adm0_%s_slope.png", ind))
}

# Additional slope maps for 1990-2019
gpkg_sf_ref <- st_read(adm0_simpl_path, quiet = TRUE) |>
  select(adm0_id, iso3, country_name)

for (ind in indicators) {
  slope_col_ref <- paste0(ind, "_slope_", first_year, "_", ref_year)

  join_df_ref <- slope_list_ref[[ind]] |>
    transmute(
      iso3,
      !!slope_col_ref := slope
    )

  gpkg_sf_ref <- gpkg_sf_ref |>
    select(-any_of(slope_col_ref)) |>
    left_join(join_df_ref, by = "iso3")
}

all_slopes_ref <- numeric(0)
for (ind in indicators) {
  all_slopes_ref <- c(
    all_slopes_ref,
    gpkg_sf_ref[[paste0(ind, "_slope_", first_year, "_", ref_year)]]
  )
}

min_s_ref <- quantile(all_slopes_ref, 0.05, na.rm = TRUE)
max_s_ref <- quantile(all_slopes_ref, 0.95, na.rm = TRUE)
abs_max_ref <- max(abs(c(min_s_ref, max_s_ref)))
accur_ref <- round((abs_max_ref * 2) / 20, digits = 5)

slope_breaks_ref <- seq(
  plyr::round_any(-abs_max_ref, accuracy = accur_ref, f = floor),
  plyr::round_any(abs_max_ref, accuracy = accur_ref, f = ceiling),
  by = plyr::round_any(abs_max_ref * 2 / 20, accuracy = 0.00001, f = floor)
)

for (ind in indicators) {
  message(sprintf("  Slope map (%d-%d): %s", first_year, ref_year, ind))

  slope_col_ref <- paste0(ind, "_slope_", first_year, "_", ref_year)

  p_ref <- f_map_slope(
    sf_in = gpkg_sf_ref,
    slope_col = slope_col_ref,
    mybreaks = slope_breaks_ref,
    pal = pal_slope,
    title_label = ind_slope_titles[[ind]]
  )

  save_map(p_ref, sprintf("figures/maps/wash_adm0_%s_slope_%d_%d.png", ind, first_year, ref_year))
  save_legend(
    p_ref,
    sprintf("figures/legends/wash_adm0_legend_%s_slope_%d_%d.png", ind, first_year, ref_year),
    width = 150,
    height = 30
  )

  message(sprintf("    Saved: wash_adm0_%s_slope_%d_%d.png", ind, first_year, ref_year))
}

message("\nDone.")
