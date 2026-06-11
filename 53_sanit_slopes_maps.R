# 53_sanit_slopes_maps.R
# ---------------------------------------------------------------------------
# Theil-Sen slope calculation and tmap maps for sanitation indicators
# at ADM1 resolution (1990-2023).
#
# Part 1 (sections 1-3): Slope calculation
# Part 2 (sections 4-6): tmap maps
#
# Inputs:
#   results_sanit/sanit_adm1_ratio_adjusted.csv
#   results_gis_sanit/sanit_adm1_values_1990_2023_wide.csv  – script 52
#   results_gis_sanit/polyg_sanit_adm1_1990_2023_geomonly_simpl.gpkg – script 05b
#   results_gis_sanit/sanit_adm1_template_5arcmin.tif      – script 05b
#
# Outputs:
#   results_gis_sanit/rast_slope_sanit_{indicator}_1990_2023.tif
#   figures_sanit/maps/sanit_{indicator}_{year}.png
#   figures_sanit/maps/sanit_{indicator}_slope.png
#   figures_sanit/legends/sanit_legend_{indicator}_{year}.png
#   figures_sanit/legends/sanit_legend_{indicator}_slope.png
# ---------------------------------------------------------------------------

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

### 0. Setup ----
years_in   <- 1990:2023
first_year <- years_in[1]
last_year  <- years_in[length(years_in)]
indicators <- c("s_bas", "s_lim", "s_unimp", "s_od")

dir.create("figures_sanit/maps",    recursive = TRUE, showWarnings = FALSE)
dir.create("figures_sanit/legends", recursive = TRUE, showWarnings = FALSE)

csv_wide_path <- sprintf("results_gis_sanit/sanit_adm1_values_%d_%d_wide.csv",
                         first_year, last_year)
gpkg_simpl_path <- sprintf("results_gis_sanit/polyg_sanit_adm1_%d_%d_geomonly_simpl.gpkg",
                           first_year, last_year)

### 1. Load data ----
message("== Loading data ==")

sanit_data <- read_csv("results_sanit/sanit_adm1_ratio_adjusted.csv",
                       show_col_types = FALSE)

r_template <- rast("results_gis_sanit/sanit_adm1_template_5arcmin.tif")

message(sprintf("  ADM1 units: %d", n_distinct(sanit_data$adm_id)))
message(sprintf("  Years: %d-%d", first_year, last_year))

### 2. Slope calculation function ----

calc_slopes <- function(indicator, data) {
  message(sprintf("\n--- Slopes: %s ---", indicator))

  vals <- data |>
    select(adm_id, year, value = all_of(indicator)) |>
    filter(!is.na(value)) |>
    distinct(adm_id, year, .keep_all = TRUE)

  slopes <- vals |>
    group_by(adm_id) |>
    mutate(time = row_number()) |>
    ungroup() |>
    nest(.by = adm_id) |>
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
    select(adm_id, slope_est = estimate, p.value)

  message(sprintf("  Slopes computed for %d regions", nrow(slopes)))

  slopes |>
    mutate(
      p.value_sig = p.value < 0.05,
      slope       = if_else(p.value_sig, slope_est, NA_real_)
    )
}

### 3. Compute slopes for all indicators ----
message("== Computing Theil-Sen slopes ==")

slope_list <- list()
for (ind in indicators) {
  slope_list[[ind]] <- calc_slopes(ind, sanit_data)
}

# Save slopes as CSV
message("\n== Saving slope results ==")
slope_df <- bind_rows(slope_list, .id = "indicator")
write_csv(slope_df, sprintf("results_sanit/sanit_slopes_adm1_%d_%d.csv",
                            first_year, last_year))
message("  Saved to results_sanit/")

# Slope rasters
message("== Creating slope rasters ==")
for (ind in indicators) {
  slopes_ind <- slope_list[[ind]]

  tmpl <- r_template
  all_ids    <- as.numeric(na.omit(unique(values(tmpl))))
  used_ids   <- as.numeric(unique(slopes_ind$adm_id))
  missing    <- all_ids[!all_ids %in% used_ids]
  if (length(missing) > 0) tmpl[tmpl %in% missing] <- NA

  r_slope <- classify(tmpl,
                      cbind(as.numeric(slopes_ind$adm_id),
                            as.numeric(slopes_ind$slope_est)),
                      others = NA)
  names(r_slope) <- paste0(ind, "_slope_", first_year, "_", last_year)

  rast_path <- sprintf("results_gis_sanit/rast_slope_sanit_%s_%d_%d.tif",
                        ind, first_year, last_year)
  writeRaster(r_slope, rast_path, gdal = "COMPRESS=LZW", overwrite = TRUE)
  message(sprintf("  Saved: %s", rast_path))
}

### 4. Define tmap plotting functions ----

adm0_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
  filter(iso_a3 != "ATA")

pal_level <- cols4all::c4a("bilbao", reverse = TRUE)
pal_slope <- cols4all::c4a("bu_br_div", reverse = FALSE)

f_map_level <- function(sf_in, data_col, mybreaks, pal, title_label) {
  tm_shape(sf_in) +
    tm_fill(
      fill       = data_col,
      fill.scale = tm_scale(
        values   = pal,
        breaks   = mybreaks,
        value.na = "white"
      ),
      fill.legend = tm_legend(
        orientation = "landscape",
        position    = tm_pos_out("center", "bottom"),
        frame       = FALSE,
        title       = title_label,
        title.size  = 1.5,
        text.size   = 1.2
      )
    ) +
    tm_shape(adm0_sf) +
    tm_borders(col = "grey35", lwd = 0.025) +
    tm_crs("+proj=robin") +
    tm_layout(main.title.position = "center", frame = FALSE)
}

f_map_slope <- function(sf_in, slope_col, mybreaks, pal, title_label) {
  tm_shape(sf_in) +
    tm_fill(
      fill       = slope_col,
      fill.scale = tm_scale(
        values   = pal,
        breaks   = mybreaks,
        midpoint = 0,
        value.na = "white"
      ),
      fill.legend = tm_legend(
        orientation = "landscape",
        position    = tm_pos_out("center", "bottom"),
        frame       = FALSE,
        title       = title_label
      )
    ) +
    tm_shape(adm0_sf) +
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
    tm_layout(legend.only = TRUE, frame = FALSE,
              outer.margins = 0, inner.margins = c(0, 0, 0, 0), font.size = 1.5)
  tmap_save(p_leg, filename = filename, width = width, height = height,
            units = units, dpi = dpi)
}

### 5. Level maps (first year and last year) ----
message("\n== Creating level maps ==")

mybreaks_level <- seq(0, 50, 5)

ind_titles <- c(
  s_bas   = "Basic sanitation service (%)",
  s_lim   = "Limited sanitation service (%)",
  s_unimp = "Unimproved sanitation (%)",
  s_od    = "Open defecation (%)"
)

gpkg_sf <- st_read(gpkg_simpl_path, quiet = TRUE) |>
  select(adm_id)

{
  full_attrs <- read_csv(csv_wide_path, show_col_types = FALSE)
  gpkg_sf <- gpkg_sf |> left_join(full_attrs, by = "adm_id")
}

# Join slope columns
for (ind in indicators) {
  slope_col  <- paste0(ind, "_slope")
  pval_col   <- paste0(ind, "_p.value")
  sig_col    <- paste0(ind, "_sig")

  join_df <- slope_list[[ind]] |>
    transmute(
      adm_id,
      !!slope_col := slope,
      !!pval_col  := p.value,
      !!sig_col   := p.value_sig
    )
  gpkg_sf <- gpkg_sf |>
    select(-any_of(c(slope_col, pval_col, sig_col))) |>
    left_join(join_df, by = "adm_id")
}

for (ind in indicators) {
  message(sprintf("  Level maps: %s", ind))

  for (yr in c(first_year, last_year)) {
    col_name <- paste0(ind, "_", yr)

    p <- f_map_level(
      sf_in       = gpkg_sf,
      data_col    = col_name,
      mybreaks    = mybreaks_level,
      pal         = pal_level,
      title_label = ind_titles[[ind]]
    )

    save_map(p,
             sprintf("figures_sanit/maps/sanit_%s_%d.png", ind, yr))
    save_legend(p,
                sprintf("figures_sanit/legends/sanit_legend_%s_%d.png", ind, yr))

    message(sprintf("    Saved: sanit_%s_%d.png", ind, yr))
  }
}

### 6. Slope maps ----
message("\n== Creating slope maps ==")

ind_slope_titles <- c(
  s_bas   = "Trend, basic sanitation service (pp per year)",
  s_lim   = "Trend, limited sanitation service (pp per year)",
  s_unimp = "Trend, unimproved sanitation (pp per year)",
  s_od    = "Trend, open defecation (pp per year)"
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
  plyr::round_any( abs_max, accuracy = accur, f = ceiling),
  by = plyr::round_any(abs_max * 2 / 20, accuracy = 0.00001, f = floor)
)

message(sprintf("  Slope range (5th-95th pct): %.4f to %.4f",
                -abs_max, abs_max))

for (ind in indicators) {
  message(sprintf("  Slope map: %s", ind))

  slope_col <- paste0(ind, "_slope")

  p <- f_map_slope(
    sf_in       = gpkg_sf,
    slope_col   = slope_col,
    mybreaks    = slope_breaks,
    pal         = pal_slope,
    title_label = ind_slope_titles[[ind]]
  )

  save_map(p,
           sprintf("figures_sanit/maps/sanit_%s_slope.png", ind))
  save_legend(p,
              sprintf("figures_sanit/legends/sanit_legend_%s_slope.png", ind),
              width = 150, height = 30)

  message(sprintf("    Saved: sanit_%s_slope.png", ind))
}

message("\nDone.")
