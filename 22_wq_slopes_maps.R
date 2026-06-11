# 22_wq_slopes_maps.R
# ---------------------------------------------------------------------------
# Theil-Sen slope calculation and tmap maps for water-quality indicators
# at ADM1 resolution (1990–2019).
#
# Part 1: Slope calculation (mblm Theil-Sen per adm_id)
# Part 2: tmap maps – level maps (first/last year) + slope maps
#
# Inputs:
#   results/wq_adm1_extracted.csv                      – script 21
#   results_gis/polyg_wq_adm1.gpkg                     – script 21
#   results_gis/wash_adm1_template_5arcmin.tif          – script 11
#
# Outputs:
#   results/wq_slopes_adm1_1990_2019.csv
#   results_gis/rast_slope_wq_{indicator}_1990_2019.tif
#   figures/maps/wq_{indicator}_{year}.png
#   figures/maps/wq_{indicator}_slope.png
#   figures/legends/wq_legend_{indicator}_{year}.png
#   figures/legends/wq_legend_{indicator}_slope.png
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
years_in   <- 1990:2019
first_year <- years_in[1]
last_year  <- years_in[length(years_in)]
wq_indicators <- c("organic", "overall", "pathogen", "salinity")

dir.create("figures/maps",    recursive = TRUE, showWarnings = FALSE)
dir.create("figures/legends", recursive = TRUE, showWarnings = FALSE)

gpkg_wq_path    <- "results_gis/polyg_wq_adm1.gpkg"
gpkg_simpl_path <- "results_gis/polyg_wash_adm1_1990_2023_geomonly_simpl.gpkg"

### 1. Load data ----
message("== Loading data ==")

wq_long <- read_csv("results/wq_adm1_extracted.csv", show_col_types = FALSE) |>
  filter(year %in% years_in)

r_template <- rast("results_gis/wash_adm1_template_5arcmin.tif")

message(sprintf("  ADM1 units: %d", n_distinct(wq_long$adm_id)))
message(sprintf("  Years: %d–%d", first_year, last_year))

### 2. Pivot wide for level maps ----
wq_wide <- wq_long |>
  pivot_wider(
    names_from  = c(indicator, year),
    values_from = wq_mean,
    names_sep   = "_"
  )

### 3. Slope calculation function ----

calc_slopes <- function(indicator, data) {
  message(sprintf("\n--- Slopes: %s ---", indicator))

  vals <- data |>
    filter(indicator == !!indicator, !is.na(wq_mean)) |>
    select(adm_id, year, value = wq_mean) |>
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

### 4. Compute slopes for all WQ indicators ----
message("== Computing Theil-Sen slopes ==")

slope_list <- list()
for (ind in wq_indicators) {
  slope_list[[ind]] <- calc_slopes(ind, wq_long)
}

# Save slopes as CSV
slope_df <- bind_rows(slope_list, .id = "indicator")
write_csv(slope_df, sprintf("results/wq_slopes_adm1_%d_%d.csv",
                            first_year, last_year))
message("  Saved slopes to results/")

# Slope rasters
message("== Creating slope rasters ==")
for (ind in wq_indicators) {
  slopes_ind <- slope_list[[ind]]

  tmpl <- r_template
  all_ids  <- as.numeric(na.omit(unique(values(tmpl))))
  used_ids <- as.numeric(unique(slopes_ind$adm_id))
  missing  <- all_ids[!all_ids %in% used_ids]
  if (length(missing) > 0) tmpl[tmpl %in% missing] <- NA

  r_slope <- classify(tmpl,
                      cbind(as.numeric(slopes_ind$adm_id),
                            as.numeric(slopes_ind$slope_est)),
                      others = NA)
  names(r_slope) <- paste0(ind, "_slope_", first_year, "_", last_year)

  rast_path <- sprintf("results_gis/rast_slope_wq_%s_%d_%d.tif",
                        ind, first_year, last_year)
  writeRaster(r_slope, rast_path, gdal = "COMPRESS=LZW", overwrite = TRUE)
  message(sprintf("  Saved: %s", rast_path))
}

### 5. Define tmap plotting functions ----

adm0_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
  filter(iso_a3 != "ATA")
#cols4all:: c4a_gui()
pal_level <- scico::scico(n = 20, palette = "roma", direction = -1, begin = 0, end = 0.6)
pal_slope <- cols4all::c4a("bu_br_div", reverse = F)

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

### 6. Level maps (first year and last year) ----
message("\n== Creating level maps ==")

# WQ class values: all indicators 1–5 (1=best, 5=worst)
ind_max_class <- c(organic = 5, overall = 5, pathogen = 5, salinity = 5)

ind_titles <- c(
  organic   = "Organic pollution (class, 5=worst)",
  overall   = "Overall water quality (class, 5=worst)",
  pathogen  = "Pathogen contamination (class, 5=worst)",
  salinity  = "Salinity (class, 5=worst)"
)

# Load simplified gpkg and join WQ data
gpkg_sf <- st_read(gpkg_simpl_path, quiet = TRUE) |>
  select(adm_id) |>
  left_join(wq_wide, by = "adm_id")

# Join slope columns
for (ind in wq_indicators) {
  slope_col <- paste0(ind, "_slope")
  pval_col  <- paste0(ind, "_p.value")
  sig_col   <- paste0(ind, "_sig")

  join_df <- slope_list[[ind]] |>
    transmute(
      adm_id,
      !!slope_col := slope,
      !!pval_col  := p.value,
      !!sig_col   := p.value_sig
    )
  gpkg_sf <- gpkg_sf |>
    left_join(join_df, by = "adm_id")
}

for (ind in wq_indicators) {
  message(sprintf("  Level maps: %s", ind))

  # Breaks adapted to each indicator's class range
  max_cls <- ind_max_class[[ind]]
  mybreaks <- seq(1, max_cls, length.out = 17)

  for (yr in c(first_year, last_year)) {
    col_name <- paste0(ind, "_", yr)

    p <- f_map_level(
      sf_in       = gpkg_sf,
      data_col    = col_name,
      mybreaks    = mybreaks,
      pal         = pal_level,
      title_label = ind_titles[[ind]]
    )

    save_map(p, sprintf("figures/maps/wq_%s_%d.png", ind, yr))
    save_legend(p, sprintf("figures/legends/wq_legend_%s_%d.png", ind, yr))
    message(sprintf("    Saved: wq_%s_%d.png", ind, yr))
  }
}

### 7. Slope maps ----
message("\n== Creating slope maps ==")

ind_slope_titles <- c(
  organic   = "Trend, organic pollution (class units/yr)",
  overall   = "Trend, overall water quality (class units/yr)",
  pathogen  = "Trend, pathogen contamination (class units/yr)",
  salinity  = "Trend, salinity (class units/yr)"
)

# Derive symmetric slope breaks from all indicators combined
all_slopes <- numeric(0)
for (ind in wq_indicators) {
  all_slopes <- c(all_slopes, gpkg_sf[[paste0(ind, "_slope")]])
}

min_s   <- quantile(all_slopes, 0.05, na.rm = TRUE)
max_s   <- quantile(all_slopes, 0.95, na.rm = TRUE)
abs_max <- max(abs(c(min_s, max_s)))
accur   <- round((abs_max * 2) / 20, digits = 5)

slope_breaks <- seq(
  plyr::round_any(-abs_max, accuracy = accur, f = floor),
  plyr::round_any( abs_max, accuracy = accur, f = ceiling),
  by = plyr::round_any(abs_max * 2 / 20, accuracy = 0.00001, f = floor)
)

message(sprintf("  Slope range (5th–95th pct): %.4f to %.4f",
                -abs_max, abs_max))

for (ind in wq_indicators) {
  message(sprintf("  Slope map: %s", ind))

  slope_col <- paste0(ind, "_slope")

  p <- f_map_slope(
    sf_in       = gpkg_sf,
    slope_col   = slope_col,
    mybreaks    = slope_breaks,
    pal         = pal_slope,
    title_label = ind_slope_titles[[ind]]
  )

  save_map(p, sprintf("figures/maps/wq_%s_slope.png", ind))
  save_legend(p, sprintf("figures/legends/wq_legend_%s_slope.png", ind),
              width = 150, height = 30)
  message(sprintf("    Saved: wq_%s_slope.png", ind))
}

message("\nDone.")
