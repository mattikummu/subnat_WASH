# 17_input_data_maps_wash.R
# ---------------------------------------------------------------------------
# Map input-data coverage for drinking-water pipeline.
#
# Produces ADM0 and ADM1 maps showing:
#   1) Number of years with observed data
#   2) Observed year range (max year - min year + 1)
#
# Inputs:
#   results/JMP_adm0_water.csv
#   results/JMP_adm0_water_historical.csv
#   results/wash_best_division_data.csv
#   results/wash_best_division.gpkg
#   results_gis/polyg_wash_adm0_geomonly_simpl.gpkg
#
# Outputs:
#   results/wash_input_data_coverage_adm0.csv
#   results/wash_input_data_coverage_adm1.csv
#   figures/maps/wash_inputdata_adm0_n_years.png
#   figures/maps/wash_inputdata_adm0_year_range.png
#   figures/maps/wash_inputdata_adm1_n_years.png
#   figures/maps/wash_inputdata_adm1_year_range.png
#   figures/legends/wash_inputdata_adm0_n_years_legend.png
#   figures/legends/wash_inputdata_adm0_year_range_legend.png
#   figures/legends/wash_inputdata_adm1_n_years_legend.png
#   figures/legends/wash_inputdata_adm1_year_range_legend.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tmap)
  library(cols4all)
  library(rnaturalearth)
})

### 0. Setup ----
adm0_geom_path <- "results_gis/polyg_wash_adm0_geomonly_simpl.gpkg"
adm0_data_path <- "results/JMP_adm0_water.csv"
adm0_hist_data_path <- "results/JMP_adm0_water_historical.csv"
adm1_data_path <- "results/wash_best_division_data.csv"
adm1_geom_path <- "results/wash_best_division.gpkg"
years_adm0 <- 1990:2024

if (!file.exists(adm0_geom_path)) {
  stop("Missing input: ", adm0_geom_path,
       "\nRun 05b_build_adm1_shared_assets.R first.")
}
if (!file.exists(adm0_data_path)) stop("Missing input: ", adm0_data_path)
if (!file.exists(adm0_hist_data_path)) stop("Missing input: ", adm0_hist_data_path)
if (!file.exists(adm1_data_path)) stop("Missing input: ", adm1_data_path)
if (!file.exists(adm1_geom_path)) stop("Missing input: ", adm1_geom_path)

dir.create("results", showWarnings = FALSE)
dir.create("figures/maps", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/legends", recursive = TRUE, showWarnings = FALSE)

tmap_mode("plot")

### 1. Helpers ----
coverage_summary <- function(data, group_cols, year_col = "year", indicator_cols) {
  data |>
    mutate(
      has_data = if_any(all_of(indicator_cols), ~ !is.na(.x))
    ) |>
    summarise(
      n_years = n_distinct(.data[[year_col]][has_data]),
      year_min = if_else(any(has_data), min(.data[[year_col]][has_data], na.rm = TRUE), NA_real_),
      year_max = if_else(any(has_data), max(.data[[year_col]][has_data], na.rm = TRUE), NA_real_),
      .by = all_of(group_cols)
    ) |>
    mutate(
      year_range = if_else(!is.na(year_min) & !is.na(year_max), year_max - year_min + 1, NA_real_)
    )
}

map_metric <- function(sf_obj, value_col, title_label, breaks = NULL, labels = NULL,
                       full_globe = FALSE) {
  world_border <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
    filter(iso_a3 != "ATA")

  vals <- sf_obj[[value_col]]
  vals <- vals[is.finite(vals)]
  mybreaks <- if (!is.null(breaks)) {
    breaks
  } else if (length(vals) > 1) {
    unique(pretty(vals, n = 6))
  } else {
    c(0, 1)
  }

  make_interval_labels <- function(breaks_in) {
    n_cls <- length(breaks_in) - 1
    out <- character(n_cls)
    for (i in seq_len(n_cls)) {
      lo <- ceiling(breaks_in[i])
      hi <- breaks_in[i + 1]
      if (is.finite(hi)) {
        hi <- floor(hi)
        out[i] <- if (lo == hi) as.character(lo) else paste0(lo, "-", hi)
      } else {
        out[i] <- paste0(lo, "+")
      }
    }
    out
  }

  legend_labels <- if (!is.null(labels)) {
    labels
  } else {
    make_interval_labels(mybreaks)
  }

  bbox_obj <- if (full_globe) sf::st_bbox(world_border) else NULL

  p <- tm_shape(sf_obj, bbox = bbox_obj) +
    tm_fill(
      fill = value_col,
      fill.scale = tm_scale(
        values = cols4all::c4a("lapaz", reverse = TRUE),
        breaks = mybreaks,
        labels = legend_labels,
        value.na = "white"
      ),
      fill.legend = tm_legend(
        orientation = "portrait",
        position = tm_pos_out("right", "center"),
        frame = FALSE,
        title = title_label
      )
    ) +
    tm_shape(world_border) +
    tm_borders(col = "grey35", lwd = 0.04) +
    tm_crs("+proj=robin") +
    tm_layout(main.title.position = "center", frame = FALSE)
}

save_map <- function(p_map, filename, width = 150, units = "mm", dpi = 450) {
  p_no_legend <- p_map + tm_layout(legend.show = FALSE)
  tmap_save(p_no_legend, filename = filename, width = width, units = units, dpi = dpi)
  message(sprintf("  Saved: %s", filename))
}

save_legend <- function(p_map, filename, width = 70, height = 70,
                        units = "mm", dpi = 600) {
  p_leg <- p_map +
    tm_layout(
      legend.only = TRUE,
      frame = FALSE,
      outer.margins = 0,
      inner.margins = c(0, 0, 0, 0),
      legend.position = tm_pos_in("center", "center")
    )
  tmap_save(p_leg, filename = filename, width = width, height = height,
            units = units, dpi = dpi)
  message(sprintf("  Saved: %s", filename))
}

### 2. ADM0 coverage ----
message("== ADM0 input-data coverage ==")

adm0_data_modern <- read_csv(adm0_data_path, show_col_types = FALSE) |>
  transmute(
    iso3 = toupper(as.character(iso3)),
    country_name = as.character(country_name),
    year = as.integer(year),
    Basic = as.numeric(Basic),
    Limited = as.numeric(Limited),
    Unimproved = as.numeric(Unimproved),
    Surface.water = as.numeric(Surface.water)
  )

adm0_data_hist <- read_csv(adm0_hist_data_path, show_col_types = FALSE) |>
  transmute(
    iso3 = toupper(as.character(iso3)),
    country_name = as.character(country_name),
    year = as.integer(year),
    Basic = as.numeric(Basic),
    Limited = NA_real_,
    Unimproved = NA_real_,
    Surface.water = as.numeric(Surface.water)
  )

adm0_data <- bind_rows(adm0_data_hist, adm0_data_modern) |>
  filter(year %in% years_adm0) |>
  arrange(iso3, year) |>
  distinct(iso3, year, .keep_all = TRUE)

adm0_cov <- coverage_summary(
  data = adm0_data,
  group_cols = c("iso3"),
  year_col = "year",
  indicator_cols = c("Basic", "Limited", "Unimproved", "Surface.water")
)

write_csv(adm0_cov, "results/wash_input_data_coverage_adm0.csv")
message("  Saved: results/wash_input_data_coverage_adm0.csv")

adm0_sf <- st_read(adm0_geom_path, quiet = TRUE) |>
  transmute(
    iso3 = toupper(as.character(iso3)),
    country_name = as.character(country_name)
  ) |>
  left_join(adm0_cov, by = "iso3")

### 3. ADM1 coverage ----
message("== ADM1 input-data coverage ==")

adm1_data <- read_csv(adm1_data_path, show_col_types = FALSE)

adm1_inds <- intersect(c("Basic", "Limited", "Unimproved", "Surface.water"), names(adm1_data))
adm1_groups <- intersect(c("iso3", "region_name", "boundary_source"), names(adm1_data))

adm1_cov <- coverage_summary(
  data = adm1_data,
  group_cols = adm1_groups,
  year_col = "year",
  indicator_cols = adm1_inds
)

write_csv(adm1_cov, "results/wash_input_data_coverage_adm1.csv")
message("  Saved: results/wash_input_data_coverage_adm1.csv")

adm1_sf <- st_read(adm1_geom_path, quiet = TRUE)

if (!"boundary_source" %in% names(adm1_sf) && "source" %in% names(adm1_sf)) {
  adm1_sf <- adm1_sf |> mutate(boundary_source = as.character(source))
}

join_keys <- intersect(c("iso3", "region_name", "boundary_source"), names(adm1_sf))

adm1_sf <- adm1_sf |>
  left_join(adm1_cov, by = join_keys)

### 4. Plot maps ----
message("== Plotting coverage maps ==")

p <- map_metric(
  sf_obj = adm0_sf,
  value_col = "n_years",
  title_label = "ADM0 years with observed water data"
)
save_map(p, "figures/maps/wash_inputdata_adm0_n_years.png")
save_legend(p, "figures/legends/wash_inputdata_adm0_n_years_legend.png")

p <- map_metric(
  sf_obj = adm0_sf,
  value_col = "year_range",
  title_label = "ADM0 observed-year range (years)"
)
save_map(p, "figures/maps/wash_inputdata_adm0_year_range.png")
save_legend(p, "figures/legends/wash_inputdata_adm0_year_range_legend.png")

p <- map_metric(
  sf_obj = adm1_sf,
  value_col = "n_years",
  title_label = "ADM1 years with observed water data",
  breaks = c(0.5, 1.5, 3.5, 5.5, 7.5, Inf),
  labels = c("1", "2-3", "4-5", "6-7", "8+"),
  full_globe = TRUE
)
save_map(p, "figures/maps/wash_inputdata_adm1_n_years.png")
save_legend(p, "figures/legends/wash_inputdata_adm1_n_years_legend.png")

p <- map_metric(
  sf_obj = adm1_sf,
  value_col = "year_range",
  title_label = "ADM1 observed-year range (years)",
  full_globe = TRUE
)
save_map(p, "figures/maps/wash_inputdata_adm1_year_range.png")
save_legend(p, "figures/legends/wash_inputdata_adm1_year_range_legend.png")

message("\nDone.")
