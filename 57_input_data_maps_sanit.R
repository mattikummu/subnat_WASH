# 57_input_data_maps_sanit.R
# ---------------------------------------------------------------------------
# Map input-data coverage for sanitation pipeline.
#
# Produces ADM0 and ADM1 maps showing:
#   1) Number of years with observed data
#   2) Observed year range (max year - min year + 1)
#
# Inputs:
#   results_sanit/JMP_adm0_sanitation.csv
#   results_sanit/sanit_best_division_data.csv
#   results/wash_best_division.gpkg
#   results_gis_sanit/polyg_sanit_adm0_geomonly_simpl.gpkg
#
# Outputs:
#   results_sanit/sanit_input_data_coverage_adm0.csv
#   results_sanit/sanit_input_data_coverage_adm1.csv
#   figures_sanit/maps/sanit_inputdata_adm0_n_years.png
#   figures_sanit/maps/sanit_inputdata_adm0_year_range.png
#   figures_sanit/maps/sanit_inputdata_adm1_n_years.png
#   figures_sanit/maps/sanit_inputdata_adm1_year_range.png
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
adm0_geom_path <- "results_gis_sanit/polyg_sanit_adm0_geomonly_simpl.gpkg"
adm0_data_path <- "results_sanit/JMP_adm0_sanitation.csv"
adm1_data_path <- "results_sanit/sanit_best_division_data.csv"
adm1_geom_path <- "results/wash_best_division.gpkg"

if (!file.exists(adm0_geom_path)) {
  stop("Missing input: ", adm0_geom_path,
       "\nRun 05b_build_adm1_shared_assets.R first.")
}
if (!file.exists(adm0_data_path)) stop("Missing input: ", adm0_data_path)
if (!file.exists(adm1_data_path)) stop("Missing input: ", adm1_data_path)
if (!file.exists(adm1_geom_path)) stop("Missing input: ", adm1_geom_path)

dir.create("results_sanit", showWarnings = FALSE)
dir.create("figures_sanit/maps", recursive = TRUE, showWarnings = FALSE)

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

map_metric <- function(sf_obj, value_col, title_label, file_out) {
  world_border <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
    filter(iso_a3 != "ATA")

  vals <- sf_obj[[value_col]]
  vals <- vals[is.finite(vals)]
  mybreaks <- if (length(vals) > 1) unique(pretty(vals, n = 10)) else c(0, 1)

  p <- tm_shape(sf_obj) +
    tm_fill(
      fill = value_col,
      fill.scale = tm_scale(
        values = cols4all::c4a("lapaz", reverse = TRUE),
        breaks = mybreaks,
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
    tm_borders(col = "grey35", lwd = 0.02) +
    tm_crs("+proj=robin") +
    tm_layout(main.title.position = "center", frame = FALSE)

  tmap_save(p, filename = file_out, width = 150, units = "mm", dpi = 450)
  message(sprintf("  Saved: %s", file_out))
}

### 2. ADM0 coverage ----
message("== ADM0 input-data coverage ==")

adm0_data <- read_csv(adm0_data_path, show_col_types = FALSE) |>
  transmute(
    iso3 = toupper(as.character(iso3)),
    country_name = as.character(country_name),
    year = as.integer(year),
    s_bas = as.numeric(s_bas),
    s_lim = as.numeric(s_lim),
    s_unimp = as.numeric(s_unimp),
    s_od = as.numeric(s_od)
  )

adm0_cov <- coverage_summary(
  data = adm0_data,
  group_cols = c("iso3"),
  year_col = "year",
  indicator_cols = c("s_bas", "s_lim", "s_unimp", "s_od")
)

write_csv(adm0_cov, "results_sanit/sanit_input_data_coverage_adm0.csv")
message("  Saved: results_sanit/sanit_input_data_coverage_adm0.csv")

adm0_sf <- st_read(adm0_geom_path, quiet = TRUE) |>
  transmute(
    iso3 = toupper(as.character(iso3)),
    country_name = as.character(country_name)
  ) |>
  left_join(adm0_cov, by = "iso3")

### 3. ADM1 coverage ----
message("== ADM1 input-data coverage ==")

adm1_data <- read_csv(adm1_data_path, show_col_types = FALSE)

adm1_inds <- intersect(c("s_bas", "s_lim", "s_unimp", "s_od"), names(adm1_data))
adm1_groups <- intersect(c("iso3", "region_name", "boundary_source"), names(adm1_data))

adm1_cov <- coverage_summary(
  data = adm1_data,
  group_cols = adm1_groups,
  year_col = "year",
  indicator_cols = adm1_inds
)

write_csv(adm1_cov, "results_sanit/sanit_input_data_coverage_adm1.csv")
message("  Saved: results_sanit/sanit_input_data_coverage_adm1.csv")

adm1_sf <- st_read(adm1_geom_path, quiet = TRUE)

if (!"boundary_source" %in% names(adm1_sf) && "source" %in% names(adm1_sf)) {
  adm1_sf <- adm1_sf |> mutate(boundary_source = as.character(source))
}

join_keys <- intersect(c("iso3", "region_name", "boundary_source"), names(adm1_sf))

adm1_sf <- adm1_sf |>
  left_join(adm1_cov, by = join_keys)

### 4. Plot maps ----
message("== Plotting coverage maps ==")

map_metric(
  sf_obj = adm0_sf,
  value_col = "n_years",
  title_label = "ADM0 years with observed sanitation data",
  file_out = "figures_sanit/maps/sanit_inputdata_adm0_n_years.png"
)

map_metric(
  sf_obj = adm0_sf,
  value_col = "year_range",
  title_label = "ADM0 observed-year range (years)",
  file_out = "figures_sanit/maps/sanit_inputdata_adm0_year_range.png"
)

map_metric(
  sf_obj = adm1_sf,
  value_col = "n_years",
  title_label = "ADM1 years with observed sanitation data",
  file_out = "figures_sanit/maps/sanit_inputdata_adm1_n_years.png"
)

map_metric(
  sf_obj = adm1_sf,
  value_col = "year_range",
  title_label = "ADM1 observed-year range (years)",
  file_out = "figures_sanit/maps/sanit_inputdata_adm1_year_range.png"
)

message("\nDone.")
