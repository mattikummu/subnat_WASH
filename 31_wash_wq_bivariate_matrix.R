# 31_wash_wq_bivariate_matrix.R
# ---------------------------------------------------------------------------
# Bivariate 3×3 matrix maps combining WASH access indicators with overall
# water quality at ADM1 resolution.
#
# Two sets of bivariate maps are produced:
#   Part A (ratio-based):
#     A1. Overall WQ  ×  Surface water (%)
#     A2. Overall WQ  ×  Unimproved (%)
#   Part B (population hotspots):
#     B1. Overall WQ  ×  log10(pop using surface water)
#     B2. Overall WQ  ×  log10(pop using unimproved water)
#
# Method:
#   1. Pool ALL years × ALL admin units to derive global percentile
#      cut-points (33rd and 67th percentiles), ignoring zero values.
#   2. Classify each unit × year into three ordinal groups (Low / Med / High).
#      Zeros are assigned to the "Low" group.
#   3. Build a bivariate 9-class matrix and map the last overlapping year
#      (2019) using a 3×3 bivariate palette.
#
# Inputs:
#   results/wq_adm1_extracted.csv                       – script 21
#   results/wash_adm1_ratio_adjusted.csv                – script 11
#   results_gis/polyg_wash_adm1_1990_2023_geomonly_simpl.gpkg – script 05b
#
# Outputs:
#   results/bivar_wq_surfwater_matrix_2019.csv
#   results/bivar_wq_unimproved_matrix_2019.csv
#   figures/maps/bivar_wq_surfwater_2019.png
#   figures/maps/bivar_wq_unimproved_2019.png
#   figures/legends/bivar_legend_wq_surfwater.png
#   figures/legends/bivar_legend_wq_unimproved.png
# ---------------------------------------------------------------------------

library(sf)
library(tmap)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

dir.create("figures/maps",    recursive = TRUE, showWarnings = FALSE)
dir.create("figures/legends", recursive = TRUE, showWarnings = FALSE)
dir.create("results",         showWarnings = FALSE)

### 0. Settings ----
map_year <- 2019
years_in <- 1990:2019   # overlapping period for both datasets

# User-defined cut-points (two values each).
# Set to NULL to use data-driven 33rd/67th percentile thresholds.
user_wq_cuts  <- c(2, 3)  # overall WQ class units on 1–5 scale (5=worst; NULL = percentile)
user_sw_cuts <- c(2, 8) # surface water %       (NULL = percentile)
user_un_cuts <- c(2, 8) # unimproved %           (NULL = percentile)
user_log_sw_cuts <- c(4.5,5.5)    # log10(pop surf. water) (NULL = percentile)
user_log_un_cuts <- c(5,6)    # log10(pop unimproved)  (NULL = percentile)

### 1. Load data ----
message("== Loading data ==")

# WQ long: adm_id, year, indicator, wq_mean
wq_long <- read_csv("results/wq_adm1_extracted.csv", show_col_types = FALSE) |>
  filter(year %in% years_in, indicator == "overall") |>
  select(adm_id, year, wq_overall = wq_mean)

# WASH ratio-adjusted: adm_id, year, pop, Basic, Limited, Unimproved, Surface.water
wash_long <- read_csv("results/wash_adm1_ratio_adjusted.csv",
                      show_col_types = FALSE) |>
  filter(year %in% years_in) |>
  select(adm_id, year, pop, Surface.water, Unimproved)

message(sprintf("  WQ records:   %d", nrow(wq_long)))
message(sprintf("  WASH records: %d", nrow(wash_long)))

# Merge
combined <- inner_join(wq_long, wash_long, by = c("adm_id", "year"))
message(sprintf("  Joined records: %d  (%d admin units)",
                nrow(combined), n_distinct(combined$adm_id)))

### 2. Global percentile thresholds (ignoring zeros) ----
message("\n== Computing percentile thresholds (zeros excluded) ==")

calc_cuts <- function(x) {
  x_nz <- x[!is.na(x) & x > 0]
  quantile(x_nz, probs = c(1/3, 2/3), na.rm = TRUE)
}

wq_cuts <- if (!is.null(user_wq_cuts)) user_wq_cuts else calc_cuts(combined$wq_overall)
sw_cuts <- if (!is.null(user_sw_cuts)) user_sw_cuts else calc_cuts(combined$Surface.water)
un_cuts <- if (!is.null(user_un_cuts)) user_un_cuts else calc_cuts(combined$Unimproved)

message(sprintf("  Overall WQ cuts:    %.3f / %.3f", wq_cuts[1], wq_cuts[2]))
message(sprintf("  Surface water cuts: %.2f%% / %.2f%%", sw_cuts[1], sw_cuts[2]))
message(sprintf("  Unimproved cuts:    %.2f%% / %.2f%%", un_cuts[1], un_cuts[2]))

### 3. Classification function ----
# 1 = Low, 2 = Medium, 3 = High
# Zero values are mapped to class 1 (Low)
classify3 <- function(x, cuts) {
  case_when(
    is.na(x)     ~ NA_integer_,
    x == 0       ~ 1L,
    x <= cuts[1] ~ 1L,
    x <= cuts[2] ~ 2L,
    x >  cuts[2] ~ 3L,
    .default     = NA_integer_
  )
}

combined <- combined |>
  mutate(
    wq_class = classify3(wq_overall,   wq_cuts),  # high class = worst WQ (5=worst)
    sw_class = classify3(Surface.water, sw_cuts),
    un_class = classify3(Unimproved,    un_cuts),
    pop_sw   = pop * Surface.water / 100,
    pop_un   = pop * Unimproved    / 100
  )

### 4. Build bivariate tables for map year ----
message(sprintf("\n== Building bivariate matrices for %d ==", map_year))

map_data <- combined |>
  filter(year == map_year)

# A. Overall WQ × Surface water
bivar_sw <- map_data |>
  filter(!is.na(wq_class), !is.na(sw_class)) |>
  mutate(
    bivar_class = paste0(wq_class, "_", sw_class),
    bivar_id    = (wq_class - 1L) * 3L + sw_class
  )

# B. Overall WQ × Unimproved
bivar_un <- map_data |>
  filter(!is.na(wq_class), !is.na(un_class)) |>
  mutate(
    bivar_class = paste0(wq_class, "_", un_class),
    bivar_id    = (wq_class - 1L) * 3L + un_class
  )

message(sprintf("  WQ × Surface water: %d units", nrow(bivar_sw)))
message(sprintf("  WQ × Unimproved:    %d units", nrow(bivar_un)))

write_csv(bivar_sw |> select(adm_id, wq_overall, Surface.water,
                             wq_class, sw_class, bivar_class),
          sprintf("results/bivar_wq_surfwater_matrix_%d.csv", map_year))
write_csv(bivar_un |> select(adm_id, wq_overall, Unimproved,
                             wq_class, un_class, bivar_class),
          sprintf("results/bivar_wq_unimproved_matrix_%d.csv", map_year))

### 5. Bivariate colour palette (3×3) ----
# Rows = WQ quality (1 Low/good → 3 High/bad)
# Cols = WASH deprivation (1 Low → 3 High)
# Corner logic:
#   1_1  Low WQ problem,  Low WASH deprivation  → light grey (best)
#   3_1  High WQ problem, Low WASH deprivation  → teal (WQ dominant)
#   1_3  Low WQ problem,  High WASH deprivation → pink/purple (WASH dominant)
#   3_3  High WQ problem, High WASH deprivation → dark (both bad)

# # Load the library
# library(pals)

# # Get all 9 hex codes for stevens.bluered
# stevens_hex <- stevens.bluered(n = 9)

# "#e8e8e8"
# "#e4acac"
# "#c85a5a"
# "#b0d5df"
# "#ad9ea5"
# "#985356"
# "#64acbe"
# "#627f8c"
# "#574249"

# bivar_palette_1 <- c(
#   "1_1" = "#E8E8E8", # low-low: light grey (best/neutral)
#   "1_2" = "#C8D8E8", # low pathogen, med surface water
#   "1_3" = "#A8C8E8", # low pathogen, high surface water: light blue
#   "2_1" = "#F0C8B0", # med pathogen, low surface water
#   "2_2" = "#C8B0C8", # med-med: middle purple-grey
#   "2_3" = "#9090C8", # med pathogen, high surface water: blue-purple
#   "3_1" = "#E89070", # high pathogen, low surface water: orange
#   "3_2" = "#C86888", # high pathogen, med surface water: red-purple
#   "3_3" = "#D7191C" # high-high: BRIGHT RED (stevens.bluered red)
# )

# bivar_palette_2 <- c(
#   "1_1" = "#F0F0F0", # Lighter grey for better contrast
#   "1_2" = "#B8D4E8",
#   "1_3" = "#6BAED6", # More saturated blue
#   "2_1" = "#FCAE91",
#   "2_2" = "#C994C7",
#   "2_3" = "#7A6FAC",
#   "3_1" = "#FB6A4A",
#   "3_2" = "#CB4C7C",
#   "3_3" = "#CB181D" # High saturation red
# )

### 5. Bivariate colour palette (3×3) ----
bivar_palette <- c(
  "1_1" = "#FCFCFC", # Lighter neutral (near white)
  "1_2" = "#DEEBF7", # Lighter blue tint
  "1_3" = "#9ECAE1", # Lighter blue (was #6BAED6)

  "2_1" = "#FEE8C8", # Lighter peach/orange tint
  "2_2" = "#D9D9D9", # Lighter neutral grey/purple
  "2_3" = "#756BB1", # [LOCKED] med WQ, high WASH

  "3_1" = "#FCBBA1", # Lighter red/salmon tint
  "3_2" = "#CA364C", # [LOCKED] high WQ, med WASH
  "3_3" = "#A50F15" # [LOCKED] high WQ, high WASH
)

class_labels_sw <- c(
  "1_1" = "Low WQ / Low surf. water",
  "1_2" = "Low WQ / Med surf. water",
  "1_3" = "Low WQ / High surf. water",
  "2_1" = "Moderate WQ / Low surf. water",
  "2_2" = "Moderate WQ / Med surf. water",
  "2_3" = "Moderate WQ / High surf. water",
  "3_1" = "High WQ / Low surf. water",
  "3_2" = "High WQ / Med surf. water",
  "3_3" = "High WQ / High surf. water"
)

class_labels_un <- c(
  "1_1" = "Low WQ / Low unimproved",
  "1_2" = "Low WQ / Med unimproved",
  "1_3" = "Low WQ / High unimproved",
  "2_1" = "Moderate WQ / Low unimproved",
  "2_2" = "Moderate WQ / Med unimproved",
  "2_3" = "Moderate WQ / High unimproved",
  "3_1" = "High WQ / Low unimproved",
  "3_2" = "High WQ / Med unimproved",
  "3_3" = "High WQ / High unimproved"
)

### 6. Load geometry and join ----
message("\n== Joining to geometry ==")

gpkg_sf <- st_read("results_gis/polyg_wash_adm1_1990_2023_geomonly_simpl.gpkg",
                   quiet = TRUE) |>
  select(adm_id)

adm0_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
  filter(iso_a3 != "ATA")

# Join bivariate classifications
map_sf_sw <- gpkg_sf |>
  left_join(
    bivar_sw |> select(adm_id, bivar_class),
    by = "adm_id"
  ) |>
  mutate(bivar_class = factor(bivar_class, levels = names(bivar_palette)))

map_sf_un <- gpkg_sf |>
  left_join(
    bivar_un |> select(adm_id, bivar_class),
    by = "adm_id"
  ) |>
  mutate(bivar_class = factor(bivar_class, levels = names(bivar_palette)))

### 7. Map function ----
make_bivar_map <- function(sf_in, palette, adm0) {
  tm_shape(sf_in) +
    tm_fill(
      fill       = "bivar_class",
      fill.scale = tm_scale_categorical(
        values   = palette,
        value.na = "#D0D0D0"
      ),
      fill.legend = tm_legend(show = FALSE)
    ) +
    tm_shape(adm0) +
    tm_borders(col = "grey35", lwd = 0.025) +
    tm_crs("+proj=robin") +
    tm_layout(frame = FALSE, legend.show = FALSE)
}

### 8. Maps ----
message("\n== Rendering bivariate maps ==")

# A. WQ × Surface water
p_sw <- make_bivar_map(map_sf_sw, bivar_palette, adm0_sf)
sw_file <- sprintf("figures/maps/bivar_wq_surfwater_%d.png", map_year)
tmap_save(p_sw, filename = sw_file, width = 80, units = "mm", dpi = 450)
message("  Saved: ", sw_file)

# B. WQ × Unimproved
p_un <- make_bivar_map(map_sf_un, bivar_palette, adm0_sf)
un_file <- sprintf("figures/maps/bivar_wq_unimproved_%d.png", map_year)
tmap_save(p_un, filename = un_file, width = 80, units = "mm", dpi = 450)
message("  Saved: ", un_file)

### 9. Bivariate legends (ggplot2 3×3 tile grid) ----
message("\n== Rendering legends ==")

make_bivar_legend <- function(palette, x_label, y_label, filename) {
  legend_df <- tibble(
    bivar_class = names(palette),
    fill_col    = unname(palette),
    row_grp     = as.integer(substr(bivar_class, 1, 1)),
    col_grp     = as.integer(substr(bivar_class, 3, 3))
  )

  row_labels <- c("1" = "Low", "2" = "Moderate", "3" = "High")
  col_labels <- c("1" = "Low", "2" = "Med", "3" = "High")

  p <- ggplot(legend_df,
              aes(x = factor(col_grp), y = factor(row_grp),
                  fill = bivar_class)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_x_discrete(name = x_label, labels = col_labels) +
    scale_y_discrete(name = y_label, labels = row_labels) +
    coord_fixed() +
    theme_minimal(base_size = 9) +
    theme(
      axis.title.x    = element_text(size = 8, hjust = 0.5),
      axis.title.y    = element_text(size = 8, hjust = 0.5),
      axis.text       = element_text(size = 7),
      panel.grid      = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(filename, p, width = 60, height = 60, units = "mm", dpi = 600)
  message("  Saved: ", filename)
}

make_bivar_legend(
  bivar_palette,
  x_label  = "Surface water use ->",
  y_label  = "Water quality risk ->",
  filename = "figures/legends/bivar_legend_wq_surfwater.png"
)

make_bivar_legend(
  bivar_palette,
  x_label  = "Unimproved water use ->",
  y_label  = "Water quality risk ->",
  filename = "figures/legends/bivar_legend_wq_unimproved.png"
)

### 10. Summary statistics ----
message(sprintf("\n=== BIVARIATE CLASS SUMMARY (%d) ===", map_year))

summarise_bivar <- function(tbl, labels) {
  tbl |>
    count(bivar_class) |>
    mutate(
      label = labels[bivar_class],
      pct   = round(n / sum(n) * 100, 1)
    ) |>
    arrange(bivar_class)
}

message("\n--- WQ × Surface water ---")
summary_sw <- summarise_bivar(bivar_sw, class_labels_sw)
summary_sw

message("\n--- WQ × Unimproved ---")
summary_un <- summarise_bivar(bivar_un, class_labels_un)
summary_un

### ========================================================================
### PART B: Population-weighted bivariate maps
### WQ × log10(pop using surface water)
### WQ × log10(pop using unimproved water)
### ========================================================================

message("\n== PART B: Population hotspot maps ==")

### 11. Compute absolute population under each WASH ladder ----
combined <- combined |>
  mutate(
    log_pop_sw = if_else(pop_sw > 0, log10(pop_sw), NA_real_),
    log_pop_un = if_else(pop_un > 0, log10(pop_un), NA_real_)
  )

### 12. Percentile thresholds on log10(pop), zeros excluded ----
message("\n== Computing log10(pop) thresholds (zeros excluded) ==")

log_sw_cuts <- if (!is.null(user_log_sw_cuts)) user_log_sw_cuts else calc_cuts(combined$log_pop_sw)
log_un_cuts <- if (!is.null(user_log_un_cuts)) user_log_un_cuts else calc_cuts(combined$log_pop_un)

message(sprintf("  log10(pop surf. water) cuts: %.3f / %.3f  (= %.0f / %.0f people)",
                log_sw_cuts[1], log_sw_cuts[2],
                10^log_sw_cuts[1], 10^log_sw_cuts[2]))
message(sprintf("  log10(pop unimproved)  cuts: %.3f / %.3f  (= %.0f / %.0f people)",
                log_un_cuts[1], log_un_cuts[2],
                10^log_un_cuts[1], 10^log_un_cuts[2]))

### 13. Classify and build bivariate tables ----
combined <- combined |>
  mutate(
    pop_sw_class = classify3(log_pop_sw, log_sw_cuts),
    pop_un_class = classify3(log_pop_un, log_un_cuts)
  )

# Zero population → class 1 (Low)
combined <- combined |>
  mutate(
    pop_sw_class = if_else(is.na(pop_sw_class) & !is.na(pop_sw) & pop_sw == 0,
                           1L, pop_sw_class),
    pop_un_class = if_else(is.na(pop_un_class) & !is.na(pop_un) & pop_un == 0,
                           1L, pop_un_class)
  )

message(sprintf("\n== Building population bivariate matrices for %d ==", map_year))

map_data_pop <- combined |>
  filter(year == map_year)

# C. WQ × log10(pop surface water)
bivar_pop_sw <- map_data_pop |>
  filter(!is.na(wq_class), !is.na(pop_sw_class)) |>
  mutate(
    bivar_class = paste0(wq_class, "_", pop_sw_class),
    bivar_id    = (wq_class - 1L) * 3L + pop_sw_class
  )

# D. WQ × log10(pop unimproved)
bivar_pop_un <- map_data_pop |>
  filter(!is.na(wq_class), !is.na(pop_un_class)) |>
  mutate(
    bivar_class = paste0(wq_class, "_", pop_un_class),
    bivar_id    = (wq_class - 1L) * 3L + pop_un_class
  )

message(sprintf("  WQ × pop surf. water: %d units", nrow(bivar_pop_sw)))
message(sprintf("  WQ × pop unimproved:  %d units", nrow(bivar_pop_un)))

write_csv(bivar_pop_sw |> select(adm_id, wq_overall, pop_sw, log_pop_sw,
                                 wq_class, pop_sw_class, bivar_class),
          sprintf("results/bivar_wq_pop_surfwater_matrix_%d.csv", map_year))
write_csv(bivar_pop_un |> select(adm_id, wq_overall, pop_un, log_pop_un,
                                 wq_class, pop_un_class, bivar_class),
          sprintf("results/bivar_wq_pop_unimproved_matrix_%d.csv", map_year))

### 14. Population bivariate labels ----
class_labels_pop_sw <- c(
  "1_1" = "Low WQ / Few pop surf. water",
  "1_2" = "Low WQ / Med pop surf. water",
  "1_3" = "Low WQ / Many pop surf. water",
  "2_1" = "Moderate WQ / Few pop surf. water",
  "2_2" = "Moderate WQ / Med pop surf. water",
  "2_3" = "Moderate WQ / Many pop surf. water",
  "3_1" = "High WQ / Few pop surf. water",
  "3_2" = "High WQ / Med pop surf. water",
  "3_3" = "High WQ / Many pop surf. water"
)

class_labels_pop_un <- c(
  "1_1" = "Low WQ / Few pop unimproved",
  "1_2" = "Low WQ / Med pop unimproved",
  "1_3" = "Low WQ / Many pop unimproved",
  "2_1" = "Moderate WQ / Few pop unimproved",
  "2_2" = "Moderate WQ / Med pop unimproved",
  "2_3" = "Moderate WQ / Many pop unimproved",
  "3_1" = "High WQ / Few pop unimproved",
  "3_2" = "High WQ / Med pop unimproved",
  "3_3" = "High WQ / Many pop unimproved"
)

### 15. Join to geometry and render population bivariate maps ----
message("\n== Rendering population bivariate maps ==")

map_sf_pop_sw <- gpkg_sf |>
  left_join(bivar_pop_sw |> select(adm_id, bivar_class), by = "adm_id") |>
  mutate(bivar_class = factor(bivar_class, levels = names(bivar_palette)))

map_sf_pop_un <- gpkg_sf |>
  left_join(bivar_pop_un |> select(adm_id, bivar_class), by = "adm_id") |>
  mutate(bivar_class = factor(bivar_class, levels = names(bivar_palette)))

# C. WQ × pop surface water
p_pop_sw <- make_bivar_map(map_sf_pop_sw, bivar_palette, adm0_sf)
pop_sw_file <- sprintf("figures/maps/bivar_wq_pop_surfwater_%d.png", map_year)
tmap_save(p_pop_sw, filename = pop_sw_file, width = 80, units = "mm", dpi = 450)
message("  Saved: ", pop_sw_file)

# D. WQ × pop unimproved
p_pop_un <- make_bivar_map(map_sf_pop_un, bivar_palette, adm0_sf)
pop_un_file <- sprintf("figures/maps/bivar_wq_pop_unimproved_%d.png", map_year)
tmap_save(p_pop_un, filename = pop_un_file, width = 80, units = "mm", dpi = 450)
message("  Saved: ", pop_un_file)

### 16. Population bivariate legends ----
make_bivar_legend(
  bivar_palette,
  x_label  = "Pop using surface water (log10) ->",
  y_label  = "Water quality risk ->",
  filename = "figures/legends/bivar_legend_wq_pop_surfwater.png"
)

make_bivar_legend(
  bivar_palette,
  x_label  = "Pop using unimproved (log10) ->",
  y_label  = "Water quality risk ->",
  filename = "figures/legends/bivar_legend_wq_pop_unimproved.png"
)

### 17. Population bivariate summary ----
message(sprintf("\n=== POPULATION BIVARIATE CLASS SUMMARY (%d) ===", map_year))

message("\n--- WQ × pop surface water ---")
summary_pop_sw <- summarise_bivar(bivar_pop_sw, class_labels_pop_sw)
summary_pop_sw

message("\n--- WQ × pop unimproved ---")
summary_pop_un <- summarise_bivar(bivar_pop_un, class_labels_pop_un)
summary_pop_un

### 18. Population per bivariate class (millions) ----
message("\n== Computing population per bivariate class ==")

pop_by_class <- function(tbl, matrix_type, labels, wq_col = "wq_class",
                         wash_col, pop_col = "pop") {
  tbl |>
    group_by(bivar_class) |>
    summarise(
      wq_class   = first(.data[[wq_col]]),
      wash_class = first(.data[[wash_col]]),
      n_units    = n(),
      pop_M      = round(sum(.data[[pop_col]], na.rm = TRUE) / 1e6, 2),
      .groups    = "drop"
    ) |>
    mutate(
      matrix_type = matrix_type,
      label       = labels[bivar_class]
    ) |>
    select(matrix_type, wq_class, wash_class, bivar_class, label, n_units, pop_M) |>
    arrange(bivar_class)
}

# Ensure all 9 classes appear even if some are empty
all_classes <- expand.grid(wq_class = 1:3, wash_class = 1:3) |>
  mutate(bivar_class = paste0(wq_class, "_", wash_class))

fill_missing <- function(df) {
  all_classes |>
    left_join(df, by = c("bivar_class", "wq_class", "wash_class")) |>
    mutate(
      n_units = replace(n_units, is.na(n_units), 0L),
      pop_M   = replace(pop_M,   is.na(pop_M),   0)
    )
}

pop_all <- bind_rows(
  pop_by_class(bivar_sw,     "wq_surfwater",     class_labels_sw,     wash_col = "sw_class",     pop_col = "pop_sw") |> fill_missing(),
  pop_by_class(bivar_un,     "wq_unimproved",    class_labels_un,     wash_col = "un_class",     pop_col = "pop_un") |> fill_missing(),
  pop_by_class(bivar_pop_sw, "wq_pop_surfwater", class_labels_pop_sw, wash_col = "pop_sw_class", pop_col = "pop_sw") |> fill_missing(),
  pop_by_class(bivar_pop_un, "wq_pop_unimproved",class_labels_pop_un, wash_col = "pop_un_class", pop_col = "pop_un") |> fill_missing()
)

pop_csv <- sprintf("results/bivar_wq_population_by_class_%d.csv", map_year)
write_csv(pop_all, pop_csv)
message("  Saved: ", pop_csv)
message(pop_all |> group_by(matrix_type) |> summarise(total_pop_M = sum(pop_M)))

### 19. Population-annotated legends (separate by result set) ----
message("\n== Rendering population-annotated legends (separate) ==")

make_bivar_legend_pop_single <- function(palette, x_label, y_label,
                                         pop_df, filename) {
  legend_df <- tibble(
    bivar_class = names(palette),
    fill_col    = unname(palette),
    row_grp     = as.integer(substr(bivar_class, 1, 1)),
    col_grp     = as.integer(substr(bivar_class, 3, 3))
  ) |>
    left_join(pop_df |> select(bivar_class, pop_M), by = "bivar_class") |>
    mutate(
      pop_M      = replace_na(pop_M, 0),
      pop_label  = format(round(pop_M, 1), nsmall = 1),
      brightness = (col2rgb(fill_col)[1,] * 299 +
                    col2rgb(fill_col)[2,] * 587 +
                    col2rgb(fill_col)[3,] * 114) / 1000,
      text_col   = if_else(brightness > 140, "grey20", "white")
    )

  row_labels <- c("1" = "Low", "2" = "Moderate", "3" = "High")
  col_labels <- c("1" = "Low", "2" = "Med", "3" = "High")

  p <- ggplot(legend_df,
              aes(x = factor(col_grp), y = factor(row_grp),
                  fill = bivar_class)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = pop_label, colour = text_col),
              size = 3.2, lineheight = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_colour_identity() +
    scale_x_discrete(name = x_label, labels = col_labels) +
    scale_y_discrete(name = y_label, labels = row_labels) +
    coord_fixed() +
    theme_minimal(base_size = 9) +
    theme(
      axis.title.x    = element_text(size = 8, hjust = 0.5),
      axis.title.y    = element_text(size = 8, hjust = 0.5),
      axis.text       = element_text(size = 7),
      panel.grid      = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(filename, p, width = 60, height = 60, units = "mm", dpi = 600)
  message("  Saved: ", filename)
}

make_bivar_legend_pop_single(
  bivar_palette,
  x_label  = "Surface water use ->",
  y_label  = "Water quality risk ->",
  pop_df   = pop_all |> filter(matrix_type == "wq_surfwater"),
  filename = "figures/legends/bivar_legend_wq_surfwater_ratio_pop.png"
)

make_bivar_legend_pop_single(
  bivar_palette,
  x_label  = "Unimproved water use ->",
  y_label  = "Water quality risk ->",
  pop_df   = pop_all |> filter(matrix_type == "wq_unimproved"),
  filename = "figures/legends/bivar_legend_wq_unimproved_ratio_pop.png"
)

make_bivar_legend_pop_single(
  bivar_palette,
  x_label  = "Pop using surface water (log10) ->",
  y_label  = "Water quality risk ->",
  pop_df   = pop_all |> filter(matrix_type == "wq_pop_surfwater"),
  filename = "figures/legends/bivar_legend_wq_pop_surfwater_pop.png"
)

make_bivar_legend_pop_single(
  bivar_palette,
  x_label  = "Pop using unimproved (log10) ->",
  y_label  = "Water quality risk ->",
  pop_df   = pop_all |> filter(matrix_type == "wq_pop_unimproved"),
  filename = "figures/legends/bivar_legend_wq_pop_unimproved_pop.png"
)

message("\nScript 31 completed successfully.")
