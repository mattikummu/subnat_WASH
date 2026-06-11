# 31b_wash_pathogen_bivariate_matrix.R
# ---------------------------------------------------------------------------
# Bivariate 3×3 matrix maps combining WASH access indicators with selected
# water quality dimensions at ADM1 resolution.
#
# Indicators processed:
#   - pathogen
#   - organic
#   - salinity
#
# For each indicator, two sets of bivariate maps are produced:
#   Part A (ratio-based):
#     A1. WQ indicator × Surface water (%)
#     A2. WQ indicator × Unimproved (%)
#   Part B (population hotspots):
#     B1. WQ indicator × log10(pop using surface water)
#     B2. WQ indicator × log10(pop using unimproved water)
# ---------------------------------------------------------------------------

library(sf)
library(tmap)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

dir.create("figures/maps", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/legends", recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

### 0. Settings ----
map_year <- 2019
years_in <- 1990:2019

# User-defined cut-points (two values each).
# Set to NULL to use data-driven 33rd/67th percentile thresholds.
user_wq_cuts <- list(
  pathogen = c(2, 3),
  organic = c(2, 3),
  salinity = c(2, 3)
)
user_sw_cuts <- c(2, 8)
user_un_cuts <- c(2, 8)
user_log_sw_cuts <- c(4.5, 5.5)
user_log_un_cuts <- c(5, 6)

indicator_cfg <- tibble::tribble(
  ~indicator, ~title, ~ylabel,
  "pathogen", "Pathogen", "Pathogen risk ->",
  "organic", "Organic", "Organic risk ->",
  "salinity", "Salinity", "Salinity risk ->"
)

### 1. Shared data and helpers ----
message("== Loading shared data ==")

wq_all <- read_csv("results/wq_adm1_extracted.csv", show_col_types = FALSE) |>
  filter(year %in% years_in)

wash_long <- read_csv("results/wash_adm1_ratio_adjusted.csv",
                      show_col_types = FALSE) |>
  filter(year %in% years_in) |>
  select(adm_id, year, pop, Surface.water, Unimproved)

message(sprintf("  WQ records:   %d", nrow(wq_all)))
message(sprintf("  WASH records: %d", nrow(wash_long)))

gpkg_sf <- st_read("results_gis/polyg_wash_adm1_1990_2023_geomonly_simpl.gpkg",
                   quiet = TRUE) |>
  select(adm_id)

adm0_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
  filter(iso_a3 != "ATA")

calc_cuts <- function(x) {
  x_nz <- x[!is.na(x) & x > 0]
  quantile(x_nz, probs = c(1/3, 2/3), na.rm = TRUE)
}

classify3 <- function(x, cuts) {
  case_when(
    is.na(x) ~ NA_integer_,
    x == 0 ~ 1L,
    x <= cuts[1] ~ 1L,
    x <= cuts[2] ~ 2L,
    x > cuts[2] ~ 3L,
    .default = NA_integer_
  )
}

make_bivar_map <- function(sf_in, palette, adm0) {
  tm_shape(sf_in) +
    tm_fill(
      fill = "bivar_class",
      fill.scale = tm_scale_categorical(
        values = palette,
        value.na = "#D0D0D0"
      ),
      fill.legend = tm_legend(show = FALSE)
    ) +
    tm_shape(adm0) +
    tm_borders(col = "grey35", lwd = 0.025) +
    tm_crs("+proj=robin") +
    tm_layout(frame = FALSE, legend.show = FALSE)
}

make_bivar_legend <- function(palette, x_label, y_label, filename) {
  legend_df <- tibble(
    bivar_class = names(palette),
    fill_col = unname(palette),
    row_grp = as.integer(substr(bivar_class, 1, 1)),
    col_grp = as.integer(substr(bivar_class, 3, 3))
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
      axis.title.x = element_text(size = 8, hjust = 0.5),
      axis.title.y = element_text(size = 8, hjust = 0.5),
      axis.text = element_text(size = 7),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(filename, p, width = 60, height = 60, units = "mm", dpi = 600)
  message("  Saved: ", filename)
}

make_bivar_legend_pop_single <- function(palette, x_label, y_label,
                                         pop_df, filename) {
  legend_df <- tibble(
    bivar_class = names(palette),
    fill_col = unname(palette),
    row_grp = as.integer(substr(bivar_class, 1, 1)),
    col_grp = as.integer(substr(bivar_class, 3, 3))
  ) |>
    left_join(pop_df |> select(bivar_class, pop_M), by = "bivar_class") |>
    mutate(
      pop_M = replace_na(pop_M, 0),
      pop_label = format(round(pop_M, 1), nsmall = 1),
      brightness = (col2rgb(fill_col)[1, ] * 299 +
                    col2rgb(fill_col)[2, ] * 587 +
                    col2rgb(fill_col)[3, ] * 114) / 1000,
      text_col = if_else(brightness > 140, "grey20", "white")
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
      axis.title.x = element_text(size = 8, hjust = 0.5),
      axis.title.y = element_text(size = 8, hjust = 0.5),
      axis.text = element_text(size = 7),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(filename, p, width = 60, height = 60, units = "mm", dpi = 600)
  message("  Saved: ", filename)
}

summarise_bivar <- function(tbl, labels) {
  tbl |>
    count(bivar_class) |>
    mutate(
      label = labels[bivar_class],
      pct = round(n / sum(n) * 100, 1)
    ) |>
    arrange(bivar_class)
}

pop_by_class <- function(tbl, matrix_type, labels, wq_col = "wq_class",
                         wash_col, pop_col = "pop") {
  tbl |>
    group_by(bivar_class) |>
    summarise(
      wq_class = first(.data[[wq_col]]),
      wash_class = first(.data[[wash_col]]),
      n_units = n(),
      pop_M = round(sum(.data[[pop_col]], na.rm = TRUE) / 1e6, 2),
      .groups = "drop"
    ) |>
    mutate(
      matrix_type = matrix_type,
      label = labels[bivar_class]
    ) |>
    select(matrix_type, wq_class, wash_class, bivar_class, label, n_units, pop_M) |>
    arrange(bivar_class)
}

all_classes <- expand.grid(wq_class = 1:3, wash_class = 1:3) |>
  mutate(bivar_class = paste0(wq_class, "_", wash_class))

fill_missing <- function(df) {
  all_classes |>
    left_join(df, by = c("bivar_class", "wq_class", "wash_class")) |>
    mutate(
      n_units = replace(n_units, is.na(n_units), 0L),
      pop_M = replace(pop_M, is.na(pop_M), 0)
    )
}

build_labels <- function(wq_name, wash_name) {
  lev_wq   <- c("Low", "Moderate", "High")
  lev_wash <- c("Low", "Med", "High")
  out <- character(9)
  nm <- character(9)
  idx <- 1
  for (wq_i in 1:3) {
    for (wash_i in 1:3) {
      nm[idx] <- paste0(wq_i, "_", wash_i)
      out[idx] <- paste0(lev_wq[wq_i], " ", wq_name, " / ", lev_wash[wash_i], " ", wash_name)
      idx <- idx + 1
    }
  }
  stats::setNames(out, nm)
}

bivar_palette <- c(
  "1_1" = "#FCFCFC",
  "1_2" = "#DEEBF7",
  "1_3" = "#9ECAE1",
  "2_1" = "#FEE8C8",
  "2_2" = "#D9D9D9",
  "2_3" = "#756BB1",
  "3_1" = "#FCBBA1",
  "3_2" = "#CA364C",
  "3_3" = "#A50F15"
)

### 2. Run workflow for each WQ indicator ----
for (i in seq_len(nrow(indicator_cfg))) {
  ind <- indicator_cfg$indicator[i]
  ind_title <- indicator_cfg$title[i]
  y_label <- indicator_cfg$ylabel[i]

  message("\n========================================================")
  message("== Indicator: ", ind)
  message("========================================================")

  wq_long <- wq_all |>
    filter(indicator == ind) |>
    select(adm_id, year, wq_value = wq_mean)

  combined <- inner_join(wq_long, wash_long, by = c("adm_id", "year"))
  message(sprintf("  WQ records:     %d", nrow(wq_long)))
  message(sprintf("  Joined records: %d (%d admin units)",
                  nrow(combined), n_distinct(combined$adm_id)))

  this_wq_cuts <- user_wq_cuts[[ind]]
  wq_cuts <- if (!is.null(this_wq_cuts)) this_wq_cuts else calc_cuts(combined$wq_value)
  sw_cuts <- if (!is.null(user_sw_cuts)) user_sw_cuts else calc_cuts(combined$Surface.water)
  un_cuts <- if (!is.null(user_un_cuts)) user_un_cuts else calc_cuts(combined$Unimproved)

  message(sprintf("  %s WQ cuts:      %.3f / %.3f", ind_title, wq_cuts[1], wq_cuts[2]))
  message(sprintf("  Surface cuts:    %.2f%% / %.2f%%", sw_cuts[1], sw_cuts[2]))
  message(sprintf("  Unimproved cuts: %.2f%% / %.2f%%", un_cuts[1], un_cuts[2]))

  combined <- combined |>
    mutate(
      wq_class = classify3(wq_value, wq_cuts),
      sw_class = classify3(Surface.water, sw_cuts),
      un_class = classify3(Unimproved, un_cuts),
      pop_sw = pop * Surface.water / 100,
      pop_un = pop * Unimproved / 100,
      log_pop_sw = if_else(pop_sw > 0, log10(pop_sw), NA_real_),
      log_pop_un = if_else(pop_un > 0, log10(pop_un), NA_real_)
    )

  map_data <- combined |>
    filter(year == map_year)

  bivar_sw <- map_data |>
    filter(!is.na(wq_class), !is.na(sw_class)) |>
    mutate(
      bivar_class = paste0(wq_class, "_", sw_class),
      bivar_id = (wq_class - 1L) * 3L + sw_class
    )

  bivar_un <- map_data |>
    filter(!is.na(wq_class), !is.na(un_class)) |>
    mutate(
      bivar_class = paste0(wq_class, "_", un_class),
      bivar_id = (wq_class - 1L) * 3L + un_class
    )

  write_csv(
    bivar_sw |>
      select(adm_id, wq_value, Surface.water, wq_class, sw_class, bivar_class),
    sprintf("results/bivar_%s_surfwater_matrix_%d.csv", ind, map_year)
  )
  write_csv(
    bivar_un |>
      select(adm_id, wq_value, Unimproved, wq_class, un_class, bivar_class),
    sprintf("results/bivar_%s_unimproved_matrix_%d.csv", ind, map_year)
  )

  map_sf_sw <- gpkg_sf |>
    left_join(bivar_sw |> select(adm_id, bivar_class), by = "adm_id") |>
    mutate(bivar_class = factor(bivar_class, levels = names(bivar_palette)))

  map_sf_un <- gpkg_sf |>
    left_join(bivar_un |> select(adm_id, bivar_class), by = "adm_id") |>
    mutate(bivar_class = factor(bivar_class, levels = names(bivar_palette)))

  p_sw <- make_bivar_map(map_sf_sw, bivar_palette, adm0_sf)
  p_un <- make_bivar_map(map_sf_un, bivar_palette, adm0_sf)

  tmap_save(p_sw,
            filename = sprintf("figures/maps/bivar_%s_surfwater_%d.png", ind, map_year),
            width = 80, units = "mm", dpi = 450)
  tmap_save(p_un,
            filename = sprintf("figures/maps/bivar_%s_unimproved_%d.png", ind, map_year),
            width = 80, units = "mm", dpi = 450)

  make_bivar_legend(
    bivar_palette,
    x_label = "Surface water use ->",
    y_label = y_label,
    filename = sprintf("figures/legends/bivar_legend_%s_surfwater.png", ind)
  )

  make_bivar_legend(
    bivar_palette,
    x_label = "Unimproved water use ->",
    y_label = y_label,
    filename = sprintf("figures/legends/bivar_legend_%s_unimproved.png", ind)
  )

  log_sw_cuts <- if (!is.null(user_log_sw_cuts)) user_log_sw_cuts else calc_cuts(combined$log_pop_sw)
  log_un_cuts <- if (!is.null(user_log_un_cuts)) user_log_un_cuts else calc_cuts(combined$log_pop_un)

  combined <- combined |>
    mutate(
      pop_sw_class = classify3(log_pop_sw, log_sw_cuts),
      pop_un_class = classify3(log_pop_un, log_un_cuts),
      pop_sw_class = if_else(is.na(pop_sw_class) & !is.na(pop_sw) & pop_sw == 0,
                             1L, pop_sw_class),
      pop_un_class = if_else(is.na(pop_un_class) & !is.na(pop_un) & pop_un == 0,
                             1L, pop_un_class)
    )

  map_data_pop <- combined |>
    filter(year == map_year)

  bivar_pop_sw <- map_data_pop |>
    filter(!is.na(wq_class), !is.na(pop_sw_class)) |>
    mutate(
      bivar_class = paste0(wq_class, "_", pop_sw_class),
      bivar_id = (wq_class - 1L) * 3L + pop_sw_class
    )

  bivar_pop_un <- map_data_pop |>
    filter(!is.na(wq_class), !is.na(pop_un_class)) |>
    mutate(
      bivar_class = paste0(wq_class, "_", pop_un_class),
      bivar_id = (wq_class - 1L) * 3L + pop_un_class
    )

  write_csv(
    bivar_pop_sw |>
      select(adm_id, wq_value, pop_sw, log_pop_sw, wq_class, pop_sw_class, bivar_class),
    sprintf("results/bivar_%s_pop_surfwater_matrix_%d.csv", ind, map_year)
  )
  write_csv(
    bivar_pop_un |>
      select(adm_id, wq_value, pop_un, log_pop_un, wq_class, pop_un_class, bivar_class),
    sprintf("results/bivar_%s_pop_unimproved_matrix_%d.csv", ind, map_year)
  )

  map_sf_pop_sw <- gpkg_sf |>
    left_join(bivar_pop_sw |> select(adm_id, bivar_class), by = "adm_id") |>
    mutate(bivar_class = factor(bivar_class, levels = names(bivar_palette)))

  map_sf_pop_un <- gpkg_sf |>
    left_join(bivar_pop_un |> select(adm_id, bivar_class), by = "adm_id") |>
    mutate(bivar_class = factor(bivar_class, levels = names(bivar_palette)))

  p_pop_sw <- make_bivar_map(map_sf_pop_sw, bivar_palette, adm0_sf)
  p_pop_un <- make_bivar_map(map_sf_pop_un, bivar_palette, adm0_sf)

  tmap_save(p_pop_sw,
            filename = sprintf("figures/maps/bivar_%s_pop_surfwater_%d.png", ind, map_year),
            width = 80, units = "mm", dpi = 450)
  tmap_save(p_pop_un,
            filename = sprintf("figures/maps/bivar_%s_pop_unimproved_%d.png", ind, map_year),
            width = 80, units = "mm", dpi = 450)

  make_bivar_legend(
    bivar_palette,
    x_label = "Pop using surface water (log10) ->",
    y_label = y_label,
    filename = sprintf("figures/legends/bivar_legend_%s_pop_surfwater.png", ind)
  )

  make_bivar_legend(
    bivar_palette,
    x_label = "Pop using unimproved (log10) ->",
    y_label = y_label,
    filename = sprintf("figures/legends/bivar_legend_%s_pop_unimproved.png", ind)
  )

  class_labels_sw <- build_labels(tolower(ind_title), "surf. water")
  class_labels_un <- build_labels(tolower(ind_title), "unimproved")
  class_labels_pop_sw <- build_labels(tolower(ind_title), "pop surf. water")
  class_labels_pop_un <- build_labels(tolower(ind_title), "pop unimproved")

  summary_sw <- summarise_bivar(bivar_sw, class_labels_sw)
  summary_un <- summarise_bivar(bivar_un, class_labels_un)
  summary_pop_sw <- summarise_bivar(bivar_pop_sw, class_labels_pop_sw)
  summary_pop_un <- summarise_bivar(bivar_pop_un, class_labels_pop_un)

  message("\n--- Ratio summary: Surface water ---")
  summary_sw
  message("\n--- Ratio summary: Unimproved ---")
  summary_un
  message("\n--- Pop summary: Surface water ---")
  summary_pop_sw
  message("\n--- Pop summary: Unimproved ---")
  summary_pop_un

  pop_all <- bind_rows(
    pop_by_class(bivar_sw, paste0(ind, "_surfwater"), class_labels_sw,
                 wash_col = "sw_class", pop_col = "pop_sw") |> fill_missing(),
    pop_by_class(bivar_un, paste0(ind, "_unimproved"), class_labels_un,
                 wash_col = "un_class", pop_col = "pop_un") |> fill_missing(),
    pop_by_class(bivar_pop_sw, paste0(ind, "_pop_surfwater"), class_labels_pop_sw,
                 wash_col = "pop_sw_class", pop_col = "pop_sw") |> fill_missing(),
    pop_by_class(bivar_pop_un, paste0(ind, "_pop_unimproved"), class_labels_pop_un,
                 wash_col = "pop_un_class", pop_col = "pop_un") |> fill_missing()
  )

  pop_csv <- sprintf("results/bivar_%s_population_by_class_%d.csv", ind, map_year)
  write_csv(pop_all, pop_csv)
  message("  Saved: ", pop_csv)

  make_bivar_legend_pop_single(
    bivar_palette,
    x_label = "Surface water use ->",
    y_label = y_label,
    pop_df = pop_all |> filter(matrix_type == paste0(ind, "_surfwater")),
    filename = sprintf("figures/legends/bivar_legend_%s_surfwater_ratio_pop.png", ind)
  )

  make_bivar_legend_pop_single(
    bivar_palette,
    x_label = "Unimproved water use ->",
    y_label = y_label,
    pop_df = pop_all |> filter(matrix_type == paste0(ind, "_unimproved")),
    filename = sprintf("figures/legends/bivar_legend_%s_unimproved_ratio_pop.png", ind)
  )

  make_bivar_legend_pop_single(
    bivar_palette,
    x_label = "Pop using surface water (log10) ->",
    y_label = y_label,
    pop_df = pop_all |> filter(matrix_type == paste0(ind, "_pop_surfwater")),
    filename = sprintf("figures/legends/bivar_legend_%s_pop_surfwater_pop.png", ind)
  )

  make_bivar_legend_pop_single(
    bivar_palette,
    x_label = "Pop using unimproved (log10) ->",
    y_label = y_label,
    pop_df = pop_all |> filter(matrix_type == paste0(ind, "_pop_unimproved")),
    filename = sprintf("figures/legends/bivar_legend_%s_pop_unimproved_pop.png", ind)
  )
}

message("\nScript 31b completed successfully for pathogen, organic, and salinity.")
