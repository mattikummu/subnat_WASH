# 14_wash_population_estimates.R
# ---------------------------------------------------------------------------
# Population on each WASH water-service ladder, 1990–2023.
# Focus: people WITHOUT basic water service (Limited + Unimproved + Surface.water).
#
# Part 1: Extract population × WASH share → pop per service level (ADM1)
# Part 2: Global / country / macro-region summaries + ggplot trend lines
# Part 3: Theil-Sen slopes for pop without basic water + trend direction maps
#
# Inputs:
#   results_gis/polyg_wash_adm1_1990_2023_geomonly_simpl.gpkg – script 05b
#   results/wash_adm1_ratio_adjusted.csv                 – script 11
#   results/global_adm1_predictors_1990_2023.csv         – script 05 (pop cache)
#   data_in/cntry_ids.csv
#
# Outputs:
#   results/wash_population_by_service_adm1.csv
#   results/wash_population_global_summary.csv
#   results/wash_population_country_summary.csv
#   results/wash_population_macroregion_summary.csv
#   results/wash_population_slopes_adm1.csv
#   results/wash_population_trend_summary_global.csv
#   results/wash_population_trend_summary_macroregion.csv
#   figures/wash_global_population_trend.png
#   figures/wash_global_headcount_trend.png
#   figures/wash_macroregion_population_trend.png
#   figures/wash_macroregion_headcount_trend.png
#   figures/maps/wash_population_{indicator}_slope.png
#   figures/legends/wash_legend_population_{indicator}_slope.png
# ---------------------------------------------------------------------------

library(terra)
library(sf)
library(tmap)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(broom)
library(mblm)
library(ggplot2)
library(ggrepel)
library(cols4all)

source("cache_prompt.R")

### 0. Setup ----
years_in   <- 1990:2023
first_year <- years_in[1]
last_year  <- years_in[length(years_in)]
indicators <- c("Basic", "Limited", "Unimproved", "Surface.water")
# Indicators of concern: everything except Basic
indicators_concern <- c("Limited", "Unimproved", "Surface.water")

dir.create("results",         showWarnings = FALSE)
dir.create("figures/maps",    recursive = TRUE, showWarnings = FALSE)
dir.create("figures/legends", recursive = TRUE, showWarnings = FALSE)

# Build country-to-region lookup from the countrycode package
cntry_info <- countrycode::codelist |>
  as_tibble() |>
  filter(!is.na(iso3c)) |>
  select(iso3 = iso3c, region_name = region) |>
  filter(!is.na(region_name)) |>
  distinct(iso3, .keep_all = TRUE)

# Create numeric region IDs from unique region names
region_names <- cntry_info |>
  distinct(region_name) |>
  arrange(region_name) |>
  mutate(region_id = row_number())

cntry_info <- cntry_info |>
  left_join(region_names, by = "region_name")

### 1. Load population from global predictor cache (script 05) ----
message("== Loading ADM1 population from global predictors (script 05) ==")

pop_long <- read_csv("results/global_adm1_predictors_1990_2023.csv",
                     show_col_types = FALSE) |>
  filter(year %in% years_in) |>
  select(adm_id, year, population = pop) |>
  distinct(adm_id, year, .keep_all = TRUE)

message(sprintf("  %d adm_id × year rows", nrow(pop_long)))

### 2. Load WASH polygons ----
message("== Loading WASH polygons ==")

gpkg_simpl_path <- sprintf("results_gis/polyg_wash_adm1_%d_%d_geomonly_simpl.gpkg",
                           first_year, last_year)

### 3. Save reused population table (for compatibility) ----
pop_cache_path <- sprintf("results/wash_population_extract_adm1_%d_%d.csv",
                          first_year, last_year)
write_csv(pop_long, pop_cache_path)
message(sprintf("  Saved: %s", pop_cache_path))

### 4. Load ratio-adjusted WASH data ----
message("== Loading WASH data ==")

wash_data <- read_csv("results/wash_adm1_ratio_adjusted.csv",
                      show_col_types = FALSE) |>
  filter(year %in% years_in) |>
  select(adm_id, iso3, year, all_of(indicators))

# Compute "not basic" = Limited + Unimproved + Surface.water
wash_data <- wash_data |>
  mutate(Not_basic = Limited + Unimproved + Surface.water)

all_indicators <- c(indicators, "Not_basic")

### 5. Join population × WASH shares ----
message("== Calculating population per service level ==")

pop_wash <- pop_long |>
  left_join(wash_data, by = c("adm_id", "year"))

# Population in each service level
for (ind in all_indicators) {
  pop_wash[[paste0("pop_", ind)]] <- pop_wash$population * (pop_wash[[ind]] / 100)
}

# Long format (one row per adm_id × year × indicator)
pop_wash_long <- pop_wash |>
  pivot_longer(
    cols      = c(all_of(all_indicators), starts_with("pop_Basic"),
                  starts_with("pop_Limited"), starts_with("pop_Unimproved"),
                  starts_with("pop_Surface.water"), starts_with("pop_Not_basic")),
    names_to  = "col",
    values_to = "val"
  ) |>
  mutate(
    type      = if_else(startsWith(col, "pop_"), "pop_count", "share_pct"),
    indicator = sub("^pop_", "", col)
  ) |>
  filter(indicator %in% all_indicators) |>
  select(-col) |>
  pivot_wider(names_from = type, values_from = val)

### 6. Global / country / macro-region summaries ----
message("== Computing summaries ==")

admID_iso3_map <- wash_data |>
  select(adm_id, iso3) |>
  distinct()

global_summary <- pop_wash_long |>
  group_by(indicator, year) |>
  summarise(
    total_population = sum(population,  na.rm = TRUE),
    total_pop_count  = sum(pop_count,   na.rm = TRUE),
    global_share_pct = total_pop_count / total_population * 100,
    .groups = "drop"
  )

country_summary <- pop_wash_long |>
  left_join(cntry_info |> select(iso3, region_id),   by = "iso3") |>
  left_join(region_names,                            by = "region_id") |>
  group_by(indicator, iso3, region_id, region_name, year) |>
  summarise(
    total_population = sum(population,  na.rm = TRUE),
    total_pop_count  = sum(pop_count,   na.rm = TRUE),
    country_share_pct = total_pop_count / total_population * 100,
    .groups = "drop"
  )

macroregion_summary <- country_summary |>
  group_by(indicator, region_id, region_name, year) |>
  summarise(
    total_population = sum(total_population, na.rm = TRUE),
    total_pop_count  = sum(total_pop_count,  na.rm = TRUE),
    region_share_pct = total_pop_count / total_population * 100,
    n_countries      = n_distinct(iso3),
    .groups = "drop"
  )

### 7. Save tabular results ----
message("== Saving tabular results ==")

write_csv(pop_wash_long,       "results/wash_population_by_service_adm1.csv")
write_csv(global_summary,      "results/wash_population_global_summary.csv")
write_csv(country_summary,     "results/wash_population_country_summary.csv")
write_csv(macroregion_summary, "results/wash_population_macroregion_summary.csv")
message("  Saved to results/")

### 8. ggplot trend figures ----
message("== Creating global trend plots ==")

ind_labels <- c(
  Basic         = "Basic",
  Limited       = "Limited",
  Unimproved    = "Unimproved",
  Surface.water = "Surface water",
  Not_basic     = "Without basic service"
)

ind_colours <- c(
  Basic         = "#2166ac",
  Limited       = "#fc8d59",
  Unimproved    = "#d73027",
  Surface.water = "#4d4d4d",
  Not_basic     = "#b2182b"
)

region_colours <- c(
  "East Asia & Pacific"        = "#7570b3",
  "Europe & Central Asia"      = "#e7298a",
  "Latin America & Caribbean"  = "#d95f02",
  "Middle East & North Africa" = "#666666",
  "North America"              = "#bebada",
  "South Asia"                 = "#66a61e",
  "Sub-Saharan Africa"         = "#a6761d"
)

# Plot 1: Global population per service level (concern indicators + Not_basic)
plot_indicators <- c(indicators_concern, "Not_basic")

p1 <- global_summary |>
  filter(indicator %in% plot_indicators) |>
  mutate(
    label = ind_labels[indicator],
    indicator = factor(indicator, levels = plot_indicators)
  ) |>
  ggplot(aes(x = year, y = total_pop_count / 1e9, colour = indicator)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_text_repel(
    data = \(d) d |> filter(year == last_year),
    aes(label = label), hjust = 0, nudge_x = 0.5,
    direction = "y", size = 3.5, segment.size = 0.2, segment.colour = "grey50"
  ) +
  scale_colour_manual(values = ind_colours) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(
    title = sprintf("Global Population Without Basic Water Service (%d\u2013%d)",
                    first_year, last_year),
    x = "Year", y = "Population (billions)"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(size = 11))

ggsave("figures/wash_global_population_trend.png",
       p1, width = 10, height = 6, dpi = 300)

# Plot 2: Global share (%) — concern indicators
p2 <- global_summary |>
  filter(indicator %in% plot_indicators) |>
  mutate(
    label = ind_labels[indicator],
    indicator = factor(indicator, levels = plot_indicators)
  ) |>
  ggplot(aes(x = year, y = global_share_pct, colour = indicator)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_text_repel(
    data = \(d) d |> filter(year == last_year),
    aes(label = label), hjust = 0, nudge_x = 0.5,
    direction = "y", size = 3.5, segment.size = 0.2, segment.colour = "grey50"
  ) +
  scale_colour_manual(values = ind_colours) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(
    title = sprintf("Global Share Without Basic Water Service (%d\u2013%d)",
                    first_year, last_year),
    x = "Year", y = "Share of population (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(size = 11))

ggsave("figures/wash_global_headcount_trend.png",
       p2, width = 10, height = 6, dpi = 300)

# Plot 3: Macro-region population by region — one panel per indicator
for (ind in plot_indicators) {
  p_reg_pop <- macroregion_summary |>
    filter(indicator == ind) |>
    ggplot(aes(x = year, y = total_pop_count / 1e6, colour = region_name)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1) +
    geom_text_repel(
      data = \(d) d |> filter(year == last_year),
      aes(label = region_name), hjust = 0, nudge_x = 1,
      direction = "y", size = 3, segment.size = 0.2, segment.colour = "grey50"
    ) +
    scale_colour_manual(values = region_colours) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
    labs(
      title = sprintf("%s: Population by Region (%d\u2013%d)",
                      ind_labels[ind], first_year, last_year),
      x = "Year", y = "Population (millions)"
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 12),
          axis.title = element_text(size = 11))

  ggsave(sprintf("figures/wash_macroregion_population_trend_%s.png", ind),
         p_reg_pop, width = 12, height = 7, dpi = 300)
}

# Plot 4: Macro-region share (%) by region — one panel per indicator
for (ind in plot_indicators) {
  p_reg_share <- macroregion_summary |>
    filter(indicator == ind) |>
    ggplot(aes(x = year, y = region_share_pct, colour = region_name)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1) +
    geom_text_repel(
      data = \(d) d |> filter(year == last_year),
      aes(label = region_name), hjust = 0, nudge_x = 1,
      direction = "y", size = 3, segment.size = 0.2, segment.colour = "grey50"
    ) +
    scale_colour_manual(values = region_colours) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
    labs(
      title = sprintf("%s: Share by Region (%d\u2013%d)",
                      ind_labels[ind], first_year, last_year),
      x = "Year", y = "Share of population (%)"
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 12),
          axis.title = element_text(size = 11))

  ggsave(sprintf("figures/wash_macroregion_headcount_trend_%s.png", ind),
         p_reg_share, width = 12, height = 7, dpi = 300)
}

# Plot 5: Stacked area — global population by Limited / Unimproved / Surface water
stack_levels <- c("Surface.water", "Unimproved", "Limited")
stack_labels <- ind_labels[stack_levels]

stack_data <- global_summary |>
  filter(indicator %in% indicators_concern) |>
  mutate(
    indicator = factor(indicator,
                       levels = stack_levels,
                       labels = stack_labels)
  )

stack_colours <- setNames(
  ind_colours[stack_levels],
  stack_labels
)

make_stack_plot <- function(df_in, end_year) {
  ggplot(df_in,
    aes(x = year, y = total_pop_count / 1e9,
        fill = indicator)) +
    geom_area(alpha = 0.85, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = stack_colours) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    labs(
      title = sprintf("Global Population by Water-Service Deficit (%d\u2013%d)",
       first_year, end_year),
      x = "Year", y = "Population (billions)",
      fill = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom",
     plot.title = element_text(face = "bold", size = 13),
     axis.title = element_text(size = 11))
}

year_ref <- 2019

p5_2023 <- make_stack_plot(stack_data, end_year = last_year)
ggsave("figures/wash_global_population_stacked_2023.png",
  p5_2023, width = 10, height = 6, dpi = 300)

# Keep legacy filename pointing to the full-period (latest year) version
ggsave("figures/wash_global_population_stacked.png",
  p5_2023, width = 10, height = 6, dpi = 300)

p5_2019 <- make_stack_plot(
  stack_data |> filter(year <= year_ref),
  end_year = year_ref
)
ggsave("figures/wash_global_population_stacked_2019.png",
  p5_2019, width = 10, height = 6, dpi = 300)

# Figure-linked table: first and last year by stacked sector
stack_first_last <- global_summary |>
  filter(indicator %in% stack_levels, year %in% c(first_year, year_ref, last_year)) |>
  mutate(pop_M = total_pop_count / 1e6) |>
  select(indicator, year, pop_M, global_share_pct) |>
  arrange(indicator, year) |>
  pivot_wider(
    names_from = year,
    values_from = c(pop_M, global_share_pct),
    names_glue = "year_{year}_{.value}"
  ) |>
  mutate(
    pop_M_first = .data[[paste0("year_", first_year, "_pop_M")]],
    pct_first   = .data[[paste0("year_", first_year, "_global_share_pct")]],
    pop_M_2019  = .data[[paste0("year_", year_ref, "_pop_M")]],
    pct_2019    = .data[[paste0("year_", year_ref, "_global_share_pct")]],
    pop_M_last  = .data[[paste0("year_", last_year, "_pop_M")]],
    pct_last    = .data[[paste0("year_", last_year, "_global_share_pct")]]
  ) |>
  select(indicator, pop_M_first, pct_first, pop_M_2019, pct_2019, pop_M_last, pct_last)

write_csv(stack_first_last, "results/fig_wash_global_population_stacked_first_last.csv")
message("  Saved: results/fig_wash_global_population_stacked_first_last.csv")

message("  ggplot figures saved to figures/")

### 9. Theil-Sen slopes for population without basic water ----
message("== Computing Theil-Sen slopes ==")

# Compute slopes for the concern indicators + Not_basic
slope_indicators <- c(indicators_concern, "Not_basic")

calc_pop_slope <- function(indicator_name) {
  message(sprintf("  %s", indicator_name))
  pop_wash_long |>
    filter(indicator == indicator_name, !is.na(pop_count)) |>
    group_by(adm_id) |>
    mutate(time = row_number()) |>
    ungroup() |>
    nest(.by = adm_id) |>
    mutate(
      model = map(data, \(d) {
        if (nrow(d) < 3) return(NULL)
        tryCatch(
          mblm(pop_count ~ time, data = d, repeated = FALSE),
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
    transmute(
      adm_id,
      indicator  = indicator_name,
      slope      = estimate,
      p.value,
      p.value_sig = p.value < 0.05,
      trend_direction = case_when(
        !p.value_sig ~ "no_change",
        slope < 0    ~ "decreasing",
        slope > 0    ~ "increasing",
        TRUE         ~ "no_change"
      )
    )
}

all_pop_slopes <- bind_rows(lapply(slope_indicators, calc_pop_slope))

# Attach iso3 / region info
all_pop_slopes_geo <- all_pop_slopes |>
  left_join(admID_iso3_map, by = "adm_id") |>
  left_join(cntry_info |> select(iso3, region_id), by = "iso3") |>
  left_join(region_names, by = "region_id")

write_csv(all_pop_slopes_geo, "results/wash_population_slopes_adm1.csv")

# Latest-year population for weighting
pop_latest <- pop_long |> filter(year == last_year)

# Global trend summary
trend_global <- all_pop_slopes |>
  left_join(pop_latest, by = "adm_id") |>
  group_by(indicator, trend_direction) |>
  summarise(
    total_population = sum(population, na.rm = TRUE),
    n_admin_units    = n(),
    .groups = "drop"
  ) |>
  group_by(indicator) |>
  mutate(pct_population = total_population / sum(total_population) * 100) |>
  ungroup()

write_csv(trend_global, "results/wash_population_trend_summary_global.csv")

# Macro-region trend summary
trend_macroregion <- all_pop_slopes_geo |>
  left_join(pop_latest, by = "adm_id") |>
  filter(!is.na(region_id)) |>
  group_by(indicator, region_id, region_name, trend_direction) |>
  summarise(
    total_population = sum(population, na.rm = TRUE),
    n_admin_units    = n(),
    .groups = "drop"
  ) |>
  group_by(indicator, region_id, region_name) |>
  mutate(pct_population = total_population / sum(total_population) * 100) |>
  ungroup()

write_csv(trend_macroregion, "results/wash_population_trend_summary_macroregion.csv")

message("  Slope results saved to results/")

### 10. Slope maps (Not_basic) ----
message("== Creating population-slope maps ==")

adm0_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
  filter(iso_a3 != "ATA")

ref_sf <- st_read(gpkg_simpl_path, quiet = TRUE) |>
  select(adm_id)

# Build slope sf for each concern indicator + Not_basic
slope_sf_list <- lapply(slope_indicators, \(ind) {
  slope_data <- all_pop_slopes |>
    filter(indicator == ind) |>
    mutate(
      slope_sig   = if_else(p.value_sig, slope, NA_real_),
      slope_log10 = sign(slope_sig) * log10(abs(slope_sig) + 1)
    ) |>
    select(adm_id, slope_sig, slope_log10, p.value, p.value_sig, trend_direction)

  ref_sf |>
    select(adm_id) |>
    left_join(slope_data, by = "adm_id")
}) |> setNames(slope_indicators)

# Symmetric breaks
all_log_slopes <- unlist(lapply(slope_sf_list, \(s) s$slope_log10))
abs_max_log    <- max(abs(quantile(all_log_slopes, c(0.05, 0.95), na.rm = TRUE)))
accur_log      <- round(abs_max_log * 2 / 20, digits = 5)

slope_breaks_log <- seq(
  plyr::round_any(-abs_max_log, accuracy = accur_log, f = floor),
  plyr::round_any( abs_max_log, accuracy = accur_log, f = ceiling),
  by = plyr::round_any(abs_max_log * 2 / 20, accuracy = 0.00001, f = floor)
)

pal_slope <- cols4all::c4a("bu_br_div")

ind_pop_slope_titles <- c(
  Limited       = "Pop. with limited service: slope (log\u2081\u2080 people/yr)",
  Unimproved    = "Pop. with unimproved source: slope (log\u2081\u2080 people/yr)",
  Surface.water = "Pop. using surface water: slope (log\u2081\u2080 people/yr)",
  Not_basic     = "Pop. without basic service: slope (log\u2081\u2080 people/yr)"
)

save_map <- function(p_map, filename, width = 80, units = "mm", dpi = 450) {
  tmap_save(p_map + tm_layout(legend.show = FALSE),
            filename = filename, width = width, units = units, dpi = dpi)
}

save_legend <- function(p_map, filename, width = 150, height = 30,
                        units = "mm", dpi = 600) {
  tmap_save(
    p_map + tm_layout(legend.only = TRUE, frame = FALSE,
                      outer.margins = 0, inner.margins = c(0, 0, 0, 0),
                      font.size = 1.5),
    filename = filename, width = width, height = height, units = units, dpi = dpi
  )
}

for (ind in slope_indicators) {
  message(sprintf("  Slope map: %s", ind))

  p <- tm_shape(slope_sf_list[[ind]]) +
    tm_fill(
      fill       = "slope_log10",
      fill.scale = tm_scale(
        values   = pal_slope,
        breaks   = slope_breaks_log,
        midpoint = 0,
        value.na = "white"
      ),
      fill.legend = tm_legend(
        orientation = "landscape",
        position    = tm_pos_out("center", "bottom"),
        frame       = FALSE,
        title       = ind_pop_slope_titles[[ind]]
      )
    ) +
    tm_shape(adm0_sf) +
    tm_borders(col = "grey35", lwd = 0.025) +
    tm_crs("+proj=robin") +
    tm_layout(main.title.position = "center", frame = FALSE)

  save_map(p,
           sprintf("figures/maps/wash_population_%s_slope.png", ind))
  save_legend(p,
              sprintf("figures/legends/wash_legend_population_%s_slope.png", ind))

  message(sprintf("    Saved: wash_population_%s_slope.png", ind))
}

### 11. Summary statistics ----
message("\n== SUMMARY STATISTICS ==")

for (ind in c("Not_basic", indicators_concern)) {
  d_first <- global_summary |> filter(indicator == ind, year == first_year)
  d_last  <- global_summary |> filter(indicator == ind, year == last_year)
  message(sprintf("%s:", ind_labels[ind]))
  message(sprintf("  %d: %.0f M (%.1f%%)",
                  first_year,
                  d_first$total_pop_count / 1e6,
                  d_first$global_share_pct))
  message(sprintf("  %d: %.0f M (%.1f%%)",
                  last_year,
                  d_last$total_pop_count / 1e6,
                  d_last$global_share_pct))
  message(sprintf("  Change: %+.0f M",
                  (d_last$total_pop_count - d_first$total_pop_count) / 1e6))
}

message("\nDone.")
