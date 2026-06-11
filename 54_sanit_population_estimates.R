# 54_sanit_population_estimates.R
# ---------------------------------------------------------------------------
# Population on each sanitation-service ladder, 1990-2023.
# Focus: people WITHOUT basic sanitation (Limited + Unimproved + Open defecation).
#
# Part 1: Extract population x sanitation share -> pop per service level (ADM1)
# Part 2: Global / country / macro-region summaries + ggplot trend lines
# Part 3: Theil-Sen slopes for pop without basic sanitation + trend maps
#
# Inputs:
#   results_gis_sanit/polyg_sanit_adm1_1990_2023_geomonly_simpl.gpkg – script 05b
#   results_sanit/sanit_adm1_ratio_adjusted.csv
#   results/global_adm1_predictors_1990_2023.csv         – script 05 (pop cache)
#   data_in/cntry_ids.csv
#
# Outputs:
#   results_sanit/sanit_population_by_service_adm1.csv
#   results_sanit/sanit_population_global_summary.csv
#   results_sanit/sanit_population_country_summary.csv
#   results_sanit/sanit_population_macroregion_summary.csv
#   results_sanit/sanit_population_slopes_adm1.csv
#   results_sanit/sanit_population_trend_summary_global.csv
#   results_sanit/sanit_population_trend_summary_macroregion.csv
#   figures_sanit/sanit_global_population_trend.png
#   figures_sanit/sanit_global_headcount_trend.png
#   figures_sanit/sanit_macroregion_population_trend_*.png
#   figures_sanit/sanit_macroregion_headcount_trend_*.png
#   figures_sanit/maps/sanit_population_{indicator}_slope.png
#   figures_sanit/legends/sanit_legend_population_{indicator}_slope.png
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
indicators <- c("s_bas", "s_lim", "s_unimp", "s_od")
# Indicators of concern: everything except Basic
indicators_concern <- c("s_lim", "s_unimp", "s_od")

dir.create("results_sanit",          showWarnings = FALSE)
dir.create("figures_sanit/maps",     recursive = TRUE, showWarnings = FALSE)
dir.create("figures_sanit/legends",  recursive = TRUE, showWarnings = FALSE)

# Build country-to-region lookup
cntry_info <- countrycode::codelist |>
  as_tibble() |>
  filter(!is.na(iso3c)) |>
  select(iso3 = iso3c, region_name = region) |>
  filter(!is.na(region_name)) |>
  distinct(iso3, .keep_all = TRUE)

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

message(sprintf("  %d adm_id x year rows", nrow(pop_long)))

### 2. Load sanitation polygons ----
message("== Loading sanitation polygons ==")

gpkg_simpl_path <- sprintf("results_gis_sanit/polyg_sanit_adm1_%d_%d_geomonly_simpl.gpkg",
                           first_year, last_year)

### 3. Save reused population table (for compatibility) ----
pop_cache_path <- sprintf("results_sanit/sanit_population_extract_adm1_%d_%d.csv",
                          first_year, last_year)
write_csv(pop_long, pop_cache_path)
message(sprintf("  Saved: %s", pop_cache_path))

### 4. Load ratio-adjusted sanitation data ----
message("== Loading sanitation data ==")

sanit_data <- read_csv("results_sanit/sanit_adm1_ratio_adjusted.csv",
                       show_col_types = FALSE) |>
  filter(year %in% years_in) |>
  select(adm_id, iso3, year, all_of(indicators))

# Compute "not basic" = s_lim + s_unimp + s_od
sanit_data <- sanit_data |>
  mutate(Not_basic = s_lim + s_unimp + s_od)

all_indicators <- c(indicators, "Not_basic")

### 5. Join population x sanitation shares ----
message("== Calculating population per service level ==")

pop_sanit <- pop_long |>
  left_join(sanit_data, by = c("adm_id", "year"))

# Population in each service level
for (ind in all_indicators) {
  pop_sanit[[paste0("pop_", ind)]] <- pop_sanit$population * (pop_sanit[[ind]] / 100)
}

# Long format
pop_sanit_long <- pop_sanit |>
  pivot_longer(
    cols      = c(all_of(all_indicators), starts_with("pop_s_bas"),
                  starts_with("pop_s_lim"), starts_with("pop_s_unimp"),
                  starts_with("pop_s_od"), starts_with("pop_Not_basic")),
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

admID_iso3_map <- sanit_data |>
  select(adm_id, iso3) |>
  distinct()

global_summary <- pop_sanit_long |>
  group_by(indicator, year) |>
  summarise(
    total_population = sum(population,  na.rm = TRUE),
    total_pop_count  = sum(pop_count,   na.rm = TRUE),
    global_share_pct = total_pop_count / total_population * 100,
    .groups = "drop"
  )

country_summary <- pop_sanit_long |>
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

write_csv(pop_sanit_long,      "results_sanit/sanit_population_by_service_adm1.csv")
write_csv(global_summary,      "results_sanit/sanit_population_global_summary.csv")
write_csv(country_summary,     "results_sanit/sanit_population_country_summary.csv")
write_csv(macroregion_summary, "results_sanit/sanit_population_macroregion_summary.csv")
message("  Saved to results_sanit/")

### 8. ggplot trend figures ----
message("== Creating global trend plots ==")

ind_labels <- c(
  s_bas     = "Basic",
  s_lim     = "Limited",
  s_unimp   = "Unimproved",
  s_od      = "Open defecation",
  Not_basic = "Without basic service"
)

ind_colours <- c(
  s_bas     = "#2166ac",
  s_lim     = "#fc8d59",
  s_unimp   = "#d73027",
  s_od      = "#4d4d4d",
  Not_basic = "#b2182b"
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

plot_indicators <- c(indicators_concern, "Not_basic")

# Plot 1: Global population per service level
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
    title = sprintf("Global Population Without Basic Sanitation (%d\u2013%d)",
                    first_year, last_year),
    x = "Year", y = "Population (billions)"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(size = 11))

ggsave("figures_sanit/sanit_global_population_trend.png",
       p1, width = 10, height = 6, dpi = 300)

# Plot 2: Global share (%)
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
    title = sprintf("Global Share Without Basic Sanitation (%d\u2013%d)",
                    first_year, last_year),
    x = "Year", y = "Share of population (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(size = 11))

ggsave("figures_sanit/sanit_global_headcount_trend.png",
       p2, width = 10, height = 6, dpi = 300)

# Plot 3: Macro-region population by region
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

  ggsave(sprintf("figures_sanit/sanit_macroregion_population_trend_%s.png", ind),
         p_reg_pop, width = 12, height = 7, dpi = 300)
}

# Plot 4: Macro-region share (%)
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

  ggsave(sprintf("figures_sanit/sanit_macroregion_headcount_trend_%s.png", ind),
         p_reg_share, width = 12, height = 7, dpi = 300)
}

# Plot 5: Stacked area
stack_data <- global_summary |>
  filter(indicator %in% indicators_concern) |>
  mutate(
    indicator = factor(indicator,
                       levels = rev(indicators_concern),
                       labels = rev(ind_labels[indicators_concern]))
  )

stack_colours <- setNames(
  ind_colours[indicators_concern],
  ind_labels[indicators_concern]
)

p5 <- ggplot(stack_data,
             aes(x = year, y = total_pop_count / 1e9,
                 fill = indicator)) +
  geom_area(alpha = 0.85) +
  scale_fill_manual(values = stack_colours) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = sprintf("Global Population by Sanitation-Service Deficit (%d\u2013%d)",
                    first_year, last_year),
    x = "Year", y = "Population (billions)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(size = 11))

ggsave("figures_sanit/sanit_global_population_stacked.png",
       p5, width = 10, height = 6, dpi = 300)

message("  ggplot figures saved to figures_sanit/")

### 9. Theil-Sen slopes for population without basic sanitation ----
message("== Computing Theil-Sen slopes ==")

slope_indicators <- c(indicators_concern, "Not_basic")

calc_pop_slope <- function(indicator_name) {
  message(sprintf("  %s", indicator_name))
  pop_sanit_long |>
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

write_csv(all_pop_slopes_geo, "results_sanit/sanit_population_slopes_adm1.csv")

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

write_csv(trend_global, "results_sanit/sanit_population_trend_summary_global.csv")

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

write_csv(trend_macroregion, "results_sanit/sanit_population_trend_summary_macroregion.csv")

message("  Slope results saved to results_sanit/")

### 10. Slope maps ----
message("== Creating population-slope maps ==")

adm0_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
  filter(iso_a3 != "ATA")

ref_sf <- st_read(gpkg_simpl_path, quiet = TRUE) |>
  select(adm_id)

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
  s_lim     = "Pop. with limited sanitation: slope (log\u2081\u2080 people/yr)",
  s_unimp   = "Pop. with unimproved sanitation: slope (log\u2081\u2080 people/yr)",
  s_od      = "Pop. practising open defecation: slope (log\u2081\u2080 people/yr)",
  Not_basic = "Pop. without basic sanitation: slope (log\u2081\u2080 people/yr)"
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
           sprintf("figures_sanit/maps/sanit_population_%s_slope.png", ind))
  save_legend(p,
              sprintf("figures_sanit/legends/sanit_legend_population_%s_slope.png", ind))

  message(sprintf("    Saved: sanit_population_%s_slope.png", ind))
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
