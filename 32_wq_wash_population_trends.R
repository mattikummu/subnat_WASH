# 32_wq_wash_population_trends.R
# ---------------------------------------------------------------------------
# Stacked area charts: population using surface water / unimproved water,
# stratified by water quality class (overall WQ and pathogen WQ), 1990–2019.
#
# Four panels arranged in a 2×2 grid:
#   [1] Overall WQ   × Surface water
#   [2] Overall WQ   × Unimproved
#   [3] Pathogen WQ  × Surface water
#   [4] Pathogen WQ  × Unimproved
#
# WQ class scale: 1–5 (1 = worst quality, 5 = best quality)
# WQ class bins used:
#   Overall WQ  : 1–2, 2–3, 3–4, 4–5
#   Pathogen WQ : 1–2, 2–3, 3–4, 4–5
#
# Population for each admin unit × year:
#   pop_sw = pop × Surface.water / 100
#   pop_un = pop × Unimproved    / 100
#
# Inputs:
#   results/wq_adm1_extracted.csv        – script 21
#   results/wash_adm1_ratio_adjusted.csv – script 11
#
# Outputs:
#   figures/wq_wash_population_trends.png   – 2×2 stacked area chart
#   results/wq_wash_population_trends.csv   – underlying trend data
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

dir.create("figures", showWarnings = FALSE)
dir.create("figures/legends", recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

### 0. Settings ----
years_in <- 1990:2019

# WQ class break points and labels — scale 1–5 (1=worst, 5=best)
overall_breaks  <- c(1, 2, 3, 4, 5)
pathogen_breaks <- c(1, 2, 3, 4, 5)

overall_labels  <- c("1\u20132", "2\u20133", "3\u20134", "4\u20135")
pathogen_labels <- c("1\u20132", "2\u20133", "3\u20134", "4\u20135")

### 1. Load data ----
message("== Loading data ==")

wq_long <- read_csv("results/wq_adm1_extracted.csv", show_col_types = FALSE) |>
  filter(year %in% years_in, indicator %in% c("overall", "pathogen")) |>
  select(adm_id, year, indicator, wq_mean)

wash_long <- read_csv("results/wash_adm1_ratio_adjusted.csv",
                      show_col_types = FALSE) |>
  filter(year %in% years_in) |>
  select(adm_id, year, pop, Surface.water, Unimproved) |>
  mutate(
    pop_sw = pop * Surface.water / 100,
    pop_un = pop * Unimproved    / 100
  )

message(sprintf("  WQ records:   %d", nrow(wq_long)))
message(sprintf("  WASH records: %d", nrow(wash_long)))

### 2. Join and classify ----
# Wide: one row per adm_id × year, with wq_overall and wq_pathogen columns
wq_wide <- wq_long |>
  pivot_wider(names_from = indicator, values_from = wq_mean,
              names_prefix = "wq_")

combined <- inner_join(wq_wide, wash_long, by = c("adm_id", "year"))
message(sprintf("  Joined records: %d  (%d admin units)",
                nrow(combined), n_distinct(combined$adm_id)))

combined <- combined |>
  mutate(
    overall_class  = cut(wq_overall,  breaks = overall_breaks,
                         labels = overall_labels,
                         include.lowest = TRUE, right = TRUE),
    pathogen_class = cut(wq_pathogen, breaks = pathogen_breaks,
                         labels = pathogen_labels,
                         include.lowest = TRUE, right = TRUE)
  )

### 3. Summarise: population (millions) by year × WQ class ----

# Helper: sum pop_sw or pop_un by year × class
agg_pop <- function(df, wq_col, pop_col, indicator_label, wash_label) {
  df |>
    filter(!is.na(.data[[wq_col]])) |>
    summarise(
      pop_M = sum(.data[[pop_col]], na.rm = TRUE) / 1e6,
      .by   = c(year, !!sym(wq_col))
    ) |>
    rename(wq_class = !!sym(wq_col)) |>
    mutate(indicator = indicator_label, wash = wash_label)
}

trend_data <- bind_rows(
  agg_pop(combined, "overall_class",  "pop_sw", "Overall WQ",  "Surface water"),
  agg_pop(combined, "overall_class",  "pop_un", "Overall WQ",  "Unimproved"),
  agg_pop(combined, "pathogen_class", "pop_sw", "Pathogen WQ", "Surface water"),
  agg_pop(combined, "pathogen_class", "pop_un", "Pathogen WQ", "Unimproved")
)

write_csv(trend_data, "results/wq_wash_population_trends.csv")
message("  Saved: results/wq_wash_population_trends.csv")

# Figure-linked table: first and last year by panel sector
first_year <- min(trend_data$year, na.rm = TRUE)
last_year  <- max(trend_data$year, na.rm = TRUE)

trend_totals <- trend_data |>
  summarise(total_pop_M = sum(pop_M, na.rm = TRUE), .by = c(indicator, wash, year))

trend_first_last <- trend_data |>
  left_join(trend_totals, by = c("indicator", "wash", "year")) |>
  mutate(pct = if_else(total_pop_M > 0, pop_M / total_pop_M * 100, NA_real_)) |>
  filter(year %in% c(first_year, last_year)) |>
  arrange(indicator, wash, wq_class, year) |>
  pivot_wider(
    names_from = year,
    values_from = c(pop_M, pct),
    names_glue = "year_{year}_{.value}"
  ) |>
  rename(
    pop_M_first = paste0("year_", first_year, "_pop_M"),
    pct_first   = paste0("year_", first_year, "_pct"),
    pop_M_last  = paste0("year_", last_year,  "_pop_M"),
    pct_last    = paste0("year_", last_year,  "_pct")
  )

write_csv(trend_first_last, "results/fig_wq_wash_population_trends_first_last.csv")
message("  Saved: results/fig_wq_wash_population_trends_first_last.csv")

### 4. Colour palettes — Roma scale, matching script 22 (c4a "roma", reverse) ----
# Blues from range c(0, 0.85), reds from full-range tail
# cols4all::c4a("roma", reverse = TRUE, n = 10) → positions 2, 4, 9, 10
wq_pal_4 <- c(
  "1\u20132" = "#badfc5", # dark blue        – class 1 (lowest contamination)
  "2\u20133" = "#cbc165", # cyan
  "3\u20134" = "#ab772b", # lighter red
  "4\u20135" = "#88321c" # dark red         – class 5 (highest contamination)
)

pal_overall  <- wq_pal_4
pal_pathogen <- wq_pal_4

### 5. Panel plot function ----
make_panel <- function(df, palette, title, show_x = TRUE, show_y = TRUE) {
  # Ensure factor order matches palette (low → high = bottom → top of stack)
  df <- df |>
    mutate(wq_class = factor(wq_class, levels = names(palette)))

  ggplot(df, aes(x = year, y = pop_M, fill = wq_class)) +
    geom_area(position = "stack", alpha = 0.92,
              colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = palette, name = NULL,
                      guide  = guide_legend(reverse = FALSE, direction = "horizontal", nrow = 1)) +
    scale_x_continuous(breaks = seq(1990, 2019, by = 5),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(
      title = title,
      x = if (show_x) "Year" else NULL,
      y = if (show_y) "Population (millions)" else NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.key.height  = unit(4, "mm"),
      legend.key.width   = unit(3, "mm"),
      legend.text        = element_text(size = 7),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title         = element_text(size = 9, face = "bold"),
      axis.text          = element_text(size = 7),
      axis.title.y       = if (show_y) element_text() else element_blank(),
      axis.text.y        = if (show_y) element_text(size = 7) else element_blank(),
      axis.ticks.y       = if (show_y) element_line() else element_blank(),
      axis.title.x       = if (show_x) element_text() else element_blank(),
      axis.text.x        = if (show_x) element_text(size = 7) else element_blank(),
      axis.ticks.x       = if (show_x) element_line() else element_blank()
    )
}

### 6. Build panels ----
message("\n== Rendering panels ==")

p1 <- make_panel(
  trend_data |> filter(indicator == "Overall WQ",  wash == "Surface water"),
  pal_overall,  "Overall WQ \u00d7 Surface water",
  show_x = TRUE, show_y = TRUE
)
p2 <- make_panel(
  trend_data |> filter(indicator == "Overall WQ",  wash == "Unimproved"),
  pal_overall,  "Overall WQ \u00d7 Unimproved",
  show_x = FALSE, show_y = TRUE
)
p3 <- make_panel(
  trend_data |> filter(indicator == "Pathogen WQ", wash == "Surface water"),
  pal_pathogen, "Pathogen WQ \u00d7 Surface water",
  show_x = TRUE, show_y = FALSE
)
p4 <- make_panel(
  trend_data |> filter(indicator == "Pathogen WQ", wash == "Unimproved"),
  pal_pathogen, "Pathogen WQ \u00d7 Unimproved",
  show_x = FALSE, show_y = FALSE
)

### 7. Combine into 2×2 and save ----
fig <- (p2 | p4) / (p1 | p3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )

legend_plot <- p1 +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )

legend_grob <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")

out_file <- "figures/wq_wash_population_trends.png"
ggsave(out_file, fig, width = 180, height = 130, units = "mm", dpi = 600)
message("  Saved: ", out_file)

legend_file <- "figures/legends/wq_wash_population_trends_legend.png"
ggsave(legend_file, legend_grob, width = 120, height = 20, units = "mm", dpi = 600)
message("  Saved: ", legend_file)

message("\nScript 32 completed successfully.")
