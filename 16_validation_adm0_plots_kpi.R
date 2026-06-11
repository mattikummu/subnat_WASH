# 16_validation_adm0_plots_kpi.R
# ---------------------------------------------------------------------------
# Plot ADM0 validation results and compute KPI metrics from script 15 outputs.
#
# Inputs:
#   results/validation_adm0_compare_common_2000_2017.csv
#   results/validation_adm0_compare_our_full.csv
#
# Outputs (tables):
#   results/validation_adm0_kpi_common_source_indicator.csv
#   results/validation_adm0_kpi_common_source_indicator_year.csv
#   results/validation_adm0_kpi_our_full_indicator.csv
#   results/validation_adm0_kpi_our_full_indicator_year.csv
#
# Outputs (figures):
#   figures/validation/validation_scatter_common.png
#   figures/validation/validation_scatter_our_full.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

dir.create("results", showWarnings = FALSE)
dir.create("figures/validation", recursive = TRUE, showWarnings = FALSE)

source_colors <- c(
  our = "#009E73",
  ihme = "#CC79A7",
  eawag = "#D55E00"
)

source_shapes <- c(
  our = 16,
  ihme = 17,
  eawag = 15
)

source_labels <- c(
  our = "our",
  ihme = "Deshpande et al (2020)",
  eawag = "Greenwood et al (2024)"
)

cmp_common <- read_csv(
  "results/validation_adm0_compare_common_2000_2017.csv",
  show_col_types = FALSE
)

cmp_our_full <- read_csv(
  "results/validation_adm0_compare_our_full.csv",
  show_col_types = FALSE
)

kpi_summarise <- function(df) {
  df |>
    summarise(
      n = sum(!is.na(value) & !is.na(jmp_value)),
      mean_diff = mean(diff, na.rm = TRUE),
      mae = mean(abs_diff, na.rm = TRUE),
      rmse = sqrt(mean(diff^2, na.rm = TRUE)),
      mape_gt1 = mean(abs(diff) / if_else(abs(jmp_value) <= 1, NA_real_, abs(jmp_value)), na.rm = TRUE) * 100,
      smape = mean(2 * abs(diff) / (abs(value) + abs(jmp_value)), na.rm = TRUE) * 100,
      pearson_r = suppressWarnings(cor(value, jmp_value, use = "complete.obs")),
      r_squared = pearson_r^2,
      .groups = "drop"
    )
}

kpi_common_src_ind <- cmp_common |>
  filter(!is.na(value), !is.na(jmp_value)) |>
  group_by(source, indicator) |>
  kpi_summarise() |>
  arrange(indicator, source)

kpi_common_src_ind_year <- cmp_common |>
  filter(!is.na(value), !is.na(jmp_value)) |>
  group_by(source, indicator, year) |>
  kpi_summarise() |>
  arrange(indicator, source, year)

kpi_our_full_ind <- cmp_our_full |>
  filter(!is.na(value), !is.na(jmp_value)) |>
  group_by(indicator) |>
  kpi_summarise() |>
  arrange(indicator)

kpi_our_full_ind_year <- cmp_our_full |>
  filter(!is.na(value), !is.na(jmp_value)) |>
  group_by(indicator, year) |>
  kpi_summarise() |>
  arrange(indicator, year)

write_csv(kpi_common_src_ind, "results/validation_adm0_kpi_common_source_indicator.csv")
write_csv(kpi_common_src_ind_year, "results/validation_adm0_kpi_common_source_indicator_year.csv")
write_csv(kpi_our_full_ind, "results/validation_adm0_kpi_our_full_indicator.csv")
write_csv(kpi_our_full_ind_year, "results/validation_adm0_kpi_our_full_indicator_year.csv")

common_clean <- cmp_common |>
  filter(!is.na(value), !is.na(jmp_value)) |>
  mutate(indicator = factor(indicator, levels = c("Unimproved", "Surface.water")))

scatter_metrics <- common_clean |>
  group_by(indicator, source) |>
  summarise(
    rmse = sqrt(mean((value - jmp_value)^2, na.rm = TRUE)),
    r_squared = suppressWarnings(cor(value, jmp_value, use = "complete.obs")^2),
    .groups = "drop"
  )

scatter_ranges <- common_clean |>
  group_by(indicator) |>
  summarise(
    x_min = min(jmp_value, na.rm = TRUE),
    x_max = max(jmp_value, na.rm = TRUE),
    y_min = min(value, na.rm = TRUE),
    y_max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

scatter_annot <- scatter_metrics |>
  left_join(scatter_ranges, by = "indicator") |>
  group_by(indicator) |>
  arrange(source, .by_group = TRUE) |>
  mutate(
    idx = row_number(),
    x_pos = x_min + 0.03 * (x_max - x_min),
    y_pos = y_max - (idx - 1) * 0.08 * (y_max - y_min),
    source_kpi = recode(source,
      ihme = "Deshpande",
      eawag = "Greenwood",
      .default = source
    ),
    label = sprintf("%s: R²=%.2f, RMSE=%.2f", source_kpi, r_squared, rmse)
  ) |>
  ungroup()

p_scatter_common <- ggplot(common_clean, aes(x = jmp_value, y = value, color = source, shape = source)) +
  geom_abline(intercept = 0, slope = 1, color = "grey45", linewidth = 0.4) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_text(
    data = scatter_annot,
    aes(x = x_pos, y = y_pos, label = label, color = source),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    show.legend = FALSE
  ) +
  facet_wrap(~ indicator, scales = "fixed") +
  coord_equal() +
  scale_color_manual(name = "Source", values = source_colors, labels = source_labels) +
  scale_shape_manual(name = "Source", values = source_shapes, labels = source_labels) +
  labs(
    x = "JMP reported value (%)",
    y = "Estimated value (%)",
    color = "Source",
    shape = "Source"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  "figures/validation/validation_scatter_common.png",
  p_scatter_common,
  width = 180,
  height = 100,
  units = "mm",
  dpi = 300
)

our_full_clean <- cmp_our_full |>
  filter(!is.na(value), !is.na(jmp_value))

scatter_our_metrics <- our_full_clean |>
  group_by(indicator) |>
  summarise(
    rmse = sqrt(mean((value - jmp_value)^2, na.rm = TRUE)),
    r_squared = suppressWarnings(cor(value, jmp_value, use = "complete.obs")^2),
    .groups = "drop"
  ) |>
  mutate(
    x_pos = -Inf,
    y_pos = Inf,
    label = sprintf("R²=%.2f, RMSE=%.2f", r_squared, rmse)
  )

p_scatter_our_full <- ggplot(our_full_clean, aes(x = jmp_value, y = value)) +
  geom_abline(intercept = 0, slope = 1, color = "grey45", linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 0.8, color = source_colors[["our"]]) +
  geom_text(
    data = scatter_our_metrics,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.35,
    size = 2.8,
    color = source_colors[["our"]]
  ) +
  facet_wrap(~ indicator, scales = "fixed") +
  coord_equal() +
  labs(
    x = "JMP reported value (%)",
    y = "Estimated value (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/validation/validation_scatter_our_full.png",
  p_scatter_our_full,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 300
)
