# 16c_validation_adm0_our_vs_eawag_plots_kpi.R
# ---------------------------------------------------------------------------
# Plot and summarise ADM0 validation for Eawag-scope countries (2020),
# analogous to script 16 structure:
#   - Sources: eawag, our
#   - Reference: reported JMP Basic+Limited
#
# Input:
#   results/validation_adm0_eawag_reported_common_2020.csv
#
# Outputs (tables):
#   results/validation_adm0_eawag_reported_kpi_common_source_indicator.csv
#   results/validation_adm0_eawag_reported_kpi_common_source_indicator_year.csv
#
# Outputs (figures):
#   figures/validation/validation_adm0_eawag_reported_scatter_common.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
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
  "results/validation_adm0_eawag_reported_common_2020.csv",
  show_col_types = FALSE
)

common_clean <- cmp_common |>
  filter(!is.na(value), !is.na(reported_value))

kpi_summarise <- function(df) {
  df |>
    summarise(
      n = sum(!is.na(value) & !is.na(reported_value)),
      mean_diff = mean(diff, na.rm = TRUE),
      mae = mean(abs_diff, na.rm = TRUE),
      rmse = sqrt(mean(diff^2, na.rm = TRUE)),
      mape_gt1 = mean(abs(diff) / if_else(abs(reported_value) <= 1, NA_real_, abs(reported_value)), na.rm = TRUE) * 100,
      smape = mean(2 * abs(diff) / (abs(value) + abs(reported_value)), na.rm = TRUE) * 100,
      pearson_r = suppressWarnings(cor(value, reported_value, use = "complete.obs")),
      r_squared = pearson_r^2,
      .groups = "drop"
    )
}

kpi_common_src_ind <- cmp_common |>
  filter(!is.na(value), !is.na(reported_value)) |>
  group_by(source, indicator) |>
  kpi_summarise() |>
  arrange(indicator, source)

kpi_common_src_ind_year <- cmp_common |>
  filter(!is.na(value), !is.na(reported_value)) |>
  group_by(source, indicator, year) |>
  kpi_summarise() |>
  arrange(indicator, source, year)

write_csv(
  kpi_common_src_ind,
  "results/validation_adm0_eawag_reported_kpi_common_source_indicator.csv"
)

write_csv(
  kpi_common_src_ind_year,
  "results/validation_adm0_eawag_reported_kpi_common_source_indicator_year.csv"
)

p_scatter_common <- ggplot(cmp_common, aes(x = reported_value, y = value, color = source, shape = source)) +
  geom_text(
    data = common_clean |>
      group_by(indicator, source) |>
      summarise(
        rmse = sqrt(mean((value - reported_value)^2, na.rm = TRUE)),
        r_squared = suppressWarnings(cor(value, reported_value, use = "complete.obs")^2),
        x_min = min(reported_value, na.rm = TRUE),
        x_max = max(reported_value, na.rm = TRUE),
        y_min = min(value, na.rm = TRUE),
        y_max = max(value, na.rm = TRUE),
        .groups = "drop"
      ) |>
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
      ungroup(),
    aes(x = x_pos, y = y_pos, label = label, color = source),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    show.legend = FALSE
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey45", linewidth = 0.4) +
  geom_point(alpha = 0.4, size = 1.1) +
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
  "figures/validation/validation_adm0_eawag_reported_scatter_common.png",
  p_scatter_common,
  width = 170,
  height = 100,
  units = "mm",
  dpi = 300
)
