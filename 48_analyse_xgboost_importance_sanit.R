# 48_analyse_xgboost_importance_sanit.R
# ─────────────────────────────────────────────────────────────────────────────
# Variable importance analysis for the trained XGBoost SANITATION models.
#
# Uses the single tuned model per indicator (xgboost_sanit_model_*.rds).
# For each indicator:
#   - Extracts gain-based importance from the underlying xgboost engine
#   - Plots individual importance bar charts
#   - Produces a combined faceted plot across all indicators
#   - Saves a summary table
#   - SHAP analysis
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidymodels)
  library(xgboost)
  library(tidyverse)
  library(vip)
  library(shapviz)
})

dir.create("figures_sanit", showWarnings = FALSE)
dir.create("results_sanit", showWarnings = FALSE)

sanit_indicators <- c("s_bas", "s_lim", "s_unimp", "s_od")

predictor_cols <- c("travel_time", "poverty", "msch", "gnic",
                    "gini", "urban", "gdi", "corruption")

predictor_labels <- c(
  travel_time = "Travel time to city",
  poverty     = "Poverty rate",
  msch        = "Mean schooling",
  gnic        = "GNI per capita",
  gini        = "Gini coefficient",
  urban       = "Urbanisation",
  gdi         = "Gender development index",
  corruption  = "Corruption index"
)

# ── 1. Extract importance from each model ────────────────────────────────────

importance_list <- map(sanit_indicators, function(ind) {
  path <- sprintf("results_sanit/xgboost_sanit_model_%s.rds", ind)
  if (!file.exists(path)) {
    message(sprintf("  ⚠ Model not found for %s, skipping.", ind))
    return(NULL)
  }

  wf <- readRDS(path)
  xgb_model <- extract_fit_engine(wf)

  imp <- xgb.importance(model = xgb_model) |>
    as_tibble() |>
    rename(variable = Feature, gain = Gain, cover = Cover, frequency = Frequency) |>
    mutate(
      indicator = ind,
      label     = predictor_labels[variable]
    ) |>
    arrange(desc(gain))

  message(sprintf("  %s: top predictor = %s (gain = %.3f)",
                  ind, imp$variable[1], imp$gain[1]))
  imp
})

importance_all <- bind_rows(importance_list)

# ── 2. Save summary table ────────────────────────────────────────────────────

write_csv(importance_all, "results_sanit/xgboost_importance_sanit.csv")
message(sprintf("Saved: results_sanit/xgboost_importance_sanit.csv (%d rows)", nrow(importance_all)))

# ── 3. Individual importance plots ───────────────────────────────────────────

walk(sanit_indicators, function(ind) {
  imp <- importance_all |> filter(indicator == ind)
  if (nrow(imp) == 0) return(invisible(NULL))

  p <- imp |>
    mutate(label = fct_reorder(label, gain)) |>
    ggplot(aes(x = gain, y = label)) +
    geom_col(fill = "#2166ac", width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", gain)), hjust = -0.1, size = 3.2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title    = sprintf("Variable importance: %s", ind),
      subtitle = "Gain = fractional contribution to model improvement",
      x        = "Gain",
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank())

  ggsave(sprintf("figures_sanit/importance_sanit_%s.png", ind), p, width = 7, height = 4.5)
})

# ── 4. Combined faceted plot ─────────────────────────────────────────────────

p_combined <- importance_all |>
  mutate(
    label     = fct_reorder(label, gain, .fun = mean),
    indicator = factor(indicator, levels = sanit_indicators)
  ) |>
  ggplot(aes(x = gain, y = label, fill = indicator)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  facet_wrap(~indicator, scales = "free_x", ncol = 2) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Variable importance by sanitation indicator",
    x     = "Gain",
    y     = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text         = element_text(face = "bold")
  )

ggsave("figures_sanit/importance_sanit_combined.png", p_combined, width = 12, height = 6)
message("Saved: figures_sanit/importance_sanit_combined.png")

# ── 5. Rank summary table ────────────────────────────────────────────────────

gain_wide <- importance_all |>
  select(variable, label, indicator, gain) |>
  pivot_wider(names_from = indicator, values_from = gain, values_fill = 0) |>
  arrange(desc(rowMeans(across(all_of(sanit_indicators)))))

message("\n── Gain by variable and indicator ──────────────────────────────────")
print(gain_wide |> select(label, all_of(sanit_indicators)), n = Inf)

# ── 6. SHAP analysis ─────────────────────────────────────────────────────────
message("\n── SHAP value analysis ──────────────────────────────────────────────")

all_data <- read_csv("results_sanit/xgb_training_data_sanit.csv", show_col_types = FALSE)

shap_importance_list <- map(sanit_indicators, function(ind) {
  path <- sprintf("results_sanit/xgboost_sanit_model_%s.rds", ind)
  if (!file.exists(path)) return(NULL)

  wf        <- readRDS(path)
  xgb_model <- extract_fit_engine(wf)

  fitted_recipe <- extract_recipe(wf)
  model_data <- all_data |>
    select(all_of(predictor_cols)) |>
    drop_na()

  baked <- bake(fitted_recipe, new_data = model_data,
                all_predictors(), composition = "matrix")

  sv <- shapviz(xgb_model, X_pred = baked, X = model_data)

  # SHAP importance bar chart
  p_bar <- sv_importance(sv, kind = "bar", show_numbers = TRUE) +
    scale_y_discrete(
      labels = function(x) coalesce(predictor_labels[x], x)
    ) +
    labs(
      title    = sprintf("SHAP importance: %s", ind),
      subtitle = "Mean |SHAP value| across all observations",
      x        = "Mean |SHAP|",
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank())

  ggsave(sprintf("figures_sanit/shap_bar_sanit_%s.png", ind), p_bar, width = 7, height = 4.5)

  # SHAP beeswarm plot
  p_bee <- sv_importance(sv, kind = "beeswarm", size = 0.8, alpha = 0.2) +
    scale_y_discrete(
      labels = function(x) coalesce(predictor_labels[x], x)
    ) +
    labs(
      title    = sprintf("SHAP summary: %s", ind),
      subtitle = "Each point = one observation; colour = feature value",
      x        = "SHAP value",
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank())

  ggsave(sprintf("figures_sanit/shap_beeswarm_sanit_%s.png", ind), p_bee, width = 7, height = 4.5)

  # Dependence plots for top 3 predictors
  top3 <- names(head(sv_importance(sv, kind = "no"), 3))

  walk(top3, function(feat) {
    p_dep <- sv_dependence(sv, v = feat, color_var = "auto") +
      labs(
        title = sprintf("SHAP dependence: %s — %s",
                        coalesce(predictor_labels[feat], feat), ind),
        x     = coalesce(predictor_labels[feat], feat),
        y     = "SHAP value"
      ) +
      theme_minimal(base_size = 12)

    ggsave(sprintf("figures_sanit/shap_dep_sanit_%s_%s.png", ind, feat), p_dep,
           width = 6, height = 5)
  })

  message(sprintf("  %s: SHAP done (top predictor = %s)", ind, top3[1]))

  shap_vals <- sv_importance(sv, kind = "no")
  tibble(
    variable      = names(shap_vals),
    mean_abs_shap = as.numeric(shap_vals)
  ) |>
    mutate(
      indicator = ind,
      label     = coalesce(predictor_labels[variable], variable)
    )
})

shap_importance_all <- bind_rows(shap_importance_list)
write_csv(shap_importance_all, "results_sanit/xgboost_shap_importance_sanit.csv")
message(sprintf("Saved: results_sanit/xgboost_shap_importance_sanit.csv (%d rows)",
                nrow(shap_importance_all)))

# Combined SHAP faceted plot
p_shap_combined <- shap_importance_all |>
  mutate(
    label     = fct_reorder(label, mean_abs_shap, .fun = mean),
    indicator = factor(indicator, levels = sanit_indicators)
  ) |>
  ggplot(aes(x = mean_abs_shap, y = label, fill = indicator)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  facet_wrap(~indicator, scales = "free_x", ncol = 2) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "SHAP importance by sanitation indicator",
    x     = "Mean |SHAP value|",
    y     = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text         = element_text(face = "bold")
  )

ggsave("figures_sanit/shap_importance_sanit_combined.png", p_shap_combined,
       width = 12, height = 6)
message("Saved: figures_sanit/shap_importance_sanit_combined.png")

message("\nDone.")
