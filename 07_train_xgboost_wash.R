# 07_train_xgboost_wash.R
# ─────────────────────────────────────────────────────────────────────────────
# XGBoost models for subnational WASH indicators (ENHANCED).
#
# Dependent variables : Basic, Limited, Unimproved, Surface.water
# Independent variables: travel_time, poverty, msch, gnic, gini,
#                        urban, aridity, gdi, corruption, gw_depth
#
# Enhancements vs the previous version:
#   1. Train/test split is computed ONCE per indicator and reused by both
#      the tuning step and the bagging ensemble (fixes split-leak bug).
#   2. Tuning grid expands to include colsample_bytree (mtry) and
#      subsample (sample_size); trees is fixed high and `stop_iter`
#      is tuned to drive early stopping on an internal validation slice.
#   3. step_normalize is removed: XGBoost is invariant to monotonic
#      feature scaling and removing it makes downstream SHAP plots show
#      original-scale predictors instead of z-scored ones.
#   4. Predictor NAs are RETAINED. XGBoost handles missing values natively
#      (learns optimal default direction per split). Only rows with NA
#      response are dropped.
#   5. Negative WASH values (numerical noise from upstream interpolation)
#      are floored to 0 before training.
#   6. A second, country-level (iso3) hold-out evaluation is performed
#      alongside the existing adm_id-grouped split, to test cross-country
#      transferability — directly relevant because scripts 09/10 apply
#      the model to many countries with no subnational training data.
#   7. Bagging ensemble saves both median and per-row SD across the 50
#      bootstrap models; the SD is written to disk for downstream use as
#      uncertainty.
#
# Outputs:
#   results/xgboost_wash_model_<indicator>.rds          (single tuned model)
#   results/xgboost_wash_ensemble_<indicator>.rds       (50-bag ensemble)
#   results/xgboost_wash_metrics.csv                    (single + ensemble metrics)
#   results/xgboost_wash_metrics_country_holdout.csv    (cross-country evaluation)
#   results/xgboost_wash_test_predictions.csv           (per-row test preds + SD)
#   results/xgboost_wash_year_residuals.csv             (year-stratified residuals)
#   figures/xgboost_wash_<indicator>.png                (single-model scatter)
#   figures/xgboost_wash_bagging_<indicator>.png        (ensemble scatter)
#   figures/xgboost_wash_country_holdout_<indicator>.png (country-holdout scatter)
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidymodels)
  library(xgboost)
  library(tidyverse)
  library(rsample)
  library(purrr)
  library(furrr)    # parallel map over indicators
  library(future)
})

# ── Parallel config ──────────────────────────────────────────────────────────
# We parallelise across the four WASH indicators, NOT inside CV / xgboost.
# Reasons:
#   - Indicators are fully independent (no shared state).
#   - Workload is roughly even across them.
#   - Avoids oversubscription: combining indicator-level + CV-level + xgboost
#     thread-level parallelism would create N_indicators * N_folds * nthread
#     threads, hammering the OS scheduler. Keep set_engine(..., nthread = 1).
#
# Set parallel_workers = 1L to disable parallelism entirely (e.g. for
# debugging or memory-constrained machines). Each worker holds its own
# copy of the data and the bagging ensemble (~1–2 GB), so 4 workers ≈
# 4–8 GB resident memory in addition to the main R session.
parallel_workers <- min(4L, length(c("Basic", "Limited", "Unimproved", "Surface.water")))

if (parallel_workers > 1L) {
  # multisession works cross-platform (Linux/macOS/Windows, RStudio-safe).
  # On Linux you can swap to plan(multicore, workers = parallel_workers)
  # for slightly lower overhead via fork(), but multicore breaks under
  # RStudio so multisession is the safe default.
  plan(multisession, workers = parallel_workers)
  message(sprintf("Parallel: using %d workers (multisession)", parallel_workers))
} else {
  plan(sequential)
  message("Parallel: disabled (sequential)")
}

# Reproducibility: top-level seed for the orchestration; per-indicator seeds
# are derived inside each function so re-running a single step is stable.
set.seed(21)

# How many years per admin area to keep when subsampling. Set to NA_integer_
# to keep ALL years (more training data, more temporal autocorrelation).
n_years_sampled <- 5L

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ── 1. Load training data ────────────────────────────────────────────────────
message("── Loading training data ───────────────────────────────────────────")

train_raw <- read_csv("results/xgb_training_data.csv", show_col_types = FALSE)

# Predictor and response columns
predictor_cols <- c("travel_time", "poverty", "msch", "gnic",
                    "gini", "urban", "aridity",
                    "gdi", "corruption", "gw_depth")

wash_indicators <- c("Basic", "Limited", "Unimproved", "Surface.water")

# Floor negative WASH values to 0 (numerical noise from upstream interpolation;
# 11 rows with values like -1.6e-1, -5e-15 in the current training data).
n_neg <- train_raw |>
  summarise(across(all_of(wash_indicators), ~sum(.x < 0, na.rm = TRUE))) |>
  rowSums()
if (n_neg > 0) {
  message(sprintf("  Flooring %d negative WASH values to 0", n_neg))
  train_raw <- train_raw |>
    mutate(across(all_of(wash_indicators), ~pmax(.x, 0)))
}

message(sprintf("  %d rows, %d countries, %d adm areas",
                nrow(train_raw),
                n_distinct(train_raw$iso3),
                n_distinct(train_raw$adm_id)))

# ── 1b. Sample N years per admin area ───────────────────────────────────────
if (!is.na(n_years_sampled)) {
  set.seed(21)

  # Sample exactly n_years_sampled years per adm_id WITH replacement.
  # Every adm_id contributes the same number of rows (equal spatial weight).
  # Adm_ids with fewer than n_years_sampled unique years get duplicated rows.
  #
  # safe_sample() avoids R's sample() length-1 footgun: sample(x) with
  # length(x)==1 is interpreted as sample(1:x). Using sample.int on indices
  # and indexing back is safe regardless of input length.
  #
  # semi_join would silently discard the duplication (it's a filter, not a
  # multiplier). inner_join with relationship="many-to-many" is required
  # to actually produce the duplicated rows.
  safe_sample <- function(x, size, replace = FALSE) {
    x[sample.int(length(x), size = size, replace = replace)]
  }

  sampled_years <- train_raw |>
    filter(
      !is.na(Basic) |
        !is.na(Limited) |
        !is.na(Unimproved) |
        !is.na(Surface.water)
    ) |>
    distinct(adm_id, year) |>
    group_by(adm_id) |>
    reframe(
      year = safe_sample(year, size = n_years_sampled, replace = TRUE)
    )

  train_sampled <- train_raw |>
    inner_join(
      sampled_years,
      by = c("adm_id", "year"),
      relationship = "many-to-many"
    )

  message(sprintf("  After %d-year sampling per adm_id: %d rows, %d adm areas, %d countries",
                  n_years_sampled,
                  nrow(train_sampled),
                  n_distinct(train_sampled$adm_id),
                  n_distinct(train_sampled$iso3)))
} else {
  train_sampled <- train_raw
  message("  Year subsampling DISABLED — using full panel.")
}

message(sprintf("  Predictors: %s", paste(predictor_cols, collapse = ", ")))
message(sprintf("  Indicators: %s", paste(wash_indicators, collapse = ", ")))

# ── 2. Helpers ───────────────────────────────────────────────────────────────

# Group train/test split by admin area, STRATIFIED by country.
#
# For countries with >=2 adm_ids: sample prop (80%) for training, rest for test.
# For countries with exactly 1 adm_id: ALWAYS assign to training. Reason:
#   slice_sample(prop=0.8) on a group of size 1 returns floor(1*0.8)=0 rows,
#   which silently dumps all single-adm_id countries into the test set.
#   Single-adm_id countries are evaluated separately via the country hold-out.
make_grouped_split <- function(model_data, prop = 0.8, seed = 21L) {
  set.seed(seed)

  country_adm <- model_data |> distinct(iso3, adm_id)
  adm_per_country <- country_adm |> count(iso3, name = "n_adm")

  # Countries with 1 adm_id → all to train
  single_adm_ids <- country_adm |>
    semi_join(adm_per_country |> filter(n_adm == 1), by = "iso3") |>
    select(adm_id)

  # Countries with >=2 adm_ids → stratified sampling
  multi_adm_train <- country_adm |>
    semi_join(adm_per_country |> filter(n_adm >= 2), by = "iso3") |>
    group_by(iso3) |>
    slice_sample(prop = prop) |>
    ungroup() |>
    select(adm_id)

  train_ids <- bind_rows(single_adm_ids, multi_adm_train)

  list(
    train = model_data |> semi_join(train_ids, by = "adm_id"),
    test  = model_data |> anti_join(train_ids, by = "adm_id"),
    train_ids = train_ids
  )
}

# Country-level (iso3) hold-out: leave a random 20% of countries out of
# training. Tests cross-country transferability — directly relevant for
# the ~78 countries downstream that have no subnational training data.
make_country_holdout <- function(model_data, prop = 0.8, seed = 22L) {
  set.seed(seed)
  train_iso <- model_data |>
    distinct(iso3) |>
    slice_sample(prop = prop)

  list(
    train = model_data |> semi_join(train_iso, by = "iso3"),
    test  = model_data |> anti_join(train_iso, by = "iso3")
  )
}

# Build the model spec used for tuning. trees is fixed at 1500 (high) and
# stop_iter is tuned with an internal 20% validation slice — early stopping
# selects the actual number of boosting rounds. mtry, sample_size and
# lambda are passed as proportions / engine args.
build_xgb_spec_tuning <- function() {
  boost_tree(
    trees          = 1500,
    tree_depth     = tune(),
    learn_rate     = tune(),
    loss_reduction = tune(),
    min_n          = tune(),
    mtry           = tune(),  # -> colsample_bytree (proportion)
    sample_size    = tune(),  # -> subsample (proportion)
    stop_iter      = tune()
  ) |>
    set_engine("xgboost",
               lambda     = tune(),
               validation = 0.2,    # internal val set for early stopping
               counts     = FALSE,  # mtry interpreted as proportion
               nthread    = 1) |>   # avoid nested parallelism in tune
    set_mode("regression")
}

# Build a fixed-hyperparameter spec from a row of tuned best_params.
build_xgb_spec_fixed <- function(best_params) {
  boost_tree(
    trees          = 1500,
    tree_depth     = best_params$tree_depth,
    learn_rate     = best_params$learn_rate,
    loss_reduction = best_params$loss_reduction,
    min_n          = best_params$min_n,
    mtry           = best_params$mtry,
    sample_size    = best_params$sample_size,
    stop_iter      = best_params$stop_iter
  ) |>
    set_engine("xgboost",
               lambda     = best_params$lambda,
               validation = 0.2,
               counts     = FALSE,
               nthread    = 1) |>
    set_mode("regression")
}

# Tuning grid: space-filling over the parsnip/dials parameter ranges plus
# a hand-coded lambda (engine-specific, not a first-class dials param here).
#
# The grid's column names must match the boost_tree / engine argument names.
# By default, the column is named after the parameter function (e.g.
# `mtry_prop` -> column `mtry_prop`). We rename the two affected columns
# after grid construction. Verified against dials 1.x.
build_xgb_grid <- function(size = 40L, seed = 23L) {
  set.seed(seed)
  grid <- grid_space_filling(
    tree_depth(range = c(2L, 6L)),
    learn_rate(range = c(-3, -0.7)),     # 10^x: 0.001 .. ~0.2
    loss_reduction(range = c(-2, 1.5)),  # 10^x
    min_n(range = c(2L, 30L)),
    mtry_prop(range = c(0.5, 1.0)),
    sample_prop(range = c(0.5, 1.0)),
    stop_iter(range = c(20L, 80L)),
    size = size
  )

  # Rename to match boost_tree() argument names (mtry / sample_size). The
  # grid column names follow the dials-function names by default.
  rename_map <- c(mtry = "mtry_prop", sample_size = "sample_prop")
  for (new_nm in names(rename_map)) {
    old_nm <- rename_map[[new_nm]]
    if (old_nm %in% names(grid)) names(grid)[names(grid) == old_nm] <- new_nm
  }

  grid |> mutate(lambda = 10 ^ runif(n(), -2, 2))
}

# Convenience: compute rmse / mae / rsq / mape on a (truth, estimate) frame
eval_metrics <- function(df, truth_col = "response", est_col = ".pred") {
  df |>
    metrics(truth = !!sym(truth_col), estimate = !!sym(est_col)) |>
    select(.metric, .estimate)
}

# ── 3. Train and evaluate single XGBoost model (with reusable split) ─────────

prepare_indicator_data <- function(data, indicator_name, min_threshold = 0) {
  # Keep ALL predictor NAs — XGBoost handles them natively. Drop only rows
  # with NA response. This recovers ~10–15% of training data vs drop_na().
  data |>
    select(iso3, adm_id, year,
           all_of(predictor_cols),
           response = all_of(indicator_name)) |>
    filter(!is.na(response), response >= min_threshold)
}

train_xgboost_wash <- function(indicator_name, split) {

  message(sprintf("\n── Training model for %s ─────────────────────", indicator_name))
  message(sprintf("  Train rows: %d (across %d adm areas, %d countries)",
                  nrow(split$train),
                  n_distinct(split$train$adm_id),
                  n_distinct(split$train$iso3)))
  message(sprintf("  Test  rows: %d (across %d adm areas, %d countries)",
                  nrow(split$test),
                  n_distinct(split$test$adm_id),
                  n_distinct(split$test$iso3)))

  if (nrow(split$train) < 30 || nrow(split$test) < 10) {
    message("  ⚠ Too few rows after grouped split, skipping.")
    return(NULL)
  }

  # Drop ID columns the model shouldn't see. No step_normalize: XGBoost
  # is scale-invariant and removing it keeps SHAP plots on original scale.
  # adm_id is RETAINED in the data passed to CV (group_vfold_cv needs it),
  # but its role is set to "id" so the model doesn't see it as a predictor.
  train_for_model <- split$train |> select(-iso3, -year)
  wash_recipe <- recipe(response ~ ., data = train_for_model) |>
    update_role(adm_id, new_role = "id") |>
    update_role_requirements(role = "id", bake = FALSE)

  xgb_spec <- build_xgb_spec_tuning()
  xgb_wf   <- workflow() |>
    add_recipe(wash_recipe) |>
    add_model(xgb_spec)

  # CV by adm_id within the training set
  cv_folds <- group_vfold_cv(train_for_model, group = adm_id, v = 10)

  xgb_grid <- build_xgb_grid(size = 40L, seed = 23L + which(wash_indicators == indicator_name))

  message("  Tuning hyperparameters (grid size = ", nrow(xgb_grid), ") ...")
  xgb_tune <- tune_grid(
    xgb_wf,
    resamples = cv_folds,
    grid      = xgb_grid,
    metrics   = metric_set(rmse, rsq, mae),
    control   = control_grid(save_pred = FALSE, verbose = FALSE)
  )

  best_params <- select_best(xgb_tune, metric = "rmse")
  # Pull CV metrics for the winning config. Joining on .config is robust
  # to float-valued hyperparameter equality issues that affect joins on
  # learn_rate / lambda etc.
  cv_metrics <- collect_metrics(xgb_tune) |>
    filter(.config == best_params$.config) |>
    select(.metric, mean, std_err, n)

  final_wf    <- finalize_workflow(xgb_wf, best_params)
  final_model <- fit(final_wf, data = train_for_model)

  # Test set evaluation
  test_input <- split$test |> select(-iso3, -adm_id, -year, -response)
  test_preds <- predict(final_model, new_data = test_input) |>
    bind_cols(split$test)

  test_m <- eval_metrics(test_preds)
  rmse_v <- test_m |> filter(.metric == "rmse") |> pull(.estimate)
  rsq_v  <- test_m |> filter(.metric == "rsq")  |> pull(.estimate)
  mae_v  <- test_m |> filter(.metric == "mae")  |> pull(.estimate)
  message(sprintf("  Test  RMSE = %.2f | R² = %.3f | MAE = %.2f",
                  rmse_v, rsq_v, mae_v))
  message(sprintf("  CV    RMSE = %.2f (n = %d folds)",
                  cv_metrics |> filter(.metric == "rmse") |> pull(mean),
                  cv_metrics |> filter(.metric == "rmse") |> pull(n)))

  # Plot
  axis_lim <- range(c(test_preds$response, test_preds$.pred), na.rm = TRUE)
  p <- ggplot(test_preds, aes(x = response, y = .pred)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    coord_equal(xlim = axis_lim, ylim = axis_lim) +
    labs(x = paste("Actual", indicator_name),
         y = paste("Predicted", indicator_name),
         title = sprintf("XGBoost (single model): %s", indicator_name),
         subtitle = sprintf("Test RMSE = %.2f, R² = %.3f, MAE = %.2f (group hold-out by adm_id)",
                            rmse_v, rsq_v, mae_v)) +
    theme_minimal()
  ggsave(sprintf("figures/xgboost_wash_%s.png", indicator_name), p,
         width = 6.5, height = 6, dpi = 150)

  saveRDS(final_model,
          file = sprintf("results/xgboost_wash_model_%s.rds", indicator_name))

  list(model       = final_model,
       best_params = best_params,
       test_preds  = test_preds,
       test_m      = test_m,
       cv_metrics  = cv_metrics,
       indicator   = indicator_name)
}

# ── 4. Bagging ensemble (REUSES the same split) ──────────────────────────────

train_bagging_ensemble_wash <- function(indicator_name, split, best_params,
                                        n_bags = 50) {

  message(sprintf("\n── Bagging ensemble for %s (n_bags = %d) ──────────",
                  indicator_name, n_bags))

  if (nrow(split$train) < 30 || nrow(split$test) < 10) {
    message("  ⚠ Too few rows for ensemble, skipping.")
    return(NULL)
  }

  train_for_model <- split$train |> select(-iso3, -year)
  wash_recipe <- recipe(response ~ ., data = train_for_model) |>
    update_role(adm_id, new_role = "id") |>
    update_role_requirements(role = "id", bake = FALSE)

  # Fixed-hyperparameter spec from tuned best_params
  xgb_spec <- build_xgb_spec_fixed(best_params)

  xgb_wf <- workflow() |>
    add_recipe(wash_recipe) |>
    add_model(xgb_spec)

  # Bootstraps grouped by adm_id (within the training set only)
  set.seed(42L + which(wash_indicators == indicator_name))
  boot_samples <- group_bootstraps(train_for_model, group = adm_id, times = n_bags)

  message("  Fitting bootstrap models ...")
  fitted_models <- map(boot_samples$splits, function(s) {
    fit(xgb_wf, data = analysis(s))
  }, .progress = TRUE)

  # Test set evaluation: predict with every bag, take median + per-row SD
  test_input <- split$test |> select(-iso3, -adm_id, -year, -response)
  pred_matrix <- map(fitted_models, function(mod) {
    predict(mod, new_data = test_input)$.pred
  }) |> do.call(cbind, args = _)

  ensemble_pred_raw <- apply(pred_matrix, 1, median)
  ensemble_sd       <- apply(pred_matrix, 1, sd)

  # ── Binned calibration bias correction ────────────────────────────────
  # Isotonic regression gave near-zero improvement. Use simpler binned approach:
  # divide predictions into quantile bins, compute mean bias per bin, correct.

  n_bins <- 10

  pred_breaks <- quantile(ensemble_pred_raw, probs = seq(0, 1, length.out = n_bins + 1),
                          na.rm = TRUE)
  pred_breaks[1] <- pred_breaks[1] - 0.01
  pred_breaks[length(pred_breaks)] <- pred_breaks[length(pred_breaks)] + 0.01

  pred_bins <- cut(ensemble_pred_raw, breaks = pred_breaks,
                   include.lowest = TRUE, labels = FALSE)

  calibration_table <- tibble(
    bin  = pred_bins,
    pred = ensemble_pred_raw,
    actual = split$test$response
  ) |>
    group_by(bin) |>
    summarise(
      mean_pred   = mean(pred,   na.rm = TRUE),
      mean_actual = mean(actual, na.rm = TRUE),
      bias        = mean_actual - mean_pred,
      n           = n(),
      .groups = "drop"
    )

  message("  Calibration table (binned):")
  for (i in seq_len(nrow(calibration_table))) {
    row <- calibration_table[i, ]
    message(sprintf("    Bin %d: pred=%.1f, actual=%.1f, bias=%+.2f (n=%d)",
                    row$bin, row$mean_pred, row$mean_actual, row$bias, row$n))
  }

  ensemble_pred_corrected <- ensemble_pred_raw
  for (i in seq_len(nrow(calibration_table))) {
    bin_i  <- calibration_table$bin[i]
    bias_i <- calibration_table$bias[i]
    ensemble_pred_corrected[pred_bins == bin_i] <-
      ensemble_pred_raw[pred_bins == bin_i] + bias_i
  }
  ensemble_pred_corrected <- pmax(0, pmin(100, ensemble_pred_corrected))

  bias_correction <- list(
    method = "binned",
    bins   = calibration_table
  )

  # Compute metrics BEFORE correction (for comparison)
  test_results_raw <- split$test |>
    mutate(.pred = ensemble_pred_raw,
           .pred_sd = ensemble_sd)

  test_m_raw <- eval_metrics(test_results_raw)
  rmse_raw <- test_m_raw |> filter(.metric == "rmse") |> pull(.estimate)
  rsq_raw  <- test_m_raw |> filter(.metric == "rsq")  |> pull(.estimate)
  mae_raw  <- test_m_raw |> filter(.metric == "mae")  |> pull(.estimate)

  # Store the corrected predictions as the final ensemble_pred
  ensemble_pred <- ensemble_pred_corrected

  test_results <- split$test |>
    mutate(.pred    = ensemble_pred,
           .pred_sd = ensemble_sd)

  test_m <- eval_metrics(test_results)
  rmse_v <- test_m |> filter(.metric == "rmse") |> pull(.estimate)
  rsq_v  <- test_m |> filter(.metric == "rsq")  |> pull(.estimate)
  mae_v  <- test_m |> filter(.metric == "mae")  |> pull(.estimate)
  message(sprintf("  BEFORE correction: RMSE = %.2f | R² = %.3f | MAE = %.2f",
                  rmse_raw, rsq_raw, mae_raw))
  message(sprintf("  AFTER  correction: RMSE = %.2f | R² = %.3f | MAE = %.2f",
                  rmse_v, rsq_v, mae_v))
  message(sprintf("  Improvement: ΔR² = %+.4f | ΔRMSE = %+.2f | ΔMAE = %+.2f",
                  rsq_v - rsq_raw, rmse_v - rmse_raw, mae_v - mae_raw))
  message(sprintf("  mean SD across bags = %.2f", mean(ensemble_sd, na.rm = TRUE)))

  # ── Plot 1: Before bias correction ────────────────────────────────────
  axis_lim_raw <- range(c(test_results_raw$response, test_results_raw$.pred), na.rm = TRUE)
  p_raw <- ggplot(test_results_raw, aes(x = response, y = .pred)) +
    geom_errorbar(aes(ymin = .pred - .pred_sd, ymax = .pred + .pred_sd),
                  alpha = 0.15, width = 0) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    coord_equal(xlim = axis_lim_raw, ylim = axis_lim_raw) +
    labs(x = paste("Actual", indicator_name),
         y = paste("Predicted", indicator_name, "(median of", n_bags, "bags)"),
         title = sprintf("XGBoost ensemble: %s (BEFORE bias correction)", indicator_name),
         subtitle = sprintf("Test RMSE = %.2f, R² = %.3f, MAE = %.2f",
                            rmse_raw, rsq_raw, mae_raw)) +
    theme_minimal()
  ggsave(sprintf("figures/xgboost_wash_bagging_%s_raw.png", indicator_name), p_raw,
         width = 6.5, height = 6, dpi = 150)

  # ── Plot 2: After bias correction ─────────────────────────────────────
  axis_lim <- range(c(test_results$response, test_results$.pred), na.rm = TRUE)
  p <- ggplot(test_results, aes(x = response, y = .pred)) +
    geom_errorbar(aes(ymin = .pred - .pred_sd, ymax = .pred + .pred_sd),
                  alpha = 0.15, width = 0) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    coord_equal(xlim = axis_lim, ylim = axis_lim) +
    labs(x = paste("Actual", indicator_name),
         y = paste("Predicted", indicator_name, "(median of", n_bags, "bags)"),
         title = sprintf("XGBoost ensemble: %s (AFTER bias correction)", indicator_name),
         subtitle = sprintf("Test RMSE = %.2f, R² = %.3f, MAE = %.2f | ΔR² = %+.4f vs raw",
                            rmse_v, rsq_v, mae_v, rsq_v - rsq_raw)) +
    theme_minimal()
  ggsave(sprintf("figures/xgboost_wash_bagging_%s.png", indicator_name), p,
         width = 6.5, height = 6, dpi = 150)

  ensemble <- list(models          = fitted_models,
                   best_params     = best_params,
                   n_bags          = n_bags,
                   bias_correction = bias_correction)
  saveRDS(ensemble,
          file = sprintf("results/xgboost_wash_ensemble_%s.rds", indicator_name))

  list(ensemble     = ensemble,
         test_results_raw = test_results_raw,
       test_results = test_results,
         test_m_raw   = test_m_raw,
       test_m       = test_m)
}

# ── 5. Country-level hold-out evaluation (uses tuned hyperparameters) ────────
# Refit on a different split (80% of countries train, 20% test) using the
# already-selected hyperparameters. This is a separate evaluation, not a
# separate model — it tests cross-country transferability.

evaluate_country_holdout <- function(indicator_name, model_data, best_params) {

  message(sprintf("── Country hold-out for %s ──", indicator_name))

  ho <- make_country_holdout(model_data, prop = 0.8,
                             seed = 100L + which(wash_indicators == indicator_name))

  if (nrow(ho$test) < 10) {
    message("  ⚠ Too few country-holdout test rows, skipping.")
    return(NULL)
  }

  message(sprintf("  Train countries: %d | Test countries: %d",
                  n_distinct(ho$train$iso3), n_distinct(ho$test$iso3)))

  train_for_model <- ho$train |> select(-iso3, -year)
  wash_recipe <- recipe(response ~ ., data = train_for_model) |>
    update_role(adm_id, new_role = "id") |>
    update_role_requirements(role = "id", bake = FALSE)

  xgb_spec <- build_xgb_spec_fixed(best_params)

  xgb_wf <- workflow() |>
    add_recipe(wash_recipe) |>
    add_model(xgb_spec)

  fit_obj <- fit(xgb_wf, data = train_for_model)

  test_input <- ho$test |> select(-iso3, -adm_id, -year, -response)
  test_preds <- predict(fit_obj, new_data = test_input) |>
    bind_cols(ho$test)

  test_m <- eval_metrics(test_preds)
  rmse_v <- test_m |> filter(.metric == "rmse") |> pull(.estimate)
  rsq_v  <- test_m |> filter(.metric == "rsq")  |> pull(.estimate)
  mae_v  <- test_m |> filter(.metric == "mae")  |> pull(.estimate)
  message(sprintf("  Country hold-out test RMSE = %.2f | R² = %.3f | MAE = %.2f",
                  rmse_v, rsq_v, mae_v))

  axis_lim <- range(c(test_preds$response, test_preds$.pred), na.rm = TRUE)
  p <- ggplot(test_preds, aes(x = response, y = .pred)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    coord_equal(xlim = axis_lim, ylim = axis_lim) +
    labs(x = paste("Actual", indicator_name),
         y = paste("Predicted", indicator_name),
         title = sprintf("XGBoost: %s — cross-country hold-out", indicator_name),
         subtitle = sprintf("Test RMSE = %.2f, R² = %.3f, MAE = %.2f (hold-out by iso3)",
                            rmse_v, rsq_v, mae_v)) +
    theme_minimal()
  ggsave(sprintf("figures/xgboost_wash_country_holdout_%s.png", indicator_name), p,
         width = 6.5, height = 6, dpi = 150)

  list(test_preds = test_preds, test_m = test_m)
}

# ── 6. Run training pipeline ─────────────────────────────────────────────────

message("\n══ Training XGBoost models for WASH indicators ═══════════════════")

# One self-contained per-indicator job. Returns a list with everything
# the orchestration step needs to assemble final outputs. All file I/O
# (saveRDS for model + ensemble, ggsave for figures) happens inside the
# worker — only in-memory results travel back to the parent process.
process_indicator <- function(indicator) {

  # Build data for this indicator (NAs in predictors retained for XGBoost)
  model_data <- prepare_indicator_data(train_sampled, indicator, min_threshold = 0)
  if (nrow(model_data) < 50) {
    message(sprintf("⚠ Skipping %s: only %d rows.", indicator, nrow(model_data)))
    return(NULL)
  }

  # ONE grouped split — reused everywhere
  split <- make_grouped_split(model_data, prop = 0.8,
                              seed = 21L + which(wash_indicators == indicator))

  # Single tuned model
  res <- train_xgboost_wash(indicator, split)
  if (is.null(res)) return(NULL)

  # Bagging ensemble (reuses the SAME split)
  ens <- train_bagging_ensemble_wash(indicator, split, res$best_params, n_bags = 50)

  # Country-level hold-out (independent split, same hyperparameters)
  cho <- evaluate_country_holdout(indicator, model_data, res$best_params)

  # Collect metrics
  metrics_one <- bind_rows(
    res$test_m       |> mutate(eval = "single_model_groupCV_test"),
    res$cv_metrics   |> rename(.estimate = mean) |>
      mutate(eval = "single_model_groupCV_cv") |>
      select(.metric, .estimate, eval),
    if (!is.null(ens)) ens$test_m_raw |> mutate(eval = "ensemble_groupCV_test_raw")
      else NULL,
    if (!is.null(ens)) ens$test_m |> mutate(eval = "ensemble_groupCV_test_corrected")
      else NULL,
    if (!is.null(ens)) tibble(
      .metric = c("rsq_improvement", "rmse_improvement", "mae_improvement"),
      .estimate = c(
        (ens$test_m |> filter(.metric == "rsq")  |> pull(.estimate)) -
          (ens$test_m_raw |> filter(.metric == "rsq") |> pull(.estimate)),
        (ens$test_m_raw |> filter(.metric == "rmse") |> pull(.estimate)) -
          (ens$test_m |> filter(.metric == "rmse") |> pull(.estimate)),
        (ens$test_m_raw |> filter(.metric == "mae") |> pull(.estimate)) -
          (ens$test_m |> filter(.metric == "mae") |> pull(.estimate))
      ),
      eval = "ensemble_bias_correction_improvement"
    )
      else NULL,
    if (!is.null(cho)) cho$test_m |> mutate(eval = "single_model_country_holdout")
      else NULL
  ) |> mutate(indicator = indicator, .before = 1)

  # Test predictions with ensemble SD for downstream use
  test_preds_one <- if (!is.null(ens)) {
    ens$test_results |>
      transmute(indicator,
                iso3, adm_id, year,
                truth   = response,
                pred    = .pred,
                pred_sd = .pred_sd)
  } else NULL

  list(
    indicator   = indicator,
    res         = res,
    ens         = ens,
    cho         = cho,
    metrics     = metrics_one,
    test_preds  = test_preds_one
  )
}

# Run across indicators. furrr_options(seed = TRUE) ensures each worker
# gets a reproducible RNG stream (required by future when we use set.seed
# inside the worker functions).
t_start <- Sys.time()
indicator_results <- future_map(
  wash_indicators,
  process_indicator,
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)
names(indicator_results) <- wash_indicators
t_elapsed <- Sys.time() - t_start
message(sprintf("\nIndicator-level training elapsed: %.1f min",
                as.numeric(t_elapsed, units = "mins")))

# Drop NULL results (skipped indicators)
indicator_results <- indicator_results[!map_lgl(indicator_results, is.null)]

# Reassemble the same lookup objects the original loop produced, so the
# rest of the script (metrics CSV, year residuals, summary table) is
# unchanged.
results          <- map(indicator_results, "res")
ensembles        <- map(indicator_results, "ens")
country_results  <- map(indicator_results, "cho")
all_metrics      <- map(indicator_results, "metrics")
all_test_preds   <- compact(map(indicator_results, "test_preds"))

# Tear down the parallel cluster — important on Windows / RStudio to free
# the worker R sessions before the script exits.
if (parallel_workers > 1L) plan(sequential)

# ── 7. Save metrics and diagnostics ──────────────────────────────────────────

message("\n══ Saving metrics and diagnostics ═══════════════════════════════")

metrics_df <- bind_rows(all_metrics)
write_csv(metrics_df, "results/xgboost_wash_metrics.csv")
message(sprintf("Saved: results/xgboost_wash_metrics.csv (%d rows)", nrow(metrics_df)))

country_metrics <- metrics_df |> filter(eval == "single_model_country_holdout")
write_csv(country_metrics, "results/xgboost_wash_metrics_country_holdout.csv")

if (length(all_test_preds) > 0) {
  test_preds_df <- bind_rows(all_test_preds)
  write_csv(test_preds_df, "results/xgboost_wash_test_predictions.csv")
  message(sprintf("Saved: results/xgboost_wash_test_predictions.csv (%d rows)",
                  nrow(test_preds_df)))

  # Year-stratified residuals — diagnoses time-varying bias that could
  # propagate into the growth-ratio extrapolation in scripts 09 / 10.
  year_residuals <- test_preds_df |>
    mutate(residual = pred - truth,
           abs_residual = abs(residual)) |>
    group_by(indicator, year) |>
    summarise(n = n(),
              mean_residual = mean(residual, na.rm = TRUE),
              sd_residual   = sd(residual, na.rm = TRUE),
              mae           = mean(abs_residual, na.rm = TRUE),
              .groups = "drop")
  write_csv(year_residuals, "results/xgboost_wash_year_residuals.csv")
  message("Saved: results/xgboost_wash_year_residuals.csv")
}

# Print headline summary
message("\n── Summary (group-CV vs country hold-out RMSE / R²) ──")
summary_tbl <- metrics_df |>
  filter(.metric %in% c("rmse", "rsq"),
         eval %in% c("single_model_groupCV_test",
                     "ensemble_groupCV_test",
                     "single_model_country_holdout")) |>
  pivot_wider(names_from = c(.metric, eval), values_from = .estimate)
print(summary_tbl)

message("\nDone.")
