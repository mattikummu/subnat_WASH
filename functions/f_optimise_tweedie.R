# Functions: f_optimise_tweedie.R
# ─────────────────────────────────────────────────────────────────────────────
# optimise_tweedie_power()
#
# Tests `reg:squarederror` against `reg:tweedie` with powers 1.1, 1.3, 1.5,
# 1.7, 1.9 using raw xgboost::xgb.train (faster than tidymodels for this
# purpose). Selects the objective that gives the lowest RMSE on a held-out
# validation slice. Also computes and reports R² for each model.
#
# Returns a list:
#   $best_power      – numeric; NA when squarederror wins
#   $use_tweedie     – logical; TRUE when any Tweedie power wins
#   $max_delta_step  – numeric; 5 for Tweedie (gradient stabilisation),
#                      0 (xgboost default) for squarederror
#   $results         – tibble; full comparison table with columns:
#                      objective, power, rsq, rmse
#   $indicator       – character; echo of the indicator_name argument
#
# Parameters:
#   train_data     – data frame with columns: response, all predictor_cols
#   indicator      – character; e.g. "Surface.water"
#   predictor_cols – character vector of predictor column names
#
# Notes:
#   - Uses a stratified grouped split (by adm_id) so the same admin area
#     does not appear in both train and validation.
#   - NAs in predictors are passed through (XGBoost handles them natively).
#   - Only rows with non-NA response are used.
# ─────────────────────────────────────────────────────────────────────────────

optimise_tweedie_power <- function(train_data, indicator, predictor_cols) {

  # ── 1. Prepare data ──────────────────────────────────────────────────────
  dat <- train_data |>
    dplyr::select(adm_id, dplyr::all_of(predictor_cols),
                  response = dplyr::all_of(indicator)) |>
    dplyr::filter(!is.na(response))

  # Stratified 80/20 split by adm_id so no adm_id straddles train/val
  set.seed(42L)
  adm_ids   <- unique(dat$adm_id)
  train_ids <- sample(adm_ids, size = floor(0.8 * length(adm_ids)))

  dat_train <- dat |> dplyr::filter( adm_id %in% train_ids)
  dat_val   <- dat |> dplyr::filter(!adm_id %in% train_ids)

  if (nrow(dat_val) < 10) {
    message(sprintf("  [optimise_tweedie] %s: too few validation rows (%d) — using squarederror",
                    indicator, nrow(dat_val)))
    return(list(
      best_power     = NA_real_,
      use_tweedie    = FALSE,
      max_delta_step = 0,
      results        = tibble::tibble(objective = "squarederror", power = NA_real_,
                                      rsq = NA_real_, rmse = NA_real_),
      indicator      = indicator
    ))
  }

  feat_cols <- predictor_cols
  dtrain <- xgboost::xgb.DMatrix(
    data  = as.matrix(dat_train[, feat_cols]),
    label = dat_train$response
  )
  dval <- xgboost::xgb.DMatrix(
    data  = as.matrix(dat_val[, feat_cols]),
    label = dat_val$response
  )

  # ── 2. Fixed common hyperparameters ────────────────────────────────────
  base_params <- list(
    eta               = 0.05,
    max_depth         = 4L,
    subsample         = 0.8,
    colsample_bytree  = 0.8,
    min_child_weight  = 5,
    lambda            = 1,
    nthread           = 1L,
    verbosity         = 0L
  )

  fit_and_eval <- function(params, nrounds = 500, early_stop = 30) {
    fit <- xgboost::xgb.train(
      params    = params,
      data      = dtrain,
      nrounds   = nrounds,
      watchlist = list(val = dval),
      early_stopping_rounds = early_stop,
      verbose   = 0
    )
    preds  <- predict(fit, dval)
    truth  <- xgboost::getinfo(dval, "label")
    rmse   <- sqrt(mean((preds - truth)^2, na.rm = TRUE))
    ss_res <- sum((preds - truth)^2, na.rm = TRUE)
    ss_tot <- sum((truth - mean(truth, na.rm = TRUE))^2, na.rm = TRUE)
    rsq    <- 1 - (ss_res / ss_tot)
    list(rmse = rmse, rsq = rsq)
  }

  # ── 3. Candidate objectives ─────────────────────────────────────────────
  candidates <- list(
    list(objective = "squarederror", power = NA_real_,
         params = c(base_params, list(objective = "reg:squarederror")))
  )

  tweedie_powers <- c(1.1, 1.3, 1.5, 1.7, 1.9)
  for (p in tweedie_powers) {
    candidates <- c(candidates, list(
      list(objective = "tweedie", power = p,
           params = c(base_params,
                      list(objective          = "reg:tweedie",
                           tweedie_variance_power = p,
                           max_delta_step     = 5)))
    ))
  }

  # ── 4. Evaluate each candidate ─────────────────────────────────────────
  results <- purrr::map_dfr(candidates, function(cand) {
    evals <- tryCatch(
      fit_and_eval(cand$params),
      error = function(e) {
        message(sprintf("  [optimise_tweedie] %s power=%.1f error: %s",
                        indicator, cand$power %||% NA, conditionMessage(e)))
        list(rmse = NA_real_, rsq = NA_real_)
      }
    )
    tibble::tibble(
      objective = cand$objective,
      power     = cand$power,
      rsq       = evals$rsq,
      rmse      = evals$rmse
    )
  })

  message(sprintf("  [optimise_tweedie] %s results:", indicator))
  for (i in seq_len(nrow(results))) {
    r <- results[i, ]
    lbl <- if (is.na(r$power)) r$objective else sprintf("tweedie(p=%.1f)", r$power)
    message(sprintf("    %-22s RMSE = %.4f  R² = %.4f", lbl, r$rmse, r$rsq))
  }

  # ── 5. Select best ──────────────────────────────────────────────────────
  best_idx <- which.min(results$rmse)
  best_row <- results[best_idx, ]

  use_tweedie    <- best_row$objective == "tweedie"
  best_power     <- if (use_tweedie) best_row$power else NA_real_
  max_delta_step <- if (use_tweedie) 5 else 0

  squarederror_rmse <- results |>
    dplyr::filter(objective == "squarederror") |>
    dplyr::pull(rmse)
  delta_rmse <- squarederror_rmse - best_row$rmse   # positive = Tweedie better

  if (use_tweedie) {
    message(sprintf("  [optimise_tweedie] %s → Tweedie(p=%.1f) wins (ΔRMSE = %+.4f vs squarederror)",
                    indicator, best_power, delta_rmse))
  } else {
    message(sprintf("  [optimise_tweedie] %s → squarederror wins", indicator))
  }

  list(
    best_power     = best_power,
    use_tweedie    = use_tweedie,
    max_delta_step = max_delta_step,
    results        = results,
    indicator      = indicator
  )
}
