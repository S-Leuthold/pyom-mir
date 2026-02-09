#' Local Run — PyOM v2 Cleaned Pipeline
#'
#' @description
#' Comprehensive local run on cleaned dataset (no litter, no BCInt, no spectral
#' outliers, no ALS). 126 configs covering 7 models × 3 transforms × 6
#' preprocessing × 1 feature selection. Designed to run unattended (~2-3 hours).
#'
#' @keywords internal
#' @name test-local


## =============================================================================
## Section 1: Setup
## =============================================================================


library(magrittr)
library(recipes)
library(horizons)
library(plsmod)
library(mixOmics)
library(pls)
library(rules)

## Work around horizons bug: sequential path never sets future::plan()
future::plan(future::multisession, workers = 10)

hd <- readRDS("data/processed/horizons_data.rds")

cat("===============================================================================\n")
cat("PyOM-MIR v2 — Full Local Run (no ALS, cleaned data)\n")
cat("===============================================================================\n\n")

cat(sprintf("Data: %d samples, %d predictors\n",
            hd$data$n_rows, hd$data$n_predictors))
cat(sprintf("Response: %d with pyc_abs\n",
            sum(!is.na(hd$data$analysis$pyc_abs))))


## =============================================================================
## Section 2: Configure
## =============================================================================
##
## 126 configs: 7 models × 3 transforms × 6 preprocessing × 1 feature selection
##
## Models: cubist (best so far), rf, xgboost, mars, elastic_net, svm_rbf, lightgbm
## Transforms: log (best so far), none, sqrt
## Preprocessing: full sweep — raw, sg, snv, deriv1, snv_deriv1, snv_deriv2
## Feature selection: pca (fast, works well enough locally)
##
## PLS excluded — engine issues through horizons/parsnip. PLS works fine
## standalone (RPD 1.49 LOO) but fails through the tidymodels stack.
## MLP excluded — slow and unlikely to beat tree/kernel models on 307 samples.


cat("\n=== Configuring ===\n")

hd |>
  horizons::configure(
    outcome           = "pyc_abs",
    models            = c("cubist", "rf", "xgboost", "mars",
                          "elastic_net", "svm_rbf", "plsr"),
    transformations   = c("log", "none", "sqrt"),
    preprocessing     = c("raw", "sg", "snv", "deriv1", "deriv2", "snv_deriv1"),
    feature_selection = c("pca", "correlation", "boruta", "none"),
    cv_folds          = 10,
    grid_size         = 20,
    bayesian_iter     = 15
  ) ->
hd

cat(sprintf("Config grid: %d configs\n", hd$config$n_configs))


## =============================================================================
## Section 3: Validate
## =============================================================================


cat("\n=== Validating ===\n")

hd |>
  horizons::validate(
    remove_outliers    = FALSE,
    spectral_method    = "mahalanobis",
    spectral_threshold = 0.975,
    response_method    = "iqr",
    response_threshold = 1.5
  ) ->
hd


## =============================================================================
## Section 4: Evaluate
## =============================================================================

#Small, local models?


cat("\n=== Evaluating ===\n")
cat("Running 126 configs — expect ~2-3 hours.\n\n")

hd |>
  horizons::evaluate(
    metric          = "rpd",
    prune           = TRUE,
    prune_threshold = 1.0,
    workers         = 10,
    seed            = 307L,
    verbose         = TRUE,
    output_dir      = "./output/pyom-initial-examination-020526"
  ) ->
result


## =============================================================================
## Section 5: Results
## =============================================================================


cat("\n===============================================================================\n")
cat("Results\n")
cat("===============================================================================\n\n")

eval_results <- result$evaluation$results

cat(sprintf("Configs: %d total\n", nrow(eval_results)))
cat(sprintf("  Succeeded: %d\n", sum(eval_results$status == "success")))
cat(sprintf("  Failed:    %d\n", sum(eval_results$status == "failed")))
cat(sprintf("  Pruned:    %d\n", sum(eval_results$status == "pruned")))

cat(sprintf("\nBest config: %s\n", result$evaluation$best_config))

successful <- eval_results |>
  dplyr::filter(status == "success") |>
  dplyr::arrange(dplyr::desc(rpd))

if (nrow(successful) > 0) {
  cat("\n=== Top 15 Configs ===\n")
  successful |>
    dplyr::select(config_id, rpd, rsq, rmse) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) |>
    print(n = 15)
}

if (sum(eval_results$status == "failed") > 0) {
  cat("\n=== Failures ===\n")
  eval_results |>
    dplyr::filter(status == "failed") |>
    dplyr::select(config_id, error_message) |>
    print(n = Inf)
}

## Save results for later analysis
saveRDS(result, "data/processed/local_eval_result.rds")

cat("\n===============================================================================\n")
cat(sprintf("Runtime: %.1f min\n", result$evaluation$runtime_secs / 60))
cat("Results saved to: data/processed/local_eval_result.rds\n")
cat("===============================================================================\n")
