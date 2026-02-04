#!/usr/bin/env Rscript
#===============================================================================
# Script:  run-pyom-analysis.R
# Purpose: Run horizons v1 evaluate() pipeline on Sybil HPC
# Author:  Sam Leuthold
# Created: 2026-02-04
# Project: PyOM-MIR SPA-C prediction from FTIR spectra
#===============================================================================
#
# Input:  data/horizons_data.rds  (404 samples, built by prepare-horizons-data.R)
# Output: results/run_TIMESTAMP/  (evaluation results, checkpoints, session info)
#
# The horizons_data object already contains averaged spectra + response.
# This script runs: configure() -> validate() -> evaluate()
#
#===============================================================================


## =============================================================================
## Section 1: Thread Control
## =============================================================================
##
## CRITICAL for Sybil shared environment: pin all libraries to 1 thread.
## Parallelism is managed by future::multisession, not by BLAS/OpenMP.


rm(list = ls())
gc(verbose = FALSE, full = TRUE, reset = TRUE)

Sys.setenv(
  OMP_NUM_THREADS        = 1,
  OPENBLAS_NUM_THREADS   = 1,
  MKL_NUM_THREADS        = 1,
  NUMEXPR_NUM_THREADS    = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  BLAS_NUM_THREADS       = 1,
  LAPACK_NUM_THREADS     = 1,
  OMP_THREAD_LIMIT       = 1
)

if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(1)
}


## =============================================================================
## Section 2: Load Packages and Configuration
## =============================================================================


suppressPackageStartupMessages({
  library(yaml)
  library(future)
  library(horizons)
})

config <- yaml::read_yaml("config-pyom-analysis.yml")

start_time <- Sys.time()
timestamp  <- format(start_time, "%Y%m%d_%H%M%S")
run_id     <- paste0("run_", timestamp)

cat("\n")
cat("===============================================================================\n")
cat("PyOM-MIR SPA-C Analysis\n")
cat("===============================================================================\n")
cat("Start time: ", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Hostname:   ", Sys.info()["nodename"], "\n")
cat("R version:  ", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat("Horizons:   ", as.character(packageVersion("horizons")), "\n")
cat("Run ID:     ", run_id, "\n")
cat("Workers:    ", config$compute$n_workers, "\n")
cat("===============================================================================\n\n")


## =============================================================================
## Section 3: Setup Parallel Backend
## =============================================================================


n_workers <- config$compute$n_workers

## Override parallelly's localhost limit
options(mc.cores = n_workers)
Sys.setenv(R_PARALLELLY_MAXWORKERS_LOCALHOST = as.character(n_workers))

future::plan(future::multisession, workers = n_workers)

cat(sprintf("Parallel backend: %d multisession workers\n", n_workers))
cat(sprintf("System cores detected: %d\n", parallel::detectCores()))
cat("\n")


## =============================================================================
## Section 4: Load Data
## =============================================================================


data_path <- config$paths$data

if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path)
}

hd <- readRDS(data_path)

cat(sprintf("Loaded horizons_data: %d samples, %d predictors\n",
            hd$data$n_rows, hd$data$n_predictors))

## Verify response variable exists
outcome <- config$response
role_map <- hd$data$role_map

if (!outcome %in% role_map$variable) {
  stop("Response variable '", outcome, "' not found in role_map")
}

n_with_response <- sum(!is.na(hd$data$analysis[[outcome]]))
cat(sprintf("Samples with response (%s): %d\n", outcome, n_with_response))
cat("\n")


## =============================================================================
## Section 5: Configure
## =============================================================================


cat("===============================================================================\n")
cat("Configuring model grid\n")
cat("===============================================================================\n\n")

hd |>
  horizons::configure(
    outcome         = outcome,
    models          = config$models$types,
    transformations = config$models$transformations,
    preprocessing   = config$models$preprocessing,
    feature_selection = config$models$feature_selection,
    cv_folds        = config$tuning$cv_folds,
    grid_size       = config$tuning$grid_size,
    bayesian_iter   = config$tuning$bayesian_iter
  ) ->
hd

cat(sprintf("\nConfiguration grid: %d configs\n", hd$config$n_configs))
cat(sprintf("  Models:            %s\n", paste(config$models$types, collapse = ", ")))
cat(sprintf("  Transformations:   %s\n", paste(config$models$transformations, collapse = ", ")))
cat(sprintf("  Preprocessing:     %s\n", paste(config$models$preprocessing, collapse = ", ")))
cat(sprintf("  Feature selection: %s\n", paste(config$models$feature_selection, collapse = ", ")))
cat(sprintf("  Tuning: %d-fold CV, grid=%d, bayes=%d\n",
            config$tuning$cv_folds, config$tuning$grid_size, config$tuning$bayesian_iter))
cat("\n")


## =============================================================================
## Section 6: Validate
## =============================================================================


cat("===============================================================================\n")
cat("Validating data\n")
cat("===============================================================================\n\n")

hd |>
  horizons::validate(
    remove_outliers    = config$validation$remove_outliers,
    spectral_method    = config$validation$spectral_method,
    spectral_threshold = config$validation$spectral_threshold,
    response_method    = config$validation$response_method,
    response_threshold = config$validation$response_threshold
  ) ->
hd

cat(sprintf("Validation passed: %s\n", hd$validation$passed))
cat("\n")


## =============================================================================
## Section 7: Evaluate
## =============================================================================


cat("===============================================================================\n")
cat("Starting evaluation\n")
cat("===============================================================================\n\n")

## Create output directory for checkpoints
results_dir <- file.path(config$paths$results, run_id)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

hd |>
  horizons::evaluate(
    metric          = config$evaluation$metric,
    prune           = config$evaluation$prune,
    prune_threshold = config$evaluation$prune_threshold,
    allow_par       = config$compute$allow_par,
    output_dir      = results_dir,
    seed            = config$evaluation$seed,
    verbose         = config$compute$verbose
  ) ->
result

cat("\n")
cat("===============================================================================\n")
cat("Evaluation complete\n")
cat("===============================================================================\n\n")


## =============================================================================
## Section 8: Save Results
## =============================================================================


## Save full result object
result_file <- file.path(results_dir, "pyom_eval_result.rds")
saveRDS(result, result_file)
cat(sprintf("Results saved: %s\n", result_file))

## Save session info
session_file <- file.path(results_dir, "session_info.txt")
sink(session_file)
cat("PyOM-MIR Analysis Session Info\n")
cat("Run ID:", run_id, "\n")
cat("Date:", format(Sys.time()), "\n\n")
sessionInfo()
sink()

## Save config copy for reproducibility
file.copy("config-pyom-analysis.yml",
          file.path(results_dir, "config-pyom-analysis.yml"))


## =============================================================================
## Section 9: Summary
## =============================================================================


end_time   <- Sys.time()
total_time <- difftime(end_time, start_time, units = "hours")

cat("\n")
cat("===============================================================================\n")
cat("PyOM-MIR Analysis Complete\n")
cat("===============================================================================\n")
cat(sprintf("Total runtime:  %.2f hours\n", as.numeric(total_time)))
cat(sprintf("Results:        %s\n", results_dir))
cat(sprintf("Best config:    %s\n", result$evaluation$best_config))
cat("===============================================================================\n")

## Shut down workers
future::plan(future::sequential)
