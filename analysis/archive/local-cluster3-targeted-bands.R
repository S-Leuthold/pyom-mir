#' Cluster 3 — Targeted PyC Bands
#'
#' @description
#' Re-run cluster 3 (no FORCE) restricted to PyC-sensitive wavenumber regions
#' instead of full spectrum. Hypothesis: removing mineralogical interference
#' from non-diagnostic bands improves prediction in this difficult cluster.
#'
#' @details
#' Targeted regions (4 bands, ~250 features vs 1699 full-spectrum):
#'   - 3000-3100 cm⁻¹: aromatic C-H stretch
#'   - 1550-1750 cm⁻¹: aromatic C=C stretch + C=O
#'   - 1200-1300 cm⁻¹: aromatic C-O stretch
#'   - 700-900 cm⁻¹:   aromatic C-H out-of-plane bending (condensed rings)

## =============================================================================
## Setup
## =============================================================================

library(dplyr, warn.conflicts = FALSE)
library(stringr)
library(cluster)
library(horizons)
library(rules)

hd_full <- readRDS("data/processed/horizons_data.rds")

## Reproduce cluster assignments from main script
has_response <- !is.na(hd_full$data$analysis$pyc_abs)
wn_cols <- hd_full$data$role_map |>
  filter(role == "predictor") |>
  pull(variable)

analysis <- hd_full$data$analysis |> filter(has_response)
spec     <- analysis |> select(all_of(wn_cols)) |> as.matrix()

pca_fit <- prcomp(spec, center = TRUE, scale. = TRUE)
cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
n_pcs   <- which(cum_var >= 0.99)[1]
scores  <- pca_fit$x[, 1:n_pcs]

spec_dist  <- dist(scores, method = "euclidean")
hclust_fit <- hclust(spec_dist, method = "ward.D2")
clusters   <- cutree(hclust_fit, k = 3)
analysis$cluster <- clusters

## Join project labels for FORCE filter
resp <- readr::read_csv("data/processed/modeling_response.csv",
                        show_col_types = FALSE)
analysis <- analysis |>
  left_join(resp |> select(Sample_ID, project),
            by = c("sample_id" = "Sample_ID"))


## =============================================================================
## Define PyC-sensitive wavenumber regions
## =============================================================================

## Extract numeric wavenumbers from column names
wn_nums <- as.numeric(gsub("^wn_", "", wn_cols))

## PyC-diagnostic bands
band_ranges <- list(
  aromatic_ch_stretch = c(3000, 3100),
  aromatic_cc_co      = c(1550, 1750),
  aromatic_co_stretch = c(1200, 1300),
  aromatic_ch_oop     = c(700, 900)
)

## Select columns falling within any band
in_band <- rep(FALSE, length(wn_nums))
for (band in band_ranges) {
  in_band <- in_band | (wn_nums >= band[1] & wn_nums <= band[2])
}

targeted_cols <- wn_cols[in_band]
cat(sprintf("Targeted bands: %d features (from %d full-spectrum)\n",
            length(targeted_cols), length(wn_cols)))

for (nm in names(band_ranges)) {
  rng <- band_ranges[[nm]]
  n <- sum(wn_nums >= rng[1] & wn_nums <= rng[2])
  cat(sprintf("  %s (%d-%d cm⁻¹): %d features\n", nm, rng[1], rng[2], n))
}


## =============================================================================
## Helper — subset horizons_data to targeted bands
## =============================================================================

subset_horizons_bands <- function(hd, sample_ids, keep_cols) {

  hd_sub <- hd

  ## Filter samples
  hd_sub$data$analysis <- hd$data$analysis |>
    filter(sample_id %in% sample_ids)
  hd_sub$data$n_rows <- nrow(hd_sub$data$analysis)

  ## Drop non-targeted predictor columns from analysis
  drop_cols <- hd$data$role_map |>
    filter(role == "predictor", !variable %in% keep_cols) |>
    pull(variable)
  hd_sub$data$analysis <- hd_sub$data$analysis |>
    select(-all_of(drop_cols))

  ## Update role_map to only include targeted predictors
  hd_sub$data$role_map <- hd$data$role_map |>
    filter(role != "predictor" | variable %in% keep_cols)

  ## Reset downstream slots
  hd_sub$config     <- NULL
  hd_sub$validation <- NULL
  hd_sub$evaluation <- NULL
  hd_sub$models     <- NULL
  hd_sub$ensemble   <- NULL
  hd_sub$artifacts  <- NULL

  hd_sub
}


## =============================================================================
## Run cluster 3 with targeted bands
## =============================================================================

future::plan(future::multisession, workers = 5)

## Cluster 3 without FORCE
cl3_data <- analysis |>
  filter(cluster == 3, project != "FORCE")
cl3_ids <- cl3_data |> pull(sample_id)

hd_cl3 <- subset_horizons_bands(hd_full, cl3_ids, targeted_cols)

n_resp <- sum(!is.na(hd_cl3$data$analysis$pyc_abs))
cat(sprintf("\nCluster 3 (targeted bands, no FORCE): %d samples (%d with response)\n",
            hd_cl3$data$n_rows, n_resp))

## Configure — same grid as main script
hd_cl3 |>
  horizons::configure(
    outcome           = "pyc_abs",
    models            = c("svm_rbf", "cubist", "elastic_net"),
    transformations   = c("log", "none"),
    preprocessing     = c("raw", "sg", "snv", "deriv1", "snv_deriv1", "deriv2"),
    feature_selection = c("pca", "correlation"),
    cv_folds          = 5,
    grid_size         = 20,
    bayesian_iter     = 15
  ) ->
hd_cl3

cat(sprintf("Configs: %d\n", hd_cl3$config$n_configs))

## Validate
hd_cl3 |>
  horizons::validate(
    remove_outliers    = FALSE,
    spectral_method    = "mahalanobis",
    spectral_threshold = 0.975,
    response_method    = "iqr",
    response_threshold = 1.5
  ) ->
hd_cl3

## Evaluate
output_dir <- "output/local-clusters-k3/cluster_3_targeted"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hd_cl3 |>
  horizons::evaluate(
    metric          = "rpd",
    prune           = TRUE,
    prune_threshold = 1.0,
    workers         = 5,
    output_dir      = output_dir,
    seed            = 307L,
    verbose         = TRUE
  ) ->
result_cl3

## Results
eval_res <- result_cl3$evaluation$results
best_id  <- result_cl3$evaluation$best_config

cat(sprintf("\n=== Cluster 3 Targeted Bands Results ===\n"))
cat(sprintf("Success: %d, Failed: %d, Pruned: %d\n",
            sum(eval_res$status == "success"),
            sum(eval_res$status == "failed"),
            sum(eval_res$status == "pruned")))
cat(sprintf("Best: %s (RPD=%.2f, R²=%.3f, RMSE=%.4f)\n",
            best_id,
            eval_res$rpd[eval_res$config_id == best_id],
            eval_res$rsq[eval_res$config_id == best_id],
            eval_res$rmse[eval_res$config_id == best_id]))

cat("\nTop 5 configs:\n")
eval_res |>
  filter(status == "success") |>
  arrange(desc(rpd)) |>
  head(5) |>
  select(config_id, rpd, rsq, rmse) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  print()

saveRDS(result_cl3, file.path(output_dir, "result.rds"))
cat(sprintf("\nSaved: %s\n", output_dir))
