## =============================================================================
## local-cluster-models.R
## Spectral clustering + per-cluster horizons evaluation
##
## Author:        Sam Leuthold
## Contact:       sam.leuthold@colostate.edu
## Last modified: 2026-02-09
##
## Description:
##   Clusters samples via Euclidean distance in PCA space (Ward.D2), then runs
##   full horizons evaluate() per cluster with SVM/Cubist/ElasticNet grid.
##   Tests whether spectrally-defined local models outperform global PyC
##   prediction.
##
## Input:
##   data/processed/horizons_data.rds  (430 rows, 307 with response)
##   data/processed/modeling_response.csv  (project labels for cluster context)
##
## Output:
##   output/local-clusters-k3/cluster_{1,2,3}/result.rds
##   output/local-clusters-k3/all_cluster_results.rds
##
## Changelog:
##   2026-02-09  Restyled from original local-cluster-models.R
##   2026-02-06  Forced k=3 (silhouette 0.275) to balance cluster heterogeneity
##               vs. sample count. k=2 silhouette better (0.41) but largest
##               cluster (n=199) mixes France/FORCE/Jens — too heterogeneous.
##   2026-02-06  Excluded FORCE from cluster 3 (n=16, SD=0.006, range
##               0.072-0.096). Response flat; held out for external validation.
## =============================================================================


pacman::p_load(dplyr, stringr, ggplot2, cluster, horizons, rules, readr,
               future, tibble)


## =============================================================================
## Step 1: Load data and extract spectral matrix
## =============================================================================

hd_full      <- readRDS("data/processed/horizons_data.rds")
has_response <- !is.na(hd_full$data$analysis$pyc_abs)

hd_full$data$role_map |>
  filter(role == "predictor") |>
  pull(variable) -> wn_cols

analysis <- hd_full$data$analysis |> filter(has_response)
spec     <- analysis |> select(all_of(wn_cols)) |> as.matrix()
pyc      <- analysis$pyc_abs


## =============================================================================
## Step 2: PCA and variance-based PC selection
## =============================================================================

pca_fit <- prcomp(spec, center = TRUE, scale. = TRUE)
cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)

# 99% variance threshold for clustering — captures mineralogical gradients
# without overfitting to noise in tail PCs

n_pcs_95 <- which(cum_var >= 0.95)[1]
n_pcs_99 <- which(cum_var >= 0.99)[1]
n_pcs    <- n_pcs_99
scores   <- pca_fit$x[, 1:n_pcs]

## =============================================================================
## Step 3: Euclidean distance clustering in PCA space
## =============================================================================

# Euclidean on raw PCA scores lets high-variance PCs (mineralogical gradients)
# dominate clustering. Mahalanobis produced lopsided clusters (282/25)
# in prior testing.

spec_dist <- dist(scores, method = "euclidean")

## =============================================================================
## Step 4: Hierarchical clustering and silhouette analysis
## =============================================================================

hclust_fit <- hclust(spec_dist, method = "ward.D2")

# k=3 chosen over k=2 despite lower silhouette (0.275 vs 0.41) because k=2's
# large cluster (n=199) mixes France/FORCE/Jens. k=3 splits into ~100
# samples/cluster with more homogeneous spectral neighborhoods.

sil_scores <- numeric(5)

for (k in 2:6) {

  clusters          <- cutree(hclust_fit, k = k)
  sil               <- silhouette(clusters, spec_dist)
  sil_scores[k - 1] <- mean(sil[, "sil_width"])

  }

k_opt <- 3

## =============================================================================
## Step 5: Assign clusters and add project labels
## =============================================================================

clusters         <- cutree(hclust_fit, k = k_opt)
analysis$cluster <- clusters

read_csv("data/processed/modeling_response.csv",
         show_col_types = FALSE) -> resp

analysis |>
  left_join(resp |> select(Sample_ID, project),
            by = c("sample_id" = "Sample_ID")) -> analysis


## =============================================================================
## Step 6: Diagnostic plots
## =============================================================================

## PCA biplot colored by cluster ----------------------------------------------

tibble(PC1     = pca_fit$x[, 1],
       PC2     = pca_fit$x[, 2],
       cluster = factor(clusters),
       pyc     = pyc) -> pca_plot_data

ggplot(pca_plot_data, aes(x     = PC1,
                          y     = PC2,
                          color = cluster)) +
  geom_point(alpha = 0.7,
             size  = 2) +
  theme_minimal()

## Silhouette scores across k values ------------------------------------------

tibble(k          = 2:6,
       silhouette = sil_scores) -> sil_data

ggplot(sil_data,
       aes(x = k,
           y = silhouette)) +
  geom_line() +
  geom_point(size = 3) +
  geom_vline(xintercept = k_opt,
             linetype   = "dashed",
             color      = "red") +
  labs(x = "k", y = "Mean Silhouette Width") +
  theme_minimal()

## =============================================================================
## Step 7: Horizons subset helper
## =============================================================================

# horizons has no native subset method. Deep-copy object, filter analysis to
# cluster sample IDs, and reset downstream slots.

subset_horizons <- function(hd, sample_ids) {

  hd_sub <- hd

  hd$data$analysis |>
    filter(sample_id %in% sample_ids) -> hd_sub$data$analysis

  hd_sub$data$n_rows <- nrow(hd_sub$data$analysis)
  hd_sub$config      <- NULL
  hd_sub$validation  <- NULL
  hd_sub$evaluation  <- NULL
  hd_sub$models      <- NULL
  hd_sub$ensemble    <- NULL
  hd_sub$artifacts   <- NULL

  hd_sub

  }

## =============================================================================
## Step 8: Per-cluster horizons pipeline
## =============================================================================

# Config grid per cluster (72 configs):
#   3 models:            svm_rbf, cubist, elastic_net
#   6 preprocessing:     raw, sg, snv, deriv1, snv_deriv1, deriv2
#   2 transforms:        log, none
#   2 feature selection: pca, correlation (boruta too slow for local runs)
#
# Deeper tuning: grid_size=20, bayesian_iter=15.

## Work around horizons bug: sequential path never sets future::plan() ---------

future::plan(future::multisession, workers = 5)

output_dir <- "output/local-clusters-k3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Run models ------------------------------------------------------------------

cluster_results <- list()

for (cl in sort(unique(clusters))) {

  ## Subset to cluster samples ------------------------------------------------

  cl_data <- analysis |> filter(cluster == cl)

  # FORCE excluded from cluster 3 — 16 samples with no PyC variance
  # (SD=0.006, range 0.072-0.096). Held out for external validation.

  if (cl == 3) {

    cl_data <- cl_data |> filter(project != "FORCE")

  }

  cl_ids <- cl_data |> pull(sample_id)
  hd_cl <- subset_horizons(hd_full, cl_ids)
  n_resp <- sum(!is.na(hd_cl$data$analysis$pyc_abs))

  # Need enough samples for train/test + 5-fold CV

  if (n_resp < 50) {

    cluster_results[[cl]] <- tibble(cluster = cl, n = n_resp,
                                    status  = "skipped")
    next
  }

  ## Configure ----------------------------------------------------------------

  configure(hd_cl,
            outcome           = "pyc_abs",
            models            = c("svm_rbf", "cubist", "elastic_net"),
            transformations   = c("log", "none"),
            preprocessing     = c("raw", "sg", "snv", "deriv1",
                                  "snv_deriv1", "deriv2"),
            feature_selection = c("pca", "correlation"),
            cv_folds          = 5,
            grid_size         = 20,
            bayesian_iter     = 15) -> hd_cl

  ## Validate -----------------------------------------------------------------

  validate(hd_cl,
           remove_outliers    = FALSE,
           spectral_method    = "mahalanobis",
           spectral_threshold = 0.975,
           response_method    = "iqr",
           response_threshold = 1.5) -> hd_cl

  ## Evaluate -----------------------------------------------------------------

  cl_output <- file.path(output_dir, sprintf("cluster_%d", cl))
  dir.create(cl_output, recursive = TRUE, showWarnings = FALSE)

  evaluate(hd_cl,
           metric          = "rpd",
           prune           = TRUE,
           prune_threshold = 1.0,
           workers         = 5,
           output_dir      = cl_output,
           seed            = 307L,
           verbose         = TRUE) -> result_cl

  ## Extract results ----------------------------------------------------------

  eval_res <- result_cl$evaluation$results
  best_id  <- result_cl$evaluation$best_config

  cluster_results[[cl]] <- eval_res |>
    mutate(cluster = cl, n_samples = n_resp)

  saveRDS(result_cl, file.path(cl_output, "result.rds"))
}


## =============================================================================
## Step 9: Collate and compare to global
## =============================================================================

## Kill the parallel backend ---------------------------------------------------

future::plan("sequential")

all_local <- bind_rows(cluster_results)

## Best model per cluster -----------------------------------------------------

all_local |>
  filter(status == "success") |>
  group_by(cluster) |>
  slice_max(rpd, n = 1) |>
  select(cluster, n_samples, config_id, rpd, rsq, rmse, ccc) |>
  mutate(across(where(is.numeric) & !c(cluster, n_samples),
                \(x) round(x, 3)))

# Global best (from overnight): SVM-RBF sg+log+PCA, RPD=1.70, R²=0.662

saveRDS(all_local, file.path(output_dir, "all_cluster_results.rds"))
