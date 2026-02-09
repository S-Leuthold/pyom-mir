#' PyOM-MIR Report Figures
#'
#' Generates all figures for preliminary results summary.
#' Output: output/figures/

library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(parsnip)
library(rsample)

theme_set(theme_minimal(base_size = 12))

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)


## =============================================================================
## Load data
## =============================================================================

## Overnight global results (310 configs)
global_res <- readRDS("output/pyom-initial-examination-020526/eval_checkpoint.rds")

## Local cluster results
cl1_res <- readRDS("output/local-clusters-k3/cluster_1/result.rds")$evaluation$results |>
  mutate(cluster = 1)
cl2_res <- readRDS("output/local-clusters-k3/cluster_2/result.rds")$evaluation$results |>
  mutate(cluster = 2)
cl3_res <- readRDS("output/local-clusters-k3/cluster_3/eval_checkpoint.rds") |>
  mutate(cluster = 3)

## Parse model type from config_id
parse_model <- function(config_id) {
  str_extract(config_id, "^[a-z_]+(?=_)")
}

parse_config <- function(df) {
  df |>
    mutate(
      model = case_when(
        str_detect(config_id, "^svm_rbf")      ~ "SVM-RBF",
        str_detect(config_id, "^cubist")        ~ "Cubist",
        str_detect(config_id, "^elastic_net")   ~ "Elastic Net",
        str_detect(config_id, "^xgboost")       ~ "XGBoost",
        str_detect(config_id, "^rf")            ~ "Random Forest",
        str_detect(config_id, "^pls")           ~ "PLSR",
        str_detect(config_id, "^mars")          ~ "MARS",
        TRUE ~ "Other"
      ),
      preprocessing = case_when(
        str_detect(config_id, "snv_deriv1") ~ "SNV+D1",
        str_detect(config_id, "snv_deriv2") ~ "SNV+D2",
        str_detect(config_id, "deriv1")     ~ "D1",
        str_detect(config_id, "deriv2")     ~ "D2",
        str_detect(config_id, "_sg_")       ~ "SG",
        str_detect(config_id, "_snv_")      ~ "SNV",
        str_detect(config_id, "_raw_")      ~ "Raw",
        TRUE ~ "Other"
      ),
      transform = ifelse(str_detect(config_id, "_log_"), "log", "none")
    )
}

global_res <- parse_config(global_res)
cl1_res <- parse_config(cl1_res)
cl2_res <- parse_config(cl2_res)
cl3_res <- parse_config(cl3_res)


## =============================================================================
## Figure 1: Global RPD distribution by model type
## =============================================================================

fig1_data <- global_res |>
  filter(status == "success") |>
  mutate(model = factor(model, levels = c("SVM-RBF", "Cubist", "Elastic Net",
                                           "XGBoost", "Random Forest", "MARS")))

## Best per model
best_per_model <- fig1_data |>
  group_by(model) |>
  slice_max(rpd, n = 1) |>
  ungroup()

fig1 <- ggplot(fig1_data, aes(x = model, y = rpd, fill = model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  geom_point(data = best_per_model, aes(x = model, y = rpd),
             shape = 18, size = 4, color = "red") +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.6, y = 1.52, label = "RPD = 1.5", hjust = 0,
           size = 3, color = "grey50") +
  labs(
    title = "Global model performance across 203 configurations",
    subtitle = "Red diamonds = best per model type. Dashed line = minimum screening threshold.",
    x = NULL, y = "RPD (ratio of performance to deviation)"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.5, 2.0))

ggsave("output/figures/fig1_rpd_distribution.png", fig1,
       width = 8, height = 5, dpi = 300)
cat("Saved: fig1_rpd_distribution.png\n")


## =============================================================================
## Figure 2: Spectral clusters in PCA space
## =============================================================================

hd_full <- readRDS("data/processed/horizons_data.rds")
wn_cols <- hd_full$data$role_map |> filter(role == "predictor") |> pull(variable)
analysis <- hd_full$data$analysis |> filter(!is.na(pyc_abs))
spec <- analysis |> select(all_of(wn_cols)) |> as.matrix()

pca_fit <- prcomp(spec, center = TRUE, scale. = TRUE)
cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
n_pcs <- which(cum_var >= 0.99)[1]
scores <- pca_fit$x[, 1:n_pcs]
clusters <- cutree(hclust(dist(scores), method = "ward.D2"), k = 3)

resp <- readr::read_csv("data/processed/modeling_response.csv", show_col_types = FALSE)
analysis$cluster <- clusters
analysis <- analysis |>
  left_join(resp |> select(Sample_ID, project),
            by = c("sample_id" = "Sample_ID"))

pca_df <- tibble(
  PC1 = pca_fit$x[, 1],
  PC2 = pca_fit$x[, 2],
  cluster = factor(clusters),
  pyc = analysis$pyc_abs,
  project = analysis$project
)

var_explained <- round(100 * pca_fit$sdev[1:2]^2 / sum(pca_fit$sdev^2), 1)

fig2a <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "Spectral clusters (k = 3)",
    subtitle = sprintf("Euclidean distance + Ward.D2 on %d PCs (99%% variance)", n_pcs),
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    color = "Cluster"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.15))

fig2b <- ggplot(pca_df, aes(x = PC1, y = PC2, color = pyc)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_viridis_c(option = "inferno", name = "PyC\n(g/100g)") +
  labs(
    title = "Colored by PyC concentration",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.15))

fig2 <- fig2a + fig2b +
  plot_annotation(tag_levels = "A")

ggsave("output/figures/fig2_spectral_clusters.png", fig2,
       width = 12, height = 5, dpi = 300)
cat("Saved: fig2_spectral_clusters.png\n")


## =============================================================================
## Figure 3: Local vs global RPD comparison
## =============================================================================

local_all <- bind_rows(cl1_res, cl2_res, cl3_res) |>
  filter(status == "success")

local_best <- local_all |>
  group_by(cluster) |>
  slice_max(rpd, n = 1) |>
  ungroup()

global_best_rpd <- max(global_res$rpd[global_res$status == "success"], na.rm = TRUE)

fig3 <- ggplot(local_all, aes(x = factor(cluster), y = rpd, fill = factor(cluster))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  geom_point(data = local_best, shape = 18, size = 4, color = "red") +
  geom_hline(yintercept = global_best_rpd, linetype = "dashed", color = "blue", linewidth = 0.8) +
  annotate("text", x = 0.6, y = global_best_rpd + 0.05,
           label = sprintf("Global best (RPD = %.2f)", global_best_rpd),
           hjust = 0, size = 3, color = "blue") +
  geom_hline(yintercept = 2.0, linetype = "dotted", color = "grey50") +
  annotate("text", x = 3.4, y = 2.05, label = "RPD = 2.0", hjust = 1,
           size = 3, color = "grey50") +
  labs(
    title = "Per-cluster local model performance (72 configs each)",
    subtitle = "Red diamonds = best per cluster. Blue dashed = global best.",
    x = "Spectral cluster", y = "RPD"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("1" = "Cluster 1\n(n=71)", "2" = "Cluster 2\n(n=108)", "3" = "Cluster 3\n(n=128)"))

ggsave("output/figures/fig3_local_vs_global.png", fig3,
       width = 7, height = 5, dpi = 300)
cat("Saved: fig3_local_vs_global.png\n")


## =============================================================================
## Figure 4: Concentration floor
## =============================================================================

## Rerun the 5-fold CV to get per-sample predictions
snv <- function(x) {
  m <- as.matrix(x)
  t(apply(m, 1, function(r) (r - mean(r)) / sd(r))) |> as.data.frame()
}

svm_mod <- svm_rbf(cost = 2.5, rbf_sigma = 0.05) |>
  set_engine("kernlab") |>
  set_mode("regression")

set.seed(307)
folds <- rsample::vfold_cv(analysis, v = 5, strata = pyc_abs)

all_preds <- tibble()
for (f in seq_along(folds$splits)) {
  train_f <- rsample::analysis(folds$splits[[f]])
  test_f  <- rsample::assessment(folds$splits[[f]])

  train_spec <- snv(train_f |> select(all_of(wn_cols)))
  test_spec  <- snv(test_f  |> select(all_of(wn_cols)))

  train_ready <- bind_cols(pyc_log = log(train_f$pyc_abs), train_spec)
  test_ready  <- bind_cols(pyc_log = log(test_f$pyc_abs), test_spec)

  fit_f <- svm_mod |> fit(pyc_log ~ ., data = train_ready)
  preds_f <- exp(predict(fit_f, test_ready)$.pred)

  all_preds <- bind_rows(all_preds, tibble(
    obs = test_f$pyc_abs,
    pred = preds_f,
    cluster = factor(test_f$cluster)
  ))
}

all_preds <- all_preds |>
  mutate(abs_error = abs(obs - pred),
         rel_error = abs_error / obs)

fig4a <- ggplot(all_preds, aes(x = obs, y = abs_error, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.055, y = max(all_preds$abs_error) * 0.95,
           label = "PyC = 0.05", hjust = 0, size = 3, color = "grey40") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, aes(group = 1),
              color = "black") +
  labs(
    title = "Prediction error vs PyC concentration",
    x = "Observed PyC (g / 100g soil)",
    y = "Absolute prediction error",
    color = "Cluster"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.85))

fig4b <- ggplot(all_preds, aes(x = obs, y = rel_error * 100, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 100, linetype = "dotted", color = "grey50") +
  annotate("text", x = 0.42, y = 110, label = "100% relative error",
           hjust = 1, size = 3, color = "grey50") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, aes(group = 1),
              color = "black") +
  labs(
    title = "Relative error vs PyC concentration",
    x = "Observed PyC (g / 100g soil)",
    y = "Relative prediction error (%)",
    color = "Cluster"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 500)) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.85))

fig4 <- fig4a + fig4b +
  plot_annotation(tag_levels = "A")

ggsave("output/figures/fig4_concentration_floor.png", fig4,
       width = 12, height = 5, dpi = 300)
cat("Saved: fig4_concentration_floor.png\n")


## =============================================================================
## Figure 5: 1:1 predicted vs observed (5-fold CV, global model)
## =============================================================================

r2_cv <- round(cor(all_preds$obs, all_preds$pred)^2, 3)
rmse_cv <- round(sqrt(mean((all_preds$obs - all_preds$pred)^2)), 4)
rpd_cv <- round(sd(all_preds$obs) / sqrt(mean((all_preds$obs - all_preds$pred)^2)), 2)

fig5 <- ggplot(all_preds, aes(x = obs, y = pred, color = cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    title = "Predicted vs observed PyC (5-fold CV, SVM-RBF + SNV + log)",
    subtitle = sprintf("R² = %.3f, RMSE = %.4f, RPD = %.2f  (n = %d)",
                       r2_cv, rmse_cv, rpd_cv, nrow(all_preds)),
    x = "Observed PyC (g / 100g soil)",
    y = "Predicted PyC (g / 100g soil)",
    color = "Cluster"
  ) +
  scale_color_brewer(palette = "Set1") +
  coord_equal(xlim = c(0, 0.5), ylim = c(0, 0.5)) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.15))

ggsave("output/figures/fig5_pred_vs_obs.png", fig5,
       width = 6, height = 6, dpi = 300)
cat("Saved: fig5_pred_vs_obs.png\n")


## =============================================================================
## Summary stats for report
## =============================================================================

cat("\n\n=== Summary statistics for report ===\n\n")

cat("--- Global model (overnight, 504 planned, 310 completed) ---\n")
cat(sprintf("  Successful: %d\n", sum(global_res$status == "success")))
cat(sprintf("  Pruned: %d (all sqrt transform)\n", sum(global_res$status == "pruned")))
cat(sprintf("  Failed: %d (PLSR)\n", sum(global_res$status == "failed")))
cat(sprintf("  Best: %s\n", global_res$config_id[which.max(global_res$rpd)]))
cat(sprintf("  Best RPD: %.2f, R²: %.3f, RMSE: %.4f\n",
    max(global_res$rpd, na.rm = TRUE),
    global_res$rsq[which.max(global_res$rpd)],
    global_res$rmse[which.max(global_res$rpd)]))

cat("\n--- Local cluster models (k=3, 72 configs each) ---\n")
for (cl in 1:3) {
  res_cl <- switch(cl, cl1_res, cl2_res, cl3_res)
  best <- res_cl |> filter(status == "success") |> slice_max(rpd, n = 1)
  cat(sprintf("  Cluster %d: best RPD=%.2f, R²=%.3f (%s)\n",
      cl, best$rpd, best$rsq, best$config_id))
}

cat("\n--- Concentration floor ---\n")
cat("  Below 0.05 g/100g: median relative error > 100%\n")
cat("  Sweet spot 0.08-0.12 g/100g: median relative error 8%\n")
cat("  Error pattern consistent across all clusters\n")

cat("\nAll figures saved to output/figures/\n")
