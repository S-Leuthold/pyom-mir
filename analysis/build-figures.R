## =============================================================================
## build-figures.R
## Generate all figures for PyOM-MIR preliminary results
##
## Author:        Sam Leuthold
## Contact:       sam.leuthold@colostate.edu
## Last modified: 2026-02-09
##
## Description:
##   Generates report figures from global and local-cluster model results.
##   Steps 1-6 are fast (load saved results, plot). Steps 7-9 refit models
##   with LOO-CV to produce per-sample predictions — these take hours.
##
## Input:
##   data/processed/horizons_data.rds
##   data/processed/modeling_response.csv
##   output/pyom-initial-examination-020526/eval_checkpoint.rds  (global)
##   output/local-clusters-k3/cluster_{1,2,3}/  (local results)
##
## Output:
##   output/figures/fig1_rpd_distribution.png
##   output/figures/fig2_spectral_clusters.png
##   output/figures/fig3_local_vs_global.png
##   output/figures/fig4_concentration_floor.png
##   output/figures/fig5_pred_vs_obs_global.png
##   output/figures/fig6_pred_vs_obs_hybrid.png
##
## Changelog:
##   2026-02-09  Merged from build-report-figures.R + build-hybrid-figure.R
##   2026-02-06  Initial versions (separate scripts)
## =============================================================================


pacman::p_load(dplyr, tidyr, stringr, ggplot2, patchwork, parsnip, rsample,
               rules, readr, prospectr)

theme_set(theme_minimal(base_size = 12))
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)


## =============================================================================
## Step 1: Load data and reproduce clusters
## =============================================================================

hd_full <- readRDS("data/processed/horizons_data.rds")

hd_full$data$role_map |>
  filter(role == "predictor") |>
  pull(variable) -> wn_cols

analysis <- hd_full$data$analysis |> filter(!is.na(pyc_abs))
spec     <- analysis |> select(all_of(wn_cols)) |> as.matrix()

pca_fit <- prcomp(spec, center = TRUE, scale. = TRUE)
cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
n_pcs   <- which(cum_var >= 0.99)[1]
scores  <- pca_fit$x[, 1:n_pcs]

clusters <- cutree(hclust(dist(scores), method = "ward.D2"), k = 3)
analysis$cluster <- clusters

read_csv("data/processed/modeling_response.csv",
         show_col_types = FALSE) -> resp

analysis <- analysis |>
  left_join(resp |> select(Sample_ID, project),
            by = c("sample_id" = "Sample_ID"))


## =============================================================================
## Step 2: Load evaluation results
## =============================================================================

global_res <- readRDS("output/pyom-initial-examination-020526/eval_checkpoint.rds")

cl1_res <- readRDS("output/local-clusters-k3/cluster_1/result.rds")$evaluation$results |>
  mutate(cluster = 1)
cl2_res <- readRDS("output/local-clusters-k3/cluster_2/result.rds")$evaluation$results |>
  mutate(cluster = 2)
cl3_res <- readRDS("output/local-clusters-k3/cluster_3/eval_checkpoint.rds") |>
  mutate(cluster = 3)

## Parse model/preprocessing/transform from config_id -------------------------

# Can't split on underscores because some model names have them (svm_rbf,
# elastic_net). Use pattern matching instead.
parse_config <- function(df) {
  df |>
    mutate(model = case_when(str_detect(config_id, "^svm_rbf")    ~ "SVM-RBF",
                             str_detect(config_id, "^cubist")      ~ "Cubist",
                             str_detect(config_id, "^elastic_net") ~ "Elastic Net",
                             str_detect(config_id, "^xgboost")     ~ "XGBoost",
                             str_detect(config_id, "^rf")          ~ "Random Forest",
                             str_detect(config_id, "^pls")         ~ "PLSR",
                             str_detect(config_id, "^mars")        ~ "MARS",
                             TRUE ~ "Other"),
           preprocessing = case_when(str_detect(config_id, "snv_deriv1") ~ "SNV+D1",
                                     str_detect(config_id, "snv_deriv2") ~ "SNV+D2",
                                     str_detect(config_id, "deriv1")     ~ "D1",
                                     str_detect(config_id, "deriv2")     ~ "D2",
                                     str_detect(config_id, "_sg_")       ~ "SG",
                                     str_detect(config_id, "_snv_")      ~ "SNV",
                                     str_detect(config_id, "_raw_")      ~ "Raw",
                                     TRUE ~ "Other"),
           transform = ifelse(str_detect(config_id, "_log_"), "log", "none"))
}

global_res <- parse_config(global_res)
cl1_res    <- parse_config(cl1_res)
cl2_res    <- parse_config(cl2_res)
cl3_res    <- parse_config(cl3_res)


## =============================================================================
## Step 3: Fig 1 — Global RPD distribution by model type
## =============================================================================

global_res |>
  filter(status == "success") |>
  mutate(model = factor(model, levels = c("SVM-RBF", "Cubist", "Elastic Net",
                                           "XGBoost", "Random Forest",
                                           "MARS"))) -> fig1_data

fig1_data |>
  group_by(model) |>
  slice_max(rpd, n = 1) |>
  ungroup() -> best_per_model

ggplot(fig1_data, aes(x = model, y = rpd, fill = model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  geom_point(data = best_per_model, aes(x = model, y = rpd),
             shape = 18, size = 4, color = "red") +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.6, y = 1.52, label = "RPD = 1.5", hjust = 0,
           size = 3, color = "grey50") +
  labs(title = "Global model performance across 203 configurations",
       subtitle = "Red diamonds = best per model type. Dashed line = minimum screening threshold.",
       x = NULL, y = "RPD (ratio of performance to deviation)") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.5, 2.0)) -> fig1

ggsave("output/figures/fig1_rpd_distribution.png", fig1,
       width = 8, height = 5, dpi = 300)


## =============================================================================
## Step 4: Fig 2 — Spectral clusters in PCA space
## =============================================================================

var_explained <- round(100 * pca_fit$sdev[1:2]^2 / sum(pca_fit$sdev^2), 1)

tibble(PC1     = pca_fit$x[, 1],
       PC2     = pca_fit$x[, 2],
       cluster = factor(clusters),
       pyc     = analysis$pyc_abs,
       project = analysis$project) -> pca_df

ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "Spectral clusters (k = 3)",
       subtitle = sprintf("Euclidean distance + Ward.D2 on %d PCs (99%% variance)", n_pcs),
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       color = "Cluster") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.15)) -> fig2a

ggplot(pca_df, aes(x = PC1, y = PC2, color = pyc)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_viridis_c(option = "inferno", name = "PyC\n(g/100g)") +
  labs(title = "Colored by PyC concentration",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.15)) -> fig2b

fig2a + fig2b + plot_annotation(tag_levels = "A") -> fig2

ggsave("output/figures/fig2_spectral_clusters.png", fig2,
       width = 12, height = 5, dpi = 300)


## =============================================================================
## Step 5: Fig 3 — Local vs global RPD comparison
## =============================================================================

bind_rows(cl1_res, cl2_res, cl3_res) |>
  filter(status == "success") -> local_all

local_all |>
  group_by(cluster) |>
  slice_max(rpd, n = 1) |>
  ungroup() -> local_best

global_best_rpd <- max(global_res$rpd[global_res$status == "success"],
                       na.rm = TRUE)

ggplot(local_all, aes(x = factor(cluster), y = rpd,
                      fill = factor(cluster))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  geom_point(data = local_best, shape = 18, size = 4, color = "red") +
  geom_hline(yintercept = global_best_rpd, linetype = "dashed",
             color = "blue", linewidth = 0.8) +
  annotate("text", x = 0.6, y = global_best_rpd + 0.05,
           label = sprintf("Global best (RPD = %.2f)", global_best_rpd),
           hjust = 0, size = 3, color = "blue") +
  geom_hline(yintercept = 2.0, linetype = "dotted", color = "grey50") +
  annotate("text", x = 3.4, y = 2.05, label = "RPD = 2.0", hjust = 1,
           size = 3, color = "grey50") +
  labs(title = "Per-cluster local model performance (72 configs each)",
       subtitle = "Red diamonds = best per cluster. Blue dashed = global best.",
       x = "Spectral cluster", y = "RPD") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("1" = "Cluster 1\n(n=71)",
                               "2" = "Cluster 2\n(n=108)",
                               "3" = "Cluster 3\n(n=128)")) -> fig3

ggsave("output/figures/fig3_local_vs_global.png", fig3,
       width = 7, height = 5, dpi = 300)


## =============================================================================
## Step 6: Fig 4 + Fig 5 — Concentration floor and global pred vs obs
## =============================================================================

## 5-fold CV refit to get per-sample predictions for the global best model
## (SVM-RBF + SNV + log). Takes a few minutes.

snv <- function(x) {
  m <- as.matrix(x)
  t(apply(m, 1, function(r) (r - mean(r)) / sd(r))) |> as.data.frame()
}

svm_mod <- svm_rbf(cost = 2.5, rbf_sigma = 0.05) |>
  set_engine("kernlab") |>
  set_mode("regression")

set.seed(307)
folds <- vfold_cv(analysis, v = 5, strata = pyc_abs)

all_preds <- tibble()
for (f in seq_along(folds$splits)) {
  train_f <- analysis(folds$splits[[f]])
  test_f  <- assessment(folds$splits[[f]])

  train_spec <- snv(train_f |> select(all_of(wn_cols)))
  test_spec  <- snv(test_f  |> select(all_of(wn_cols)))

  train_ready <- bind_cols(pyc_log = log(train_f$pyc_abs), train_spec)
  test_ready  <- bind_cols(pyc_log = log(test_f$pyc_abs), test_spec)

  fit_f <- svm_mod |> fit(pyc_log ~ ., data = train_ready)
  preds_f <- exp(predict(fit_f, test_ready)$.pred)

  all_preds <- bind_rows(all_preds, tibble(obs     = test_f$pyc_abs,
                                            pred    = preds_f,
                                            cluster = factor(test_f$cluster)))
}

all_preds <- all_preds |>
  mutate(abs_error = abs(obs - pred),
         rel_error = abs_error / obs)

## Fig 4a: absolute error vs concentration ------------------------------------

ggplot(all_preds, aes(x = obs, y = abs_error, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.055, y = max(all_preds$abs_error) * 0.95,
           label = "PyC = 0.05", hjust = 0, size = 3, color = "grey40") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8,
              aes(group = 1), color = "black") +
  labs(title = "Prediction error vs PyC concentration",
       x = "Observed PyC (g / 100g soil)",
       y = "Absolute prediction error",
       color = "Cluster") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.85)) -> fig4a

## Fig 4b: relative error vs concentration ------------------------------------

ggplot(all_preds, aes(x = obs, y = rel_error * 100, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 100, linetype = "dotted", color = "grey50") +
  annotate("text", x = 0.42, y = 110, label = "100% relative error",
           hjust = 1, size = 3, color = "grey50") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8,
              aes(group = 1), color = "black") +
  labs(title = "Relative error vs PyC concentration",
       x = "Observed PyC (g / 100g soil)",
       y = "Relative prediction error (%)",
       color = "Cluster") +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 500)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.85)) -> fig4b

fig4a + fig4b + plot_annotation(tag_levels = "A") -> fig4

ggsave("output/figures/fig4_concentration_floor.png", fig4,
       width = 12, height = 5, dpi = 300)

## Fig 5: global pred vs obs --------------------------------------------------

r2_cv   <- round(cor(all_preds$obs, all_preds$pred)^2, 3)
rmse_cv <- round(sqrt(mean((all_preds$obs - all_preds$pred)^2)), 4)
rpd_cv  <- round(sd(all_preds$obs) / sqrt(mean((all_preds$obs - all_preds$pred)^2)), 2)

ggplot(all_preds, aes(x = obs, y = pred, color = cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50") +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "Predicted vs observed PyC (5-fold CV, SVM-RBF + SNV + log)",
       subtitle = sprintf("R\u00b2 = %.3f, RMSE = %.4f, RPD = %.2f  (n = %d)",
                          r2_cv, rmse_cv, rpd_cv, nrow(all_preds)),
       x = "Observed PyC (g / 100g soil)",
       y = "Predicted PyC (g / 100g soil)",
       color = "Cluster") +
  scale_color_brewer(palette = "Set1") +
  coord_equal(xlim = c(0, 0.5), ylim = c(0, 0.5)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.15)) -> fig5

ggsave("output/figures/fig5_pred_vs_obs_global.png", fig5,
       width = 6, height = 6, dpi = 300)


## =============================================================================
## Step 7: Preprocessing and LOO-CV helpers
## =============================================================================
##
## Everything below here runs LOO-CV per cluster to build the hybrid figure.
## This takes hours. Run Steps 1-6 independently if you only need the fast
## figures.

sg_smooth <- function(x) {
  savitzkyGolay(as.matrix(x), m = 0, p = 2, w = 11) |> as.data.frame()
}

deriv1 <- function(x) {
  savitzkyGolay(as.matrix(x), m = 1, p = 2, w = 11) |> as.data.frame()
}

snv_deriv1 <- function(x) {
  deriv1(snv(x))
}

cor_select <- function(train_x, train_y, n = 100) {
  cors <- abs(cor(as.matrix(train_x), train_y))
  top <- order(cors, decreasing = TRUE)[1:min(n, ncol(train_x))]
  names(train_x)[top]
}

pca_reduce <- function(train_x, test_x) {
  pc <- prcomp(as.matrix(train_x), center = TRUE, scale. = TRUE)
  cv <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
  n_pc <- which(cv >= 0.99)[1]
  train_pc <- as.data.frame(pc$x[, 1:n_pc])
  test_pc  <- as.data.frame(predict(pc, as.matrix(test_x))[, 1:n_pc,
                                                             drop = FALSE])
  list(train = train_pc, test = test_pc)
}

loo_cv <- function(df, wn_cols, preprocess_fn, feature_fn, model_fn) {
  n <- nrow(df)
  preds <- numeric(n)

  for (i in seq_len(n)) {
    train <- df[-i, ]
    test  <- df[i, , drop = FALSE]

    train_x <- preprocess_fn(train |> select(all_of(wn_cols)))
    test_x  <- preprocess_fn(test  |> select(all_of(wn_cols)))

    feat    <- feature_fn(train_x, train, test_x)
    train_x <- feat$train
    test_x  <- feat$test

    train_ready <- bind_cols(pyc_log = log(train$pyc_abs), train_x)
    test_ready  <- bind_cols(pyc_log = log(test$pyc_abs), test_x)

    fit <- model_fn(train_ready)
    preds[i] <- exp(predict(fit, test_ready)$.pred)
  }

  tibble(obs = df$pyc_abs, pred = preds)
}


## =============================================================================
## Step 8: Per-cluster LOO-CV
## =============================================================================

## Cluster 1: Cubist + SNV+D1 + log + correlation -----------------------------

loo_cv(df            = analysis |> filter(cluster == 1),
       wn_cols       = wn_cols,
       preprocess_fn = snv_deriv1,
       feature_fn    = function(train_x, train, test_x) {
         keep <- cor_select(train_x, log(train$pyc_abs), n = 100)
         list(train = train_x[, keep], test = test_x[, keep])
       },
       model_fn      = function(d) {
         cubist_rules(committees = 92, neighbors = 1, max_rules = 396) |>
           set_engine("Cubist") |>
           set_mode("regression") |>
           fit(pyc_log ~ ., data = d)
       }) |>
  mutate(cluster = 1L) -> preds_cl1

## Cluster 2: Elastic Net + D1 + log + correlation ----------------------------

loo_cv(df            = analysis |> filter(cluster == 2),
       wn_cols       = wn_cols,
       preprocess_fn = deriv1,
       feature_fn    = function(train_x, train, test_x) {
         keep <- cor_select(train_x, log(train$pyc_abs), n = 100)
         list(train = train_x[, keep], test = test_x[, keep])
       },
       model_fn      = function(d) {
         linear_reg(penalty = 6.16e-05, mixture = 0.2) |>
           set_engine("glmnet") |>
           set_mode("regression") |>
           fit(pyc_log ~ ., data = d)
       }) |>
  mutate(cluster = 2L) -> preds_cl2

## Cluster 3: Global SVM-RBF + SG + log + PCA (LOO on full dataset) -----------

# LOO on the full dataset; only cluster 3 predictions used in the hybrid plot
n_all     <- nrow(analysis)
preds_all <- numeric(n_all)

for (i in seq_len(n_all)) {
  train <- analysis[-i, ]
  test  <- analysis[i, , drop = FALSE]

  train_x <- sg_smooth(train |> select(all_of(wn_cols)))
  test_x  <- sg_smooth(test  |> select(all_of(wn_cols)))

  pc <- pca_reduce(train_x, test_x)

  train_ready <- bind_cols(pyc_log = log(train$pyc_abs), pc$train)
  test_ready  <- bind_cols(pyc_log = log(test$pyc_abs), pc$test)

  fit <- svm_rbf(cost = 2.56, rbf_sigma = 0.048) |>
    set_engine("kernlab") |>
    set_mode("regression") |>
    fit(pyc_log ~ ., data = train_ready)

  preds_all[i] <- exp(predict(fit, test_ready)$.pred)
}

# Extract cluster 3 predictions
cl3_idx <- analysis$cluster == 3

tibble(obs     = analysis$pyc_abs[cl3_idx],
       pred    = preds_all[cl3_idx],
       cluster = 3L) -> preds_cl3


## =============================================================================
## Step 9: Fig 6 — Hybrid pred vs obs
## =============================================================================

bind_rows(preds_cl1, preds_cl2, preds_cl3) |>
  mutate(cluster = factor(cluster)) -> pooled

rpd_pooled  <- sd(pooled$obs) / sqrt(mean((pooled$obs - pooled$pred)^2))
r2_pooled   <- cor(pooled$obs, pooled$pred)^2
rmse_pooled <- sqrt(mean((pooled$obs - pooled$pred)^2))
ccc_num     <- 2 * cor(pooled$obs, pooled$pred) * sd(pooled$obs) * sd(pooled$pred)
ccc_den     <- var(pooled$obs) + var(pooled$pred) +
  (mean(pooled$obs) - mean(pooled$pred))^2
ccc_pooled  <- ccc_num / ccc_den

## Per-cluster breakdown -------------------------------------------------------

pooled |>
  group_by(cluster) |>
  summarise(n    = n(),
            rpd  = round(sd(obs) / sqrt(mean((obs - pred)^2)), 2),
            r2   = round(cor(obs, pred)^2, 3),
            rmse = round(sqrt(mean((obs - pred)^2)), 4),
            .groups = "drop")

## Hybrid 1:1 plot -------------------------------------------------------------

ggplot(pooled, aes(x = obs, y = pred, color = cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50") +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "Hybrid approach: local models (clusters 1-2) + global model (cluster 3)",
       subtitle = sprintf("Pooled: R\u00b2 = %.3f, RMSE = %.4f, RPD = %.2f, CCC = %.3f  (n = %d)",
                          r2_pooled, rmse_pooled, rpd_pooled, ccc_pooled,
                          nrow(pooled)),
       x = "Observed PyC (g / 100g soil)",
       y = "Predicted PyC (g / 100g soil)",
       color = "Cluster") +
  scale_color_brewer(palette = "Set1",
                     labels = c("1" = "Cl 1 \u2014 local Cubist",
                                "2" = "Cl 2 \u2014 local Elastic Net",
                                "3" = "Cl 3 \u2014 global SVM-RBF")) +
  coord_equal(xlim = c(0, 0.5), ylim = c(0, 0.5)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.18)) -> fig6

ggsave("output/figures/fig6_pred_vs_obs_hybrid.png", fig6,
       width = 7, height = 7, dpi = 300)
