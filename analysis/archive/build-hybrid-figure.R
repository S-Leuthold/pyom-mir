#' Hybrid 1:1 figure: local models for clusters 1-2, global for cluster 3
#' LOO-CV for stable estimates with small cluster sizes.
#' Tuned hyperparameters pulled from horizons grid search results.

library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(parsnip)
library(rules)

hd_full <- readRDS("data/processed/horizons_data.rds")
wn_cols <- hd_full$data$role_map |> filter(role == "predictor") |> pull(variable)
analysis <- hd_full$data$analysis |> filter(!is.na(pyc_abs))
spec <- analysis |> select(all_of(wn_cols)) |> as.matrix()

## Reproduce clusters
pca_fit <- prcomp(spec, center = TRUE, scale. = TRUE)
cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
n_pcs <- which(cum_var >= 0.99)[1]
scores <- pca_fit$x[, 1:n_pcs]
clusters <- cutree(hclust(dist(scores), method = "ward.D2"), k = 3)
analysis$cluster <- clusters


## =============================================================================
## Preprocessing helpers
## =============================================================================

snv <- function(x) {
  m <- as.matrix(x)
  t(apply(m, 1, function(r) (r - mean(r)) / sd(r))) |> as.data.frame()
}

sg_smooth <- function(x) {
  prospectr::savitzkyGolay(as.matrix(x), m = 0, p = 2, w = 11) |> as.data.frame()
}

deriv1 <- function(x) {
  prospectr::savitzkyGolay(as.matrix(x), m = 1, p = 2, w = 11) |> as.data.frame()
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
  test_pc <- as.data.frame(predict(pc, as.matrix(test_x))[, 1:n_pc, drop = FALSE])
  list(train = train_pc, test = test_pc)
}


## =============================================================================
## LOO helper
## =============================================================================

loo_cv <- function(df, wn_cols, preprocess_fn, feature_fn, model_fn, label) {
  n <- nrow(df)
  preds <- numeric(n)

  t0 <- Sys.time()
  for (i in seq_len(n)) {
    train <- df[-i, ]
    test  <- df[i, , drop = FALSE]

    train_x <- preprocess_fn(train |> select(all_of(wn_cols)))
    test_x  <- preprocess_fn(test  |> select(all_of(wn_cols)))

    feat <- feature_fn(train_x, train, test_x)
    train_x <- feat$train
    test_x  <- feat$test

    train_ready <- bind_cols(pyc_log = log(train$pyc_abs), train_x)
    test_ready  <- bind_cols(pyc_log = log(test$pyc_abs), test_x)

    fit <- model_fn(train_ready)
    preds[i] <- exp(predict(fit, test_ready)$.pred)

    if (i %% 25 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("  [%d/%d] %.1f min\n", i, n, elapsed))
    }
  }

  obs <- df$pyc_abs
  rmse <- sqrt(mean((obs - preds)^2))
  rpd <- sd(obs) / rmse
  r2 <- cor(obs, preds)^2

  cat(sprintf("  %s LOO: RPD=%.2f, R²=%.3f, RMSE=%.4f (n=%d)\n",
      label, rpd, r2, rmse, n))

  tibble(obs = obs, pred = preds)
}


## =============================================================================
## Cluster 1: Cubist + SNV+D1 + log + correlation (LOO)
## =============================================================================

cat("=== Cluster 1: Cubist + SNV+D1 + log + correlation (LOO) ===\n")
cl1 <- analysis |> filter(cluster == 1)

preds_cl1 <- loo_cv(
  df = cl1, wn_cols = wn_cols,
  preprocess_fn = snv_deriv1,
  feature_fn = function(train_x, train, test_x) {
    keep <- cor_select(train_x, log(train$pyc_abs), n = 100)
    list(train = train_x[, keep], test = test_x[, keep])
  },
  model_fn = function(d) {
    cubist_rules(committees = 92, neighbors = 1, max_rules = 396) |>
      set_engine("Cubist") |>
      set_mode("regression") |>
      fit(pyc_log ~ ., data = d)
  },
  label = "Cluster 1"
) |> mutate(cluster = 1L)


## =============================================================================
## Cluster 2: Elastic Net + D1 + log + correlation (LOO)
## =============================================================================

cat("=== Cluster 2: Elastic Net + D1 + log + correlation (LOO) ===\n")
cl2 <- analysis |> filter(cluster == 2)

preds_cl2 <- loo_cv(
  df = cl2, wn_cols = wn_cols,
  preprocess_fn = deriv1,
  feature_fn = function(train_x, train, test_x) {
    keep <- cor_select(train_x, log(train$pyc_abs), n = 100)
    list(train = train_x[, keep], test = test_x[, keep])
  },
  model_fn = function(d) {
    linear_reg(penalty = 6.16e-05, mixture = 0.2) |>
      set_engine("glmnet") |>
      set_mode("regression") |>
      fit(pyc_log ~ ., data = d)
  },
  label = "Cluster 2"
) |> mutate(cluster = 2L)


## =============================================================================
## Cluster 3: Global SVM-RBF + SG + log + PCA (LOO on full dataset, extract cl3)
## =============================================================================

cat("=== Cluster 3: Global SVM-RBF + SG + log + PCA (LOO) ===\n")

## LOO on full dataset, only keep cluster 3 predictions
n_all <- nrow(analysis)
preds_all <- numeric(n_all)

t0 <- Sys.time()
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

  if (i %% 50 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  [%d/%d] %.1f min\n", i, n_all, elapsed))
  }
}

## Extract cluster 3 predictions
cl3_idx <- analysis$cluster == 3
preds_cl3 <- tibble(
  obs = analysis$pyc_abs[cl3_idx],
  pred = preds_all[cl3_idx],
  cluster = 3L
)

obs3 <- preds_cl3$obs
pred3 <- preds_cl3$pred
cat(sprintf("  Cluster 3 LOO: RPD=%.2f, R²=%.3f, RMSE=%.4f (n=%d)\n",
    sd(obs3) / sqrt(mean((obs3 - pred3)^2)),
    cor(obs3, pred3)^2,
    sqrt(mean((obs3 - pred3)^2)),
    nrow(preds_cl3)))

## Also report global LOO metrics
cat(sprintf("  Global LOO:    RPD=%.2f, R²=%.3f, RMSE=%.4f (n=%d)\n",
    sd(analysis$pyc_abs) / sqrt(mean((analysis$pyc_abs - preds_all)^2)),
    cor(analysis$pyc_abs, preds_all)^2,
    sqrt(mean((analysis$pyc_abs - preds_all)^2)),
    n_all))


## =============================================================================
## Pool and plot
## =============================================================================

cat("\n=== Pooled hybrid results ===\n")

pooled <- bind_rows(preds_cl1, preds_cl2, preds_cl3) |>
  mutate(cluster = factor(cluster))

rpd_pooled <- sd(pooled$obs) / sqrt(mean((pooled$obs - pooled$pred)^2))
r2_pooled <- cor(pooled$obs, pooled$pred)^2
rmse_pooled <- sqrt(mean((pooled$obs - pooled$pred)^2))
ccc_num <- 2 * cor(pooled$obs, pooled$pred) * sd(pooled$obs) * sd(pooled$pred)
ccc_den <- var(pooled$obs) + var(pooled$pred) + (mean(pooled$obs) - mean(pooled$pred))^2
ccc_pooled <- ccc_num / ccc_den

cat(sprintf("Total samples: %d\n", nrow(pooled)))
cat(sprintf("RPD:  %.2f\n", rpd_pooled))
cat(sprintf("R²:   %.3f\n", r2_pooled))
cat(sprintf("RMSE: %.4f\n", rmse_pooled))
cat(sprintf("CCC:  %.3f\n", ccc_pooled))

cat("\nPer-cluster breakdown:\n")
pooled |>
  group_by(cluster) |>
  summarise(
    n = n(),
    rpd = round(sd(obs) / sqrt(mean((obs - pred)^2)), 2),
    r2 = round(cor(obs, pred)^2, 3),
    rmse = round(sqrt(mean((obs - pred)^2)), 4),
    .groups = "drop"
  ) |>
  print()

fig5 <- ggplot(pooled, aes(x = obs, y = pred, color = cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    title = "Hybrid approach: local models (clusters 1-2) + global model (cluster 3)",
    subtitle = sprintf("Pooled: R\u00b2 = %.3f, RMSE = %.4f, RPD = %.2f, CCC = %.3f  (n = %d)",
                       r2_pooled, rmse_pooled, rpd_pooled, ccc_pooled, nrow(pooled)),
    x = "Observed PyC (g / 100g soil)",
    y = "Predicted PyC (g / 100g soil)",
    color = "Cluster"
  ) +
  scale_color_brewer(palette = "Set1",
    labels = c("1" = "Cl 1 \u2014 local Cubist",
               "2" = "Cl 2 \u2014 local Elastic Net",
               "3" = "Cl 3 \u2014 global SVM-RBF")) +
  coord_equal(xlim = c(0, 0.5), ylim = c(0, 0.5)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.18))

ggsave("output/figures/fig5_pred_vs_obs.png", fig5,
       width = 7, height = 7, dpi = 300)
cat("\nSaved: fig5_pred_vs_obs.png\n")
