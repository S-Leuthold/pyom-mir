#' Quick test: SVM-RBF on cluster 3 targeted bands
#' Single model, no horizons, no parallelism. Safe to run alongside other jobs.

library(dplyr, warn.conflicts = FALSE)
library(parsnip)
library(rsample)
library(yardstick)

hd_full <- readRDS("data/processed/horizons_data.rds")

## Reproduce clustering
wn_cols <- hd_full$data$role_map |>
  filter(role == "predictor") |> pull(variable)
analysis <- hd_full$data$analysis |> filter(!is.na(pyc_abs))
spec <- analysis |> select(all_of(wn_cols)) |> as.matrix()

pca_fit <- prcomp(spec, center = TRUE, scale. = TRUE)
cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
scores  <- pca_fit$x[, 1:which(cum_var >= 0.99)[1]]

clusters <- cutree(hclust(dist(scores), method = "ward.D2"), k = 3)

## Cluster 3, no FORCE
resp <- readr::read_csv("data/processed/modeling_response.csv",
                        show_col_types = FALSE)
analysis$cluster <- clusters
analysis <- analysis |>
  left_join(resp |> select(Sample_ID, project),
            by = c("sample_id" = "Sample_ID"))

cl3 <- analysis |> filter(cluster == 3, project != "FORCE")
cat(sprintf("Cluster 3 (no FORCE): %d samples\n", nrow(cl3)))

## Targeted PyC bands
wn_nums <- as.numeric(gsub("^wn_", "", wn_cols))
in_band <- (wn_nums >= 3000 & wn_nums <= 3100) |
           (wn_nums >= 1550 & wn_nums <= 1750) |
           (wn_nums >= 1200 & wn_nums <= 1300) |
           (wn_nums >= 700  & wn_nums <= 900)
targeted_cols <- wn_cols[in_band]
cat(sprintf("Targeted features: %d\n", length(targeted_cols)))

## Build modeling frame (keep all wn cols for full-spectrum comparison)
all_wn <- wn_cols[wn_cols %in% names(cl3)]
model_df <- cl3 |>
  select(pyc_abs, all_of(all_wn))

## Train/test split (80/20, stratified)
set.seed(307)
split <- initial_split(model_df, prop = 0.8, strata = pyc_abs)
train <- training(split)
test  <- testing(split)
cat(sprintf("Train: %d, Test: %d\n", nrow(train), nrow(test)))

## Log-transform response
train$pyc_log <- log(train$pyc_abs)
test$pyc_log  <- log(test$pyc_abs)

## SNV preprocess (row-normalize)
snv <- function(x) {
  m <- as.matrix(x)
  t(apply(m, 1, function(r) (r - mean(r)) / sd(r))) |> as.data.frame()
}

train_spec <- snv(train |> select(all_of(targeted_cols)))
test_spec  <- snv(test  |> select(all_of(targeted_cols)))

train_ready <- bind_cols(pyc_log = train$pyc_log, train_spec)
test_ready  <- bind_cols(pyc_log = test$pyc_log, test_spec)

## SVM-RBF — use hyperparams close to global best (cost=2.5, sigma=0.05)
svm_spec <- svm_rbf(cost = 2.5, rbf_sigma = 0.05) |>
  set_engine("kernlab") |>
  set_mode("regression")

cat("Fitting SVM-RBF (SNV + log + targeted bands)...\n")
fit <- svm_spec |> fit(pyc_log ~ ., data = train_ready)

## Predict + back-transform
preds_log <- predict(fit, test_ready)$.pred
preds     <- exp(preds_log)
obs       <- test$pyc_abs

## Metrics
results <- tibble(obs = obs, pred = preds) |>
  mutate(residual = obs - pred)

rmse_val <- sqrt(mean(results$residual^2))
r2_val   <- cor(results$obs, results$pred)^2
sd_obs   <- sd(results$obs)
rpd_val  <- sd_obs / rmse_val

cat(sprintf("\n=== Results: SVM-RBF + SNV + log + targeted bands ===\n"))
cat(sprintf("RPD:  %.2f\n", rpd_val))
cat(sprintf("R²:   %.3f\n", r2_val))
cat(sprintf("RMSE: %.4f\n", rmse_val))
cat(sprintf("SD:   %.4f\n", sd_obs))
cat(sprintf("n_test: %d\n", nrow(test)))

## Quick comparison: same model on full spectrum
cat("\n--- Full spectrum comparison ---\n")
train_full <- snv(train |> select(all_of(all_wn)))
test_full  <- snv(test  |> select(all_of(all_wn)))

train_full_ready <- bind_cols(pyc_log = train$pyc_log, train_full)
test_full_ready  <- bind_cols(pyc_log = test$pyc_log, test_full)

cat("Fitting SVM-RBF (SNV + log + full spectrum)...\n")
fit_full <- svm_spec |> fit(pyc_log ~ ., data = train_full_ready)

preds_full <- exp(predict(fit_full, test_full_ready)$.pred)
res_full   <- tibble(obs = obs, pred = preds_full) |>
  mutate(residual = obs - pred)
rmse_full <- sqrt(mean(res_full$residual^2))
rpd_full  <- sd_obs / rmse_full
r2_full   <- cor(res_full$obs, res_full$pred)^2

cat(sprintf("RPD:  %.2f\n", rpd_full))
cat(sprintf("R²:   %.3f\n", r2_full))
cat(sprintf("RMSE: %.4f\n", rmse_full))

cat(sprintf("\nDelta RPD (targeted - full): %+.2f\n", rpd_val - rpd_full))
