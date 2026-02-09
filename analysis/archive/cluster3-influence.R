#' Cluster 3 — LOO Influence Analysis
#'
#' For each sample in cluster 3, remove it from training, refit SVM-RBF,
#' predict the held-out test set, and measure RPD. Samples that improve RPD
#' when removed are hurting the model.

library(dplyr, warn.conflicts = FALSE)
library(parsnip)
library(rsample)

hd_full <- readRDS("data/processed/horizons_data.rds")
wn_cols <- hd_full$data$role_map |> filter(role == "predictor") |> pull(variable)
analysis <- hd_full$data$analysis |> filter(!is.na(pyc_abs))
spec <- analysis |> select(all_of(wn_cols)) |> as.matrix()

pca_fit <- prcomp(spec, center = TRUE, scale. = TRUE)
cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
scores  <- pca_fit$x[, 1:which(cum_var >= 0.99)[1]]
clusters <- cutree(hclust(dist(scores), method = "ward.D2"), k = 3)
analysis$cluster <- clusters

resp <- readr::read_csv("data/processed/modeling_response.csv", show_col_types = FALSE)
analysis <- analysis |>
  left_join(resp |> select(Sample_ID, project), by = c("sample_id" = "Sample_ID"))

## Cluster 3 — all samples including FORCE
cl3 <- analysis |> filter(cluster == 3)
cat(sprintf("Cluster 3: %d samples\n", nrow(cl3)))


## =============================================================================
## Setup: fixed train/test split, SNV preprocessing, log transform
## =============================================================================

model_df <- cl3 |> select(sample_id, project, pyc_abs, all_of(wn_cols))

set.seed(307)
split <- initial_split(model_df, prop = 0.8, strata = pyc_abs)
train <- training(split)
test  <- testing(split)

cat(sprintf("Train: %d, Test: %d\n", nrow(train), nrow(test)))

snv <- function(x) {
  m <- as.matrix(x)
  t(apply(m, 1, function(r) (r - mean(r)) / sd(r))) |> as.data.frame()
}

## Preprocess test set once (it stays fixed)
test_spec <- snv(test |> select(all_of(wn_cols)))
test_ready <- bind_cols(pyc_log = log(test$pyc_abs), test_spec)
test_obs <- test$pyc_abs

svm_spec <- svm_rbf(cost = 2.5, rbf_sigma = 0.05) |>
  set_engine("kernlab") |>
  set_mode("regression")


## =============================================================================
## Baseline: fit on full training set
## =============================================================================

train_spec <- snv(train |> select(all_of(wn_cols)))
train_ready <- bind_cols(pyc_log = log(train$pyc_abs), train_spec)

fit_base <- svm_spec |> fit(pyc_log ~ ., data = train_ready)
preds_base <- exp(predict(fit_base, test_ready)$.pred)
rmse_base <- sqrt(mean((test_obs - preds_base)^2))
rpd_base <- sd(test_obs) / rmse_base

cat(sprintf("\nBaseline RPD: %.3f (RMSE: %.4f)\n", rpd_base, rmse_base))


## =============================================================================
## LOO influence: remove each training sample, refit, measure test RPD
## =============================================================================

cat(sprintf("\nRunning LOO influence on %d training samples...\n", nrow(train)))

results <- tibble(
  sample_id = train$sample_id,
  project   = train$project,
  pyc_abs   = train$pyc_abs,
  rpd_without = NA_real_,
  rpd_delta   = NA_real_
)

t0 <- Sys.time()
for (i in seq_len(nrow(train))) {

  if (i %% 10 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta <- elapsed / i * (nrow(train) - i)
    cat(sprintf("  [%d/%d] %.1f min elapsed, ~%.1f min remaining\n",
                i, nrow(train), elapsed, eta))
  }

  ## Remove sample i
  train_i <- train_ready[-i, ]

  ## Refit
  fit_i <- tryCatch(
    svm_spec |> fit(pyc_log ~ ., data = train_i),
    error = function(e) NULL
  )
  if (is.null(fit_i)) next

  ## Predict test set
  preds_i <- exp(predict(fit_i, test_ready)$.pred)
  rmse_i <- sqrt(mean((test_obs - preds_i)^2))
  rpd_i <- sd(test_obs) / rmse_i

  results$rpd_without[i] <- rpd_i
  results$rpd_delta[i] <- rpd_i - rpd_base
}

elapsed_total <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nDone in %.1f minutes\n", elapsed_total))


## =============================================================================
## Results
## =============================================================================

cat("\n=== Most harmful samples (RPD improves when removed) ===\n")
results |>
  filter(!is.na(rpd_delta)) |>
  arrange(desc(rpd_delta)) |>
  head(20) |>
  mutate(across(c(pyc_abs, rpd_without, rpd_delta), ~ round(.x, 4))) |>
  print(n = 20)

cat("\n=== Most helpful samples (RPD drops when removed) ===\n")
results |>
  filter(!is.na(rpd_delta)) |>
  arrange(rpd_delta) |>
  head(10) |>
  mutate(across(c(pyc_abs, rpd_without, rpd_delta), ~ round(.x, 4))) |>
  print(n = 10)

cat("\n=== Influence by project ===\n")
results |>
  filter(!is.na(rpd_delta)) |>
  group_by(project) |>
  summarise(
    n = n(),
    mean_delta = round(mean(rpd_delta), 4),
    max_delta  = round(max(rpd_delta), 4),
    n_harmful  = sum(rpd_delta > 0),
    .groups = "drop"
  ) |>
  arrange(desc(mean_delta)) |>
  print()

## Save for later
saveRDS(results, "output/local-clusters-k3/cluster3_influence.rds")
cat("\nSaved: output/local-clusters-k3/cluster3_influence.rds\n")
