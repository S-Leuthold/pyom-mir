#' Sanity Check — Bypass horizons, fit models directly
#'
#' @description
#' Quick standalone test to rule out horizons pipeline bugs. Loads the same
#' data, extracts spectra + response manually, fits PLS and RF with a simple
#' 70/30 train/test split. If this gives RPD ~1.05 too, it's the data.
#' If it gives something better, horizons is broken.
#'
#' @keywords internal
#' @name sanity-check


## =============================================================================
## Section 1: Load and extract data
## =============================================================================


hd <- readRDS("data/processed/horizons_data.rds")

## Get predictor columns (wavenumbers)
wn_cols <- hd$data$role_map |>
  dplyr::filter(role == "predictor") |>
  dplyr::pull(variable)

## Build a clean data frame: response + predictors
df <- hd$data$analysis[, c("pyc_abs", wn_cols)]

## Drop rows with NA response
df <- df[!is.na(df$pyc_abs), ]

cat(sprintf("Samples with response: %d\n", nrow(df)))
cat(sprintf("Predictors: %d\n", length(wn_cols)))
cat(sprintf("Response range: %.3f — %.3f (mean %.3f, median %.3f)\n",
            min(df$pyc_abs), max(df$pyc_abs),
            mean(df$pyc_abs), median(df$pyc_abs)))


## =============================================================================
## Section 2: Train/test split
## =============================================================================


set.seed(307)
n <- nrow(df)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))
test_idx  <- setdiff(seq_len(n), train_idx)

train <- df[train_idx, ]
test  <- df[test_idx, ]

cat(sprintf("\nTrain: %d, Test: %d\n", nrow(train), nrow(test)))

## Helper: compute RPD, R², RMSE
eval_metrics <- function(observed, predicted) {
  rmse <- sqrt(mean((observed - predicted)^2))
  ss_res <- sum((observed - predicted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  rsq <- 1 - ss_res / ss_tot
  rpd <- sd(observed) / rmse
  c(RPD = rpd, Rsq = rsq, RMSE = rmse)
}


## =============================================================================
## Section 3: Baseline — predict the mean
## =============================================================================


cat("\n=== Baseline: predict train mean ===\n")
pred_mean <- rep(mean(train$pyc_abs), nrow(test))
print(round(eval_metrics(test$pyc_abs, pred_mean), 3))


## =============================================================================
## Section 4: PCA + linear model (simplest possible)
## =============================================================================


cat("\n=== PCA (10 components) + linear model ===\n")

X_train <- as.matrix(train[, wn_cols])
X_test  <- as.matrix(test[, wn_cols])

pca_fit <- prcomp(X_train, center = TRUE, scale. = FALSE)

## Project both sets
pc_train <- pca_fit$x[, 1:10]
pc_test  <- predict(pca_fit, X_test)[, 1:10]

## Simple linear regression on PCs
pc_df_train <- data.frame(pyc_abs = train$pyc_abs, pc_train)
pc_df_test  <- data.frame(pyc_abs = test$pyc_abs, pc_test)

lm_fit <- lm(pyc_abs ~ ., data = pc_df_train)
pred_lm <- predict(lm_fit, pc_df_test)

print(round(eval_metrics(test$pyc_abs, pred_lm), 3))

## Show what PCs correlate with
pc_cors <- cor(pca_fit$x[, 1:10], train$pyc_abs)
cat("\nPC correlations with pyc_abs:\n")
for (i in 1:10) {
  cat(sprintf("  PC%d: r = %.3f (%.1f%% var)\n", i, pc_cors[i],
              100 * summary(pca_fit)$importance[2, i]))
}


## =============================================================================
## Section 5: PLS regression (direct, no parsnip)
## =============================================================================


cat("\n=== PLS regression (10 components) ===\n")

library(pls)

pls_fit <- plsr(pyc_abs ~ ., data = data.frame(pyc_abs = train$pyc_abs, X_train),
                ncomp = 10, validation = "CV", segments = 5)

## Find optimal ncomp via CV
cv_rmse <- RMSEP(pls_fit, estimate = "CV")
best_ncomp <- which.min(cv_rmse$val[1, 1, -1])
cat(sprintf("Best ncomp by CV: %d\n", best_ncomp))

pred_pls <- predict(pls_fit, newdata = data.frame(X_test), ncomp = best_ncomp)[, , 1]
cat("PLS metrics:\n")
print(round(eval_metrics(test$pyc_abs, pred_pls), 3))


## =============================================================================
## Section 6: Random Forest (ranger, no feature selection)
## =============================================================================


cat("\n=== Random Forest (ranger, all features) ===\n")

rf_fit <- ranger::ranger(pyc_abs ~ ., data = train, num.trees = 500,
                         mtry = floor(sqrt(length(wn_cols))),
                         seed = 307)

pred_rf <- predict(rf_fit, data = test)$predictions
cat("RF metrics:\n")
print(round(eval_metrics(test$pyc_abs, pred_rf), 3))


## =============================================================================
## Section 7: Summary
## =============================================================================


cat("\n===============================================================================\n")
cat("Sanity Check Summary\n")
cat("===============================================================================\n")
cat("If RPD > 1.2 here but ~1.05 in horizons → horizons bug\n")
cat("If RPD ~1.05 here too → spectra genuinely lack signal\n")
cat("===============================================================================\n")
