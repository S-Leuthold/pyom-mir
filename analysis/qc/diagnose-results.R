## =============================================================================
## diagnose-results.R
## Interactive diagnostics for HPC evaluate() results
##
## Author:        Sam Leuthold
## Contact:       sam.leuthold@colostate.edu
## Last modified: 2026-02-09
##
## Description:
##   Post-hoc diagnostic script for the Sybil HPC evaluate() run. Run
##   section-by-section in RStudio to investigate model performance, response
##   distribution, data quality, and spectral structure.
##
## Input:
##   output/models/run_20260204_121214/pyom_eval_result.rds
##   data/processed/horizons_data.rds
##   data/processed/modeling_response.csv
##
## Output:
##   Interactive plots and tables (no files written)
##
## Changelog:
##   2026-02-09  Restyled from original diagnose-results.R
##   2026-02-05  Initial version for investigating RPD 1.36 vs expected ~2.0
## =============================================================================


pacman::p_load(tidyverse, horizons)


## =============================================================================
## Step 1: Load results
## =============================================================================

result <- readRDS("output/models/run_20260204_121214/pyom_eval_result.rds")
hd     <- readRDS("data/processed/horizons_data.rds")
resp   <- read_csv("data/processed/modeling_response.csv",
                   show_col_types = FALSE)

eval_results <- result$evaluation$results


## =============================================================================
## Step 2: Run overview
## =============================================================================

## Status counts ---------------------------------------------------------------

eval_results |>
  count(status)

## Train/test split and runtime ------------------------------------------------

result$evaluation$n_train
result$evaluation$n_test
result$evaluation$runtime_secs / 60
result$evaluation$best_config


## =============================================================================
## Step 3: Failure analysis
## =============================================================================

## Failures by model type ------------------------------------------------------

eval_results |>
  filter(status == "failed") |>
  mutate(model = str_extract(config_id, "^[^_]+")) |>
  count(model, name = "n_failed")

## Distinct error messages -----------------------------------------------------

eval_results |>
  filter(status == "failed") |>
  distinct(error_message) |>
  pull(error_message) |>
  head(3)


## =============================================================================
## Step 4: Model performance landscape
## =============================================================================

## Parse config components from the config table, not the config_id string
## (some model names have underscores: svm_rbf, elastic_net)

config_lookup <- result$config$configs |>
  select(config_id, model, preprocessing, transformation, feature_selection)

successful <- eval_results |>
  filter(status == "success") |>
  left_join(config_lookup, by = "config_id")

## Top 10 by RPD ---------------------------------------------------------------

successful |>
  arrange(desc(rpd)) |>
  select(config_id, rpd, rsq, rmse, ccc, mae) |>
  head(10) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

## RPD distribution ------------------------------------------------------------

ggplot(successful, aes(x = rpd)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 1.36, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 2.0, color = "darkgreen", linetype = "dashed") +
  annotate("text", x = 1.36, y = Inf, label = "Best (1.36)",
           vjust = 2, hjust = -0.1, color = "red", size = 3) +
  annotate("text", x = 2.0, y = Inf, label = "Target (2.0)",
           vjust = 2, hjust = -0.1, color = "darkgreen", size = 3) +
  labs(title = "RPD Distribution Across Successful Configs",
       x = "RPD", y = "Count") +
  theme_minimal()

## Heatmap: preprocessing x model, best RPD per combo -------------------------

successful |>
  group_by(model, preprocessing) |>
  summarise(best_rpd = max(rpd, na.rm = TRUE),
            .groups = "drop") -> heatmap_data

ggplot(heatmap_data, aes(x = preprocessing, y = model, fill = best_rpd)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(best_rpd, 2)), size = 3) +
  scale_fill_viridis_c(option = "C", direction = 1) +
  labs(title = "Best RPD by Model x Preprocessing",
       fill = "RPD") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Transformation effect -------------------------------------------------------

successful |>
  group_by(transformation) |>
  summarise(n       = n(),
            rpd_med = median(rpd, na.rm = TRUE),
            rpd_max = max(rpd, na.rm = TRUE),
            rsq_med = median(rsq, na.rm = TRUE),
            .groups = "drop")

## Feature selection effect ----------------------------------------------------

successful |>
  group_by(feature_selection) |>
  summarise(n       = n(),
            rpd_med = median(rpd, na.rm = TRUE),
            rpd_max = max(rpd, na.rm = TRUE),
            .groups = "drop")


## =============================================================================
## Step 5: Best model diagnostics
## =============================================================================

## The evaluate() result stores the train/test split but not individual
## predictions. Need to reconstruct best model to get pred vs obs.

best_id  <- result$evaluation$best_config
best_row <- eval_results |> filter(config_id == best_id)

best_row |>
  select(config_id, rpd, rsq, rmse, ccc, mae)

best_row$best_params[[1]]

## Check what's available in the result object
has_cv_preds <- !is.null(result$artifacts$cv_preds$path)

## Split info
if (!is.null(result$evaluation$split)) {
  split_obj <- result$evaluation$split
  attr(split_obj, "strata")
}


## =============================================================================
## Step 6: Response variable analysis
## =============================================================================

outcome <- hd$data$analysis$spac_g_kg_c
outcome_clean <- outcome[!is.na(outcome)]

## Summary statistics ----------------------------------------------------------

summary(outcome_clean)
sd(outcome_clean)

# Skewness and kurtosis (manual — no extra package)
skew <- mean(((outcome_clean - mean(outcome_clean)) / sd(outcome_clean))^3)
kurt <- mean(((outcome_clean - mean(outcome_clean)) / sd(outcome_clean))^4) - 3
cv_pct <- 100 * sd(outcome_clean) / mean(outcome_clean)

## Histogram -------------------------------------------------------------------

tibble(spac = outcome_clean) |>
  ggplot(aes(x = spac)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  geom_vline(xintercept = median(outcome_clean), color = "red",
             linetype = "dashed") +
  labs(title = "Distribution of SPA-C (g/kg C)",
       subtitle = sprintf("n = %d, median = %.1f, skew = %.2f",
                          length(outcome_clean), median(outcome_clean), skew),
       x = "SPA-C (g/kg C)", y = "Count") +
  theme_minimal()

## By project ------------------------------------------------------------------

hd$data$analysis |>
  filter(!is.na(spac_g_kg_c)) |>
  left_join(resp |> select(Sample_ID, project),
            by = c("sample_id" = "Sample_ID")) -> analysis

analysis |>
  group_by(project) |>
  summarise(n    = n(),
            mean = round(mean(spac_g_kg_c), 1),
            sd   = round(sd(spac_g_kg_c), 1),
            min  = round(min(spac_g_kg_c), 1),
            max  = round(max(spac_g_kg_c), 1),
            .groups = "drop")

ggplot(analysis, aes(x = project, y = spac_g_kg_c, fill = project)) +
  geom_boxplot(outlier.shape = 21) +
  labs(title = "SPA-C Distribution by Project",
       y = "SPA-C (g/kg C)") +
  theme_minimal() +
  theme(legend.position = "none")

## Extreme values (> 3xIQR fence) ---------------------------------------------

q1 <- quantile(outcome_clean, 0.25)
q3 <- quantile(outcome_clean, 0.75)
upper_fence <- q3 + 3 * (q3 - q1)

analysis |>
  filter(spac_g_kg_c > upper_fence | spac_g_kg_c == 0) |>
  select(sample_id, project, spac_g_kg_c) |>
  arrange(desc(spac_g_kg_c))


## =============================================================================
## Step 7: Data quality audit
## =============================================================================

## Known anomalies — check which are still in modeling data --------------------

anomaly_ids <- c("AUS_54_Soil",                    # zero SPAC
                 "7_HTC_5_15_Soil",                # 309 g/kg C at 1.13% C
                 "Poudre_Watershed_Hist_BC_Soil",  # triple duplicate
                 "2_LTC_5_15_Soil",                # mislabeled depth
                 "793_796_P1U4_Litter_Litter",     # litter >400
                 "729_732_P3U2_Litter_Litter")     # litter >400

analysis |>
  filter(sample_id %in% anomaly_ids) |>
  select(sample_id, project, spac_g_kg_c)

## Replicate completeness (from deduplicate-ftir.R manifest) -------------------

manifest_path <- "data/processed/ftir_manifest.csv"

if (file.exists(manifest_path)) {
  manifest <- read_csv(manifest_path, show_col_types = FALSE)

  manifest |> count(flag, name = "n_stems")

  # Which incomplete stems are in the modeling data?
  manifest |>
    filter(flag == "incomplete") |>
    mutate(in_model = stem %in% analysis$sample_id) |>
    select(stem, n_replicates, in_model)
}


## =============================================================================
## Step 8: Spectral overview
## =============================================================================

## Extract spectral matrix -----------------------------------------------------

hd$data$role_map |>
  filter(role == "predictor") |>
  pull(variable) -> wn_cols

wn_values <- as.numeric(str_remove(wn_cols, "^wn_"))

analysis |>
  select(all_of(wn_cols)) |>
  as.matrix() -> spec_matrix

## Mean spectrum +/- SD --------------------------------------------------------

spec_mean <- colMeans(spec_matrix, na.rm = TRUE)
spec_sd   <- apply(spec_matrix, 2, sd, na.rm = TRUE)

tibble(wavenumber = wn_values,
       mean       = spec_mean,
       sd         = spec_sd,
       upper      = spec_mean + spec_sd,
       lower      = spec_mean - spec_sd) -> spec_summary

ggplot(spec_summary, aes(x = wavenumber)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = mean), color = "steelblue") +
  annotate("rect", xmin = 700, xmax = 900, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  annotate("text", x = 800, y = Inf, label = "Aromatic C-H\n700-900",
           vjust = 2, size = 2.5, color = "red") +
  annotate("rect", xmin = 1550, xmax = 1650, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  annotate("text", x = 1600, y = Inf, label = "Aromatic C=C\n~1600",
           vjust = 2, size = 2.5, color = "red") +
  scale_x_reverse() +
  labs(title = "Mean MIR Spectrum +/- 1 SD",
       subtitle = sprintf("n = %d samples, %d wavenumbers",
                          nrow(spec_matrix), length(wn_cols)),
       x = expression(Wavenumber~(cm^{-1})),
       y = "Absorbance") +
  theme_minimal()

## PCA colored by SPAC ---------------------------------------------------------

pca_fit <- prcomp(spec_matrix, center = TRUE, scale. = TRUE)

tibble(PC1     = pca_fit$x[, 1],
       PC2     = pca_fit$x[, 2],
       PC3     = pca_fit$x[, 3],
       spac    = analysis$spac_g_kg_c,
       project = analysis$project) -> pca_scores

var_explained <- summary(pca_fit)$importance[2, 1:3] * 100

ggplot(pca_scores, aes(x = PC1, y = PC2, color = spac)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_viridis_c(option = "C") +
  labs(title = "PCA of Spectra - Colored by SPA-C",
       subtitle = sprintf("PC1: %.1f%%, PC2: %.1f%% variance explained",
                          var_explained[1], var_explained[2]),
       color = "SPA-C\n(g/kg C)") +
  theme_minimal()

## PCA colored by project ------------------------------------------------------

ggplot(pca_scores, aes(x = PC1, y = PC2, color = project)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA of Spectra - Colored by Project",
       subtitle = sprintf("PC1: %.1f%%, PC2: %.1f%% variance explained",
                          var_explained[1], var_explained[2])) +
  theme_minimal()

## PC vs SPAC correlations -----------------------------------------------------

# Is spectral variation aligned with PyC?
tibble(pc  = paste0("PC", 1:3),
       r   = c(cor(pca_scores$PC1, pca_scores$spac, use = "complete.obs"),
               cor(pca_scores$PC2, pca_scores$spac, use = "complete.obs"),
               cor(pca_scores$PC3, pca_scores$spac, use = "complete.obs")))

ggplot(pca_scores, aes(x = PC1, y = spac)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "PC1 vs SPA-C",
       x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
       y = "SPA-C (g/kg C)") +
  theme_minimal()
