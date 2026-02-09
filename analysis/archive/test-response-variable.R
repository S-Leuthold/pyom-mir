#' Test Response Variable Options
#'
#' @description
#' Compare spectral correlations with different response formulations:
#'   1. spac_g_kg_c (current — PyC as fraction of total C)
#'   2. pyc_abs (absolute PyC — g PyC per 100g soil)
#'   3. sqrt and log transforms of each
#'
#' Uses ALS baseline-corrected spectra (our best preprocessing so far).
#'
#' @keywords internal
#' @name test-response-variable


## =============================================================================
## Section 1: Setup
## =============================================================================


library(tidyverse)
library(prospectr)
library(baseline)

hd   <- readRDS("data/processed/horizons_data.rds")
resp <- readr::read_csv("data/processed/modeling_response.csv",
                        show_col_types = FALSE)

wn_cols <- hd$data$role_map |>
  filter(role == "predictor") |>
  pull(variable)

wn_values <- as.numeric(str_remove(wn_cols, "^wn_"))

analysis <- hd$data$analysis |>
  filter(!is.na(spac_g_kg_c)) |>
  left_join(resp |> select(Sample_ID, pct_c, project),
            by = c("sample_id" = "Sample_ID"))

spec_matrix <- analysis |>
  select(all_of(wn_cols)) |>
  as.matrix()


## =============================================================================
## Section 2: Derive Response Variables
## =============================================================================


responses <- tibble(
  sample_id   = analysis$sample_id,
  project     = analysis$project,
  pct_c       = analysis$pct_c,
  spac_ratio  = analysis$spac_g_kg_c,
  pyc_abs     = analysis$pct_c * analysis$spac_g_kg_c / 1000
)

cat("=== Response Variable Summary ===\n\n")

cat("spac_ratio (g PyC / kg C) — current response:\n")
cat(sprintf("  n: %d, mean: %.2f, median: %.2f, SD: %.2f, range: %.2f-%.2f\n",
            sum(!is.na(responses$spac_ratio)),
            mean(responses$spac_ratio, na.rm = TRUE),
            median(responses$spac_ratio, na.rm = TRUE),
            sd(responses$spac_ratio, na.rm = TRUE),
            min(responses$spac_ratio, na.rm = TRUE),
            max(responses$spac_ratio, na.rm = TRUE)))

cat("\npyc_abs (g PyC / 100g soil) — absolute concentration:\n")
cat(sprintf("  n: %d, mean: %.4f, median: %.4f, SD: %.4f, range: %.4f-%.4f\n",
            sum(!is.na(responses$pyc_abs)),
            mean(responses$pyc_abs, na.rm = TRUE),
            median(responses$pyc_abs, na.rm = TRUE),
            sd(responses$pyc_abs, na.rm = TRUE),
            min(responses$pyc_abs, na.rm = TRUE),
            max(responses$pyc_abs, na.rm = TRUE)))

cat("\npct_c (total C %):\n")
cat(sprintf("  n: %d, mean: %.2f, median: %.2f, SD: %.2f, range: %.2f-%.2f\n",
            sum(!is.na(responses$pct_c)),
            mean(responses$pct_c, na.rm = TRUE),
            median(responses$pct_c, na.rm = TRUE),
            sd(responses$pct_c, na.rm = TRUE),
            min(responses$pct_c, na.rm = TRUE),
            max(responses$pct_c, na.rm = TRUE)))


## =============================================================================
## Section 3: ALS Baseline Correction
## =============================================================================


wn_ascending <- rev(wn_values)
spec_ascending <- spec_matrix[, rev(seq_len(ncol(spec_matrix)))]

cat("\nApplying ALS baseline correction...")
bc_als <- baseline::baseline(spec_ascending, method = "als",
                             lambda = 6, p = 0.05)
spec_als <- baseline::getCorrected(bc_als)
spec_als <- spec_als[, rev(seq_len(ncol(spec_als)))]
cat(" done\n")


## =============================================================================
## Section 4: Correlation Spectra — Ratio vs Absolute
## =============================================================================


## Filter to samples with both response variables
has_both <- !is.na(responses$spac_ratio) & !is.na(responses$pyc_abs)
cat(sprintf("\nSamples with both responses: %d\n", sum(has_both)))

## Raw spectra correlations
cor_ratio_raw <- apply(spec_matrix[has_both, ], 2, cor,
                       y = responses$spac_ratio[has_both],
                       use = "complete.obs")
cor_abs_raw   <- apply(spec_matrix[has_both, ], 2, cor,
                       y = responses$pyc_abs[has_both],
                       use = "complete.obs")

## ALS-corrected correlations
cor_ratio_als <- apply(spec_als[has_both, ], 2, cor,
                       y = responses$spac_ratio[has_both],
                       use = "complete.obs")
cor_abs_als   <- apply(spec_als[has_both, ], 2, cor,
                       y = responses$pyc_abs[has_both],
                       use = "complete.obs")

cor_all <- tibble(
  wavenumber  = rep(wn_values, 4),
  correlation = c(cor_ratio_raw, cor_abs_raw, cor_ratio_als, cor_abs_als),
  response    = rep(c("Ratio (g/kg C)", "Absolute (g/100g soil)",
                      "Ratio (g/kg C)", "Absolute (g/100g soil)"),
                    each = length(wn_values)),
  spectra     = rep(c("Raw", "Raw", "ALS-corrected", "ALS-corrected"),
                    each = length(wn_values))
) |>
  mutate(
    label = paste(spectra, "—", response),
    label = factor(label, levels = c(
      "Raw — Ratio (g/kg C)", "Raw — Absolute (g/100g soil)",
      "ALS-corrected — Ratio (g/kg C)", "ALS-corrected — Absolute (g/100g soil)"
    ))
  )

## Faceted comparison
ggplot(cor_all, aes(x = wavenumber, y = correlation)) +
  geom_line(linewidth = 0.4, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("rect", xmin = 700, xmax = 900, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  annotate("rect", xmin = 1550, xmax = 1650, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  scale_x_reverse() +
  facet_wrap(~ label, ncol = 1) +
  labs(title = "Correlation Spectra: Ratio vs Absolute Response",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Pearson r") +
  theme_minimal()

## Summary table
cat("\n=== Max |r| by Response × Spectra ===\n")
cor_all |>
  group_by(label) |>
  summarise(
    max_abs_r = round(max(abs(correlation), na.rm = TRUE), 3),
    best_wn   = wavenumber[which.max(abs(correlation))],
    .groups   = "drop"
  ) |>
  print()


## =============================================================================
## Section 5: PCA vs Response Variables
## =============================================================================


## PCA on ALS-corrected spectra
nonzero_var <- apply(spec_als[has_both, ], 2, var) > 0
pca_als <- prcomp(spec_als[has_both, nonzero_var], center = TRUE, scale. = TRUE)

var_expl <- summary(pca_als)$importance[2, 1:5] * 100

pca_scores <- tibble(
  PC1        = pca_als$x[, 1],
  PC2        = pca_als$x[, 2],
  PC3        = pca_als$x[, 3],
  PC4        = pca_als$x[, 4],
  PC5        = pca_als$x[, 5],
  spac_ratio = responses$spac_ratio[has_both],
  pyc_abs    = responses$pyc_abs[has_both],
  pct_c      = responses$pct_c[has_both],
  project    = responses$project[has_both]
)

## PC correlations with each response
cat("\n=== PC Correlations (ALS-corrected spectra) ===\n")
cat(sprintf("%-6s  %12s  %12s  %12s\n", "PC", "Ratio", "Absolute", "Total C"))
for (i in 1:5) {
  pc_col <- paste0("PC", i)
  r_ratio <- cor(pca_scores[[pc_col]], pca_scores$spac_ratio, use = "complete.obs")
  r_abs   <- cor(pca_scores[[pc_col]], pca_scores$pyc_abs, use = "complete.obs")
  r_c     <- cor(pca_scores[[pc_col]], pca_scores$pct_c, use = "complete.obs")
  cat(sprintf("PC%-4d  %12.3f  %12.3f  %12.3f  (%.1f%% var)\n",
              i, r_ratio, r_abs, r_c, var_expl[i]))
}

## Scatter: PC1 vs each response
p1 <- ggplot(pca_scores, aes(x = PC1, y = spac_ratio, color = project)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(title = "PC1 vs Ratio (g PyC / kg C)",
       x = sprintf("PC1 (%.1f%%)", var_expl[1]),
       y = "SPAC (g/kg C)") +
  theme_minimal()

p2 <- ggplot(pca_scores, aes(x = PC1, y = pyc_abs, color = project)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(title = "PC1 vs Absolute (g PyC / 100g soil)",
       x = sprintf("PC1 (%.1f%%)", var_expl[1]),
       y = "PyC (g/100g soil)") +
  theme_minimal()

print(p1)
print(p2)

## Also check: does absolute PyC correlate with total C?
cat(sprintf("\nCorrelation between responses:\n"))
cat(sprintf("  Ratio vs Absolute:  r = %.3f\n",
            cor(pca_scores$spac_ratio, pca_scores$pyc_abs, use = "complete.obs")))
cat(sprintf("  Ratio vs Total C:   r = %.3f\n",
            cor(pca_scores$spac_ratio, pca_scores$pct_c, use = "complete.obs")))
cat(sprintf("  Absolute vs Total C: r = %.3f\n",
            cor(pca_scores$pyc_abs, pca_scores$pct_c, use = "complete.obs")))


## =============================================================================
## Section 6: Distribution Comparison
## =============================================================================


resp_long <- responses |>
  filter(has_both) |>
  select(spac_ratio, pyc_abs) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  mutate(variable = case_when(
    variable == "spac_ratio" ~ "Ratio (g/kg C)",
    variable == "pyc_abs"    ~ "Absolute (g/100g soil)"
  ))

ggplot(resp_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~ variable, scales = "free", ncol = 1) +
  labs(title = "Response Variable Distributions",
       x = "Value", y = "Count") +
  theme_minimal()
