#' Test Baseline Correction Effect on PyC Signal
#'
#' @description
#' Applies several baseline correction methods to the existing spectra and
#' compares correlation spectra before/after. Goal: does baseline correction
#' expose hidden PyC signal?
#'
#' @keywords internal
#' @name test-baseline-correction


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
  left_join(resp |> select(Sample_ID, project),
            by = c("sample_id" = "Sample_ID"))

spec_matrix <- analysis |>
  select(all_of(wn_cols)) |>
  as.matrix()

spac <- analysis$spac_g_kg_c

cat(sprintf("Samples: %d, Wavenumbers: %d\n", nrow(spec_matrix), ncol(spec_matrix)))


## =============================================================================
## Section 2: Apply Baseline Corrections
## =============================================================================


## Wavenumbers need to be in ascending order for prospectr
wn_ascending <- rev(wn_values)
spec_ascending <- spec_matrix[, rev(seq_len(ncol(spec_matrix)))]

cat("\nApplying baseline corrections...\n")

## Method 1: Convex hull (rubber band) — what horizons uses
cat("  Convex hull (prospectr::baseline)...")
bc_hull <- prospectr::baseline(X = spec_ascending, wav = wn_ascending)
## Reverse back to descending order
bc_hull <- bc_hull[, rev(seq_len(ncol(bc_hull)))]
cat(" done\n")

## Method 2: Rolling ball
cat("  Rolling ball...")
bc_rb <- baseline::baseline(spec_ascending, method = "rollingBall",
                            wm = 200, ws = 200)
bc_rb_corrected <- baseline::getCorrected(bc_rb)
bc_rb_corrected <- bc_rb_corrected[, rev(seq_len(ncol(bc_rb_corrected)))]
cat(" done\n")

## Method 3: ALS (Asymmetric Least Squares)
cat("  ALS...")
bc_als <- baseline::baseline(spec_ascending, method = "als",
                             lambda = 6, p = 0.05)
bc_als_corrected <- baseline::getCorrected(bc_als)
bc_als_corrected <- bc_als_corrected[, rev(seq_len(ncol(bc_als_corrected)))]
cat(" done\n")


## =============================================================================
## Section 3: Correlation Spectra — Before vs After
## =============================================================================


## Compute correlation with SPAC for each method
cor_raw  <- apply(spec_matrix, 2, cor, y = spac, use = "complete.obs")
cor_hull <- apply(bc_hull, 2, cor, y = spac, use = "complete.obs")
cor_rb   <- apply(bc_rb_corrected, 2, cor, y = spac, use = "complete.obs")
cor_als  <- apply(bc_als_corrected, 2, cor, y = spac, use = "complete.obs")

cor_compare <- tibble(
  wavenumber = rep(wn_values, 4),
  correlation = c(cor_raw, cor_hull, cor_rb, cor_als),
  method = rep(c("Raw", "Convex Hull", "Rolling Ball", "ALS"),
               each = length(wn_values))
) |>
  mutate(method = factor(method, levels = c("Raw", "Convex Hull",
                                            "Rolling Ball", "ALS")))

ggplot(cor_compare, aes(x = wavenumber, y = correlation, color = method)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("rect", xmin = 700, xmax = 900, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.05) +
  annotate("rect", xmin = 1550, xmax = 1650, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.05) +
  scale_x_reverse() +
  labs(title = "Correlation Spectrum: Raw vs Baseline-Corrected",
       subtitle = "Does baseline correction expose PyC signal?",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Pearson r with SPAC") +
  theme_minimal()

## Faceted version for clarity
ggplot(cor_compare, aes(x = wavenumber, y = correlation)) +
  geom_line(linewidth = 0.4, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("rect", xmin = 700, xmax = 900, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  annotate("rect", xmin = 1550, xmax = 1650, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  scale_x_reverse() +
  facet_wrap(~ method, ncol = 1) +
  labs(title = "Correlation Spectrum by Baseline Method",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Pearson r with SPAC") +
  theme_minimal()

## Report max correlations per method
cat("\n=== Max |r| by Method ===\n")
cor_compare |>
  group_by(method) |>
  summarise(
    max_abs_r = max(abs(correlation), na.rm = TRUE),
    best_wn   = wavenumber[which.max(abs(correlation))],
    .groups   = "drop"
  ) |>
  print()


## =============================================================================
## Section 4: Difference Spectra — Before vs After
## =============================================================================


spac_quartile <- cut(spac,
                     breaks = quantile(spac, probs = c(0, 0.25, 0.75, 1)),
                     labels = c("Q1", "mid", "Q4"),
                     include.lowest = TRUE)

q4_idx <- spac_quartile == "Q4"
q1_idx <- spac_quartile == "Q1"

diff_raw  <- colMeans(spec_matrix[q4_idx, ]) - colMeans(spec_matrix[q1_idx, ])
diff_hull <- colMeans(bc_hull[q4_idx, ]) - colMeans(bc_hull[q1_idx, ])
diff_rb   <- colMeans(bc_rb_corrected[q4_idx, ]) - colMeans(bc_rb_corrected[q1_idx, ])
diff_als  <- colMeans(bc_als_corrected[q4_idx, ]) - colMeans(bc_als_corrected[q1_idx, ])

diff_compare <- tibble(
  wavenumber = rep(wn_values, 4),
  diff = c(diff_raw, diff_hull, diff_rb, diff_als),
  method = rep(c("Raw", "Convex Hull", "Rolling Ball", "ALS"),
               each = length(wn_values))
) |>
  mutate(method = factor(method, levels = c("Raw", "Convex Hull",
                                            "Rolling Ball", "ALS")))

ggplot(diff_compare, aes(x = wavenumber, y = diff)) +
  geom_line(linewidth = 0.4, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("rect", xmin = 700, xmax = 900, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  annotate("rect", xmin = 1550, xmax = 1650, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  scale_x_reverse() +
  facet_wrap(~ method, ncol = 1, scales = "free_y") +
  labs(title = "Difference Spectrum (Q4 - Q1) by Baseline Method",
       subtitle = "Positive = more absorbance in high-SPAC samples",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Absorbance Difference") +
  theme_minimal()


## =============================================================================
## Section 5: Example Spectra — Before vs After
## =============================================================================
##
## Show a few individual spectra before/after correction so we can see
## what the baseline removal actually does.


example_idx <- c(
  which.min(spac),                          # lowest SPAC
  which.min(abs(spac - median(spac))),      # median SPAC
  which.max(spac)                           # highest SPAC
)

example_labels <- sprintf("%s (SPAC=%.0f)",
                          analysis$sample_id[example_idx],
                          spac[example_idx])

example_raw <- spec_matrix[example_idx, ]
example_hull <- bc_hull[example_idx, ]

example_df <- bind_rows(
  tibble(
    sample = rep(example_labels, each = length(wn_values)),
    wavenumber = rep(wn_values, length(example_idx)),
    absorbance = as.vector(t(example_raw)),
    type = "Raw"
  ),
  tibble(
    sample = rep(example_labels, each = length(wn_values)),
    wavenumber = rep(wn_values, length(example_idx)),
    absorbance = as.vector(t(example_hull)),
    type = "Convex Hull Corrected"
  )
)

ggplot(example_df, aes(x = wavenumber, y = absorbance, color = type)) +
  geom_line(linewidth = 0.4) +
  scale_x_reverse() +
  facet_wrap(~ sample, ncol = 1, scales = "free_y") +
  labs(title = "Example Spectra: Raw vs Baseline-Corrected",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Absorbance") +
  theme_minimal() +
  theme(legend.position = "top")


## =============================================================================
## Section 6: PCA After Baseline Correction
## =============================================================================
##
## Does baseline correction change the PCA structure?
## Are projects still clustering, or does the correction remove that?


## Drop zero-variance columns (created by baseline correction at edges)
nonzero_var <- apply(bc_hull, 2, var) > 0
bc_hull_pca <- bc_hull[, nonzero_var]
cat(sprintf("Dropped %d zero-variance columns for PCA\n",
            sum(!nonzero_var)))

pca_hull <- prcomp(bc_hull_pca, center = TRUE, scale. = TRUE)

pca_scores_hull <- tibble(
  PC1     = pca_hull$x[, 1],
  PC2     = pca_hull$x[, 2],
  spac    = spac,
  project = analysis$project
)

var_expl_hull <- summary(pca_hull)$importance[2, 1:3] * 100

## PCA by project — after correction
ggplot(pca_scores_hull, aes(x = PC1, y = PC2, color = project)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA After Convex Hull Baseline Correction — by Project",
       subtitle = sprintf("PC1: %.1f%%, PC2: %.1f%%",
                          var_expl_hull[1], var_expl_hull[2])) +
  theme_minimal()

## PCA by SPAC — after correction
ggplot(pca_scores_hull, aes(x = PC1, y = PC2, color = spac)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_viridis_c(option = "C") +
  labs(title = "PCA After Convex Hull Baseline Correction — by SPAC",
       subtitle = sprintf("PC1: %.1f%%, PC2: %.1f%%",
                          var_expl_hull[1], var_expl_hull[2]),
       color = "SPA-C\n(g/kg C)") +
  theme_minimal()

cat(sprintf("\nPCA-SPAC correlations after convex hull correction:\n"))
cat(sprintf("  PC1 vs SPAC: r = %.3f (was %.3f)\n",
            cor(pca_scores_hull$PC1, spac), 0.170))
cat(sprintf("  PC2 vs SPAC: r = %.3f (was %.3f)\n",
            cor(pca_scores_hull$PC2, spac), -0.165))
