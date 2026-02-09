#' Assess Spectral Quality and PyC Signal
#'
#' @description
#' Visual assessment of MIR spectra to understand where (if anywhere) the PyC
#' signal lives, whether baseline differences are driving project clustering,
#' and what preprocessing might help.
#'
#' @keywords internal
#' @name assess-spectra


## =============================================================================
## Section 1: Setup
## =============================================================================


library(tidyverse)

hd   <- readRDS("data/processed/horizons_data.rds")
resp <- readr::read_csv("data/processed/modeling_response.csv",
                        show_col_types = FALSE)

## Extract spectral matrix and metadata
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
proj <- analysis$project

## SPAC quartiles for grouping
spac_quartile <- cut(spac,
                     breaks = quantile(spac, probs = c(0, 0.25, 0.5, 0.75, 1)),
                     labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"),
                     include.lowest = TRUE)

cat(sprintf("Samples: %d, Wavenumbers: %d\n", nrow(spec_matrix), ncol(spec_matrix)))
cat(sprintf("SPAC range: %.1f - %.1f\n", min(spac), max(spac)))


## =============================================================================
## Section 2: Mean Spectra by Project
## =============================================================================
##
## Are there visible baseline differences between projects?


project_means <- analysis |>
  select(project, all_of(wn_cols)) |>
  group_by(project) |>
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(-project, names_to = "wn_col", values_to = "absorbance") |>
  mutate(wavenumber = as.numeric(str_remove(wn_col, "^wn_")))

ggplot(project_means, aes(x = wavenumber, y = absorbance, color = project)) +
  geom_line(linewidth = 0.6) +
  scale_x_reverse() +
  labs(title = "Mean Spectrum by Project",
       subtitle = "Baseline differences between projects?",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Absorbance") +
  theme_minimal()


## =============================================================================
## Section 3: Spaghetti Plot by Project (subsample)
## =============================================================================
##
## Individual spectra colored by project — shows spread within and between.


set.seed(42)
n_per_project <- 15

subsample_idx <- analysis |>
  mutate(row_idx = row_number()) |>
  group_by(project) |>
  slice_sample(n = n_per_project) |>
  pull(row_idx)

spec_long_sub <- tibble(
  sample_id = analysis$sample_id[subsample_idx],
  project   = proj[subsample_idx]
) |>
  bind_cols(as_tibble(spec_matrix[subsample_idx, ])) |>
  pivot_longer(-c(sample_id, project),
               names_to = "wn_col", values_to = "absorbance") |>
  mutate(wavenumber = as.numeric(str_remove(wn_col, "^wn_")))

ggplot(spec_long_sub, aes(x = wavenumber, y = absorbance,
                          group = sample_id, color = project)) +
  geom_line(alpha = 0.4, linewidth = 0.3) +
  scale_x_reverse() +
  labs(title = "Individual Spectra by Project (15 per project)",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Absorbance") +
  theme_minimal()


## =============================================================================
## Section 4: Mean Spectra by SPAC Quartile
## =============================================================================
##
## Does high vs low PyC look different spectrally?


quartile_means <- tibble(quartile = spac_quartile) |>
  bind_cols(as_tibble(spec_matrix)) |>
  group_by(quartile) |>
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(-quartile, names_to = "wn_col", values_to = "absorbance") |>
  mutate(wavenumber = as.numeric(str_remove(wn_col, "^wn_")))

ggplot(quartile_means, aes(x = wavenumber, y = absorbance, color = quartile)) +
  geom_line(linewidth = 0.6) +
  scale_x_reverse() +
  scale_color_viridis_d(option = "C") +
  labs(title = "Mean Spectrum by SPAC Quartile",
       subtitle = "Q1 = lowest PyC, Q4 = highest PyC",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Absorbance") +
  theme_minimal()


## =============================================================================
## Section 5: Difference Spectrum (High - Low SPAC)
## =============================================================================
##
## Where is the PyC signal? Subtract mean of bottom quartile from top quartile.


mean_q4 <- colMeans(spec_matrix[spac_quartile == "Q4 (high)", ], na.rm = TRUE)
mean_q1 <- colMeans(spec_matrix[spac_quartile == "Q1 (low)", ], na.rm = TRUE)
diff_spec <- mean_q4 - mean_q1

diff_df <- tibble(
  wavenumber = wn_values,
  diff       = diff_spec
)

ggplot(diff_df, aes(x = wavenumber, y = diff)) +
  geom_line(color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  ## PyC regions
  annotate("rect", xmin = 700, xmax = 900, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  annotate("rect", xmin = 1550, xmax = 1650, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  annotate("text", x = 800, y = max(diff_spec) * 0.9,
           label = "Aromatic C-H", size = 2.5, color = "red") +
  annotate("text", x = 1600, y = max(diff_spec) * 0.9,
           label = "Aromatic C=C", size = 2.5, color = "red") +
  scale_x_reverse() +
  labs(title = "Difference Spectrum: Q4 (high SPAC) minus Q1 (low SPAC)",
       subtitle = "Positive = more absorbance in high-PyC samples",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Absorbance Difference") +
  theme_minimal()


## =============================================================================
## Section 6: Correlation Spectrum
## =============================================================================
##
## Pearson correlation of each wavenumber with SPAC.
## Shows which spectral regions carry predictive information.


wn_cor <- apply(spec_matrix, 2, cor, y = spac, use = "complete.obs")

cor_df <- tibble(
  wavenumber   = wn_values,
  correlation  = wn_cor
)

ggplot(cor_df, aes(x = wavenumber, y = correlation)) +
  geom_line(color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  ## PyC regions
  annotate("rect", xmin = 700, xmax = 900, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  annotate("rect", xmin = 1550, xmax = 1650, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) +
  scale_x_reverse() +
  labs(title = "Correlation Spectrum: r(wavenumber, SPAC)",
       subtitle = "Which spectral regions predict pyrogenic carbon?",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Pearson r") +
  theme_minimal()

## Report strongest correlations
cat("\n=== Strongest Correlations with SPAC ===\n")
cor_df |>
  arrange(desc(abs(correlation))) |>
  head(20) |>
  print()

cat(sprintf("\nMax |r|: %.3f at %.0f cm-1\n",
            max(abs(wn_cor)), wn_values[which.max(abs(wn_cor))]))


## =============================================================================
## Section 7: Zoom into PyC Regions
## =============================================================================
##
## 700-900 cm-1 (aromatic C-H out-of-plane bending)
## 1500-1700 cm-1 (aromatic C=C stretching)


spec_long_quartile <- tibble(
  sample_id = analysis$sample_id,
  quartile  = spac_quartile
) |>
  bind_cols(as_tibble(spec_matrix)) |>
  pivot_longer(-c(sample_id, quartile),
               names_to = "wn_col", values_to = "absorbance") |>
  mutate(wavenumber = as.numeric(str_remove(wn_col, "^wn_")))

## Aromatic C-H region
ggplot(spec_long_quartile |> filter(wavenumber >= 700, wavenumber <= 900),
       aes(x = wavenumber, y = absorbance, group = sample_id, color = quartile)) +
  geom_line(alpha = 0.15, linewidth = 0.3) +
  scale_x_reverse() +
  scale_color_viridis_d(option = "C") +
  labs(title = "Aromatic C-H Region (700-900 cm-1) by SPAC Quartile",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Absorbance") +
  theme_minimal()

## Aromatic C=C region
ggplot(spec_long_quartile |> filter(wavenumber >= 1500, wavenumber <= 1700),
       aes(x = wavenumber, y = absorbance, group = sample_id, color = quartile)) +
  geom_line(alpha = 0.15, linewidth = 0.3) +
  scale_x_reverse() +
  scale_color_viridis_d(option = "C") +
  labs(title = "Aromatic C=C Region (1500-1700 cm-1) by SPAC Quartile",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Absorbance") +
  theme_minimal()


## =============================================================================
## Section 8: Baseline Shape Assessment
## =============================================================================
##
## Compare baseline shape across projects.
## Take the minimum absorbance in each of several windows to trace the baseline.


baseline_windows <- list(
  c(650, 700), c(900, 1000), c(1200, 1300),
  c(1800, 1900), c(2400, 2500), c(3800, 3950)
)

get_baseline_point <- function(spec_row, wn_vals, window) {
  idx <- wn_vals >= window[1] & wn_vals <= window[2]
  if (sum(idx) == 0) return(NA_real_)
  min(spec_row[idx], na.rm = TRUE)
}

baseline_points <- map_dfr(seq_len(nrow(spec_matrix)), function(i) {
  pts <- map_dbl(baseline_windows, ~ get_baseline_point(spec_matrix[i, ], wn_values, .x))
  tibble(
    sample_id  = analysis$sample_id[i],
    project    = proj[i],
    spac       = spac[i],
    wn         = map_dbl(baseline_windows, ~ mean(.x)),
    baseline   = pts
  )
})

ggplot(baseline_points, aes(x = wn, y = baseline, group = sample_id, color = project)) +
  geom_line(alpha = 0.3, linewidth = 0.3) +
  scale_x_reverse() +
  labs(title = "Approximate Baseline Shape by Project",
       subtitle = "Min absorbance in 6 spectral windows — traces the baseline",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Baseline Absorbance") +
  theme_minimal()

## Mean baseline by project
baseline_points |>
  group_by(project, wn) |>
  summarise(mean_baseline = mean(baseline, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = wn, y = mean_baseline, color = project)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_reverse() +
  labs(title = "Mean Baseline Shape by Project",
       subtitle = "Do projects have systematically different baselines?",
       x = expression(Wavenumber~(cm^{-1})),
       y = "Mean Baseline Absorbance") +
  theme_minimal()
