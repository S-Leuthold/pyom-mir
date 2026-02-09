# build-figures.R

**File:** `analysis/build-figures.R`
**Added:** 2026-02-06
**Updated:** 2026-02-09 (merged from build-report-figures.R + build-hybrid-figure.R)

## Purpose

`build-figures.R` generates all figures for the PyOM-MIR preliminary results. It reads saved model evaluation results (from the global HPC run and per-cluster local runs), reproduces the spectral clustering, and produces six figures covering model performance, spectral structure, and prediction diagnostics.

The script has two runtime profiles. Steps 1–6 are fast: they load pre-computed results and plot. Steps 7–9 refit models with leave-one-out cross-validation to produce per-sample predictions for the hybrid figure — this takes hours. You can run the fast section independently.

## Mechanics

1. **Load and cluster** (Step 1) — reads the horizons data object, reproduces the PCA + Ward.D2 clustering (same pipeline as `local-cluster-models.R`), and joins project labels. This must exactly reproduce the cluster assignments from the modeling run — same 99% variance threshold, same Euclidean + Ward.D2, same k=3.

2. **Load evaluation results** (Step 2) — reads the global overnight checkpoint (`eval_checkpoint.rds`, 310 configs) and three per-cluster result objects. Cluster 3 uses `eval_checkpoint.rds` instead of `result.rds` because the cluster 3 run didn't fully complete. `parse_config()` extracts model type, preprocessing, and transform from config IDs using pattern matching (can't split on underscores because model names like `svm_rbf` contain them).

3. **Fig 1: Global RPD distribution** (Step 3) — boxplot + jitter of RPD by model type across 203 successful global configs. Red diamonds mark the best per model type. Dashed line at RPD = 1.5 (minimum screening threshold). Shows that SVM-RBF and Cubist dominate, with XGBoost and Random Forest clustered lower.

4. **Fig 2: Spectral clusters** (Step 4) — two-panel figure (patchwork). Panel A: PCA biplot colored by cluster assignment. Panel B: same biplot colored by PyC concentration. Shows whether the spectral clusters (defined by distance in PC space) correspond to PyC gradients or are orthogonal to them.

5. **Fig 3: Local vs global RPD** (Step 5) — boxplot + jitter of per-cluster RPD distributions with the global best RPD as a blue dashed reference line. Shows cluster 2 (Australia + Jens) crossing the RPD 2.0 threshold, cluster 1 performing near global, and cluster 3 underperforming.

6. **Fig 4: Concentration floor** (Step 6) — the key diagnostic. A 5-fold CV refit of the global best model (SVM-RBF + SNV + log) produces per-sample predictions. Panel A: absolute error vs observed PyC. Panel B: relative error vs observed PyC. Both show a floor around PyC = 0.05 g/100g soil where relative error explodes above 100%. This is the practical detection limit — the model can't distinguish PyC signals below this concentration from spectral noise.

   The 5-fold CV uses `rsample::vfold_cv()` with stratification on pyc_abs. The SNV preprocessing is applied manually (not through horizons) because horizons doesn't expose a predict pathway. The SVM hyperparameters (cost=2.5, rbf_sigma=0.05) are hardcoded from the best global config.

7. **Fig 5: Global pred vs obs** (Step 6 continued) — 1:1 plot from the same 5-fold CV predictions. Reports R², RMSE, RPD in the subtitle. Colored by cluster to show spatial patterns in prediction quality.

8. **Preprocessing helpers** (Step 7) — manual implementations of the spectral preprocessing that horizons applies internally. These exist because horizons has no `predict()` function, so producing per-sample predictions requires reconstructing the pipeline outside horizons:
   - `snv()` — standard normal variate (row-wise centering and scaling)
   - `sg_smooth()` — Savitzky-Golay smoothing (order 0, polynomial 2, window 11)
   - `deriv1()` — first derivative via SG (order 1)
   - `snv_deriv1()` — SNV then first derivative (composition)
   - `cor_select()` — top-N predictors by absolute correlation with response
   - `pca_reduce()` — PCA to 99% variance, applied to train, projected onto test

   `loo_cv()` is a generic LOO helper: takes a data frame, wavenumber columns, and function arguments for preprocessing, feature selection, and model fitting. Returns a tibble of observed vs predicted values.

9. **Per-cluster LOO-CV** (Step 8) — three LOO-CV runs using the best model/preprocessing/feature_selection combo identified by `evaluate()` for each cluster:
   - Cluster 1: Cubist + SNV+D1 + log + correlation (committees=92, neighbors=1, max_rules=396)
   - Cluster 2: Elastic Net + D1 + log + correlation (penalty=6.16e-05, mixture=0.2)
   - Cluster 3: Global SVM-RBF + SG + log + PCA (cost=2.56, rbf_sigma=0.048) — uses the global model, not a local one, because cluster 3 performed worse locally

   Cluster 3 runs LOO on the full dataset (all 307 samples) and extracts only the cluster 3 predictions. This matches how the hybrid approach would work in practice: assign an incoming sample to a cluster, then route it to the appropriate model.

10. **Fig 6: Hybrid pred vs obs** (Step 9) — pools predictions from all three clusters (local for 1-2, global for 3) and plots 1:1 with per-cluster coloring. Reports pooled metrics (RPD, R², RMSE, CCC) and per-cluster breakdown. This is the main result figure: it shows the hybrid approach's combined performance.

## Key decisions

- **Reproduce clustering rather than loading cluster labels** — the clustering is reproduced from scratch rather than loading saved labels. This ensures the figure script is self-contained and doesn't depend on intermediate files that might not exist. The trade-off is that any change to `prepare-data.R` could change the clustering, making previously saved model results inconsistent with the figures.

- **LOO-CV over k-fold for the hybrid figure** — leave-one-out gives one prediction per sample with no randomness in the split, making the results deterministic and maximizing the number of predictions. The downside is runtime: 307 × 3 model fits for cluster 3.

- **Hardcoded hyperparameters** — the SVM, Cubist, and Elastic Net hyperparameters are hardcoded from the horizons evaluation results. If the models are re-evaluated with different tuning, these values must be manually updated.

- **5-fold CV for figs 4-5, LOO for fig 6** — different CV strategies for different figures. 5-fold is faster and sufficient for the concentration floor analysis (we care about the pattern, not exact per-sample values). LOO is necessary for the hybrid figure because pooling per-sample predictions from three different models requires every sample to have exactly one prediction.

## Dependencies

| Package | Used for |
|---------|----------|
| dplyr, tidyr, stringr | Data manipulation |
| ggplot2 | All figures |
| patchwork | Multi-panel layouts (figs 2, 4) |
| parsnip | Model specs (SVM, Cubist, Elastic Net) |
| rsample | `vfold_cv()`, `analysis()`, `assessment()` for CV |
| rules | Required for Cubist models via parsnip |
| readr | Reading modeling_response.csv |
| prospectr | `savitzkyGolay()` for SG smoothing and derivatives |

## Known limitations

1. **Cluster 3 uses eval_checkpoint.rds** — the cluster 3 run didn't complete. The checkpoint contains partial results. If the run is completed later, the path should change to `result.rds` like clusters 1 and 2.

2. **parse_config() is regex-based** — it pattern-matches model names from config_id strings. If horizons changes the config_id format, this breaks silently (model defaults to "Other").

3. **No horizons predict()** — the entire LOO-CV section (Steps 7-9) exists because horizons doesn't export a predict function. The preprocessing helpers manually replicate what horizons does internally. If horizons' internal preprocessing changes (e.g., different SG window size), these helpers will produce different results than what was evaluated.

4. **Runtime** — Steps 7-9 take hours. Cluster 3's LOO-CV runs 307 model fits (full dataset LOO). There's no checkpointing — if interrupted, the entire LOO run must restart.

5. **Memory for LOO** — each LOO iteration fits a full SVM or Cubist model. 307 iterations × model objects that aren't cleaned up can accumulate. R's garbage collector handles this reasonably for parsnip models, but it's worth monitoring.

6. **Figure filenames are hardcoded** — all six PNGs are saved to `output/figures/` with fixed names. Re-running the script overwrites previous figures without warning.
