# local-cluster-models.R

**File:** `analysis/local-cluster-models.R`
**Added:** 2026-02-06
**Updated:** 2026-02-09 (restyled)

## Purpose

`local-cluster-models.R` tests the hypothesis that spectrally-defined local models outperform a single global model for PyC prediction. The idea is that soil spectra span diverse mineralogies and organic matter contexts — a single model tries to learn one mapping from spectral features to PyC across all of them. Clustering samples into spectrally homogeneous neighborhoods and fitting per-cluster models could reduce this heterogeneity and improve predictions.

The script clusters the spectral matrix using PCA + Euclidean distance + Ward.D2 hierarchical clustering, then runs the full horizons `configure()` → `validate()` → `evaluate()` pipeline independently per cluster. Each cluster gets 72 model configurations (3 models × 6 preprocessing × 2 transforms × 2 feature selection methods) with deeper tuning (grid_size=20, bayesian_iter=15) than the global run.

This is a long-running script — expect several hours per cluster on a local machine with 5 workers.

## Mechanics

1. **Load and extract** (Step 1) — reads the horizons data object, filters to samples with response values, extracts the spectral matrix and wavenumber column names from the role map.

2. **PCA** (Step 2) — `prcomp()` on the centered and scaled spectral matrix. Cumulative variance is computed to determine how many PCs to retain. The 99% variance threshold is used for clustering — this captures the major mineralogical gradients without overfitting to noise in tail PCs.

3. **Euclidean distance** (Step 3) — distance matrix computed on raw PCA scores. Euclidean was chosen over Mahalanobis because Mahalanobis over-weights low-variance components and produced a lopsided split (282/25) in testing. With Euclidean, the high-variance PCs (which correspond to mineralogical gradients) naturally dominate the distance metric.

4. **Hierarchical clustering** (Step 4) — Ward.D2 linkage on the distance matrix. Silhouette scores are computed for k=2 through k=6. k=3 was chosen despite having a lower silhouette (0.275) than k=2 (0.41) because k=2's largest cluster (n=199) mixes France, FORCE, and Jens project samples — too heterogeneous to model as a unit. k=3 splits into ~100 samples per cluster with more homogeneous spectral neighborhoods.

5. **Assign clusters** (Step 5) — cluster labels are added to the analysis data frame. Project labels are joined from modeling_response.csv for diagnostic context. The project composition of each cluster tells you whether the clustering is picking up geographic/methodological differences or genuine spectral structure.

6. **Diagnostic plots** (Step 6) — PCA biplot colored by cluster (visual check that clusters separate in PC space) and silhouette score vs k plot (shows the quality/k trade-off). These are interactive — meant to be inspected in RStudio, not saved.

7. **Subset helper** (Step 7) — `subset_horizons()` is a workaround for horizons lacking a native subset method. It deep-copies the object, filters the analysis data frame to the target sample IDs, resets `n_rows`, and NULLs out all downstream slots (config, validation, evaluation, models, ensemble, artifacts). This gives a clean object that can be configured and evaluated independently.

8. **Per-cluster pipeline** (Step 8) — iterates over clusters. For each:
   - Subsets the horizons object to cluster sample IDs
   - For cluster 3, excludes FORCE project samples (16 samples with essentially no PyC variance — SD=0.006, range 0.072-0.096 — held out for external validation)
   - Skips clusters with fewer than 50 response samples (not enough for reliable train/test + 5-fold CV)
   - Runs `configure()` with the 72-config grid
   - Runs `validate()` with outlier detection disabled (the global pipeline already cleaned outliers)
   - Runs `evaluate()` with RPD metric, pruning at RPD < 1.0, 5 parallel workers, and checkpointing to the output directory
   - Saves per-cluster result objects and collects evaluation results

   The `future::plan(multisession, workers = 5)` at the top works around a horizons bug where the sequential code path never sets up a plan. This means all clusters run with 5 workers for the inner CV parallelism.

9. **Collate** (Step 9) — binds per-cluster results into a single data frame, extracts the best model per cluster, and saves the collated results. The global best RPD (1.70, SVM-RBF) is noted as a comment for comparison.

## Key decisions

- **k=3 over k=2** — silhouette says k=2, domain knowledge says k=3. The k=2 split creates one large heterogeneous cluster mixing three projects with very different soil contexts. k=3 produces clusters that are more interpretable and more likely to benefit from local modeling. This is a case where statistical optimality and modeling utility diverge.

- **FORCE exclusion from cluster 3** — FORCE samples have fine spectra but a flat response (PyC range: 0.072–0.096). Including them in a local model where the response has no variance is pointless. They're better used as an external validation set for the global model.

- **Euclidean over Mahalanobis** — Mahalanobis equalizes variance across PCs, which sounds desirable but in practice over-weights noise. The mineralogical gradients that matter for clustering are in the high-variance PCs. Euclidean preserves that natural weighting.

- **validate(remove_outliers = FALSE)** — outlier removal was already done in `prepare-data.R`. Running it again per-cluster would be overly aggressive (the thresholds would tighten with smaller n, potentially removing legitimate edge-of-cluster samples).

- **grid_size=20, bayesian_iter=15** — deeper tuning than the global run (which used grid_size=10, bayesian_iter=10) because per-cluster sample sizes are smaller and hyperparameter sensitivity is higher.

## Dependencies

| Package | Used for |
|---------|----------|
| dplyr | Data manipulation |
| stringr | String operations |
| ggplot2 | Diagnostic plots |
| cluster | `silhouette()` for cluster quality |
| horizons | `configure()`, `validate()`, `evaluate()` |
| rules | Required for Cubist models via parsnip |
| readr | Reading modeling_response.csv |
| future | Parallel worker setup |
| tibble | Data frame construction |

## Known limitations

1. **Cluster assignments are not deterministic across data changes** — adding or removing samples from `prepare-data.R` will change the PCA scores, which changes the distance matrix, which changes cluster membership. The k=3 structure and project composition could shift entirely.

2. **No cross-cluster validation** — each cluster is evaluated independently. The final comparison (cluster RPDs vs global RPD) doesn't account for the information leakage of choosing k after seeing the data.

3. **subset_horizons() is fragile** — it manually NULLs out slots by name. If horizons adds new slots, this function won't know about them and the leftover state could confuse downstream verbs.

4. **Memory** — each cluster's `evaluate()` call runs 72 configs × 5-fold CV × grid search + Bayesian tuning. With 5 workers, peak memory is substantial. The global run OOM'd at ~310/504 configs; per-cluster runs are smaller but still heavy.

5. **The cluster 3 FORCE exclusion is hardcoded** — it checks `if (cl == 3)` and filters by project name. If re-clustering changes which cluster FORCE lands in, this logic won't follow.

6. **future::plan() side effect** — the script sets `future::plan(multisession, workers = 5)` and resets to sequential at the end, but if the script is interrupted mid-run, the plan persists in the session. Unlike horizons' `evaluate()`, which uses `on.exit()` to restore the previous plan, this script doesn't protect against interruption.
