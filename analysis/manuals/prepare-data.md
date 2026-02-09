# prepare-data.R

**File:** `analysis/prepare-data.R`
**Added:** 2026-02-04
**Updated:** 2026-02-09 (consolidated from two scripts)

## Purpose

`prepare-data.R` is the entry point for the PyOM-MIR pipeline. It takes the raw HyPy response spreadsheet (436 samples across 7 projects) and the FTIR OPUS spectral files (1,696 scans), cleans and filters the response data, derives the modeling target variable (`pyc_abs`), builds a horizons data object from the spectra, and joins the two. Everything downstream — clustering, modeling, figures — reads from its outputs.

The script was consolidated from two earlier scripts (`prepare-response.R` and `prepare-horizons-data.R`) that were separated during rapid iteration but shared too much context to justify the split.

## Mechanics

1. **Read xlsx** (Step 1) — the HyPy spreadsheet has a non-standard layout: row 1 is annotation, row 2 is blank, row 3 is the actual header. `skip = 2` handles this. `clean_names()` standardizes column names. The `ftir_1` column is renamed but not parsed yet.

2. **Drop header row** — one metadata row from the xlsx sometimes comes through as data (where `project == "Project"`). Filtered out explicitly.

3. **Parse numeric columns** — the xlsx stores everything as text. `as.numeric()` with `suppressWarnings()` coerces the three key columns (pct_c, spac_g_kg_c, spac_toc_pct). NAs from non-numeric strings are expected and handled downstream.

4. **Back-calculate SPAC** — Australia and France projects report SPAC as a percentage of TOC (`spac_toc_pct`) rather than the standard g/kg C. The conversion is `spac_g_kg_c = spac_toc_pct * 10`. The `case_when` prefers the direct value when available.

5. **Extract FTIR stem** — the `ftir_1` column has the full FTIR filename with replicate suffix (e.g., `SampleName-S2_A4.0`). The stem is extracted by stripping `-S[1-4]...` and `_Whole Soil` suffixes. This stem must match what `horizons::parse_ids()` extracts from the OPUS filenames — that's the join key.

6. **Drop standards and references** (Step 2) — three non-soil materials are removed: NDN_GroundBulk (instrument check), BPCAStd (20% biochar standard), CIG_Biochar (reference biochar). These have extreme SPAC values that would dominate any model.

7. **Flag SPAC > TOC** — physically impossible (pyrogenic carbon can't exceed total carbon). Flagged and excluded in the next filter.

8. **Keep complete cases** — only samples with non-NA SPAC, non-NA FTIR stem, and no QC flags proceed.

9. **Litter drop** (Step 3) — litter samples have fundamentally different FTIR signatures from soil (cellulose/lignin peaks dominate over mineral peaks) and 10x higher PyC concentrations. Including them collapsed RPD from 1.26 to 0.48 in testing.

10. **BCInt drop** — biochar amendment plots from the BCInt project have artificially elevated PyC from experimental additions. Different land use context from the fire-affected soils in other projects. Residuals were +0.17 to +0.29 when included.

11. **Named drops** — five individual samples with data errors or uncertain response values: a zero-SPAC sample (AUS_54_Soil) that broke log transforms, a suspicious ratio (7_HTC_5_15_Soil), an extreme outlier (Swiss_F_12_Soil), conflicting measurements (Poudre_Watershed), and a mislabeled depth (2_LTC_5_15_Soil).

12. **Spectral outlier drops** — 18 samples flagged by Mahalanobis distance on PC1:10 exceeding the chi-squared threshold (df=10, p=0.975). These are deep subsoils, corrupted scans, or atypical mineralogy that create high leverage in models. The Mahalanobis values range from 22 to 282 and are documented inline.

13. **Derive pyc_abs** (Step 4) — the modeling target is absolute pyrogenic carbon: `pyc_abs = pct_c * spac_g_kg_c / 1000`, giving g PyC per 100g soil. This is easier to predict from MIR than the ratio (g/kg C) because spectra see total C directly (PC1 correlates with total C at r = 0.54, so pct_c is partially encoded in the spectral signal).

14. **Trim response outliers** — samples beyond mean + 2 SD of `pyc_abs` are removed. The remaining high-PyC samples (e.g., Swiss_F_6_Soil at 2.43) dominate the loss function. Trimming keeps the distribution learnable.

15. **Write response CSV** (Step 5) — `modeling_response.csv` with columns Sample_ID, spac_g_kg_c, pyc_abs, pct_c, project. Deduplicated on Sample_ID (first occurrence kept). This CSV is the join source for the horizons object and is also used by `local-cluster-models.R` for project labels.

16. **Read OPUS files** (Step 6) — `horizons::spectra()` reads all 1,696 OPUS files from `data/raw/Final FTIR/`. `standardize()` resamples to 2 cm⁻¹ resolution and trims to 600–4000 cm⁻¹.

17. **Parse sample IDs** — `horizons::parse_ids()` extracts sample stems from OPUS filenames. The standard naming convention is `{stem}-S{1-4}_{position}.0`, but FORCE project samples use `{stem}_{position}.0` without the `-S` replicate marker. `too_few = "keep_original"` handles the fallback, then FORCE stems are cleaned up by stripping trailing position codes (`_A8`, `_B3`) and duplicate-scan suffixes (` 2`).

18. **Average replicates** — `horizons::average()` groups by `sample_id` and averages spectral replicates. `quality_check = TRUE` with `correlation_threshold = 0.996` flags replicates with low inter-scan agreement. Poorly correlating replicates are averaged anyway but the flag is stored.

19. **Join response** — `horizons::add_response()` joins `modeling_response.csv` to the horizons object by sample_id = Sample_ID. Samples without a response match remain in the object (they're excluded at modeling time by the `has_response` filter).

20. **Save** — `horizons_data.rds`, the single input for all downstream scripts.

## Key decisions

- **pyc_abs over spac_g_kg_c** — the absolute concentration (g/100g soil) was chosen over the ratio (g SPAC / kg C) because MIR spectra partially encode total carbon. The ratio strips that information out, forcing the model to predict a signal it can't see. RPD improved from ~1.2 to ~1.5 when switching.

- **ALS baseline correction disabled** — tested and removed. SNV and SG derivatives through the horizons preprocessing pipeline handle baseline variation more effectively. RPD was 1.46 with ALS vs 1.51 without.

- **Mean + 2SD trim** — standard outlier rule on the response. More aggressive than 3SD but justified because the response distribution is right-skewed and a few extreme values (Swiss forest soils with high organic C) dominated model training.

- **Mahalanobis threshold at chi-sq(0.975, df=10)** — conservative but necessary. The 18 dropped samples are genuine outliers (deep subsoils hitting bedrock, corrupted OPUS files, atypical soil mineralogy). The threshold was validated by inspecting each sample's spectral profile against the population.

## Dependencies

| Package | Used for |
|---------|----------|
| tidyverse | Core data manipulation and piping |
| readxl | Reading the HyPy xlsx |
| janitor | `clean_names()` for column standardization |
| horizons | `spectra()`, `standardize()`, `parse_ids()`, `average()`, `add_response()` |
| stringr | Regex operations for stem extraction and cleanup |

## Known limitations

1. **Hardcoded file paths** — the xlsx filename includes a date stamp (`2025_12_09`). If the source data is updated, the path must change. Same for the OPUS directory.

2. **FORCE naming convention is fragile** — the cleanup regex (`_[A-H][0-9]+$` and ` [0-9]+$`) is tuned to the specific FORCE naming pattern. A new project with different position codes would need its own cleanup rule.

3. **Deduplication is first-occurrence** — when multiple rows share a `ftir_stem`, only the first is kept. No attempt to choose the "best" measurement or average conflicting values.

4. **The 2SD trim is global** — it doesn't account for project-level differences in PyC distribution. A legitimate high-PyC sample from a fire-affected project could be trimmed if it exceeds the global cutoff.

5. **No provenance tracking** — the script doesn't record which samples were dropped or why in the output files. The decisions are documented in comments and the changelog, but the output CSV has no audit trail. You'd need to re-run with breakpoints to reconstruct the filtering steps.
