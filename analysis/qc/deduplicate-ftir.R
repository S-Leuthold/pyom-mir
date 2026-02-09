#' Deduplicate FTIR Files
#'
#' @description
#' Audits the raw FTIR OPUS file inventory. Parses sample stems from filenames,
#' counts replicates per stem, flags anomalies, and produces a manifest CSV.
#'
#' @details
#' Input:  data/raw/Final FTIR/*.0  (OPUS files)
#' Output: data/processed/ftir_manifest.csv
#'
#' Two naming conventions:
#'   Standard: {stem}-S{1-4}_{position}.0   → stem extracted via -S split
#'   FORCE:    {stem}_{position}.0           → stem extracted via position strip
#'
#' This is an audit script. It reports issues but does not delete files.
#'
#' @keywords internal
#' @name deduplicate-ftir


## =============================================================================
## Section 1: Setup
## =============================================================================


library(tidyverse)

ftir_dir <- "data/raw/Final FTIR"


## =============================================================================
## Section 2: Inventory OPUS Files
## =============================================================================


## List only .0 files (OPUS spectra, not companion files)
all_files <- list.files(ftir_dir, pattern = "\\.0$", full.names = FALSE)

cat(sprintf("Total OPUS files (.0): %d\n", length(all_files)))
cat(sprintf("Total items in folder: %d\n",
            length(list.files(ftir_dir, all.files = FALSE))))


## =============================================================================
## Section 3: Parse Sample Stems
## =============================================================================
##
## Mirrors the logic in prepare-horizons-data.R:
##   1. Standard files: strip "-S[1-4]_..." suffix to get stem
##   2. FORCE files: strip trailing "_[A-H]\d+" position code
##   3. Handle " 2" suffix on duplicate scans
##   4. Trim whitespace


inventory <- tibble(filename = all_files) |>
  mutate(
    ## Classify pattern
    pattern = case_when(
      str_detect(filename, "-S[1-4]_") ~ "standard",
      TRUE ~ "force"
    ),

    ## Extract stem
    stem = case_when(
      ## Standard: everything before -S (trim leading space before -S too)
      pattern == "standard" ~ str_remove(filename, "\\s*-S[1-4]_.*$"),
      ## FORCE: strip .0 extension, then strip position code
      TRUE ~ filename |>
        str_remove("\\.0$") |>
        str_remove("_[A-H][0-9]+$") |>
        str_remove(" [0-9]+$")
    ),

    ## Clean up
    stem = str_trim(stem),

    ## Extract replicate number (standard only)
    replicate = case_when(
      pattern == "standard" ~ str_extract(filename, "-S([1-4])", group = 1),
      TRUE ~ NA_character_
    ),

    ## Flag " 2" duplicate scans
    is_rescan = str_detect(filename, " [0-9]+\\.")
  )


## =============================================================================
## Section 4: Replicate Counts
## =============================================================================


stem_summary <- inventory |>
  group_by(stem, pattern) |>
  summarise(
    n_replicates = n(),
    n_rescans    = sum(is_rescan),
    files        = paste(filename, collapse = "; "),
    .groups      = "drop"
  ) |>
  arrange(stem)

cat(sprintf("\nUnique sample stems: %d\n", nrow(stem_summary)))


## ---------------------------------------------------------------------------
## Replicate distribution
## ---------------------------------------------------------------------------

cat("\n=== Replicate Counts ===\n")
stem_summary |>
  count(n_replicates, name = "n_stems") |>
  print()


## ---------------------------------------------------------------------------
## Flag anomalies
## ---------------------------------------------------------------------------

anomalies <- stem_summary |>
  filter(n_replicates != 4) |>
  select(stem, pattern, n_replicates, n_rescans, files)

cat(sprintf("\n=== Anomalous Stems (n != 4 replicates): %d ===\n", nrow(anomalies)))

if (nrow(anomalies) > 0) {
  anomalies |>
    select(stem, n_replicates, n_rescans, pattern) |>
    print(n = Inf)
}


## =============================================================================
## Section 5: Check for Duplicate Stems Across Patterns
## =============================================================================
##
## A stem should appear in only one naming pattern. If it shows up in both
## standard and FORCE, something is wrong.


cross_pattern <- stem_summary |>
  group_by(stem) |>
  filter(n() > 1)

if (nrow(cross_pattern) > 0) {
  cat("\n=== WARNING: Stems appearing in multiple naming patterns ===\n")
  print(cross_pattern)
} else {
  cat("\nNo cross-pattern duplicates found.\n")
}


## =============================================================================
## Section 6: Rescan Audit
## =============================================================================


rescans <- inventory |>
  filter(is_rescan)

cat(sprintf("\n=== Rescans (' 2' suffix files): %d ===\n", nrow(rescans)))

if (nrow(rescans) > 0) {
  rescans |>
    select(filename, stem, pattern) |>
    print(n = Inf)
}


## =============================================================================
## Section 7: Save Manifest
## =============================================================================


manifest <- stem_summary |>
  select(stem, pattern, n_replicates, n_rescans, files) |>
  mutate(
    flag = case_when(
      n_replicates < 4 ~ "incomplete",
      n_replicates > 4 ~ "excess",
      n_rescans > 0    ~ "has_rescans",
      TRUE             ~ "ok"
    )
  )

readr::write_csv(manifest, "data/processed/ftir_manifest.csv")
cat(sprintf("\nManifest saved: data/processed/ftir_manifest.csv (%d stems)\n",
            nrow(manifest)))

## Summary
cat("\n=== Final Summary ===\n")
cat(sprintf("  Total OPUS files:     %d\n", length(all_files)))
cat(sprintf("  Unique stems:         %d\n", nrow(manifest)))
cat(sprintf("  Complete (4 scans):   %d\n", sum(manifest$flag == "ok")))
cat(sprintf("  Incomplete (<4):      %d\n", sum(manifest$flag == "incomplete")))
cat(sprintf("  Excess (>4):          %d\n", sum(manifest$flag == "excess")))
cat(sprintf("  Has rescans:          %d\n", sum(manifest$flag == "has_rescans")))
