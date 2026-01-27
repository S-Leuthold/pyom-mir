#' Build Sample-FTIR Lookup Table
#'
#' @description
#' Creates a canonical mapping between samples and FTIR spectra files.
#'
#' @details
#' Output: data/processed/sample_ftir_lookup.csv
#'
#' ID scheme:
#'   pyom-0001    = sample (links to SPAC value)
#'   pyom-0001-s1 = spectrum replicate 1 (derived at runtime, not stored)
#'
#' Matching is performed in order:
#'   1. Explicit FTIR ID from cleaned data (if present)
#'   2. Project-specific transforms (Australia, Jens, France)
#'   3. Generic normalized matching (fallback)
#'
#' @keywords internal
#' @name build-sample-lookup


## =============================================================================
## Section 1: Setup
## =============================================================================


library(tidyverse)


## =============================================================================
## Section 2: Helper Functions
## =============================================================================


## ---------------------------------------------------------------------------
## normalize — Standardize string for fuzzy matching
## ---------------------------------------------------------------------------

#' Normalize string for matching
#'
#' @description
#' Lowercases, replaces whitespace/hyphens with underscores, removes _soil suffix.
#'
#' @param x Character. String to normalize.
#'
#' @return Character. Normalized string.
#'
#' @noRd
normalize <- function(x) {

  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("\\s+", "_") %>%
    stringr::str_replace_all("-", "_") %>%
    stringr::str_remove("_soil$") %>%
    stringr::str_remove("_+$") %>%
    stringr::str_trim() -> result

  result

}


## ---------------------------------------------------------------------------
## find_ftir_match — Project-specific matching logic
## ---------------------------------------------------------------------------

#' Find FTIR stem for a sample
#'
#' @description
#' Attempts to match a sample to its FTIR spectrum files using project-specific
#' naming conventions, falling back to normalized string matching.
#'
#' @param sample_name Character. Original sample name.
#' @param project Character. Project identifier.
#' @param category_1 Character. Category 1 value (unused, kept for extensibility).
#' @param ftir_stems Character vector. Available FTIR stem names.
#'
#' @return List with `stem` (matched FTIR stem or NA) and `method` (how matched).
#'
#' @noRd
find_ftir_match <- function(sample_name, project, category_1, ftir_stems) {

  ## Handle NA/NULL inputs ----------------------------------------------------

  if (is.na(sample_name) || is.null(sample_name) || sample_name == "") {

    return(list(stem = NA_character_, method = "unmatched"))

  }

  ## Australia: sample "1" → "AUS_1_Soil" -------------------------------------

  if (!is.na(project) && project == "Australia") {

    candidate <- paste0("AUS_", sample_name, "_Soil")

    if (candidate %in% ftir_stems) {

      return(list(stem = candidate, method = "auto_australia"))

    }

  }

  ## Jens: "1 HBC 5-15" → "1_HBC_5_15_Soil" -----------------------------------

  if (!is.na(project) && project == "Jens") {

    sample_name %>%
      stringr::str_replace_all(" ", "_") %>%
      stringr::str_replace_all("-", "_") %>%
      paste0("_Soil") -> candidate

    if (candidate %in% ftir_stems) {

      return(list(stem = candidate, method = "auto_jens"))

    }

  }

  ## France: "Can_F_EUP1_0-15" → "Can_F_EUP_1_0_15_Soil" ----------------------

  if (!is.na(project) && project == "France") {

    sample_name %>%
      stringr::str_replace("([A-Za-z])(\\d)", "\\1_\\2") %>%
      stringr::str_replace_all("-", "_") %>%
      paste0("_Soil") -> candidate

    if (candidate %in% ftir_stems) {

      return(list(stem = candidate, method = "auto_france"))

    }

  }

  ## Generic normalized matching (fallback) -----------------------------------

  sample_norm <- normalize(sample_name)

  for (stem in ftir_stems) {

    if (normalize(stem) == sample_norm) {

      return(list(stem = stem, method = "auto_normalized"))

    }

  }

  ## No match found -----------------------------------------------------------

  list(stem = NA_character_, method = "unmatched")

}


## =============================================================================
## Section 3: Load Data
## =============================================================================


## ---------------------------------------------------------------------------
## Load cleaned sample data
## ---------------------------------------------------------------------------

samples <- readr::read_csv(
  "data/processed/hypy_clean.csv",
  show_col_types = FALSE
)


## ---------------------------------------------------------------------------
## Build FTIR stem inventory from file system
## ---------------------------------------------------------------------------

ftir_files <- list.files("data/raw/Final FTIR", full.names = FALSE)

ftir_files %>%
  stringr::str_remove("-S[1-4]_[A-H]\\d+\\.0$") %>%
  stringr::str_remove("_[A-H]\\d+\\.0$") %>%
  stringr::str_trim() %>%
  unique() %>%
  sort() -> ftir_stems


## =============================================================================
## Section 4: Match Samples to FTIR
## =============================================================================


## ---------------------------------------------------------------------------
## Process each sample
## ---------------------------------------------------------------------------

results <- vector("list", nrow(samples))

for (i in seq_len(nrow(samples))) {

  row     <- samples[i, ]
  ftir_id <- as.character(row$ftir_id)

  ## Check if explicit FTIR ID exists ---------------------------------------

  has_explicit <- !is.na(ftir_id) &&
                  ftir_id != "NA" &&
                  ftir_id != "" &&
                  nchar(trimws(ftir_id)) > 0

  if (has_explicit) {

    results[[i]] <- list(stem = ftir_id, method = "explicit")

  } else {

    results[[i]] <- find_ftir_match(
      sample_name = as.character(row$sample_name_orig),
      project     = as.character(row$project),
      category_1  = as.character(row$category_1),
      ftir_stems  = ftir_stems
    )

  }

}


## ---------------------------------------------------------------------------
## Extract results into vectors
## ---------------------------------------------------------------------------

ftir_stems_out <- character(length(results))
match_methods  <- character(length(results))

for (i in seq_along(results)) {

  r <- results[[i]]

  ftir_stems_out[i] <- if (is.null(r$stem) || is.na(r$stem)) NA_character_ else r$stem
  match_methods[i]  <- if (is.null(r$method) || length(r$method) == 0) "error" else r$method

}


## ---------------------------------------------------------------------------
## Build lookup table
## ---------------------------------------------------------------------------

samples %>%
  dplyr::select(sample_id, project, sample_name_orig) %>%
  dplyr::mutate(
    ftir_stem    = ftir_stems_out,
    match_method = match_methods
  ) -> lookup


## Add canonical pyom_id ------------------------------------------------------

lookup %>%
  dplyr::mutate(
    pyom_id  = sprintf("pyom-%04d", dplyr::row_number()),
    verified = FALSE
  ) %>%
  dplyr::select(
    pyom_id,
    project,
    sample_id,
    sample_name_orig,
    ftir_stem,
    match_method,
    verified
  ) -> lookup


## =============================================================================
## Section 5: Summary and Output
## =============================================================================


## ---------------------------------------------------------------------------
## Print matching summary
## ---------------------------------------------------------------------------

cat("\n=== Matching Summary ===\n")
lookup %>%
  dplyr::count(match_method) %>%
  print()

cat("\n=== By Project ===\n")
lookup %>%
  dplyr::group_by(project) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    matched   = sum(!is.na(ftir_stem)),
    unmatched = sum(is.na(ftir_stem)),
    .groups   = "drop"
  ) %>%
  print()


## ---------------------------------------------------------------------------
## Report orphan FTIR files
## ---------------------------------------------------------------------------

lookup %>%
  dplyr::filter(!is.na(ftir_stem)) %>%
  dplyr::pull(ftir_stem) %>%
  unique() -> matched_stems

orphan_ftir <- setdiff(ftir_stems, matched_stems)

cat(sprintf("\nOrphan FTIR files (no sample match): %d\n", length(orphan_ftir)))

if (length(orphan_ftir) > 0 && length(orphan_ftir) <= 20) {

  cat(" ", paste(orphan_ftir[1:min(10, length(orphan_ftir))], collapse = "\n  "), "\n")

}


## ---------------------------------------------------------------------------
## Write output
## ---------------------------------------------------------------------------

readr::write_csv(lookup, "data/processed/sample_ftir_lookup.csv")
cat("\nSaved to: data/processed/sample_ftir_lookup.csv\n")
