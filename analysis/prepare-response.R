#' Prepare Modeling Response Data
#'
#' @description
#' Reads Michelle's curated HyPy spreadsheet, cleans it, and produces the
#' response CSV that horizons' add_response() can join against FTIR spectra.
#'
#' @details
#' Input:  data/raw/Engine_hypydata_Combined_2025_12_09.xlsx
#' Output: data/processed/modeling_response.csv
#'
#' This script replaces the previous 3-script pipeline:
#'   - clean-hypy-data.R
#'   - build-sample-lookup.R
#'   - build-modeling-response.R
#'
#' Michelle's Feb 2026 data upload curated the dataset to 436 samples, all
#' with explicit FTIR IDs. The complex auto-matching logic is no longer needed.
#'
#' @keywords internal
#' @name prepare-response


## =============================================================================
## Section 1: Setup
## =============================================================================


library(tidyverse)
library(readxl)


## =============================================================================
## Section 2: Load and Standardize
## =============================================================================


## ---------------------------------------------------------------------------
## Read xlsx
## ---------------------------------------------------------------------------
##
## The xlsx has a weird header structure:
##   Row 1: NA | "All in final, 2nd check" | NA | ...  (annotation)
##   Row 2: blank
##   Row 3: Project | sample name | Category 1 | ...   (actual header)
##   Row 4+: data
##
## skip = 2 uses row 3 as column names.

readxl::read_xlsx(
  "data/raw/Engine_hypydata_Combined_2025_12_09.xlsx",
  skip = 2
) %>%
  janitor::clean_names() %>%
  dplyr::rename(
    sample_name  = sample_name,
    pct_c        = percent_c,
    spac_g_kg_c  = g_spa_c_kg_c,
    spac_toc_pct = spac_toc_percent,
    ftir_1       = ftir_spectra_1
  ) -> df


## ---------------------------------------------------------------------------
## Drop header row that got read as data
## ---------------------------------------------------------------------------

df <- df[df$project != "Project" & !is.na(df$project), ]

cat(sprintf("Loaded %d samples from xlsx\n", nrow(df)))


## =============================================================================
## Section 3: Parse Numeric Columns
## =============================================================================
##
## The xlsx stores everything as text. Coerce numeric columns.


df %>%
  dplyr::mutate(
    pct_c        = suppressWarnings(as.numeric(pct_c)),
    spac_g_kg_c  = suppressWarnings(as.numeric(spac_g_kg_c)),
    spac_toc_pct = suppressWarnings(as.numeric(spac_toc_pct))
  ) -> df


## =============================================================================
## Section 4: Back-Calculate SPAC
## =============================================================================
##
## Australia and France don't have direct g SPA-C/kg C values.
## They have SPAC/TOC % which converts: spac_g_kg_c = spac_toc_pct * 10


df %>%
  dplyr::mutate(
    spac_g_kg_c = dplyr::case_when(
      !is.na(spac_g_kg_c)  ~ spac_g_kg_c,
      !is.na(spac_toc_pct) ~ round(spac_toc_pct * 10, 2),
      TRUE ~ NA_real_
    )
  ) -> df


## =============================================================================
## Section 5: Extract FTIR Stem
## =============================================================================
##
## The ftir_1 column contains the FTIR filename with replicate suffix.
## Strip "-S[1-4]..." to get the sample stem that matches parse_ids() output.


df %>%
  dplyr::mutate(
    ftir_stem = ftir_1 %>%
      stringr::str_remove("-S[1-4].*$") %>%
      stringr::str_trim() %>%
      dplyr::if_else(. == "" | is.na(.), NA_character_, .)
  ) -> df


## =============================================================================
## Section 6: QC and Filter
## =============================================================================


## Drop standards and non-soil reference materials --------------------------
##
## NDN_GroundBulk = instrument check standard (no HyPy, but catch it here)
## BPCAStd        = 20% biochar standard (SPAC ~850 g/kg C)
## CIG_Biochar    = reference biochar

df %>%
  dplyr::filter(
    !stringr::str_detect(ftir_stem, "(?i)NDN|BPCAStd|CIG_Biochar")
  ) -> df

cat(sprintf("  After dropping standards: %d\n", nrow(df)))


## QC flags -----------------------------------------------------------------

df %>%
  dplyr::mutate(
    qc_flag = dplyr::case_when(
      spac_toc_pct > 100 ~ "spac_exceeds_toc",
      TRUE ~ NA_character_
    )
  ) -> df


## Print QC summary ---------------------------------------------------------

cat("\n=== QC Summary ===\n")
cat(sprintf("  Total samples: %d\n", nrow(df)))
cat(sprintf("  Has SPAC g/kg C: %d\n", sum(!is.na(df$spac_g_kg_c))))
cat(sprintf("  Has FTIR stem: %d\n", sum(!is.na(df$ftir_stem))))
cat(sprintf("  QC flagged: %d\n", sum(!is.na(df$qc_flag))))


## Filter to usable samples ------------------------------------------------

df %>%
  dplyr::filter(
    !is.na(spac_g_kg_c),
    !is.na(ftir_stem),
    is.na(qc_flag)
  ) -> df_clean

cat(sprintf("  After filtering: %d\n", nrow(df_clean)))


## =============================================================================
## Section 7: Summary by Project
## =============================================================================


cat("\n=== By Project ===\n")
df_clean %>%
  dplyr::count(project) %>%
  print()


## =============================================================================
## Section 8: Output
## =============================================================================
##
## Sample_ID must match the stem that parse_ids() extracts from OPUS filenames.
## Both strip "-S[1-4]" and trailing position codes.


df_clean %>%
  dplyr::transmute(
    Sample_ID   = ftir_stem,
    spac_g_kg_c = spac_g_kg_c,
    pct_c       = pct_c,
    project     = project
  ) -> modeling_response


## Check for duplicates -----------------------------------------------------

dupes <- modeling_response$Sample_ID[duplicated(modeling_response$Sample_ID)]

if (length(dupes) > 0) {

  cat(sprintf("\n⚠ %d duplicate Sample_IDs found:\n", length(dupes)))
  cat(paste("  ", unique(dupes), collapse = "\n"), "\n")
  cat("Keeping first occurrence.\n")

  modeling_response <- modeling_response[!duplicated(modeling_response$Sample_ID), ]

}


## Write output -------------------------------------------------------------

readr::write_csv(modeling_response, "data/processed/modeling_response.csv")
cat(sprintf("\nSaved %d samples to: data/processed/modeling_response.csv\n", nrow(modeling_response)))
