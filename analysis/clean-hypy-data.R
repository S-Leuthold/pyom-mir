#' Clean HyPy Data
#'
#' @description
#' Standardizes and cleans the raw HyPy dataset for PyOM-MIR modeling.
#'
#' @details
#' Input:  data/raw/Engine_hypydata_Combined_2025_12_09(Sheet1).csv
#' Output: data/processed/hypy_clean.csv
#'
#' Processing steps:
#'
#' 1. Standardize column names
#' 2. Parse depth from various formats (project-specific)
#' 3. Extract site/treatment/replicate where possible
#' 4. Back-calculate g SPA-C/kg C from SPAC/TOC % where missing
#' 5. Flag QC issues (outliers >100%, BDL, duplicates)
#' 6. Standardize FTIR IDs for matching
#'
#' @keywords internal
#' @name clean-hypy-data


## =============================================================================
## Section 1: Setup
## =============================================================================


library(tidyverse)
library(janitor)


## =============================================================================
## Section 2: Helper Functions
## =============================================================================


## ---------------------------------------------------------------------------
## parse_depth — Extract top/bottom depth from string
## ---------------------------------------------------------------------------

#' Parse depth string into numeric top/bottom values
#'
#' @description
#' Handles various depth formats: "0-5", "15-30", "0-5 cm", etc.
#'
#' @param x Character. Depth string to parse.
#'
#' @return Numeric vector of length 2: c(top, bottom). Returns c(NA, NA) if
#'   parsing fails.
#'
#' @noRd
parse_depth <- function(x) {


  if (is.na(x) || x == "" || tolower(x) == "na") {

    return(c(NA_real_, NA_real_))

  }

  ## Remove "cm" suffix -------------------------------------------------------


  x <- stringr::str_remove(x, "\\s*cm$")

  ## Handle range format: "0-5", "15-30", etc. --------------------------------

  if (stringr::str_detect(x, "^\\d+[_-]\\d+$")) {

    parts <- stringr::str_split(x, "[_-]")[[1]]
    return(as.numeric(parts))

  }

  ## Handle single value — can't determine range ------------------------------

  if (stringr::str_detect(x, "^\\d+$")) {

    return(c(NA_real_, NA_real_))

  }

  c(NA_real_, NA_real_)

}


## ---------------------------------------------------------------------------
## standardize_ftir_id — Clean FTIR column values
## ---------------------------------------------------------------------------

#' Standardize FTIR ID strings
#'
#' @description
#' Removes replicate suffixes (-S1, -S2, etc.) from FTIR column values.
#'
#' @param x Character vector. FTIR ID strings.
#'
#' @return Character vector with suffixes removed, empty strings converted to NA.
#'
#' @noRd
standardize_ftir_id <- function(x) {

  x %>%
    stringr::str_remove("-S[1-4].*$") %>%
    stringr::str_trim() %>%
    dplyr::if_else(. == "" | is.na(.), NA_character_, .) -> result

  result

}


## =============================================================================
## Section 3: Load and Rename
## =============================================================================


## ---------------------------------------------------------------------------
## Load raw data
## ---------------------------------------------------------------------------

raw <- readr::read_csv(
  "data/raw/Engine_hypydata_Combined_2025_12_09(Sheet1).csv",
  show_col_types = FALSE
)


## ---------------------------------------------------------------------------
## Standardize column names
## ---------------------------------------------------------------------------

raw %>%
  janitor::clean_names() %>%
  dplyr::rename(
    sample_name_orig = sample_name,
    pct_c            = percent_c,
    spac_g_kg_dry    = g_spa_c_kg_dry_mass,
    bpca_g_kg_dry    = g_bpca_c_kg_dry_mass,
    spac_g_kg_c_orig = g_spa_c_kg_c,
    spac_toc_pct_orig = spac_toc_percent,
    pyc_toc          = py_c_toc,
    ftir_1           = ftir_spectra_1,
    ftir_2           = ftir_spectra_2,
    ftir_3           = ftir_spectra_3,
    ftir_4           = ftir_spectra_4
  ) -> df


## =============================================================================
## Section 4: Parse Depth
## =============================================================================
##
## Depth appears in different columns per project:
##   BCInt:     category_3 (e.g., "0-5", "15-30") or category_6 (e.g., "0-5 cm")
##   France:    category_3 (e.g., "0-10", "0-15", "0-20")
##   Jens:      category_5 (e.g., "0-5", "15-30", "30-50")
##   Australia: category_5 (e.g., "5", "25" — single values, unusable)
##   FORCE:     not clearly specified
##


df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    depth_source = dplyr::case_when(
      project == "BCInt" & !is.na(category_6) & stringr::str_detect(category_6, "\\d") ~ category_6,
      project == "BCInt" & !is.na(category_3) & stringr::str_detect(category_3, "^\\d") ~ category_3,
      project == "France" & !is.na(category_3) ~ category_3,
      project == "Jens" & !is.na(category_5) ~ category_5,
      TRUE ~ NA_character_
    ),
    depth_parsed = list(parse_depth(depth_source)),
    depth_top_cm = depth_parsed[[1]],
    depth_bot_cm = depth_parsed[[2]]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-depth_parsed, -depth_source) -> df


## =============================================================================
## Section 5: Extract Site, Treatment, Material Type
## =============================================================================


df %>%
  dplyr::mutate(

    ## Material type (BCInt only has this clearly) ----------------------------

    material_type = dplyr::case_when(
      project == "BCInt" & tolower(category_1) %in% c("soil", "char", "litter", "sediment bank") ~ tolower(category_1),
      TRUE ~ "soil"
    ),

    ## Site extraction (project-specific) -------------------------------------

    site = dplyr::case_when(
      project == "BCInt"     ~ category_2,
      project == "FORCE"     ~ category_1,
      project == "France"    ~ category_2,
      project == "Jens"      ~ category_2,
      project == "Australia" ~ category_3,
      TRUE ~ NA_character_
    ),

    ## Treatment (where applicable) -------------------------------------------

    treatment = dplyr::case_when(
      project == "BCInt" ~ category_4,
      project == "FORCE" ~ category_2,
      project == "Jens"  ~ category_3,
      TRUE ~ NA_character_
    ),

    ## Replicate --------------------------------------------------------------

    replicate = dplyr::case_when(
      project == "BCInt"  ~ category_5,
      project == "FORCE"  ~ category_4,
      project == "France" ~ category_4,
      TRUE ~ NA_character_
    )

  ) -> df


## =============================================================================
## Section 6: Back-Calculate SPAC
## =============================================================================
##
## g SPA-C/kg C = SPAC/TOC % * 10
##


df %>%
  dplyr::mutate(

    ## Parse SPAC/TOC % (handle "<DL" values) ---------------------------------

    spac_toc_pct_numeric = dplyr::case_when(
      stringr::str_detect(spac_toc_pct_orig, "^<") ~ NA_real_,
      tolower(spac_toc_pct_orig) == "na"           ~ NA_real_,
      TRUE ~ suppressWarnings(as.numeric(spac_toc_pct_orig))
    ),

    ## Back-calculate where original is missing -------------------------------

    spac_g_kg_c = dplyr::case_when(
      !is.na(spac_g_kg_c_orig)      ~ spac_g_kg_c_orig,
      !is.na(spac_toc_pct_numeric)  ~ round(spac_toc_pct_numeric * 10, 2),
      TRUE ~ NA_real_
    )

  ) -> df


## =============================================================================
## Section 7: QC Flags
## =============================================================================


## Identify duplicate sample names --------------------------------------------

df %>%
  dplyr::count(sample_name_orig) %>%
  dplyr::filter(n > 1) %>%
  dplyr::pull(sample_name_orig) -> dupe_names


## Apply QC flags -------------------------------------------------------------

df %>%
  dplyr::mutate(
    qc_flag = dplyr::case_when(
      stringr::str_detect(spac_toc_pct_orig, "^<") ~ "below_detection_limit",
      spac_toc_pct_numeric > 100                   ~ "spac_exceeds_toc",
      sample_name_orig %in% dupe_names             ~ "duplicate_name",
      TRUE ~ NA_character_
    )
  ) -> df


## Print QC summary -----------------------------------------------------------

cat("\n=== QC Flag Summary ===\n")
df %>%
  dplyr::count(qc_flag) %>%
  print()


## =============================================================================
## Section 8: Create Unique Sample ID
## =============================================================================


df %>%
  dplyr::group_by(sample_name_orig) %>%
  dplyr::mutate(
    dupe_n  = dplyr::row_number(),
    n_dupes = dplyr::n(),
    sample_id = dplyr::if_else(
      n_dupes > 1,
      paste0(sample_name_orig, "_", dupe_n),
      sample_name_orig
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-dupe_n, -n_dupes) -> df


## =============================================================================
## Section 9: Standardize FTIR IDs
## =============================================================================


df %>%
  dplyr::mutate(
    ftir_id = dplyr::coalesce(
      standardize_ftir_id(ftir_1),
      standardize_ftir_id(ftir_2),
      standardize_ftir_id(ftir_3),
      standardize_ftir_id(ftir_4)
    )
  ) -> df


## =============================================================================
## Section 10: Select Final Columns
## =============================================================================


df %>%
  dplyr::select(
    sample_id,
    project,
    sample_name_orig,
    site,
    depth_top_cm,
    depth_bot_cm,
    treatment,
    replicate,
    material_type,
    pct_c,
    spac_g_kg_c,
    spac_toc_pct_orig,
    spac_g_kg_c_orig,
    qc_flag,
    ftir_id,
    dplyr::starts_with("category_"),
    spac_g_kg_dry,
    bpca_g_kg_dry,
    pyc_toc,
    ftir_1, ftir_2, ftir_3, ftir_4
  ) -> df_clean


## =============================================================================
## Section 11: Summary and Output
## =============================================================================


## Print summary --------------------------------------------------------------

cat("\n=== Dataset Summary ===\n")
cat(sprintf("Total samples: %d\n", nrow(df_clean)))
cat(sprintf("Samples with SPAC g/kg C: %d\n", sum(!is.na(df_clean$spac_g_kg_c))))
cat(sprintf("Samples with FTIR ID: %d\n", sum(!is.na(df_clean$ftir_id))))
cat(sprintf("Samples with QC flags: %d\n", sum(!is.na(df_clean$qc_flag))))
cat(sprintf(
  "Clean samples (no QC flag + has SPAC + has FTIR): %d\n",
  sum(is.na(df_clean$qc_flag) & !is.na(df_clean$spac_g_kg_c) & !is.na(df_clean$ftir_id))
))

cat("\n=== By Project ===\n")
df_clean %>%
  dplyr::group_by(project) %>%
  dplyr::summarise(
    n        = dplyr::n(),
    has_spac = sum(!is.na(spac_g_kg_c)),
    has_ftir = sum(!is.na(ftir_id)),
    flagged  = sum(!is.na(qc_flag)),
    .groups  = "drop"
  ) %>%
  print()


## Write output ---------------------------------------------------------------

readr::write_csv(df_clean, "data/processed/hypy_clean.csv")
cat("\nSaved to: data/processed/hypy_clean.csv\n")
