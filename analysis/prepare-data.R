## =============================================================================
## prepare-data.R
## Prepare PyOM-MIR modeling data
##
## Author:        Sam Leuthold
## Contact:       sam.leuthold@colostate.edu
## Last modified: 2026-02-09
##
## Description:
##   Cleans the curated HyPy response data (436 samples), applies QC filters
##   and sample drops identified during diagnostic work, derives absolute PyC
##   concentrations, then builds a horizons data object from FTIR OPUS files
##   and joins the response.
##
## Input:
##   data/raw/Engine_hypydata_Combined_2025_12_09.xlsx
##   data/raw/Final FTIR/  (1,696 OPUS files)
##
## Output:
##   data/processed/modeling_response.csv
##   data/processed/horizons_data.rds
##
## Changelog:
##   2026-02-09  Consolidated from prepare-response.R + prepare-horizons-data.R
##   2026-02-06  Added spectral outlier drops (Mahalanobis), pyc_abs response
##   2026-02-05  Litter drop (RPD 0.48 -> 1.26), BCInt drop (residuals +0.17
##               to +0.29), disabled ALS (RPD 1.46 vs 1.51 without)
##   2026-02-04  Initial version (prepare-response.R + prepare-horizons-data.R)
## =============================================================================


pacman::p_load(tidyverse, readxl, janitor, horizons, stringr)


## =============================================================================
## Step 1: Response data
## =============================================================================


## Read xlsx -------------------------------------------------------------------

# The xlsx has a weird header structure: row 1 is an annotation row,
# row 2 is blank, row 3 is the actual header. skip = 2 uses row 3.

read_xlsx("data/raw/Engine_hypydata_Combined_2025_12_09.xlsx",
          skip = 2) %>%
  clean_names() %>%
  rename(sample_name  = sample_name,
         pct_c        = percent_c,
         spac_g_kg_c  = g_spa_c_kg_c,
         spac_toc_pct = spac_toc_percent,
         ftir_1       = ftir_spectra_1) -> df

## Drop header row that got read as data ---------------------------------------

df <- df[df$project != "Project" & !is.na(df$project), ]


## Parse numeric columns ------------------------------------------------------

# The xlsx stores everything as text

df %>%
  mutate(pct_c        = suppressWarnings(as.numeric(pct_c)),
         spac_g_kg_c  = suppressWarnings(as.numeric(spac_g_kg_c)),
         spac_toc_pct = suppressWarnings(as.numeric(spac_toc_pct))) -> df


## Back-calculate SPAC ---------------------------------------------------------

# Australia and France don't have direct g SPA-C/kg C values.
# They have SPAC/TOC % which converts: spac_g_kg_c = spac_toc_pct * 10

df %>%
  mutate(spac_g_kg_c = case_when(!is.na(spac_g_kg_c)  ~ spac_g_kg_c,
                                 !is.na(spac_toc_pct) ~ round(spac_toc_pct * 10, 2),
                                 TRUE ~ NA_real_)) -> df


## Extract FTIR stem -----------------------------------------------------------

# The ftir_1 column has the FTIR filename with replicate suffix.
# Strip "-S[1-4]..." to get the sample stem that matches parse_ids() output.

df %>%
  mutate(ftir_stem = ftir_1 %>%
           str_remove("-S[1-4].*$") %>%
           str_remove("_Whole Soil$") %>%
           str_trim() %>%
           if_else(. == "" | is.na(.), NA_character_, .)) -> df


## =============================================================================
## Step 2: QC and filter
## =============================================================================


## Drop standards and non-soil reference materials -----------------------------

# NDN_GroundBulk = instrument check standard
# BPCAStd        = 20% biochar standard (SPAC ~850 g/kg C)
# CIG_Biochar    = reference biochar

df %>%
  filter(!str_detect(ftir_stem, "(?i)NDN|BPCAStd|CIG_Biochar")) -> df

## Flag SPAC > TOC (physically impossible) -------------------------------------

df %>%
  mutate(qc_flag = case_when(
           spac_toc_pct > 100 ~ "spac_exceeds_toc",
           TRUE ~ NA_character_)) -> df

## Keep only samples with SPAC, FTIR, and no QC flags --------------------------

df %>%
  filter(!is.na(spac_g_kg_c),
         !is.na(ftir_stem),
         is.na(qc_flag)) -> df_clean

## =============================================================================
## Step 3: Sample drops
## =============================================================================

## Litter — different matrix, different spectral properties --------------------

# Different FTIR signatures from soil. 10x higher PyC concentrations.

df_clean <- df_clean |>
  filter(!str_detect(ftir_stem, "(?i)Litter"))

## BCInt — biochar amendment experimental plots --------------------------------

# Artificially elevated PyC from biochar amendments. Different land use
# context from the natural/fire-affected soils in other projects.

df_clean <- df_clean |>
  filter(project != "BCInt")

## Named drops — data errors or uncertain response -----------------------------

drop_ids <- c("AUS_54_Soil",                    # SPAC = 0.0 g/kg C (probable data error)
              "7_HTC_5_15_Soil",                # SPAC 309 at 1.13% C (suspicious ratio)
              "Swiss_F_12_Soil",                # Above 3xIQR fence (SPAC = 314)
              "Poudre_Watershed_Hist_BC_Soil",  # 3 conflicting HyPy measurements
              "2_LTC_5_15_Soil"                 # Mislabeled depth
              )

## Spectral outliers — Mahalanobis distance on PC1:10 > chi-sq(0.975, df=10)

# Unusual spectral signatures (deep subsoils, corrupted scans, atypical
# mineralogy) that create high leverage in models.

spectral_outlier_ids <- c("Swiss_F_51_Soil",                # Mahalanobis = 282 (extreme outlier)
                          "AUS_75_Soil",                    # Mahal = 68
                          "AUS_59_Soil",                    # Mahal = 58
                          "AUS_39_Soil",                    # Mahal = 49
                          "Swiss_F_20_Soil",                # Mahal = 46
                          "6_LTC_30_50_Soil",               # Mahal = 43 (deep subsoil)
                          "7_HTC_50_75_Soil",               # Mahal = 40 (deep subsoil)
                          "6_LTC_50_75_Soil",               # Mahal = 37 (deep subsoil)
                          "AUS_5_Soil",                     # Mahal = 33
                          "AUS_70_Soil",                    # Mahal = 27
                          "AUS_8_Soil",                     # Mahal = 26
                          "Swiss_F_23_Soil",                # Mahal = 25
                          "AUS_47_Soil",                    # Mahal = 24
                          "1_HTC_30_50_Soil",               # Mahal = 23 (deep subsoil)
                          "AUS_7_Soil",                     # Mahal = 22
                          "AUS_85_Soil",                    # Mahal = 22
                          "7_HTC_30_50_Soil",               # Mahal = 22 (deep subsoil)
                          "AUS_19_Soil"                     # Mahal = 22
                          )

df_clean <- df_clean |>
  filter(!(ftir_stem %in% c(drop_ids, spectral_outlier_ids)))

## =============================================================================
## Step 4: Derive absolute PyC response
## =============================================================================

# pyc_abs = g PyC per 100g soil = pct_c * spac_g_kg_c / 1000
# Easier to predict from MIR than the ratio (g/kg C) because spectra see
# total C directly (PC1 vs total C: r = 0.54).

df_clean <- df_clean |>
  mutate(pyc_abs = pct_c * spac_g_kg_c / 1000)

## Trim response outliers beyond mean + 2 SD -----------------------------------

# Remaining high-PyC samples (e.g., Swiss_F_6_Soil at 2.43) dominate the
# loss function. Trimming to mean + 2 SD keeps the distribution learnable.

pyc_cutoff <- mean(df_clean$pyc_abs, na.rm = TRUE) +
  2 * sd(df_clean$pyc_abs, na.rm = TRUE)

df_clean <- df_clean |>
  filter(pyc_abs <= pyc_cutoff)

## =============================================================================
## Step 5: Write response CSV
## =============================================================================

# Sample_ID must match the stem that parse_ids() extracts from OPUS filenames.
# Both strip "-S[1-4]" and trailing position codes.

df_clean %>%
  transmute(Sample_ID   = ftir_stem,
            spac_g_kg_c = spac_g_kg_c,
            pyc_abs     = pyc_abs,
            pct_c       = pct_c,
            project     = project) -> modeling_response

## Deduplicate — keep first occurrence if any stems collide --------------------

modeling_response <- modeling_response[!duplicated(modeling_response$Sample_ID), ]

write_csv(modeling_response, "data/processed/modeling_response.csv")

## =============================================================================
## Step 6: Build horizons data object
## =============================================================================

## Read and standardize spectra ------------------------------------------------

spectra("data/raw/Final FTIR",
        type = "opus") |>
  standardize(resample = 2,
              trim     = c(600, 4000)) -> hd

# ALS baseline correction tested and disabled — SNV and SG through horizons
# handle baseline variation more effectively.

## Parse sample IDs ------------------------------------------------------------

# Two FTIR naming conventions:
#   Standard: {stem}-S{1-4}_{position}.0  -> parse_ids extracts the stem
#   FORCE:    {stem}_{position}.0         -> falls through, cleaned up below

parse_ids(hd,
          patterns = c(sampleid  = "[^-]+",
                       "-S",
                       replicate = "[1-4]",
                       "_.*"),
          too_few  = "keep_original") -> hd

# Clean up FORCE sample IDs — strip trailing position codes (_A8, _B3, etc.)
# and duplicate-scan suffixes (" 2")

hd$data$analysis$sample_id |>
  str_remove("_[A-H][0-9]+$") |>
  str_remove(" [0-9]+$") |>
  str_trim() -> hd$data$analysis$sample_id


## Average replicates ----------------------------------------------------------

average(hd,
        by                    = "sample_id",
        quality_check         = TRUE,
        correlation_threshold = 0.996) -> hd


## Join response ---------------------------------------------------------------

add_response(hd,
             source   = "data/processed/modeling_response.csv",
             variable = "pyc_abs",
             by       = c("sample_id" = "Sample_ID")) -> hd


## Save ------------------------------------------------------------------------

saveRDS(hd, "data/processed/horizons_data.rds")
