#' Prepare Horizons Data Object
#'
#' @description
#' Runs the horizons v1 pipeline on PyOM FTIR spectra: reads OPUS files,
#' standardizes, extracts sample IDs, averages replicates, and joins with
#' response data. Produces a horizons_data RDS ready for configure() + evaluate().
#'
#' @details
#' Input:
#'   - data/raw/Final FTIR/ (1,696 OPUS files)
#'   - data/processed/modeling_response.csv (from prepare-response.R)
#'
#' Output:
#'   - data/processed/horizons_data.rds
#'
#' Two FTIR naming conventions exist in this dataset:
#'
#'   Standard (most samples):  {stem}-S{1-4}_{position}.0
#'     e.g., 1_4_5_S1U1_0_5_Soil-S1_A11.0
#'     stem = 1_4_5_S1U1_0_5_Soil
#'
#'   FORCE (Akron/Hoytville/KBS): {stem}_{position}.0
#'     e.g., AK_307_1_A8.0
#'     stem = AK_307_1
#'
#' parse_ids() handles the standard pattern. FORCE files fall through to
#' "keep_original" and get cleaned up with a manual position-code strip.
#'
#' @keywords internal
#' @name prepare-horizons-data


## =============================================================================
## Section 1: Setup
## =============================================================================


library(horizons)


## =============================================================================
## Section 2: Read and Standardize Spectra
## =============================================================================


horizons::spectra("data/raw/Final FTIR", type = "opus") |>
  horizons::standardize(
    resample = 2,
    trim     = c(600, 4000)
  ) ->
hd

cat(sprintf("\nLoaded %d spectra, %d wavenumbers\n",
            hd$data$n_rows, hd$data$n_predictors))


## =============================================================================
## Section 3: Parse Sample IDs
## =============================================================================
##
## Standard filenames: {stem}-S{replicate}_{position}
##   parse_ids extracts the stem as sample_id
##
## FORCE filenames: {stem}_{position}  (no -S[1-4])
##   parse_ids keeps original → manual cleanup below


horizons::parse_ids(
  hd,
  patterns = c(
    sampleid  = "[^-]+",
    "-S",
    replicate = "[1-4]",
    "_.*"
  ),
  too_few = "keep_original"
) ->
hd


## ---------------------------------------------------------------------------
## Clean up FORCE sample IDs
## ---------------------------------------------------------------------------
##
## Files that didn't match the -S[1-4] pattern still have trailing position
## codes: "AK_307_1_A8" needs to become "AK_307_1"
##
## Strategy: if sample_id still contains an _[A-H]\d+$ suffix (position code),
## strip it. This is safe because real sample stems never end in this pattern.

hd$data$analysis$sample_id <- hd$data$analysis$sample_id |>
  stringr::str_remove("_[A-H][0-9]+$")

## Handle " 2" suffix on duplicate OPUS scans (e.g., "AK_205_1_A6.0 2")
## and trim any trailing whitespace from OPUS filenames
hd$data$analysis$sample_id <- hd$data$analysis$sample_id |>
  stringr::str_remove(" [0-9]+$") |>
  stringr::str_trim()

cat(sprintf("Unique sample IDs after parsing: %d\n",
            length(unique(hd$data$analysis$sample_id))))


## =============================================================================
## Section 4: Average Replicates
## =============================================================================


horizons::average(
  hd,
  by                    = "sample_id",
  quality_check         = TRUE,
  correlation_threshold = 0.996
) ->
hd

cat(sprintf("After averaging: %d unique samples\n", hd$data$n_rows))


## =============================================================================
## Section 5: Join Response Data
## =============================================================================


horizons::add_response(
  hd,
  source   = "data/processed/modeling_response.csv",
  variable = "spac_g_kg_c",
  by       = c("sample_id" = "Sample_ID")
) ->
hd


## ---------------------------------------------------------------------------
## Report match rate
## ---------------------------------------------------------------------------

outcome_col  <- hd$data$role_map$variable[hd$data$role_map$role == "response"]
n_matched    <- sum(!is.na(hd$data$analysis[[outcome_col]]))
n_total      <- hd$data$n_rows

cat(sprintf("\nResponse join: %d/%d samples matched (%.1f%%)\n",
            n_matched, n_total, 100 * n_matched / n_total))
cat(sprintf("Unmatched samples: %d (standards, chars, or orphan spectra)\n",
            n_total - n_matched))


## =============================================================================
## Section 6: Save
## =============================================================================


saveRDS(hd, "data/processed/horizons_data.rds")
cat(sprintf("\nSaved horizons_data to: data/processed/horizons_data.rds\n"))
cat(sprintf("Object size: %.1f MB\n", object.size(hd) / 1e6))
