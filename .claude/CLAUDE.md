# PyOM-MIR

Predicting pyrogenic carbon (PyC) from mid-infrared spectra using the horizons package.

## Project status

Back-burner after Feb 9, 2026 cleanup session. Preliminary results delivered to Francesca. Next meeting Feb 16.

## Pipeline

```
prepare-data.R  →  (hpc/ or local-cluster-models.R)  →  diagnose-results.R  →  build-figures.R
```

## R script style

These are interactive analysis scripts, not package code. Written for someone running line-by-line in RStudio.

**Assignment:**
- Right assignment (`-> name`) for multi-line piped expressions
- Left assignment (`name <-`) for single-line assignments

**Formatting:**
- Function name and first argument on the same line
- Subsequent arguments aligned under the first argument
- Closing paren on the last argument line (not its own line)
- Pipe operator on the same line as the closing paren

```r
# yes
read_xlsx("data/raw/file.xlsx",
          skip = 2) %>%
  clean_names() %>%
  rename(pct_c = percent_c,
         spac  = g_spa_c) -> df

# no
readxl::read_xlsx(
  "data/raw/file.xlsx",
  skip = 2
) %>%
  janitor::clean_names() %>%
  dplyr::rename(
    pct_c = percent_c,
    spac  = g_spa_c
  ) -> df
```

**Packages:**
- Load all packages at the top with `pacman::p_load()`
- No namespace prefixes in the script body — just bare function names

**Comments:**
- Keep section headers with full-width separators
- Comments explain "why" not "what"
- No roxygen headers (not a package)

**Output:**
- No cat/sprintf messaging — the person running it can inspect objects
- No progress bars or status updates
