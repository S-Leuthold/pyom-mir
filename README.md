# PyOM-MIR

Predicting stable polycyclic aromatic carbon (SPAC) from mid-infrared spectra.

## Overview

Pyrogenic carbon (PyC) is a long-lived pool of soil carbon generated from biomass combustion during wildfire, prescribed fire, or biochar application. Hydrogen pyrolysis (HyPy) provides precise PyC measurements but is slow and expensive. This project evaluates a range of model configurations and performs predictive analysis of SPAC from mid-infrared (MIR) spectra, enabling rapid, non-destructive estimation of PyC at scale.

Spectral processing uses the [horizons](https://github.com/S-Leuthold/horizons) package.

## Data

Data are not included in this repository. See [`data/README.md`](data/README.md) for structure and sources.

## Requirements

- R >= 4.x
- Dependencies managed with [renv](https://rstudio.github.io/renv/)

```r
install.packages("renv")
renv::restore()
```

## Reproduction

Run scripts in `analysis/` in numerical/alphabetical order. Output is written to `data/processed/`.

## Funding

Supported by the NSF Regional Innovation Engines program through the [Colorado-Wyoming Climate Resilience Engine](https://www.co-wyengine.org/).

## License

[MIT](LICENSE)
