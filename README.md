# PyOM-MIR

Predicting pyrogenic carbon (PyC) from mid-infrared spectra.

## Overview

Pyrogenic carbon is a long-lived pool of soil carbon generated from biomass combustion during wildfire, prescribed fire, or biochar application. Hydrogen pyrolysis (HyPy) provides precise PyC measurements but is slow and expensive. This project builds predictive models for PyC from mid-infrared (MIR) spectra, enabling rapid, non-destructive estimation at scale.

The dataset spans 436 samples across 7 projects, reduced to 307 after QC filtering. The best global model (SVM-RBF + SNV + log) achieves RPD 1.70. A hybrid approach using spectrally-defined local models pushes cluster-level performance above RPD 2.0.

Spectral processing uses the [horizons](https://github.com/S-Leuthold/horizons) package.

## Pipeline

```
prepare-data.R ──> local-cluster-models.R ──> build-figures.R
       │                                            │
       ▼                                            ▼
  horizons_data.rds                         output/figures/
  modeling_response.csv
```

| Script | Purpose | Manual |
|--------|---------|--------|
| [`prepare-data.R`](analysis/prepare-data.R) | Clean response data, build horizons object | [manual](analysis/manuals/prepare-data.md) |
| [`local-cluster-models.R`](analysis/local-cluster-models.R) | PCA clustering + per-cluster model evaluation | [manual](analysis/manuals/local-cluster-models.md) |
| [`build-figures.R`](analysis/build-figures.R) | Generate all report figures | [manual](analysis/manuals/build-figures.md) |

The `analysis/hpc/` directory contains scripts for running the global model evaluation on a Linux cluster.

## Directory structure

```
analysis/
├── prepare-data.R              Core pipeline scripts
├── local-cluster-models.R
├── build-figures.R
├── hpc/                        HPC batch scripts for Sybil
├── qc/                         QC and diagnostic scripts
│   ├── assess-spectra.R
│   ├── deduplicate-ftir.R
│   └── diagnose-results.R
├── manuals/                    Script documentation
└── archive/                    Superseded scripts
data/                           Not tracked (see data/README.md)
output/                         Not tracked
```

## Data

Data are not included in this repository. See [`data/README.md`](data/README.md) for structure and sources.

## Requirements

- R >= 4.x
- Dependencies managed with [renv](https://rstudio.github.io/renv/)

```r
install.packages("renv")
renv::restore()
```

## Funding

Supported by the NSF Regional Innovation Engines program through the [Colorado-Wyoming Climate Resilience Engine](https://www.co-wyengine.org/).

## License

[MIT](LICENSE)
