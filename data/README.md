# Data

Data are not included in this repository but will be made available upon publication.

## Structure

```
data/
├── raw/
│   ├── Engine_hypydata_Combined_2025_12_09.xlsx
│   └── Final FTIR/  (1,696 OPUS files)
└── processed/
    ├── modeling_response.csv
    └── horizons_data.rds
```

## Raw data

| File | Description | Source |
|------|-------------|--------|
| `Engine_hypydata_Combined_2025_12_09.xlsx` | HyPy SPAC measurements (436 samples, 7 projects) | Michelle Haddix |
| `Final FTIR/*.0` | MIR spectra (1,696 OPUS files, multiple replicates per sample) | CSU SoIL |

## Processed data

| File | Description |
|------|-------------|
| `modeling_response.csv` | Cleaned response data with PyC concentrations (307 samples after QC) |
| `horizons_data.rds` | Horizons data object with averaged spectra + joined response |

## Access

Contact sam.leuthold@colostate.edu for data access prior to publication.
