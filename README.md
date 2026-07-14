# ITGAX and AXL in tumour dendritic cells

Analysis of ITGAX and AXL co-expression in tumour-infiltrating dendritic cells,
using publicly available single-cell RNA-sequencing data.

## Summary

Dendritic cells from colorectal tumours (4,751 cells, 62 donors: 31 male,
31 female; CELLxGENE Discover census) were analysed at two levels:

- **Per-cell** — Spearman correlation between ITGAX and AXL across individual
  cells, with equal-n downsampling to control for unequal cell numbers.
- **Per-donor** — each patient treated as the unit of observation, avoiding the
  pseudoreplication that arises from the non-independence of cells within a
  donor. Includes pseudobulk AXL comparison, leave-one-donor-out sensitivity
  analysis, and within-donor correlation coefficients.

## Repository structure

```
.
├── README.md
├── requirements.txt
├── run_all.py                    # pipeline order
├── src/
│   ├── fetch_data.py             # retrieve DCs from CELLxGENE census
│   ├── analysis_per_cell.py      # single-cell level analysis
│   ├── analysis_per_donor.py     # donor level analysis
│   └── figures.py                # figure generation
├── results/                      # analysis outputs (CSV)
└── figures/                      # generated figures
```

## Data

Data are retrieved programmatically from the CELLxGENE Discover census
(version 2025-11-08) by `src/fetch_data.py`. No raw data are stored in this
repository.

| Tissue | Dataset ID | Cells | Donors |
|--------|------------|-------|--------|
| Colorectal | `16023185-de21-4c0d-a9c8-73abdd52d142` | 4,751 DCs | 31 M / 31 F, single study |
| Gastric | `e6aaf5a4-16e9-4ea6-9733-4eafd4e473d3` | 1,323 DCs | 18 M / 10 F, multi-study |

Census expression values are **raw counts** and are log1p-transformed within
the analysis scripts.

## Reproducing the analysis

```bash
pip install -r requirements.txt
python src/fetch_data.py           # retrieve data
python src/analysis_per_cell.py    # per-cell statistics
python src/analysis_per_donor.py   # per-donor statistics
python src/figures.py              # figures
```

`run_all.py` documents the intended order.

## Notes on the analysis

- Cells within a donor are not independent; per-cell P values are therefore
  anticonservative and are reported alongside the donor-level analysis, which
  treats each patient as the unit of observation.
- Within-donor correlations are restricted to donors contributing at least 10
  cells of the relevant subtype with non-zero variance in both genes.
- Paths at the top of each script point to the analysis environment and should
  be edited to match your setup.

## Requirements

Python 3.11; see `requirements.txt`.
