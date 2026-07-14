#!/usr/bin/env python3
"""
fetch_data.py - Retrieve tumour dendritic cells from the CELLxGENE census.

Downloads dendritic cells from colorectal (and optionally lung and gastric)
tumour datasets, restricted to the genes of interest, and writes AnnData
objects for downstream analysis.

Datasets (CELLxGENE Discover census, version 2025-11-08):
    colorectal : 16023185-de21-4c0d-a9c8-73abdd52d142
                 4,751 DCs, 62 donors (31 male / 31 female), single study
    lung       : LuCA core atlas (multi-study; see LUCA_DATASET below)
    gastric    : e6aaf5a4-16e9-4ea6-9733-4eafd4e473d3 (multi-study)

Note: census X values are RAW COUNTS and must be log1p-transformed before use.

Outputs (OUTDIR):
    CRC_DCs.h5ad, LuCA_core_DCs.h5ad, Gastric_DCs.h5ad

Usage:
    python fetch_data.py
"""
import os
import numpy as np
import pandas as pd
import anndata as ad_mod
import cellxgene_census

CENSUS_VERSION = "2025-11-08"
OUTDIR = "/data/external"

GENES = ["AXL", "ITGAX", "ITGAE", "CLEC9A", "XCR1"]

DC_TYPES = [
    "dendritic cell",
    "conventional dendritic cell",
    "plasmacytoid dendritic cell",
    "mature conventional dendritic cell",
    "CD1c-positive myeloid dendritic cell",
    "CD141-positive myeloid dendritic cell",
    "dendritic cell, human",
    "plasmacytoid dendritic cell, human",
]

DATASETS = {
    "CRC": {
        "dataset_id": "16023185-de21-4c0d-a9c8-73abdd52d142",
        "disease": "colon adenocarcinoma",
        "outfile": "CRC_DCs.h5ad",
    },
    "Gastric": {
        "dataset_id": "e6aaf5a4-16e9-4ea6-9733-4eafd4e473d3",
        "disease": "gastric cancer",
        "outfile": "Gastric_DCs.h5ad",
    },
}


def fetch(census, name, cfg):
    """Pull DCs for one dataset and write to disk."""
    types = ",".join(f"'{t}'" for t in DC_TYPES)
    obs_filter = (f"dataset_id == '{cfg['dataset_id']}' and "
                  f"disease == '{cfg['disease']}' and "
                  f"cell_type in [{types}]")
    adata = cellxgene_census.get_anndata(
        census, "homo_sapiens",
        obs_value_filter=obs_filter,
        var_value_filter=f"feature_name in {GENES}",
    )
    print(f"\n[{name}] cells: {adata.n_obs}")
    print(f"  genes: {adata.var['feature_name'].tolist()}")
    print(f"  sex: {adata.obs['sex'].value_counts().to_dict()}")
    print(f"  donors by sex: {adata.obs.groupby('sex')['donor_id'].nunique().to_dict()}")

    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    clean = ad_mod.AnnData(
        X=np.asarray(X),
        obs=pd.DataFrame(adata.obs).reset_index(drop=True),
        var=pd.DataFrame(adata.var).reset_index(drop=True),
    )
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, cfg["outfile"])
    clean.write(path)
    print(f"  saved {path}")
    return clean


def main():
    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)
    try:
        for name, cfg in DATASETS.items():
            fetch(census, name, cfg)
    finally:
        census.close()
    print("\ndone.")


if __name__ == "__main__":
    main()
