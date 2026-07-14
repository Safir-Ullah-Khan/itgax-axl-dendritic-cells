#!/usr/bin/env python3
"""
analysis_per_cell.py - Single-cell level analysis of ITGAX and AXL in tumour DCs.

Computes, at single-cell resolution:
  1. Spearman correlation between ITGAX and AXL, per sex and DC subtype
  2. The same correlations after downsampling both sexes to equal cell number,
     to confirm that differences are not driven by unequal sampling

Note: cells within a donor are not independent, so the P values reported here
are anticonservative. See analysis_per_donor.py for the donor-level analysis
that treats each patient as the unit of observation.

INPUT:
    /data/external/CRC_DCs.h5ad   (census raw counts; log1p applied here)
OUTPUT:
    results/per_cell_correlations.csv

Usage:
    python analysis_per_cell.py
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

DATA = "/data/external/CRC_DCs.h5ad"
OUTDIR = "results"
SEED = 0


def load(path):
    """Load DC object, log-transform, and label cDC / pDC."""
    ad = sc.read_h5ad(path)

    def gv(g):
        i = list(ad.var["feature_name"]).index(g)
        x = ad.X[:, i]
        return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()

    d = pd.DataFrame({
        "donor": ad.obs["donor_id"].astype(str).values,
        "sex": ad.obs["sex"].astype(str).values,
        "cell_type": ad.obs["cell_type"].astype(str).values,
        "AXL": gv("AXL"),
        "ITGAX": gv("ITGAX"),
    })
    d["DC"] = d["cell_type"].apply(lambda c: "pDC" if "plasmacytoid" in c else "cDC")
    d = d[d["sex"].isin(["male", "female"])].copy()
    # census X is raw counts
    d["AXL_log"] = np.log1p(d["AXL"])
    d["ITGAX_log"] = np.log1p(d["ITGAX"])
    return d


def correlations(d):
    """Spearman rho between ITGAX and AXL, per DC subtype and sex."""
    rows = []
    for dc in ["cDC", "pDC"]:
        for sex in ["male", "female"]:
            s = d[(d.DC == dc) & (d.sex == sex)]
            if len(s) < 30 or s["ITGAX_log"].std() == 0 or s["AXL_log"].std() == 0:
                continue
            rho, p = spearmanr(s["ITGAX_log"], s["AXL_log"])
            rows.append({"DC": dc, "sex": sex, "n_cells": len(s),
                         "rho": rho, "pval": p})
    return pd.DataFrame(rows)


def equal_n_correlations(d, seed=SEED):
    """Recompute correlations after downsampling both sexes to equal n."""
    rows = []
    for dc in ["cDC", "pDC"]:
        m = d[(d.DC == dc) & (d.sex == "male")]
        f = d[(d.DC == dc) & (d.sex == "female")]
        n = min(len(m), len(f))
        if n < 30:
            continue
        ms = m.sample(n, random_state=seed)
        fs = f.sample(n, random_state=seed)
        rm = spearmanr(ms["ITGAX_log"], ms["AXL_log"])[0]
        rf = spearmanr(fs["ITGAX_log"], fs["AXL_log"])[0]
        rows.append({"DC": dc, "n_each": n, "rho_male": rm,
                     "rho_female": rf, "difference": rm - rf})
    return pd.DataFrame(rows)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    d = load(DATA)

    print("=" * 62)
    print("PER-CELL: ITGAX-AXL correlation by sex")
    print("=" * 62)
    corr = correlations(d)
    print(corr.to_string(index=False))

    print("\n" + "=" * 62)
    print("PER-CELL: equal-n downsampled correlations")
    print("=" * 62)
    eq = equal_n_correlations(d)
    print(eq.to_string(index=False))

    corr.to_csv(os.path.join(OUTDIR, "per_cell_correlations.csv"), index=False)
    eq.to_csv(os.path.join(OUTDIR, "per_cell_correlations_equal_n.csv"), index=False)
    print(f"\nsaved -> {OUTDIR}/per_cell_correlations.csv")
    print(f"saved -> {OUTDIR}/per_cell_correlations_equal_n.csv")


if __name__ == "__main__":
    main()
