#!/usr/bin/env python3
"""
analysis_per_donor.py - Donor-level analysis of ITGAX and AXL in tumour DCs.

Treats each patient as the unit of observation, avoiding the pseudoreplication
that arises from the non-independence of cells sampled from the same donor.

Performs:
  1. AXL expression by sex, using per-donor pseudobulk (mean log-normalized
     expression per donor), compared by Mann-Whitney U
  2. Leave-one-donor-out sensitivity analysis of (1)
  3. Within-donor ITGAX-AXL Spearman correlation, tested against zero by
     one-sample Wilcoxon signed-rank, and compared between sexes by
     Mann-Whitney U

Donors contributing fewer than MIN_CELLS cells of a subtype, or with zero
variance in either gene, are excluded from the correlation analysis.

INPUT:
    /data/external/CRC_DCs.h5ad
OUTPUTS:
    results/per_donor_axl.csv          - per-donor pseudobulk AXL
    results/per_donor_rho.csv          - within-donor ITGAX-AXL rho
    results/per_donor_axl_stats.csv    - sex comparison summary
    results/leave_one_out.csv          - sensitivity analysis

Usage:
    python analysis_per_donor.py
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, mannwhitneyu, wilcoxon

DATA = "/data/external/CRC_DCs.h5ad"
OUTDIR = "results"
MIN_CELLS = 10


def load(path):
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
    d["AXL_log"] = np.log1p(d["AXL"])
    d["ITGAX_log"] = np.log1p(d["ITGAX"])
    return d


def pseudobulk_axl(d):
    """Mean log-normalized AXL per donor, per DC subtype."""
    rows = []
    for dc in ["cDC", "pDC"]:
        dd = d[d.DC == dc]
        ps = (dd.groupby(["donor", "sex"], observed=True)
                .agg(AXL=("AXL_log", "mean"), n_cells=("AXL_log", "size"))
                .reset_index())
        ps["DC"] = dc
        rows.append(ps)
    return pd.concat(rows, ignore_index=True)


def axl_by_sex(ps):
    """Compare per-donor AXL between sexes."""
    rows = []
    for dc in ["cDC", "pDC"]:
        p = ps[ps.DC == dc]
        m = p[p.sex == "male"]["AXL"].values
        f = p[p.sex == "female"]["AXL"].values
        if len(m) < 3 or len(f) < 3:
            continue
        stat, pval = mannwhitneyu(m, f)
        rows.append({"DC": dc, "n_male": len(m), "n_female": len(f),
                     "mean_male": m.mean(), "mean_female": f.mean(),
                     "median_male": np.median(m), "median_female": np.median(f),
                     "pval": pval})
    return pd.DataFrame(rows)


def leave_one_out(ps):
    """Does any single donor drive the sex comparison?"""
    rows = []
    for dc in ["cDC", "pDC"]:
        p = ps[ps.DC == dc]
        base = mannwhitneyu(p[p.sex == "male"]["AXL"].values,
                            p[p.sex == "female"]["AXL"].values)[1]
        for donor in p["donor"].unique():
            keep = p[p.donor != donor]
            m = keep[keep.sex == "male"]["AXL"].values
            f = keep[keep.sex == "female"]["AXL"].values
            if len(m) >= 3 and len(f) >= 3:
                pval = mannwhitneyu(m, f)[1]
                rows.append({"DC": dc, "donor_removed": donor,
                             "sex_removed": p[p.donor == donor]["sex"].iloc[0],
                             "p_baseline": base, "p_without": pval,
                             "delta": pval - base})
    return pd.DataFrame(rows)


def within_donor_rho(d, min_cells=MIN_CELLS):
    """Spearman rho between ITGAX and AXL computed inside each donor."""
    rows = []
    for dc in ["cDC", "pDC"]:
        dd = d[d.DC == dc]
        for (donor, sex), g in dd.groupby(["donor", "sex"], observed=True):
            if len(g) < min_cells:
                continue
            if g["ITGAX_log"].std() == 0 or g["AXL_log"].std() == 0:
                continue
            rho, p = spearmanr(g["ITGAX_log"], g["AXL_log"])
            if np.isnan(rho):
                continue
            rows.append({"DC": dc, "donor": donor, "sex": sex,
                         "rho": rho, "n_cells": len(g)})
    return pd.DataFrame(rows)


def rho_stats(r):
    """Is the correlation different from zero? Does it differ by sex?"""
    rows = []
    for dc in ["cDC", "pDC"]:
        sub = r[r.DC == dc]
        for lab in ["male", "female", "all"]:
            v = sub["rho"].values if lab == "all" else sub[sub.sex == lab]["rho"].values
            if len(v) < 5:
                continue
            stat, p = wilcoxon(v)
            rows.append({"DC": dc, "group": lab, "n_donors": len(v),
                         "n_positive": int((v > 0).sum()),
                         "median_rho": float(np.median(v)),
                         "wilcoxon_p": p})
        m = sub[sub.sex == "male"]["rho"].values
        f = sub[sub.sex == "female"]["rho"].values
        if len(m) >= 3 and len(f) >= 3:
            stat, p = mannwhitneyu(m, f)
            rows.append({"DC": dc, "group": "male_vs_female", "n_donors": len(m) + len(f),
                         "n_positive": np.nan,
                         "median_rho": float(np.median(m) - np.median(f)),
                         "wilcoxon_p": p})
    return pd.DataFrame(rows)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    d = load(DATA)

    ps = pseudobulk_axl(d)
    ps.to_csv(os.path.join(OUTDIR, "per_donor_axl.csv"), index=False)

    print("=" * 62)
    print("PER-DONOR: AXL expression by sex (pseudobulk)")
    print("=" * 62)
    stats = axl_by_sex(ps)
    print(stats.to_string(index=False))
    stats.to_csv(os.path.join(OUTDIR, "per_donor_axl_stats.csv"), index=False)

    print("\n" + "=" * 62)
    print("PER-DONOR: leave-one-donor-out sensitivity")
    print("=" * 62)
    loo = leave_one_out(ps)
    for dc in ["cDC", "pDC"]:
        s = loo[loo.DC == dc]
        if len(s):
            print(f"  {dc}: baseline p={s['p_baseline'].iloc[0]:.4f}, "
                  f"leave-one-out range {s['p_without'].min():.4f}-{s['p_without'].max():.4f}")
    loo.to_csv(os.path.join(OUTDIR, "leave_one_out.csv"), index=False)

    print("\n" + "=" * 62)
    print("PER-DONOR: within-donor ITGAX-AXL correlation")
    print("=" * 62)
    r = within_donor_rho(d)
    r.to_csv(os.path.join(OUTDIR, "per_donor_rho.csv"), index=False)
    rs = rho_stats(r)
    print(rs.to_string(index=False))
    rs.to_csv(os.path.join(OUTDIR, "per_donor_rho_stats.csv"), index=False)

    print(f"\nsaved results to {OUTDIR}/")


if __name__ == "__main__":
    main()
