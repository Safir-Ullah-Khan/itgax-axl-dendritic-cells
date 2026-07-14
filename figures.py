#!/usr/bin/env python3
"""
figures.py - Figures for the ITGAX / AXL dendritic cell analysis.

Generates two figures:

  per-donor figure (4 panels):
    a  AXL in cDCs by sex (per-donor pseudobulk means)
    b  AXL in pDCs by sex
    c  within-donor ITGAX-AXL Spearman rho, each point one donor
    d  single-cell ITGAX vs AXL scatter in pDCs

  per-cell figure (3 panels):
    a  AXL in cDCs by sex
    b  AXL in pDCs by sex
    c  single-cell ITGAX vs AXL scatter in pDCs

INPUTS:
    /data/external/CRC_DCs.h5ad
    results/per_donor_rho.csv   (produced by analysis_per_donor.py)
OUTPUTS:
    figures/CRC_AXL_ITGAX_perdonor.{png,pdf}
    figures/CRC_AXL_ITGAX_percell.{png,pdf}

Usage:
    python figures.py
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr

DATA = "/data/external/CRC_DCs.h5ad"
RHO_CSV = "results/per_donor_rho.csv"
FIGDIR = "figures"

BLUE, PINK = "#3B6FB6", "#C74B8B"
SCAT_M, SCAT_F = "#4A6FA5", "#B5638F"

mpl.rcParams.update({
    "font.size": 11, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.2, "xtick.major.width": 1.2, "ytick.major.width": 1.2,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "axes.labelsize": 12,
    "pdf.fonttype": 42, "figure.dpi": 120,
})


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
        "AXL": gv("AXL"), "ITGAX": gv("ITGAX"),
    })
    d["DC"] = d["cell_type"].apply(lambda c: "pDC" if "plasmacytoid" in c else "cDC")
    d = d[d["sex"].isin(["male", "female"])].copy()
    d["AXL_log"] = np.log1p(d["AXL"])
    d["ITGAX_log"] = np.log1p(d["ITGAX"])
    return d


def _plabel(ax, letter):
    ax.text(-0.26, 1.10, letter, transform=ax.transAxes, fontsize=17,
            fontweight="bold", va="top", ha="left")


def _axl_bar(ax, d, dc, letter, title):
    _plabel(ax, letter)
    dd = d[d.DC == dc]
    ps = dd.groupby(["donor", "sex"], observed=True)["AXL_log"].mean().reset_index()
    m = ps[ps.sex == "male"]["AXL_log"].values
    f = ps[ps.sex == "female"]["AXL_log"].values
    means = [m.mean(), f.mean()]
    ax.bar([0, 1], means, width=0.6, color=[BLUE, PINK], alpha=0.85,
           edgecolor="black", linewidth=1.2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Male\n(n={len(m)})", f"Female\n(n={len(f)})"],
                       fontweight="bold", fontsize=12)
    ax.set_ylabel("AXL expression\n(log-normalized, per donor)",
                  fontweight="bold", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=14)
    ax.set_ylim(0, max(means) * 1.28)


def _percell_scatter(ax, d, letter, title):
    _plabel(ax, letter)
    dd = d[d.DC == "pDC"]
    rhos = {}
    for sex, col, lab in [("male", SCAT_M, "Male"), ("female", SCAT_F, "Female")]:
        s = dd[dd.sex == sex]
        rr, _ = spearmanr(s["ITGAX_log"], s["AXL_log"])
        rhos[lab] = (rr, col)
        ax.scatter(s["ITGAX_log"], s["AXL_log"], s=9, alpha=0.28, c=col,
                   edgecolors="none", rasterized=True)
        z = np.polyfit(s["ITGAX_log"], s["AXL_log"], 1)
        xs = np.linspace(s["ITGAX_log"].min(), s["ITGAX_log"].max(), 50)
        ax.plot(xs, np.polyval(z, xs), c=col, lw=2.6)
    ax.set_xlabel("ITGAX expression (log)", fontweight="bold", fontsize=12)
    ax.set_ylabel("AXL expression (log)", fontweight="bold", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=14)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.28)
    ypos = 0.96
    for lab, (rr, col) in rhos.items():
        ax.text(0.09, ypos, f"{lab}   \u03c1 = {rr:.2f}", transform=ax.transAxes,
                fontsize=11.5, fontweight="bold", color=col, va="top", ha="left")
        ypos -= 0.085
    for i, (lab, (rr, col)) in enumerate(rhos.items()):
        ax.scatter([0.045], [0.945 - i * 0.085], s=42, c=col, transform=ax.transAxes,
                   clip_on=False, edgecolors="none")


def _perdonor_rho(ax, rho_df, letter):
    _plabel(ax, letter)
    r = rho_df[rho_df.DC == "pDC"]
    rng = np.random.default_rng(1)
    for i, (sex, col) in enumerate([("male", BLUE), ("female", PINK)]):
        v = r[r.sex == sex]["rho"].values
        ax.scatter(rng.normal(i, 0.09, len(v)), v, s=42, c=col, alpha=0.75,
                   edgecolors="white", linewidth=0.6, zorder=3)
        ax.hlines(np.median(v), i - 0.26, i + 0.26, colors="black", lw=2.6, zorder=4)
    ax.axhline(0, ls="--", c="#888", lw=1.2, zorder=1)
    ax.set_xticks([0, 1])
    m_n = (r.sex == "male").sum()
    f_n = (r.sex == "female").sum()
    ax.set_xticklabels([f"Male\n(n={m_n})", f"Female\n(n={f_n})"],
                       fontweight="bold", fontsize=12)
    ax.set_ylabel("Within-donor ITGAX\u2013AXL\nSpearman \u03c1",
                  fontweight="bold", fontsize=12)
    ax.set_title("Per-donor correlation (pDCs)", fontsize=13, fontweight="bold", pad=14)
    ax.set_xlim(-0.6, 1.6)
    npos = int((r["rho"] > 0).sum())
    ntot = len(r)
    med = r["rho"].median()
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + (ymax - ymin) * 0.30)
    ax.text(0.97, 0.97, f"{npos}/{ntot} donors +\nmedian \u03c1 = {med:+.2f}\nP < 0.0001",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            fontweight="bold", linespacing=1.5)


def figure_per_donor(d, rho_df):
    fig = plt.figure(figsize=(15, 5.2))
    gs = GridSpec(1, 4, figure=fig, wspace=0.50, left=0.055, right=0.985,
                  top=0.80, bottom=0.17, width_ratios=[1, 1, 1.35, 1.2])
    _axl_bar(fig.add_subplot(gs[0, 0]), d, "cDC", "a", "AXL in cDCs")
    _axl_bar(fig.add_subplot(gs[0, 1]), d, "pDC", "b", "AXL in pDCs")
    _perdonor_rho(fig.add_subplot(gs[0, 2]), rho_df, "c")
    _percell_scatter(fig.add_subplot(gs[0, 3]), d, "d", "Single-cell view (pDCs)")
    fig.suptitle("AXL expression and the ITGAX\u2013AXL correlation in "
                 "colorectal tumour dendritic cells",
                 fontsize=14, fontweight="bold", y=0.95)
    _save(fig, "CRC_AXL_ITGAX_perdonor")


def figure_per_cell(d):
    fig = plt.figure(figsize=(12, 5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.38, left=0.08, right=0.98,
                  top=0.84, bottom=0.16)
    _axl_bar(fig.add_subplot(gs[0, 0]), d, "cDC", "a", "AXL in cDCs")
    _axl_bar(fig.add_subplot(gs[0, 1]), d, "pDC", "b", "AXL in pDCs")
    _percell_scatter(fig.add_subplot(gs[0, 2]), d, "c",
                     "ITGAX\u2013AXL correlation in pDCs")
    fig.suptitle("ITGAX and AXL expression in colorectal tumour dendritic cells",
                 fontsize=14, fontweight="bold", y=0.97)
    _save(fig, "CRC_AXL_ITGAX_percell")


def _save(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    out = os.path.join(FIGDIR, name)
    fig.savefig(out + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}.png / .pdf")


def main():
    d = load(DATA)
    rho_df = pd.read_csv(RHO_CSV)
    figure_per_donor(d, rho_df)
    figure_per_cell(d)


if __name__ == "__main__":
    main()
