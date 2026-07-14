#!/usr/bin/env python3
"""
run_all.py - Pipeline order for the ITGAX / AXL dendritic cell analysis.

Each script is self-contained; edit the paths at the top of each to match your
environment. Run this file to print the intended order.
"""
STEPS = [
    ("src/fetch_data.py",         "Retrieve tumour DCs from the CELLxGENE census"),
    ("src/analysis_per_cell.py",  "Per-cell ITGAX-AXL correlations (+ equal-n check)"),
    ("src/analysis_per_donor.py", "Per-donor AXL by sex, leave-one-out, within-donor rho"),
    ("src/figures.py",            "Generate per-donor and per-cell figures"),
]

if __name__ == "__main__":
    print(f"{'Script':<32}Description")
    print("-" * 84)
    for s, d in STEPS:
        print(f"{s:<32}{d}")
    print("\nRun in order, e.g.:  python src/fetch_data.py")
