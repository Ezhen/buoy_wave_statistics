"""
Stage 4b - Distribution fitting for storm climatology characterization.

Fits Rayleigh, Weibull, and log-normal to the RAW (cleaned) level series -
not the differenced residual. These distributions describe the physical
Hs climatology (always positive, right-skewed); the residual from Stage 3
is a different object used only for the whiteness check.

Usage:
    python 06_distribution_fit.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from utils import default_paths

CANDIDATES = {
    "rayleigh": stats.rayleigh,
    "weibull": stats.weibull_min,
    "lognormal": stats.lognorm,
}


def fit_and_test(data, name, dist):
    params = dist.fit(data)
    ks_stat, ks_p = stats.kstest(data, dist.name, args=params)
    try:
        ad_result = stats.anderson(data, dist=dist.name if dist.name in
                                    ("norm", "expon", "logistic", "gumbel", "gumbel_l", "gumbel_r")
                                    else "norm")  # AD only supports a few named dists directly
        ad_stat = ad_result.statistic
    except Exception:
        ad_stat = None
    return {"name": name, "params": params, "ks_stat": ks_stat, "ks_pvalue": ks_p, "ad_stat": ad_stat}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    args = parser.parse_args()

    in_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    s = pd.read_csv(in_path, index_col=0, parse_dates=True)[args.var].dropna()
    data = s.values
    data = data[data > 0]  # distributions below require strictly positive support

    out_dir = default_paths("06_distribution_fit")

    results = []
    for name, dist in CANDIDATES.items():
        r = fit_and_test(data, name, dist)
        results.append(r)
        print(f"{name:10s}: KS stat={r['ks_stat']:.4f}  p={r['ks_pvalue']:.4f}")

    best = min(results, key=lambda r: r["ks_stat"])
    print(f"\nBest KS fit (lowest statistic): {best['name']}")
    print("Note: with a few months of data, treat the p-value as a rough guide - "
          "compare the Q-Q plots below visually before committing to one.")

    # --- Q-Q plots, one panel per candidate ---
    fig, axes = plt.subplots(1, len(CANDIDATES), figsize=(5 * len(CANDIDATES), 5))
    for ax, r in zip(axes, results):
        dist = CANDIDATES[r["name"]]
        stats.probplot(data, dist=dist, sparams=r["params"], plot=ax)
        ax.set_title(f"{r['name']} (KS p={r['ks_pvalue']:.3f})")
    fig.suptitle(f"{args.buoy} — {args.var} Q-Q plots vs candidate distributions")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_qq_plots.png", dpi=150)
    plt.close(fig)

    # --- Histogram with fitted PDFs overlaid ---
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(data, bins=40, density=True, alpha=0.4, color="gray", label="observed")
    x = np.linspace(data.min(), data.max(), 300)
    for r in results:
        dist = CANDIDATES[r["name"]]
        ax.plot(x, dist.pdf(x, *r["params"]), label=r["name"], lw=2)
    ax.set_xlabel(args.var)
    ax.set_ylabel("density")
    ax.set_title(f"{args.buoy} — {args.var} distribution fits")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_fitted_pdfs.png", dpi=150)
    plt.close(fig)

    summary = pd.DataFrame([
        {"distribution": r["name"], "ks_stat": r["ks_stat"], "ks_pvalue": r["ks_pvalue"]}
        for r in results
    ])
    summary.to_csv(out_dir / f"{args.buoy}_{args.var}_fit_summary.csv", index=False)

    print(f"\nSaved plots + summary to {out_dir}")


if __name__ == "__main__":
    main()
