"""
Confidence intervals on Stage 11's pairwise buoy correlations.

The classic Fisher z CI assumes each sample contributes independent
information. With Hs autocorrelated at every buoy (Ljung-Box fails
everywhere), the raw N overstates how much independent information a
correlation pair actually has - so a plain Fisher z CI here would look
artificially tight.

Correction used: a simplified AR(1)-based effective sample size
    N_eff = N * (1 - r1*r2) / (1 + r1*r2)
where r1, r2 are the two series' lag-1 autocorrelations (from Stage 11b).
This is the Dawdy-Matalas approximation for two AR(1)-like series - a
simplification (it only uses lag-1, not the full ACF), stated explicitly
so it isn't over-trusted as more rigorous than it is. It needs only
each buoy's lag1_acf, which Stage 11b already computes, rather than the
full cross-series Bartlett formula which would need the complete ACF at
matching lags for every pair.

Usage:
    python 12b_correlation_confidence.py --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from utils import default_paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", default="VHM0")
    args = parser.parse_args()

    corr_path = Path("pipeline_out/11_spatial_statistics") / f"{args.var}_pairwise_correlation.csv"
    if not corr_path.exists():
        print(f"No Stage 11 output found at {corr_path} - run 11_spatial_statistics.py first.")
        return

    corr = pd.read_csv(corr_path, index_col=0)
    buoys = corr.columns.tolist()

    lag1 = {}
    dep_dir = Path("pipeline_out/11b_dependence_structure")
    for buoy in buoys:
        dep_path = dep_dir / f"{buoy}_{args.var}_dependence_summary.json"
        if dep_path.exists():
            with open(dep_path) as f:
                lag1[buoy] = json.load(f)["lag1_acf"]

    missing = [b for b in buoys if b not in lag1]
    if missing:
        print(f"WARNING: {len(missing)} buoy(s) missing Stage 11b output - their pairs "
              f"will be skipped: {missing}")

    # Approximate N per pair from Stage 11's min_periods behavior: use the
    # smaller of the two buoys' Stage 0 sample counts as a stand-in for the
    # pair's overlap size (exact overlap count isn't saved by Stage 11).
    n_samples = {}
    for buoy in buoys:
        load_path = Path("pipeline_out/01_load_clean") / f"{buoy}_{args.var}_load_summary.json"
        if load_path.exists():
            with open(load_path) as f:
                n_samples[buoy] = json.load(f)["n_samples_regularized"]

    rows = []
    for i, b1 in enumerate(buoys):
        for b2 in buoys[i + 1:]:
            r = corr.loc[b1, b2]
            if pd.isna(r) or b1 not in lag1 or b2 not in lag1:
                continue
            r1, r2 = lag1[b1], lag1[b2]
            N = min(n_samples.get(b1, 2850), n_samples.get(b2, 2850))
            n_eff = N * (1 - r1 * r2) / (1 + r1 * r2)
            n_eff = max(4.0, n_eff)  # Fisher z needs n_eff > 3

            r_clamped = np.clip(r, -0.9999, 0.9999)
            z = np.arctanh(r_clamped)
            se = 1.0 / np.sqrt(n_eff - 3)
            z_lo, z_hi = z - 1.96 * se, z + 1.96 * se
            r_lo, r_hi = np.tanh(z_lo), np.tanh(z_hi)

            rows.append({
                "buoy1": b1, "buoy2": b2, "correlation": r,
                "n_raw": N, "n_effective": round(n_eff, 1),
                "ci_low": r_lo, "ci_high": r_hi,
                "ci_width": r_hi - r_lo,
            })

    df = pd.DataFrame(rows)
    out_dir = default_paths("12b_correlation_confidence")
    df.to_csv(out_dir / f"{args.var}_correlation_ci.csv", index=False)

    print(f"--- Correlation CIs, {len(df)} pairs ---")
    print(f"Mean N_raw: {df['n_raw'].mean():.0f}, mean N_effective: {df['n_effective'].mean():.0f} "
          f"({100 * (1 - df['n_effective'].mean() / df['n_raw'].mean()):.0f}% reduction from "
          f"raw N, on average, due to autocorrelation)")
    print(f"Mean CI width: {df['ci_width'].mean():.3f}")

    widest = df.loc[df["ci_width"].idxmax()]
    narrowest = df.loc[df["ci_width"].idxmin()]
    print(f"\nWidest CI: {widest['buoy1']} vs {widest['buoy2']}, "
          f"r={widest['correlation']:.3f} [{widest['ci_low']:.3f}, {widest['ci_high']:.3f}]")
    print(f"Narrowest CI: {narrowest['buoy1']} vs {narrowest['buoy2']}, "
          f"r={narrowest['correlation']:.3f} [{narrowest['ci_low']:.3f}, {narrowest['ci_high']:.3f}]")

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
