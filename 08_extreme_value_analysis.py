"""
Extreme value analysis: Peaks-over-Threshold (POT) + Generalized Pareto fit.

Different question from Stage 06's bulk Weibull fit: not "what's the
everyday Hs distribution," but "how large can storm peaks get" - the
standard approach for engineering design values (e.g. "100-year Hs").

Runs on the RAW cleaned level series (Stage 0 output), same reasoning as
Stage 06: these are physical extremes of Hs itself, not of a residual.

CAUTION: the record here is a few months. Return levels for periods much
longer than the record are extrapolations with large uncertainty - the
script flags any return period beyond ~2x the record length rather than
presenting them as equally trustworthy.

Usage:
    python 08_extreme_value_analysis.py --buoy WesthinderBuoy --var VHM0 \
        --min-separation-hours 48 --threshold-percentile 95
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import find_peaks

from utils import default_paths


def decluster_peaks(s: pd.Series, threshold: float, min_separation_hours: float) -> pd.Series:
    """Keep only the max value within each min_separation_hours window above
    threshold, so each storm contributes one independent peak."""
    dt_hours = (s.index[1] - s.index[0]).total_seconds() / 3600.0
    min_distance_samples = max(1, int(round(min_separation_hours / dt_hours)))

    vals = s.values.copy()
    vals[np.isnan(vals)] = -np.inf
    peak_idx, _ = find_peaks(vals, height=threshold, distance=min_distance_samples)
    return pd.Series(s.values[peak_idx], index=s.index[peak_idx], name="peak")


def mean_residual_life(data: np.ndarray, thresholds: np.ndarray):
    means, ses = [], []
    for u in thresholds:
        excess = data[data > u] - u
        if len(excess) < 5:
            means.append(np.nan)
            ses.append(np.nan)
            continue
        means.append(excess.mean())
        ses.append(excess.std(ddof=1) / np.sqrt(len(excess)))
    return np.array(means), np.array(ses)


def return_level(u, sigma, xi, lam_per_year, return_period_years):
    """x_N solving: expected exceedances of x_N in T years = 1, under a GPD
    excess model with exceedance rate lam_per_year above threshold u."""
    m = lam_per_year * return_period_years
    if abs(xi) < 1e-8:
        return u + sigma * np.log(m)
    return u + (sigma / xi) * (m ** xi - 1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--min-separation-hours", default=48.0, type=float,
                         help="minimum time between independent storm peaks")
    parser.add_argument("--threshold-percentile", default=95.0, type=float,
                         help="percentile of the raw series used as the GPD threshold")
    parser.add_argument("--return-periods-years", default="1,5,10,25,50",
                         help="comma-separated return periods to report, in years")
    args = parser.parse_args()

    in_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    s = pd.read_csv(in_path, index_col=0, parse_dates=True)[args.var].dropna()

    record_years = (s.index[-1] - s.index[0]).total_seconds() / (3600 * 24 * 365.25)
    print(f"Record length: {record_years:.3f} years ({len(s)} samples)")

    out_dir = default_paths("08_extreme_value_analysis")

    # --- Threshold selection diagnostic: mean residual life plot ---
    candidate_thresholds = np.linspace(
        np.percentile(s.values, 80), np.percentile(s.values, 99), 25
    )
    mrl_means, mrl_ses = mean_residual_life(s.values, candidate_thresholds)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(candidate_thresholds, mrl_means, marker="o", ms=3)
    ax.fill_between(candidate_thresholds, mrl_means - 1.96 * mrl_ses,
                     mrl_means + 1.96 * mrl_ses, alpha=0.2)
    ax.set_xlabel("Threshold u (m)")
    ax.set_ylabel("Mean excess above u")
    ax.set_title(f"{args.buoy} — mean residual life plot\n"
                 f"(pick u where the curve is roughly linear/stable)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_mean_residual_life.png", dpi=150)
    plt.close(fig)

    threshold = float(np.percentile(s.values, args.threshold_percentile))
    print(f"Using threshold u = {threshold:.3f} m ({args.threshold_percentile}th percentile) "
          f"— check *_mean_residual_life.png to confirm this sits in a stable region")

    # --- Decluster to independent storm peaks ---
    peaks = decluster_peaks(s, threshold, args.min_separation_hours)
    print(f"{len(peaks)} independent storm peaks above threshold "
          f"(min separation {args.min_separation_hours}h)")

    if len(peaks) < 10:
        print("WARNING: fewer than 10 exceedances - GPD fit will be unstable. "
              "Consider lowering --threshold-percentile or --min-separation-hours.")
        if len(peaks) < 3:
            print("Too few peaks to fit. Stopping.")
            return

    excess = peaks.values - threshold

    # --- Fit GPD (loc fixed at 0: modeling excess above threshold) ---
    shape, loc, scale = stats.genpareto.fit(excess, floc=0)
    print(f"GPD fit: shape (xi)={shape:.4f}, scale (sigma)={scale:.4f}")
    if shape < -0.05:
        print("xi < 0 -> bounded upper tail (finite maximum Hs implied)")
    elif shape > 0.05:
        print("xi > 0 -> heavy tail (no finite upper bound implied)")
    else:
        print("xi ~ 0 -> roughly exponential tail")

    ks_stat, ks_p = stats.kstest(excess, "genpareto", args=(shape, loc, scale))
    print(f"KS on exceedances vs fitted GPD: stat={ks_stat:.4f}, p={ks_p:.4f}")

    # --- Return levels ---
    lam_per_year = len(peaks) / record_years
    return_periods = [float(x) for x in args.return_periods_years.split(",")]
    rl_rows = []
    print("\nReturn levels:")
    for rp in return_periods:
        rl = return_level(threshold, scale, shape, lam_per_year, rp)
        illustrative = rp > 2 * record_years
        rl_rows.append({
            "return_period_years": rp,
            "return_level_m": rl,
            "illustrative_only": illustrative,
        })
        flag = "  <-- beyond ~2x record length, illustrative only" if illustrative else ""
        print(f"  {rp:6.1f}-year: Hs = {rl:.2f} m{flag}")

    rl_df = pd.DataFrame(rl_rows)
    rl_df.to_csv(out_dir / f"{args.buoy}_{args.var}_return_levels.csv", index=False)

    import json
    with open(out_dir / f"{args.buoy}_{args.var}_eva_summary.json", "w") as f:
        json.dump({
            "record_years": record_years,
            "threshold": threshold,
            "n_peaks": int(len(peaks)),
            "gpd_shape": float(shape),
            "gpd_scale": float(scale),
            "ks_pvalue": float(ks_p),
            "lam_per_year": float(lam_per_year),
            "fit_reliable": bool(len(peaks) >= 10),
        }, f, indent=2)

    # --- Plot: declustered peaks over the series ---
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(s.index, s.values, lw=0.4, color="steelblue", label=args.var)
    ax.scatter(peaks.index, peaks.values, color="crimson", zorder=3, s=25,
               label="independent storm peaks")
    ax.axhline(threshold, color="gray", ls="--", label=f"threshold ({threshold:.2f} m)")
    ax.legend()
    ax.set_title(f"{args.buoy} — declustered storm peaks above threshold")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_declustered_peaks.png", dpi=150)
    plt.close(fig)

    # --- Plot: GPD fit vs. observed exceedances ---
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.linspace(0, excess.max() * 1.2, 200)
    ax.hist(excess, bins=min(15, max(5, len(excess) // 3)), density=True,
            alpha=0.4, color="gray", label="observed excesses")
    ax.plot(x, stats.genpareto.pdf(x, shape, loc, scale), color="firebrick", lw=2,
            label="fitted GPD")
    ax.set_xlabel("Excess above threshold (m)")
    ax.set_ylabel("density")
    ax.legend()
    ax.set_title(f"{args.buoy} — GPD fit to exceedances (KS p={ks_p:.3f})")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_gpd_fit.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved diagnostics + return levels to {out_dir}")


if __name__ == "__main__":
    main()
