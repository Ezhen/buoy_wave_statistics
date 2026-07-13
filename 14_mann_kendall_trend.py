"""
Mann-Kendall trend test, with the Hamed-Rao variance correction for
autocorrelation - NOT the textbook version.

Why the correction is mandatory here, not optional: Ljung-Box confirms
strong serial correlation at every buoy in this network, persisting to
at least 24-48h lags at 30-min sampling. That's exactly the condition
known to inflate plain Mann-Kendall's false-positive rate on trend
detection - running the textbook test here would be a real regression
in rigor relative to everything else in this pipeline (the KS caveat,
the Fisher z effective-N correction, the block bootstrap block length -
all exist for the same underlying reason).

Runs on ANNUAL aggregates, not raw high-frequency data - two reasons:
  1. Computational: Mann-Kendall's S statistic is O(N^2) pairwise
     comparisons. At 527,620 raw samples (Westhinder), that's ~1.4e11
     pairs - infeasible. At ~36 annual points, trivial.
  2. Statistical: "is the wave climate trending" is a question about
     annual/seasonal-scale structure, not 30-min noise. This is also
     standard practice for climate/environmental trend detection (WMO
     guidelines use annual or seasonal aggregates, not raw
     high-frequency series).

Two separate trend tests per buoy: annual MEAN Hs (typical wave
climate) and annual 95th-percentile Hs (storm intensity) - these can
trend differently and answer different questions.

Only meaningful for buoys with a long enough record - default minimum
10 years, but treat results under ~20-30 years with real caution; this
network's 6 long-record buoys (30+ years) are the ones this was really
built for.

Usage:
    python 14_mann_kendall_trend.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from utils import default_paths


def mann_kendall_s(x: np.ndarray):
    """Classic MK S statistic. O(n^2) - fine at annual-aggregate sample
    sizes (tens of points)."""
    n = len(x)
    s = 0
    for i in range(n - 1):
        s += np.sum(np.sign(x[i + 1:] - x[i]))
    return s


def hamed_rao_variance_correction(x: np.ndarray):
    """Effective-sample-size correction factor n/n* from Hamed & Rao
    (1998), using the autocorrelation of the RANKS of x (not the raw
    values - standard for this correction), including only lags whose
    autocorrelation is statistically significant at the 95% level.
    Returns the correction factor (>= 1) to multiply the naive
    Var(S) by."""
    n = len(x)
    ranks = pd.Series(x).rank().values
    ranks_centered = ranks - ranks.mean()
    denom = np.sum(ranks_centered ** 2)

    # significance threshold for autocorrelation of ranks under IID null
    sig_threshold = (-1 + 1.96 * np.sqrt(n - 2)) / (n - 1)

    correction_sum = 0.0
    for k in range(1, n - 1):
        num = np.sum(ranks_centered[:n - k] * ranks_centered[k:])
        rho_k = num / denom if denom > 0 else 0.0
        if rho_k > sig_threshold:  # only significant positive autocorrelation counted
            correction_sum += (n - k) * (n - k - 1) * (n - k - 2) * rho_k

    n_star_ratio = 1 + (2.0 / (n * (n - 1) * (n - 2))) * correction_sum
    return max(1.0, n_star_ratio)  # correction can't reduce variance below the naive value


def sens_slope(t_years: np.ndarray, x: np.ndarray):
    """Median of all pairwise slopes - robust trend magnitude estimate,
    in units of x per year."""
    n = len(x)
    slopes = []
    for i in range(n - 1):
        dt = t_years[i + 1:] - t_years[i]
        dx = x[i + 1:] - x[i]
        valid = dt != 0
        slopes.extend((dx[valid] / dt[valid]).tolist())
    return float(np.median(slopes)) if slopes else float("nan")


def mann_kendall_hamed_rao(t_years: np.ndarray, x: np.ndarray, alpha: float = 0.05):
    n = len(x)
    if n < 8:
        return {"error": f"only {n} points - too few for a meaningful trend test (need >= 8)"}

    s = mann_kendall_s(x)
    var_s_naive = n * (n - 1) * (2 * n + 5) / 18.0

    correction = hamed_rao_variance_correction(x)
    var_s_corrected = var_s_naive * correction

    if s > 0:
        z = (s - 1) / np.sqrt(var_s_corrected)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s_corrected)
    else:
        z = 0.0

    p_value = 2 * (1 - norm.cdf(abs(z)))
    slope = sens_slope(t_years, x)

    if p_value < alpha:
        trend = "increasing" if s > 0 else "decreasing"
    else:
        trend = "no significant trend"

    return {
        "n": n, "S": int(s), "var_s_naive": float(var_s_naive),
        "hamed_rao_correction_factor": float(correction),
        "var_s_corrected": float(var_s_corrected),
        "z": float(z), "p_value": float(p_value),
        "sens_slope_per_year": slope, "trend": trend,
        "significant_at_alpha": bool(p_value < alpha),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--min-years", default=10.0, type=float,
                         help="refuse to run below this record length - trend "
                              "detection needs real multi-year coverage")
    parser.add_argument("--alpha", default=0.05, type=float)
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"

    hs = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var].dropna()

    record_years = None
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            record_years = json.load(f).get("record_years")
    if record_years is None:
        record_years = (hs.index[-1] - hs.index[0]).days / 365.25

    print(f"--- {args.buoy} / {args.var} Mann-Kendall trend test (Hamed-Rao corrected) ---")
    print(f"Record length: {record_years:.1f} years")
    if record_years < args.min_years:
        print(f"Record is under --min-years={args.min_years} - refusing to run. "
              f"Trend detection on a short record is unreliable regardless of "
              f"the autocorrelation correction; this isn't a sample-size problem "
              f"the correction can fix.")
        return
    if record_years < 20:
        print(f"NOTE: {record_years:.1f} years clears --min-years but is still "
              f"short by climate-trend-detection convention (10-20+ years is "
              f"typical, 20-30+ preferred). Treat any result here with real caution.")

    annual = hs.groupby(hs.index.year)
    annual_mean = annual.mean()
    annual_p95 = annual.quantile(0.95)
    annual_n = annual.count()

    # Flag any year with unusually sparse coverage (year got hit hard by a gap)
    expected_per_year = 365.25 * 24 / (
        (hs.index[1] - hs.index[0]).total_seconds() / 3600 if len(hs) > 1 else 0.5)
    sparse_years = annual_n[annual_n < 0.5 * expected_per_year].index.tolist()
    if sparse_years:
        print(f"\nNOTE: {len(sparse_years)} year(s) have <50% expected coverage "
              f"(gap-heavy): {sparse_years}. Their annual mean/p95 may be noisy "
              f"or unrepresentative of that full year.")

    t_years = annual_mean.index.values.astype(float)

    print(f"\n[Annual mean Hs - typical wave climate]")
    result_mean = mann_kendall_hamed_rao(t_years, annual_mean.values, args.alpha)
    if "error" in result_mean:
        print(f"  {result_mean['error']}")
    else:
        print(f"  n={result_mean['n']} years, S={result_mean['S']}, "
              f"Hamed-Rao correction factor={result_mean['hamed_rao_correction_factor']:.2f}")
        print(f"  Z={result_mean['z']:.3f}, p={result_mean['p_value']:.4f} "
              f"-> {result_mean['trend']}")
        print(f"  Sen's slope: {result_mean['sens_slope_per_year']:+.4f} m/year")

    print(f"\n[Annual p95 Hs - storm intensity]")
    result_p95 = mann_kendall_hamed_rao(t_years, annual_p95.values, args.alpha)
    if "error" in result_p95:
        print(f"  {result_p95['error']}")
    else:
        print(f"  n={result_p95['n']} years, S={result_p95['S']}, "
              f"Hamed-Rao correction factor={result_p95['hamed_rao_correction_factor']:.2f}")
        print(f"  Z={result_p95['z']:.3f}, p={result_p95['p_value']:.4f} "
              f"-> {result_p95['trend']}")
        print(f"  Sen's slope: {result_p95['sens_slope_per_year']:+.4f} m/year")

    out_dir = default_paths("14_mann_kendall_trend")
    with open(out_dir / f"{args.buoy}_{args.var}_mann_kendall.json", "w") as f:
        json.dump({
            "record_years": record_years,
            "sparse_years": sparse_years,
            "annual_mean_trend": result_mean,
            "annual_p95_trend": result_p95,
        }, f, indent=2)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    for ax, series, result, label in [
        (axes[0], annual_mean, result_mean, "Annual mean Hs"),
        (axes[1], annual_p95, result_p95, "Annual p95 Hs (storm intensity)"),
    ]:
        ax.plot(series.index, series.values, marker="o", ms=4, lw=1, color="steelblue")
        if "error" not in result:
            slope = result["sens_slope_per_year"]
            intercept = np.median(series.values) - slope * np.median(t_years)
            trend_line = intercept + slope * t_years
            color = "firebrick" if result["significant_at_alpha"] else "gray"
            ls = "-" if result["significant_at_alpha"] else "--"
            ax.plot(t_years, trend_line, color=color, ls=ls,
                     label=f"Sen's slope: {slope:+.4f} m/yr (p={result['p_value']:.3f})")
            ax.legend()
        ax.set_ylabel(label)
        ax.set_xlabel("Year")
        ax.grid(alpha=0.3)
    fig.suptitle(f"{args.buoy} — {args.var} Mann-Kendall trend (Hamed-Rao corrected)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_trend.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
