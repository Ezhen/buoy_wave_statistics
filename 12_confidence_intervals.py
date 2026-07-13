"""
Statistical confidence: report uncertainty, not just point estimates.

Every series in this pipeline is confirmed autocorrelated (Ljung-Box
fails at all 19 buoys), so an ordinary/IID bootstrap here would
understate uncertainty. This uses:

  - MOVING BLOCK BOOTSTRAP for Hs mean/quantile CIs, with block length
    taken from Stage 11b's per-buoy integral timescale (falls back to a
    sqrt(n) heuristic with a loud warning if 11b hasn't been run - that
    fallback is a much weaker justification and should be treated as
    such)
  - ORDINARY bootstrap for the GPD shape parameter (xi), resampling
    Stage 08's DECLUSTERED storm peaks directly - ordinary bootstrap is
    reasonable here specifically because declustering already targets
    quasi-independence between events, unlike raw Hs
  - The same block-bootstrap resamples used for Hs quantiles are reused
    to build a CI band on the fitted Weibull PDF from Stage 06

Usage:
    python 12_confidence_intervals.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.special import gamma as gamma_fn
from scipy.optimize import brentq

from utils import default_paths


def fast_weibull_moment_fit(x: np.ndarray):
    """Method-of-moments Weibull fit (loc fixed at 0) via the
    coefficient-of-variation relation, instead of full MLE.

    Why: scipy.stats.weibull_min.fit() runs a numerical optimizer that
    evaluates the log-likelihood over the WHOLE array at every iteration
    - at ~2s per fit on a 500k-sample array, 1000 bootstrap refits costs
    30-40 minutes, scaling worse on multi-year buoys with 1M+ samples.
    This needs only mean/std (one O(n) pass, numpy-vectorized) plus a
    scalar root-find independent of n - milliseconds instead of seconds.

    Simplification stated explicitly: fixes loc=0 (reasonable for a
    physically-bounded-below quantity like Hs) and doesn't estimate loc
    per resample the way full MLE would. Fine for building a CI *band*
    (visual spread), not intended to replace Stage 06's actual point
    estimate, which still uses full MLE.
    """
    x = x[x > 0]
    mean, std = x.mean(), x.std()
    if mean <= 0 or std <= 0:
        raise ValueError("degenerate resample (zero mean or variance)")
    cv = std / mean

    def cv_gap(k):
        return np.sqrt(gamma_fn(1 + 2 / k) / gamma_fn(1 + 1 / k) ** 2 - 1) - cv

    shape = brentq(cv_gap, 0.05, 50.0)
    scale = mean / gamma_fn(1 + 1 / shape)
    return shape, 0.0, scale


def moving_block_bootstrap_resample(x: np.ndarray, block_length: int, rng: np.random.Generator):
    n = len(x)
    block_length = max(1, min(block_length, n))
    n_blocks = int(np.ceil(n / block_length))
    starts = rng.integers(0, n - block_length + 1, size=n_blocks)
    out = np.concatenate([x[s:s + block_length] for s in starts])[:n]
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--n-bootstrap", default=None, type=int,
                         help="default auto-scales down for large records (bootstrap CI "
                              "precision barely improves past a few hundred resamples "
                              "once the underlying sample size itself is huge) - pass "
                              "explicitly to override")
    parser.add_argument("--quantiles", default="50,90,95",
                         help="comma-separated percentiles of Hs to CI, besides the mean")
    parser.add_argument("--random-state", default=0, type=int)
    args = parser.parse_args()

    rng = np.random.default_rng(args.random_state)
    out_dir = default_paths("12_confidence_intervals")

    # ---------- Hs mean/quantile block bootstrap ----------
    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    hs = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var].dropna().values

    if args.n_bootstrap is None:
        if len(hs) >= 500_000:
            args.n_bootstrap = 200
        elif len(hs) >= 50_000:
            args.n_bootstrap = 500
        else:
            args.n_bootstrap = 1000
        print(f"Auto-scaled --n-bootstrap to {args.n_bootstrap} for n={len(hs)} samples "
              f"(pass --n-bootstrap explicitly to override)")

    dep_path = Path("pipeline_out/11b_dependence_structure") / f"{args.buoy}_{args.var}_dependence_summary.json"
    if dep_path.exists():
        with open(dep_path) as f:
            dep = json.load(f)
        block_length = dep["suggested_block_length_samples"]
        block_source = "Stage 11b integral timescale"
        if dep.get("hit_max_lag_ceiling"):
            print("WARNING: Stage 11b's block length is itself a lower bound (hit its "
                  "search ceiling) - treat these CIs as possibly too narrow.")
    else:
        block_length = max(1, int(np.sqrt(len(hs))))
        block_source = "sqrt(n) fallback - Stage 11b not found, this is a WEAK default"
        print(f"WARNING: no Stage 11b output for {args.buoy} - using a sqrt(n) block "
              f"length ({block_length} samples) instead of a persistence-based one. "
              f"Run 11b_dependence_structure.py first for a real justification.")

    print(f"--- {args.buoy} / {args.var} confidence intervals ---")
    print(f"Block length: {block_length} samples (source: {block_source})")
    print(f"Bootstrap resamples: {args.n_bootstrap}")

    quantile_list = [float(q) for q in args.quantiles.split(",")]
    boot_means = np.empty(args.n_bootstrap)
    boot_quantiles = {q: np.empty(args.n_bootstrap) for q in quantile_list}
    boot_weibull_curves = []

    x_grid = np.linspace(max(0.01, hs.min()), hs.max() * 1.1, 200)
    weibull_fit_ok = 0

    for i in range(args.n_bootstrap):
        resample = moving_block_bootstrap_resample(hs, block_length, rng)
        boot_means[i] = resample.mean()
        for q in quantile_list:
            boot_quantiles[q][i] = np.percentile(resample, q)
        try:
            shape, loc, scale = fast_weibull_moment_fit(resample)
            boot_weibull_curves.append(stats.weibull_min.pdf(x_grid, shape, loc, scale))
            weibull_fit_ok += 1
        except Exception:
            pass

    def ci(arr, alpha=0.05):
        return float(np.percentile(arr, 100 * alpha / 2)), float(np.percentile(arr, 100 * (1 - alpha / 2)))

    point_mean = float(hs.mean())
    mean_ci = ci(boot_means)
    print(f"\nHs mean: {point_mean:.3f} m, 95% block-bootstrap CI [{mean_ci[0]:.3f}, {mean_ci[1]:.3f}]")

    quantile_results = {}
    for q in quantile_list:
        point_q = float(np.percentile(hs, q))
        q_ci = ci(boot_quantiles[q])
        quantile_results[q] = {"point": point_q, "ci_low": q_ci[0], "ci_high": q_ci[1]}
        print(f"Hs p{q:g}: {point_q:.3f} m, 95% CI [{q_ci[0]:.3f}, {q_ci[1]:.3f}]")

    # ---------- Weibull PDF CI band ----------
    weibull_band = None
    if weibull_fit_ok > args.n_bootstrap * 0.5:
        curves = np.array(boot_weibull_curves)
        band_low = np.percentile(curves, 2.5, axis=0)
        band_high = np.percentile(curves, 97.5, axis=0)
        band_mid = np.percentile(curves, 50, axis=0)
        weibull_band = {"x": x_grid, "low": band_low, "high": band_high, "mid": band_mid}
        print(f"\nWeibull PDF CI band computed from {weibull_fit_ok}/{args.n_bootstrap} "
              f"successful bootstrap refits.")
    else:
        print(f"\nWeibull refit succeeded on only {weibull_fit_ok}/{args.n_bootstrap} "
              f"bootstrap resamples - skipping the PDF CI band, too unstable to show.")

    # ---------- GPD xi bootstrap (ordinary, on declustered peaks) ----------
    peaks_path = Path("pipeline_out/08_extreme_value_analysis") / f"{args.buoy}_{args.var}_declustered_peaks.csv"
    eva_summary_path = Path("pipeline_out/08_extreme_value_analysis") / f"{args.buoy}_{args.var}_eva_summary.json"
    xi_ci = None
    if peaks_path.exists() and eva_summary_path.exists():
        peaks = pd.read_csv(peaks_path, index_col=0)["peak"].values
        with open(eva_summary_path) as f:
            eva = json.load(f)
        threshold = eva["threshold"]
        excess = peaks - threshold

        if len(excess) < 5:
            print(f"\nOnly {len(excess)} storm peaks - skipping GPD xi bootstrap, "
                  f"not enough data to resample meaningfully.")
        else:
            boot_xi = []
            n_fail = 0
            for _ in range(args.n_bootstrap):
                resample = rng.choice(excess, size=len(excess), replace=True)
                try:
                    shape, _, _ = stats.genpareto.fit(resample, floc=0)
                    boot_xi.append(shape)
                except Exception:
                    n_fail += 1
            boot_xi = np.array(boot_xi)
            if len(boot_xi) > args.n_bootstrap * 0.5:
                xi_point = eva["gpd_shape"]
                xi_ci_bounds = ci(boot_xi)
                xi_ci = {"point": xi_point, "ci_low": xi_ci_bounds[0], "ci_high": xi_ci_bounds[1],
                         "n_peaks": len(excess), "n_bootstrap_failed": n_fail}
                print(f"\nGPD xi: {xi_point:.3f}, 95% bootstrap CI "
                      f"[{xi_ci_bounds[0]:.3f}, {xi_ci_bounds[1]:.3f}] "
                      f"(from {len(excess)} storm peaks, {n_fail} bootstrap refits failed)")
                if len(excess) < 10:
                    print("NOTE: fewer than 10 peaks feeding this bootstrap - the CI "
                          "width itself is not very trustworthy at this sample size, "
                          "even though the bootstrap machinery ran fine.")
            else:
                print(f"\nGPD xi bootstrap: {n_fail}/{args.n_bootstrap} refits failed - "
                      f"too unstable to report a CI.")
    else:
        if peaks_path.exists() and not eva_summary_path.exists():
            print(f"\nStage 08's declustered peaks exist for {args.buoy} but its GPD fit "
                  f"didn't succeed (too few peaks) - skipping GPD xi CI.")
        else:
            print(f"\nNo Stage 08 output found for {args.buoy} - skipping GPD xi CI. "
                  f"Run 08_extreme_value_analysis.py first.")

    # ---------- Save ----------
    summary = {
        "block_length_samples": block_length,
        "block_length_source": block_source,
        "n_bootstrap": args.n_bootstrap,
        "hs_mean_point": point_mean,
        "hs_mean_ci_low": mean_ci[0],
        "hs_mean_ci_high": mean_ci[1],
        "hs_quantiles": quantile_results,
        "gpd_xi_ci": xi_ci,
    }
    with open(out_dir / f"{args.buoy}_{args.var}_confidence_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ---------- Plots ----------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(boot_means, bins=40, alpha=0.6, color="steelblue")
    ax.axvline(point_mean, color="black", lw=2, label="point estimate")
    ax.axvline(mean_ci[0], color="firebrick", ls="--", label="95% CI")
    ax.axvline(mean_ci[1], color="firebrick", ls="--")
    ax.set_xlabel("bootstrap Hs mean (m)")
    ax.set_title(f"{args.buoy} — block bootstrap distribution of Hs mean\n"
                 f"(block length {block_length} samples, {block_source})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_mean_bootstrap.png", dpi=150)
    plt.close(fig)

    if weibull_band is not None:
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.fill_between(weibull_band["x"], weibull_band["low"], weibull_band["high"],
                         alpha=0.3, color="steelblue", label="95% CI band")
        ax.plot(weibull_band["x"], weibull_band["mid"], color="steelblue", lw=2, label="median refit")
        ax.hist(hs, bins=40, density=True, alpha=0.3, color="gray", label="observed")
        ax.set_xlabel(args.var)
        ax.set_ylabel("density")
        ax.set_title(f"{args.buoy} — Weibull PDF with block-bootstrap CI band")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / f"{args.buoy}_{args.var}_weibull_ci_band.png", dpi=150)
        plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
