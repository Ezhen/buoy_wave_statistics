"""
Stability analysis: do the pipeline's conclusions hold up under
perturbation, or are they fragile to which sub-window or which single
storm happened to be in the record?

Three checks:
  A. Moving-window stability - recompute mean Hs and best-fit distribution
     across non-overlapping sub-windows. Does the "best distribution"
     verdict from Stage 06 hold in every window, or only on average?
  B. Jackknife / drop-biggest-storm sensitivity - refit the Weibull
     distribution (Stage 06) and the GPD shape parameter (Stage 08) with
     the single largest storm removed. Compares the shift against Stage
     12's bootstrap CI width where available: a drop-one shift comparable
     to the CI width is the SAME instability signal Stage 12 already
     found, not a new independent scare.
  C. Regime fraction stability - block bootstrap the Stage 10 regime
     label SEQUENCE itself (not re-fitting the GMM each time - that
     would be its own, more expensive question). This quantifies
     sampling uncertainty in the time-fraction estimates given the
     existing regime assignment, not uncertainty in the assignment
     itself; stated explicitly since it's a real simplification.

Usage:
    python 13_stability_analysis.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from utils import default_paths, resolve_block_length

CANDIDATES = {
    "rayleigh": stats.rayleigh,
    "weibull": stats.weibull_min,
    "lognormal": stats.lognorm,
}


def best_fit(data):
    best_name, best_ks = None, np.inf
    for name, dist in CANDIDATES.items():
        try:
            params = dist.fit(data)
            ks_stat, _ = stats.kstest(data, dist.name, args=params)
            if ks_stat < best_ks:
                best_name, best_ks = name, ks_stat
        except Exception:
            continue
    return best_name, best_ks


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
    parser.add_argument("--n-windows", default=4, type=int)
    parser.add_argument("--n-bootstrap", default=1000, type=int)
    parser.add_argument("--random-state", default=0, type=int)
    args = parser.parse_args()

    rng = np.random.default_rng(args.random_state)
    out_dir = default_paths("13_stability_analysis")

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    hs = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var].dropna()

    print(f"--- {args.buoy} / {args.var} stability analysis ---")

    # ============ A. Moving-window stability ============
    print(f"\n[A] Moving-window stability ({args.n_windows} windows)")
    n = len(hs)
    edges = np.linspace(0, n, args.n_windows + 1).astype(int)
    window_rows = []
    for i in range(args.n_windows):
        window = hs.values[edges[i]:edges[i + 1]]
        if len(window) < 30:
            print(f"  window {i+1}: only {len(window)} samples - skipped, too few to fit")
            continue
        name, ks = best_fit(window[window > 0])
        window_rows.append({
            "window": i + 1, "n": len(window),
            "start": str(hs.index[edges[i]]), "end": str(hs.index[edges[i + 1] - 1]),
            "mean_hs": float(window.mean()), "best_distribution": name, "ks_stat": ks,
        })
        print(f"  window {i+1}: n={len(window):4d}  mean={window.mean():.3f}  "
              f"best={name} (KS={ks:.4f})")

    window_df = pd.DataFrame(window_rows)
    window_df.to_csv(out_dir / f"{args.buoy}_{args.var}_window_stability.csv", index=False)

    overall_best = best_fit(hs.values[hs.values > 0])[0]
    dist_agreement = (window_df["best_distribution"] == overall_best).mean() if len(window_df) else np.nan
    mean_range = window_df["mean_hs"].max() - window_df["mean_hs"].min() if len(window_df) else np.nan
    if len(window_df):
        print(f"\n  Overall best-fit (Stage 06 style, full record): {overall_best}")
        print(f"  {dist_agreement*100:.0f}% of windows agree with the overall best-fit distribution")
        print(f"  Mean Hs range across windows: {mean_range:.3f} m")
        if dist_agreement < 0.75:
            print(f"  NOTE: best-fit distribution is NOT stable across the record - "
                  f"the Stage 06 verdict may depend on which part of the record dominates, "
                  f"not a property that holds throughout.")

    # ============ B. Jackknife / drop-biggest-storm ============
    print(f"\n[B] Drop-biggest-storm sensitivity")
    hs_pos = hs.values[hs.values > 0]
    max_idx = np.argmax(hs_pos)
    hs_dropped = np.delete(hs_pos, max_idx)

    shape_full, _, scale_full = stats.weibull_min.fit(hs_pos)
    shape_drop, _, scale_drop = stats.weibull_min.fit(hs_dropped)
    shape_pct_change = 100 * (shape_drop - shape_full) / shape_full
    scale_pct_change = 100 * (scale_drop - scale_full) / scale_full
    print(f"  Weibull shape: {shape_full:.4f} -> {shape_drop:.4f} ({shape_pct_change:+.1f}%)")
    print(f"  Weibull scale: {scale_full:.4f} -> {scale_drop:.4f} ({scale_pct_change:+.1f}%)")
    if abs(shape_pct_change) > 10 or abs(scale_pct_change) > 10:
        print(f"  NOTE: >10% shift from dropping a single storm - the Weibull fit is "
              f"somewhat driven by one event, not purely the bulk distribution.")

    peaks_path = Path("pipeline_out/08_extreme_value_analysis") / f"{args.buoy}_{args.var}_declustered_peaks.csv"
    eva_path = Path("pipeline_out/08_extreme_value_analysis") / f"{args.buoy}_{args.var}_eva_summary.json"
    xi_jackknife = None
    if peaks_path.exists() and eva_path.exists():
        peaks = pd.read_csv(peaks_path, index_col=0)["peak"].values
        with open(eva_path) as f:
            eva = json.load(f)
        threshold = eva["threshold"]
        excess = peaks - threshold
        if len(excess) >= 5:
            excess_dropped = np.delete(excess, np.argmax(excess))
            try:
                xi_full = eva["gpd_shape"]
                xi_drop, _, _ = stats.genpareto.fit(excess_dropped, floc=0)
                xi_shift = xi_drop - xi_full
                print(f"  GPD xi: {xi_full:.4f} -> {xi_drop:.4f} (shift {xi_shift:+.4f}, "
                      f"from {len(excess)} to {len(excess_dropped)} peaks)")

                ci_path = Path("pipeline_out/12_confidence_intervals") / f"{args.buoy}_{args.var}_confidence_summary.json"
                if ci_path.exists():
                    with open(ci_path) as f:
                        ci_summary = json.load(f)
                    xi_ci = ci_summary.get("gpd_xi_ci")
                    if xi_ci:
                        ci_width = xi_ci["ci_high"] - xi_ci["ci_low"]
                        relation = ("comparable to" if abs(xi_shift) > 0.3 * ci_width
                                    else "small relative to")
                        print(f"  For context, Stage 12's bootstrap CI width on xi is "
                              f"{ci_width:.3f} - {relation} this drop-one shift.")
                        if abs(xi_shift) > 0.3 * ci_width:
                            print(f"  This is the SAME instability Stage 12 already flagged, "
                                  f"not a new independent issue.")
                xi_jackknife = {"xi_full": xi_full, "xi_dropped": float(xi_drop), "shift": float(xi_shift),
                                 "n_peaks_full": len(excess)}
            except Exception as e:
                print(f"  GPD refit on dropped-peak set failed: {e}")
        else:
            print(f"  Only {len(excess)} storm peaks - skipping GPD jackknife, too few to drop one meaningfully.")
    else:
        print("  No Stage 08 output found - skipping GPD jackknife.")

    # ============ C. Regime fraction bootstrap ============
    print("\n[C] Regime fraction stability")
    regime_path = Path("pipeline_out/10_regime_identification") / f"{args.buoy}_{args.var}_regime_labels.csv"
    regime_fraction_ci = None
    if regime_path.exists():
        regimes = pd.read_csv(regime_path, index_col=0)["regime"].dropna().values.astype(int)
        n_regimes = regimes.max() + 1
        block_length, block_source, hit_ceiling = resolve_block_length(args.buoy, args.var, len(regimes))
        print(f"  Block length: {block_length} samples (source: {block_source})")

        boot_fractions = np.zeros((args.n_bootstrap, n_regimes))
        for i in range(args.n_bootstrap):
            resample = moving_block_bootstrap_resample(regimes, block_length, rng)
            for r in range(n_regimes):
                boot_fractions[i, r] = np.mean(resample == r)

        point_fractions = np.array([np.mean(regimes == r) for r in range(n_regimes)])
        regime_fraction_ci = []
        for r in range(n_regimes):
            lo, hi = np.percentile(boot_fractions[:, r], [2.5, 97.5])
            regime_fraction_ci.append({"regime_id": r, "point_fraction": float(point_fractions[r]),
                                        "ci_low": float(lo), "ci_high": float(hi)})
            print(f"  regime {r}: {point_fractions[r]*100:.1f}% of time, "
                  f"95% CI [{lo*100:.1f}%, {hi*100:.1f}%]")

        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(n_regimes)
        err = [(point_fractions[r] - regime_fraction_ci[r]["ci_low"]) * 100 for r in range(n_regimes)]
        ax.bar(x, point_fractions * 100, yerr=err, capsize=5, color="steelblue", alpha=0.7)
        ax.set_xticks(x)
        ax.set_xlabel("regime id (0=calmest)")
        ax.set_ylabel("% of time")
        ax.set_title(f"{args.buoy} — regime time-fraction with block-bootstrap CI")
        fig.tight_layout()
        fig.savefig(out_dir / f"{args.buoy}_{args.var}_regime_fraction_ci.png", dpi=150)
        plt.close(fig)
    else:
        print("  No Stage 10 regime labels found - skipping. Run 10_regime_identification.py first.")

    # ============ Save summary ============
    summary = {
        "window_stability": {
            "overall_best_distribution": overall_best,
            "fraction_windows_agreeing": None if pd.isna(dist_agreement) else float(dist_agreement),
            "mean_hs_range_across_windows": None if pd.isna(mean_range) else float(mean_range),
        },
        "weibull_jackknife": {
            "shape_full": float(shape_full), "shape_dropped": float(shape_drop),
            "shape_pct_change": float(shape_pct_change),
            "scale_full": float(scale_full), "scale_dropped": float(scale_drop),
            "scale_pct_change": float(scale_pct_change),
        },
        "gpd_xi_jackknife": xi_jackknife,
        "regime_fraction_ci": regime_fraction_ci,
    }
    with open(out_dir / f"{args.buoy}_{args.var}_stability_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
