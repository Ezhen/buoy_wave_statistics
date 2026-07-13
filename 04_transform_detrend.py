"""
Stage 3 - Detrend (differencing), on top of the variance-stabilized,
tidally-detided series from Stage 2b.

If pipeline_out/03b_tidal_notch/<buoy>_<var>_detided_boxcox.csv exists, this
script differences THAT (Box-Cox already applied, M2 tide already removed).
Otherwise it falls back to doing its own Box-Cox on the raw level series with
a warning - only appropriate if Stage 02's periodogram showed no meaningful
tidal peak for this buoy/variable.

Usage:
    python 04_transform_detrend.py --buoy WesthinderBuoy --var VHM0 --diff-order 1
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import boxcox
from statsmodels.tsa.stattools import adfuller

from utils import default_paths, longest_contiguous_segment


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--diff-order", default=1, type=int,
                         help="1 or 2; rarely go beyond 2")
    args = parser.parse_args()

    detided_path = Path("pipeline_out/03b_tidal_notch") / f"{args.buoy}_{args.var}_detided_boxcox.csv"
    raw_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"

    if detided_path.exists():
        s_transformed_full = pd.read_csv(detided_path, index_col=0, parse_dates=True).iloc[:, 0]
        lam = None
        print(f"Using detided series from Stage 2b: {detided_path}")
        s_dropna = pd.read_csv(raw_path, index_col=0, parse_dates=True)[args.var].dropna()  # for plotting only
    else:
        print("WARNING: no Stage 2b output found - Box-Cox transforming the raw series "
              "directly. Only fine if Stage 02's periodogram showed no significant "
              "tidal peak; otherwise run 03b_tidal_notch.py first.")
        s = pd.read_csv(raw_path, index_col=0, parse_dates=True)[args.var]
        s_dropna = s.dropna()
        shifted = s_dropna + 1e-6 if (s_dropna <= 0).any() else s_dropna
        transformed_vals, lam = boxcox(shifted.values)
        s_transformed_full = pd.Series(transformed_vals, index=s_dropna.index, name=f"{args.var}_boxcox")
        print(f"Box-Cox lambda = {lam:.4f}")

    # Differencing and the ADF check that follows are lag-based - restrict
    # to the longest contiguous segment so a splice across a long gap
    # doesn't get treated as temporally adjacent by adfuller()'s own
    # internal lag construction (diff() itself is positionally safe once
    # NaN gaps are preserved - this guards adfuller() specifically).
    s_transformed, seg_meta = longest_contiguous_segment(s_transformed_full)
    if seg_meta["n_segments"] > 1:
        print(f"Record has {seg_meta['n_segments']} contiguous segments - using the "
              f"longest ({seg_meta['pct_of_valid_used']}% of valid samples, "
              f"{seg_meta['segment_start']} to {seg_meta['segment_end']}) for "
              f"differencing/ADF, not the full gap-spliced record.")

    # --- Detrend via differencing on the TRANSFORMED (and, if available, detided) series ---
    s_diff = s_transformed.copy()
    for order in range(args.diff_order):
        s_diff = s_diff.diff().dropna()

    adf_after = adfuller(s_diff)
    print(f"ADF after Box-Cox + order-{args.diff_order} differencing: "
          f"stat={adf_after[0]:.3f}, p={adf_after[1]:.4f} "
          f"-> {'stationary' if adf_after[1] < 0.05 else 'still not stationary'}")

    out_dir = default_paths("04_transform_detrend")
    s_transformed.to_csv(out_dir / f"{args.buoy}_{args.var}_boxcox.csv", header=[s_transformed.name])
    s_diff.to_csv(out_dir / f"{args.buoy}_{args.var}_residual.csv", header=["residual"])

    import json
    with open(out_dir / f"{args.buoy}_{args.var}_detrend_summary.json", "w") as f:
        json.dump({
            "used_detided_input": bool(detided_path.exists()),
            "diff_order": args.diff_order,
            "adf_stat_after": float(adf_after[0]),
            "adf_pvalue_after": float(adf_after[1]),
            "n_gap_segments": seg_meta["n_segments"],
            "longest_segment_pct_of_valid": seg_meta["pct_of_valid_used"],
        }, f, indent=2)

    lam_label = f"lambda={lam:.3f}" if lam is not None else "from Stage 2b (detided)"
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=False)
    axes[0].plot(s_dropna.index, s_dropna.values, lw=0.6)
    axes[0].set_title(f"{args.var} — raw (level)")
    axes[1].plot(s_transformed.index, s_transformed.values, lw=0.6, color="darkorange")
    axes[1].set_title(f"{args.var} — after Box-Cox ({lam_label})")
    axes[2].plot(s_diff.index, s_diff.values, lw=0.6, color="firebrick")
    axes[2].set_title(f"{args.var} — after Box-Cox + order-{args.diff_order} differencing (residual)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_transform_stages.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved residual + plots to {out_dir}")


if __name__ == "__main__":
    main()
