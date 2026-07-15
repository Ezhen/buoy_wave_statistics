"""
Stage A - persistence baseline: Hs(t+h) = Hs(t), at a defined set of
operationally meaningful horizons. This is the number every later
forecasting model (ARMA, ARMA-GARCH, ARMAX) has to beat - computed
first, before anything fancier, using the shared harness in
forecast_utils.py that those models will reuse.

IMPORTANT: loads Stage 0's FULL REGULARIZED grid (gaps preserved as
NaN), NOT the dropna'd series - the backtest harness converts horizons
between hours and samples assuming uniform sample spacing
(horizon_samples = horizon_hours / dt_hours). Using a dropna'd,
gap-collapsed series here would silently make every horizon wrong,
the same category of mistake fixed across 02/03/03b/04/11b earlier -
this stage doesn't need lag-based gap segmentation the way those did
(a persistence forecast only needs one recent valid point, not a whole
contiguous stretch), but it does need the real time axis intact.

This stage only establishes the baseline error curve - a skill score
isn't computable yet (nothing to compare against until a real model
exists in a later stage).

Usage:
    python 17_forecast_baseline.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import default_paths
from forecast_utils import persistence_forecast, rolling_origin_backtest, summarize_backtest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--horizons-hours", default="1,3,6,12,24",
                         help="comma-separated forecast horizons in hours")
    parser.add_argument("--origin-step-hours", default=6.0, type=float,
                         help="spacing between backtest origins, in hours - "
                              "doesn't need to be small for validity, just "
                              "reduces redundant computation")
    parser.add_argument("--min-history-hours", default=24.0, type=float)
    parser.add_argument("--max-lookback-hours", default=3.0, type=float,
                         help="refuse to persist from a value staler than this")
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"

    # Full regularized grid, gaps preserved as NaN - NOT dropna'd (see module docstring)
    hs_full = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var]

    dt_hours = 0.5
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            dt_hours = json.load(f).get("sampling_interval_hours", 0.5)

    horizons_hours = [float(h) for h in args.horizons_hours.split(",")]
    horizons_samples = [max(1, round(h / dt_hours)) for h in horizons_hours]
    origin_step = max(1, round(args.origin_step_hours / dt_hours))
    min_history = max(1, round(args.min_history_hours / dt_hours))
    max_lookback = max(1, round(args.max_lookback_hours / dt_hours))

    print(f"--- {args.buoy} / {args.var} persistence baseline ---")
    print(f"Sampling interval: {dt_hours}h, n={len(hs_full)} (full grid, gaps as NaN)")
    print(f"Horizons: {horizons_hours} hours -> {horizons_samples} samples")
    print(f"Origin step: {args.origin_step_hours}h ({origin_step} samples), "
          f"min history: {args.min_history_hours}h, max lookback: {args.max_lookback_hours}h")

    results = rolling_origin_backtest(
        hs_full, persistence_forecast, horizons_samples,
        origin_step=origin_step, min_history=min_history, max_lookback_samples=max_lookback,
        history_window=max_lookback + 5)

    if len(results) == 0:
        print("No valid (origin, horizon) pairs produced - check the series has "
              "enough valid data given the gap structure.")
        return

    summary = summarize_backtest(results, dt_hours)
    print(f"\nUsed {results['origin_idx'].nunique()} distinct origins, "
          f"{len(results)} total (origin, horizon) scored pairs.")
    print("\nPer-horizon error (this IS the number every later model has to beat):")
    print(summary[["horizon_hours", "n", "rmse", "mae"]].to_string(index=False))

    out_dir = default_paths("17_forecast_baseline")
    summary.to_csv(out_dir / f"{args.buoy}_{args.var}_persistence_summary.csv", index=False)
    results.to_csv(out_dir / f"{args.buoy}_{args.var}_persistence_results.csv", index=False)

    with open(out_dir / f"{args.buoy}_{args.var}_persistence_meta.json", "w") as f:
        json.dump({
            "horizons_hours": horizons_hours,
            "origin_step_hours": args.origin_step_hours,
            "min_history_hours": args.min_history_hours,
            "max_lookback_hours": args.max_lookback_hours,
            "n_origins": int(results["origin_idx"].nunique()),
            "n_scored_pairs": len(results),
        }, f, indent=2)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(summary["horizon_hours"], summary["rmse"], marker="o", label="RMSE")
    ax.plot(summary["horizon_hours"], summary["mae"], marker="s", label="MAE")
    ax.set_xlabel("forecast horizon (hours)")
    ax.set_ylabel(f"{args.var} error (m)")
    ax.set_title(f"{args.buoy} — persistence baseline error vs. horizon")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_error_vs_horizon.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
