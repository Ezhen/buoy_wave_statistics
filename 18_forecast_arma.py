"""
Stage C - ARMA forecasting: does modeling autocorrelation structure
beyond "last observed value" actually help, or was persistence already
capturing most of the achievable short-horizon skill? Forecasts RAW Hs
directly (same target as Stage A's persistence baseline), not the
detided/Box-Cox series - keeps this a fair, direct comparison without
needing to invert transforms and reconstruct the tidal component
(which Stage 03b doesn't currently save the coefficients for). Since
the M2 tide is itself deterministic and predictable at 30-min sampling,
an ARMA model with enough AR terms picking up on it is a legitimate
source of short-horizon skill here, not something to filter out first.

Order selection: small AIC grid search over a modest range, on a
representative training chunk (not the whole multi-year record - full
MLE fitting is expensive, and order shouldn't need the full record to
identify). Differencing order d defaults to 0 ONLY (not 0 and 1) - not
assumed from Stage 04's d=1, which was tuned for the detided series;
originally searched both, but a real finding on A2Buoy showed AIC-based
d selection can actively hurt forecast skill even though it improves
in-sample fit (see --d-range help for the specific numbers). d=1 can
still be requested explicitly via --d-range if there's a specific
reason to.

Performance: fit ONCE, not at every backtest origin - see
forecast_utils.make_arma_forecast_fn for why (repeating the persistence
baseline's O(n^2) mistake at MLE-refit cost would be far worse). Each
origin only updates the model's state on a BOUNDED recent window via
.apply(refit=False), reusing the fixed coefficients from the one fit.

Runs on the longest contiguous segment only (ARMA fitting needs
genuinely contiguous data, same reasoning as the lag-based
characterization stages earlier in this pipeline) - reports what
fraction of the record that represents.

Usage:
    python 18_forecast_arma.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.arima.model import ARIMA

from utils import default_paths, longest_contiguous_segment
from forecast_utils import (rolling_origin_backtest, summarize_backtest,
                             skill_score, make_arma_forecast_fn, persistence_forecast,
                             select_order)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--horizons-hours", default="1,3,6,12,24")
    parser.add_argument("--origin-step-hours", default=6.0, type=float)
    parser.add_argument("--min-history-hours", default=24.0, type=float)
    parser.add_argument("--order-search-samples", default=4000, type=int,
                         help="how many samples of the training data to use for "
                              "the (cheap) AIC order search - doesn't need the "
                              "full record")
    parser.add_argument("--fit-samples", default=8000, type=int,
                         help="how many samples to fit the final model on, once "
                              "order is chosen")
    parser.add_argument("--state-window-hours", default=200.0, type=float,
                         help="bounded recent-history window used to update model "
                              "state at each backtest origin - large enough to "
                              "cover typical ARMA memory without O(n) cost per origin")
    parser.add_argument("--d-range", default="0", help="differencing orders to "
                         "search over. DEFAULT CHANGED TO d=0 ONLY after a real "
                         "finding on A2Buoy: AIC-based order selection picked d=1, "
                         "which gave NEGATIVE skill vs. persistence at every "
                         "horizon beyond 1h (-0.007 to -0.019) - a d=1 model "
                         "structurally discards mean-reversion and degenerates "
                         "toward persistence at long horizons, which AIC (an "
                         "in-sample fit measure) has no way to penalize for a "
                         "forecasting use case. Forcing d=0 on the same buoy "
                         "flipped every horizon positive (0.074 to 0.216) and "
                         "matched the pattern already seen at 2 other buoys. "
                         "Pass --d-range 0,1 explicitly to let AIC choose again "
                         "if you have a specific reason to.")
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"

    hs_full = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var]
    dt_hours = 0.5
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            dt_hours = json.load(f).get("sampling_interval_hours", 0.5)

    print(f"--- {args.buoy} / {args.var} ARMA forecasting ---")

    segment, seg_meta = longest_contiguous_segment(hs_full)
    print(f"Longest contiguous segment: {len(segment)} samples "
          f"({seg_meta['pct_of_valid_used']}% of valid data), "
          f"{seg_meta['segment_start']} to {seg_meta['segment_end']}")
    if seg_meta["n_segments"] > 1:
        print(f"(record has {seg_meta['n_segments']} segments total - "
              f"ARMA fit and backtest both restricted to the longest one)")

    if len(segment) < args.fit_samples + 1000:
        print(f"Longest segment too short for a meaningful fit+backtest split "
              f"(need > {args.fit_samples + 1000}, have {len(segment)}) - stopping.")
        return

    # --- Order selection on a modest training chunk ---
    order_train = segment.iloc[:args.order_search_samples]
    d_range = [int(d) for d in args.d_range.split(",")]
    print(f"\nSearching ARMA order on first {len(order_train)} samples "
          f"(p,q in 0..3, d in {d_range})...")
    order, aic = select_order(order_train, p_range=range(0, 4), q_range=range(0, 4), d_range=d_range)
    if order is None:
        print("Order search failed to converge on any candidate - stopping.")
        return
    print(f"Selected order: {order} (AIC={aic:.1f})")

    # --- Fit once on a larger chunk, fixed params from here on ---
    fit_train = segment.iloc[:args.fit_samples]
    print(f"Fitting final model on first {len(fit_train)} samples...")
    fitted = ARIMA(fit_train.values, order=order).fit()
    print(f"Fit AIC={fitted.aic:.1f}")

    max_horizon_samples = max(round(float(h) / dt_hours) for h in args.horizons_hours.split(","))
    freq_str = f"{int(round(dt_hours * 60))}min"
    arma_fn = make_arma_forecast_fn(fitted, max_horizon_samples, freq_str=freq_str)

    horizons_hours = [float(h) for h in args.horizons_hours.split(",")]
    horizons_samples = [max(1, round(h / dt_hours)) for h in horizons_hours]
    origin_step = max(1, round(args.origin_step_hours / dt_hours))
    min_history = max(args.fit_samples, round(args.min_history_hours / dt_hours))
    state_window = max(1, round(args.state_window_hours / dt_hours))

    print(f"\nBacktesting on the segment (origin step {args.origin_step_hours}h, "
          f"state window {args.state_window_hours}h)...")
    results = rolling_origin_backtest(
        segment, arma_fn, horizons_samples,
        origin_step=origin_step, min_history=min_history,
        max_lookback_samples=state_window, history_window=state_window)

    if len(results) == 0:
        print("No valid forecasts produced - stopping.")
        return

    summary = summarize_backtest(results, dt_hours)
    print(f"\nUsed {results['origin_idx'].nunique()} distinct origins, "
          f"{len(results)} scored (origin, horizon) pairs.")
    print("\nARMA per-horizon error:")
    print(summary[["horizon_hours", "n", "rmse", "mae"]].to_string(index=False))

    # --- Fair comparison: local persistence baseline on the SAME segment
    #     and origins ARMA was actually tested on, not Stage A's full-record
    #     baseline. Found necessary after the first real run: ARMA's segment
    #     covered only 6.9% of Westhinder's record (36,170 samples, one
    #     ~2-year window) while Stage A's baseline spans the full 36 years -
    #     comparing across mismatched coverage isn't a valid skill score,
    #     since a segment could be atypically easy or hard to forecast
    #     relative to the rest of the record. ---
    print("\nComputing LOCAL persistence baseline on the same segment/origins "
          "(for a fair comparison - NOT reusing Stage A's full-record numbers)...")
    local_baseline_results = rolling_origin_backtest(
        segment, persistence_forecast, horizons_samples,
        origin_step=origin_step, min_history=min_history,
        max_lookback_samples=state_window)
    local_baseline_summary = summarize_backtest(local_baseline_results, dt_hours)

    out_dir = default_paths("18_forecast_arma")
    if len(local_baseline_summary) > 0:
        sk = skill_score(summary, local_baseline_summary)
        print("\nSkill score vs. LOCAL persistence, same segment "
              "(positive = ARMA better, 0 = tied, negative = worse):")
        print(sk.to_string(index=False))
        sk.to_csv(out_dir / f"{args.buoy}_{args.var}_skill_vs_local_persistence.csv", index=False)
        local_baseline_summary.to_csv(
            out_dir / f"{args.buoy}_{args.var}_local_persistence_summary.csv", index=False)
    else:
        print("Local persistence baseline produced no results - skipping skill score.")
        sk = None

    # Also show the full-record baseline for context, explicitly labeled as
    # NOT a fair comparison (different coverage) - useful to see the gap
    # between "this segment" and "the whole record", not to compute skill from.
    baseline_path = Path("pipeline_out/17_forecast_baseline") / f"{args.buoy}_{args.var}_persistence_summary.csv"
    if baseline_path.exists():
        full_record_baseline = pd.read_csv(baseline_path)
        print(f"\n(For reference only - NOT a fair comparison, different coverage: "
              f"Stage A's FULL-RECORD persistence RMSE at each horizon)")
        print(full_record_baseline[["horizon_hours", "rmse"]].to_string(index=False))

    summary.to_csv(out_dir / f"{args.buoy}_{args.var}_arma_summary.csv", index=False)
    with open(out_dir / f"{args.buoy}_{args.var}_arma_meta.json", "w") as f:
        json.dump({
            "order": list(order), "aic": float(aic),
            "fit_samples": args.fit_samples,
            "segment_pct_of_valid": seg_meta["pct_of_valid_used"],
            "n_origins": int(results["origin_idx"].nunique()),
        }, f, indent=2)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(summary["horizon_hours"], summary["rmse"], marker="o", label=f"ARMA{order}")
    if len(local_baseline_summary) > 0:
        ax.plot(local_baseline_summary["horizon_hours"], local_baseline_summary["rmse"],
                marker="s", label="persistence (same segment - fair comparison)", ls="--")
    if baseline_path.exists():
        ax.plot(full_record_baseline["horizon_hours"], full_record_baseline["rmse"],
                marker="^", label="persistence (full record - context only)", ls=":", alpha=0.6)
    ax.set_xlabel("forecast horizon (hours)")
    ax.set_ylabel(f"{args.var} RMSE (m)")
    ax.set_title(f"{args.buoy} — ARMA{order} vs. persistence "
                 f"(segment = {seg_meta['pct_of_valid_used']}% of record)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_arma_vs_persistence.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
