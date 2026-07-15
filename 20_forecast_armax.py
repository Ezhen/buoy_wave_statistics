"""
Stage E - ARMAX: adds lagged ERA5 wind speed as an exogenous regressor
to the ARMA mean model, using the lag Priority 5 already found
(0-3h, wind leads Hs, up to R^2=0.71-0.77 at Westhinder) rather than
guessing one.

CRITICAL DESIGN ISSUE, worth being explicit about rather than getting
wrong silently: for genuine multi-step-ahead forecasting, SARIMAX needs
FUTURE exog values at the forecast horizon - but in a real deployment
we don't know future wind either. Using the true future ERA5 wind in a
backtest would be look-ahead bias (cheating - the model would appear to
work because it's secretly being handed the answer). Handled honestly:
- The lagged exog series exog(t) = wind(t - lag_hours) is fully known
  up to time (origin + lag_samples), since that only needs wind
  observed at or before "now" (origin).
- Beyond that horizon, the model does NOT get real future wind - the
  last genuinely-known wind observation is persisted forward instead
  (the same honest simplification Stage A used for Hs itself). This
  means ARMAX's advantage over plain ARMA should be concentrated at
  short horizons (within the lag window) and should shrink or vanish
  at longer horizons, where it degrades toward using stale wind info -
  a real, testable prediction, not just a caveat.

Compares against BOTH the local persistence baseline (same segment) AND
Stage C's plain-ARMA result (same segment, same order search) - the
second comparison isolates whether wind adds anything BEYOND what
ARMA's own autoregressive structure already captures, which is the
real question here.

Usage:
    python 20_forecast_armax.py --buoy WesthinderBuoy --var VHM0 --exog-lag-hours 3.0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.statespace.sarimax import SARIMAX
import xarray as xr

from utils import default_paths, longest_contiguous_segment, get_scalar_latlon, load_era5_for_buoy
from forecast_utils import (select_order, persistence_forecast, rolling_origin_backtest,
                             summarize_backtest, skill_score)


def select_order_exog(train_endog: np.ndarray, train_exog: np.ndarray, p_range, q_range, d_range):
    best_aic = np.inf
    best_order = None
    for d in d_range:
        for p in p_range:
            for q in q_range:
                if p == 0 and q == 0:
                    continue
                try:
                    fit = SARIMAX(train_endog, exog=train_exog, order=(p, d, q),
                                   enforce_stationarity=False, enforce_invertibility=False).fit(disp=False)
                    if fit.aic < best_aic:
                        best_aic = fit.aic
                        best_order = (p, d, q)
                except Exception:
                    continue
    return best_order, best_aic


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data_multiyear", type=Path)
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--era5-dir", default="meteo_era5")
    parser.add_argument("--exog-lag-hours", default=3.0, type=float,
                         help="lag to apply to wind_speed before using it as exog - "
                              "default matches Priority 5's finding for Westhinder; "
                              "check per-buoy if running elsewhere")
    parser.add_argument("--horizons-hours", default="1,3,6,12,24")
    parser.add_argument("--origin-step-hours", default=6.0, type=float)
    parser.add_argument("--min-history-hours", default=24.0, type=float)
    parser.add_argument("--order-search-samples", default=4000, type=int)
    parser.add_argument("--fit-samples", default=8000, type=int)
    parser.add_argument("--state-window-hours", default=200.0, type=float)
    parser.add_argument("--d-range", default="0")
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"

    hs_full = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var]
    dt_hours = 0.5
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            dt_hours = json.load(f).get("sampling_interval_hours", 0.5)

    print(f"--- {args.buoy} / {args.var} ARMAX forecasting (exog: wind_speed, "
          f"lag={args.exog_lag_hours}h) ---")

    nc_path = args.data_dir / f"{args.buoy}.nc"
    with xr.open_dataset(nc_path) as ds:
        lat, lon = get_scalar_latlon(ds)
    era5 = load_era5_for_buoy(lat, lon, args.era5_dir)
    print(f"ERA5: {len(era5)} 3-hourly samples")

    # Interpolate wind onto Hs's native grid, then apply the lag
    wind_on_hs_grid = era5["wind_speed"].reindex(
        era5["wind_speed"].index.union(hs_full.index)).interpolate(method="time").reindex(hs_full.index)
    lag_samples = max(1, round(args.exog_lag_hours / dt_hours))
    exog_lagged = wind_on_hs_grid.shift(lag_samples)
    print(f"Lag: {args.exog_lag_hours}h = {lag_samples} samples")

    segment, seg_meta = longest_contiguous_segment(hs_full)
    exog_segment = exog_lagged.reindex(segment.index)
    print(f"Longest contiguous segment: {len(segment)} samples "
          f"({seg_meta['pct_of_valid_used']}% of valid data)")
    print(f"Exog coverage in segment: {exog_segment.notna().sum()}/{len(exog_segment)}")

    if len(segment) < args.fit_samples + 1000:
        print(f"Segment too short - stopping.")
        return

    # --- Order selection + fit, same discipline as Stage C ---
    d_range = [int(d) for d in args.d_range.split(",")]
    order_train_y = segment.iloc[:args.order_search_samples].values
    order_train_x = exog_segment.iloc[:args.order_search_samples].values
    valid_mask = ~(np.isnan(order_train_y) | np.isnan(order_train_x))
    print(f"\nSearching ARMAX order on first {valid_mask.sum()} valid samples...")
    order, aic = select_order_exog(order_train_y[valid_mask], order_train_x[valid_mask],
                                     p_range=range(0, 4), q_range=range(0, 4), d_range=d_range)
    if order is None:
        print("Order search failed - stopping.")
        return
    print(f"Selected order: {order} (AIC={aic:.1f})")

    fit_y = segment.iloc[:args.fit_samples].values
    fit_x = exog_segment.iloc[:args.fit_samples].values
    fit_valid = ~(np.isnan(fit_y) | np.isnan(fit_x))
    if fit_valid.mean() < 0.9:
        print(f"WARNING: only {fit_valid.mean()*100:.0f}% of the fit window has valid "
              f"exog - ERA5 coverage gap within the buoy's own longest segment.")
    print(f"Fitting ARMAX on {fit_valid.sum()} valid samples...")
    armax_fitted = SARIMAX(fit_y[fit_valid], exog=fit_x[fit_valid], order=order,
                            enforce_stationarity=False, enforce_invertibility=False).fit(disp=False)

    # --- Backtest: honest exog handling per the module docstring ---
    horizons_hours = [float(h) for h in args.horizons_hours.split(",")]
    horizons_samples = [max(1, round(h / dt_hours)) for h in horizons_hours]
    max_h = max(horizons_samples)
    origin_step = max(1, round(args.origin_step_hours / dt_hours))
    min_history = max(args.fit_samples, round(args.min_history_hours / dt_hours))
    state_window = max(1, round(args.state_window_hours / dt_hours))
    freq_str = f"{int(round(dt_hours * 60))}min"

    print(f"\nBacktesting (origin step {args.origin_step_hours}h, "
          f"state window {args.state_window_hours}h)...")
    print(f"Exog known for h <= {lag_samples} samples ({args.exog_lag_hours}h); "
          f"beyond that, last known wind value is persisted forward (not real "
          f"future wind - see module docstring).")

    n = len(segment)
    values = segment.values
    wind_raw = wind_on_hs_grid.reindex(segment.index).values  # un-lagged, for persistence beyond the lag horizon
    rows = []

    for origin in range(min_history, n, origin_step):
        y_window = segment.iloc[max(0, origin + 1 - state_window):origin + 1]
        x_window = exog_segment.iloc[max(0, origin + 1 - state_window):origin + 1]
        window_valid = y_window.notna() & x_window.notna()
        if window_valid.sum() < 5:
            continue
        y_w = y_window[window_valid]
        x_w = x_window[window_valid]
        if y_w.index.freq is None:
            try:
                y_w = y_w.asfreq(freq_str)
                x_w = x_w.reindex(y_w.index)
            except Exception:
                continue
        if x_w.isna().any():
            continue

        # Build the forecast exog array: known-lagged values where available,
        # persisted last-known wind beyond that.
        last_known_wind = wind_raw[origin] if not np.isnan(wind_raw[origin]) else None
        if last_known_wind is None:
            continue
        forecast_exog = np.empty(max_h)
        for h in range(1, max_h + 1):
            src_idx = origin + h - lag_samples
            if src_idx <= origin and 0 <= src_idx < n and not np.isnan(wind_raw[src_idx]):
                forecast_exog[h - 1] = wind_raw[src_idx]
            else:
                forecast_exog[h - 1] = last_known_wind

        try:
            applied = armax_fitted.apply(y_w.values, exog=x_w.values, refit=False)
            point_fc = applied.forecast(steps=max_h, exog=forecast_exog)
        except Exception:
            continue

        for h in horizons_samples:
            target_idx = origin + h
            if target_idx >= n:
                continue
            actual = values[target_idx]
            if np.isnan(actual):
                continue
            predicted = float(point_fc[h - 1])
            rows.append({"origin_idx": origin, "horizon_samples": h,
                          "actual": actual, "predicted": predicted,
                          "error": actual - predicted})

    results = pd.DataFrame(rows)
    if len(results) == 0:
        print("No valid forecasts produced - stopping.")
        return
    summary = summarize_backtest(results, dt_hours)

    print(f"\nUsed {results['origin_idx'].nunique()} distinct origins, "
          f"{len(results)} scored pairs.")
    print("\nARMAX per-horizon error:")
    print(summary[["horizon_hours", "n", "rmse", "mae"]].to_string(index=False))

    out_dir = default_paths("20_forecast_armax")

    # Fair local persistence baseline, same segment/origins pattern as Stage C
    local_baseline_results = rolling_origin_backtest(
        segment, persistence_forecast, horizons_samples,
        origin_step=origin_step, min_history=min_history, max_lookback_samples=state_window)
    local_baseline_summary = summarize_backtest(local_baseline_results, dt_hours)
    if len(local_baseline_summary) > 0:
        sk_persist = skill_score(summary, local_baseline_summary)
        print("\nSkill score vs. LOCAL persistence:")
        print(sk_persist.to_string(index=False))
        sk_persist.to_csv(out_dir / f"{args.buoy}_{args.var}_skill_vs_persistence.csv", index=False)

    # Comparison against Stage C's own plain-ARMA result, if available - the
    # real question (does wind help BEYOND ARMA's own structure)
    arma_path = Path("pipeline_out/18_forecast_arma") / f"{args.buoy}_{args.var}_arma_summary.csv"
    if arma_path.exists():
        arma_summary = pd.read_csv(arma_path)
        sk_arma = skill_score(summary, arma_summary)
        print("\nSkill score vs. Stage C's plain ARMA (positive = wind helps "
              "beyond ARMA's own structure, negative = wind doesn't add anything):")
        print(sk_arma.to_string(index=False))
        sk_arma.to_csv(out_dir / f"{args.buoy}_{args.var}_skill_vs_arma.csv", index=False)
    else:
        print(f"\nNo Stage C result found at {arma_path} - run 18_forecast_arma.py "
              f"first for the ARMA-vs-ARMAX comparison (the more important one here).")

    summary.to_csv(out_dir / f"{args.buoy}_{args.var}_armax_summary.csv", index=False)
    with open(out_dir / f"{args.buoy}_{args.var}_armax_meta.json", "w") as f:
        json.dump({"order": list(order), "aic": float(aic),
                    "exog_lag_hours": args.exog_lag_hours,
                    "segment_pct_of_valid": seg_meta["pct_of_valid_used"],
                    "n_origins": int(results["origin_idx"].nunique())}, f, indent=2)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(summary["horizon_hours"], summary["rmse"], marker="o", label=f"ARMAX{order}")
    if arma_path.exists():
        ax.plot(arma_summary["horizon_hours"], arma_summary["rmse"], marker="^",
                label="ARMA (Stage C, no wind)", ls=":")
    if len(local_baseline_summary) > 0:
        ax.plot(local_baseline_summary["horizon_hours"], local_baseline_summary["rmse"],
                marker="s", label="persistence", ls="--")
    ax.axvline(args.exog_lag_hours, color="gray", ls=":", alpha=0.5,
               label=f"exog lag horizon ({args.exog_lag_hours}h)")
    ax.set_xlabel("forecast horizon (hours)")
    ax.set_ylabel(f"{args.var} RMSE (m)")
    ax.set_title(f"{args.buoy} — ARMAX vs. ARMA vs. persistence")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_armax_comparison.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
