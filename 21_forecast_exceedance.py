"""
Stage F - exceedance forecasting: reframes the question as "will Hs
exceed threshold X at some point in the next H hours" - binary/
probabilistic, not point forecasting. Possibly more operationally
relevant than point forecasts given the DEME-style monitoring context
mentioned early in this project (vessel/operation Hs limits are
threshold-based, not point-estimate-based). Connects naturally to
industrial-monitoring ideas (CUSUM/EWMA) raised in an earlier external
review and parked at the time - this is where they'd actually plug in,
as a monitoring layer on top of an exceedance forecast, not a
standalone addition.

Label definition matters and is stated explicitly: exceedance is
"Hs exceeds the threshold at ANY point within the next H hours" (a max-
based window label), not just "at exactly H hours ahead" - this is the
operationally useful framing ("will conditions get bad at some point
in this window"), not a single-point check.

Model: logistic regression on 3 features computed STRICTLY from data at
or before the origin (no look-ahead) - current Hs level, recent trend
(3h), recent volatility (24h rolling std, motivated by Stage 07's
universal ARCH finding - volatility itself should be informative about
exceedance risk, not just level and trend). Fit ONCE on a training
window, applied at each backtest origin - unlike Stages C/D/E, this
doesn't need a fit-once-apply-cheaply state-management pattern, since
logistic regression prediction is O(1) (no Kalman-filter-style state).

Evaluated via Brier score and ROC-AUC (not RMSE - a probabilistic
classification problem needs classification metrics), compared against
a naive constant-probability baseline (historical base rate) via a
Brier Skill Score, directly analogous to the point-forecast skill score
used in Stages C/D/E.

Usage:
    python 21_forecast_exceedance.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss

from utils import default_paths, longest_contiguous_segment, get_scalar_latlon, load_era5_for_buoy


def build_features(values: np.ndarray, origin: int, trend_lag_samples: int, vol_window_samples: int,
                    wind_values: np.ndarray = None):
    """Features from data STRICTLY at or before origin - no look-ahead.
    Unlike Stage E (ARMAX), exceedance classification only ever needs
    wind AT OR BEFORE the origin as a "current conditions" snapshot
    feature - it's not forecasting the exog variable forward the way
    SARIMAX needed future exog, so there's no honest-persistence-beyond-
    the-lag complexity here. If wind_values is provided, uses
    wind_values[origin] directly (the most recent available ERA5
    observation, already interpolated onto the Hs grid)."""
    current = values[origin]
    if np.isnan(current):
        return None
    trend_src = origin - trend_lag_samples
    if trend_src >= 0 and not np.isnan(values[trend_src]):
        trend = current - values[trend_src]
    else:
        trend = 0.0
    vol_start = max(0, origin - vol_window_samples + 1)
    window_vals = values[vol_start:origin + 1]
    window_vals = window_vals[~np.isnan(window_vals)]
    vol = float(np.std(window_vals)) if len(window_vals) > 5 else 0.0
    feat = [float(current), float(trend), vol]
    if wind_values is not None:
        w = wind_values[origin]
        if np.isnan(w):
            return None
        feat.append(float(w))
    return feat


def build_label(values: np.ndarray, origin: int, horizon_samples: int, threshold: float):
    """1 if Hs exceeds threshold at ANY point in (origin, origin+horizon] - a
    window/max-based label, the operationally useful framing."""
    target_window = values[origin + 1:origin + 1 + horizon_samples]
    target_window = target_window[~np.isnan(target_window)]
    if len(target_window) == 0:
        return None
    return int(np.any(target_window > threshold))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--thresholds-percentile", default="90,95",
                         help="comma-separated percentiles of the buoy's own Hs "
                              "distribution to use as exceedance thresholds")
    parser.add_argument("--horizons-hours", default="6,12,24",
                         help="exceedance window lengths - operationally meaningful "
                              "planning horizons, not the same short horizons used "
                              "for point forecasting")
    parser.add_argument("--trend-lag-hours", default=3.0, type=float)
    parser.add_argument("--vol-window-hours", default=24.0, type=float)
    parser.add_argument("--train-samples", default=8000, type=int)
    parser.add_argument("--origin-step-hours", default=6.0, type=float)
    parser.add_argument("--use-wind", action="store_true",
                         help="add current ERA5 wind_speed as a 4th feature - "
                              "compares against the Hs-only result if available")
    parser.add_argument("--data-dir", default="data_multiyear", type=Path)
    parser.add_argument("--era5-dir", default="meteo_era5")
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"

    hs_full = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var]
    dt_hours = 0.5
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            dt_hours = json.load(f).get("sampling_interval_hours", 0.5)

    print(f"--- {args.buoy} / {args.var} exceedance forecasting ---")

    segment, seg_meta = longest_contiguous_segment(hs_full)
    print(f"Longest contiguous segment: {len(segment)} samples "
          f"({seg_meta['pct_of_valid_used']}% of valid data)")

    thresholds_pct = [float(p) for p in args.thresholds_percentile.split(",")]
    thresholds = {p: float(np.percentile(segment.dropna().values, p)) for p in thresholds_pct}
    print(f"Thresholds: " + ", ".join(f"p{p:.0f}={t:.2f}m" for p, t in thresholds.items()))

    horizons_hours = [float(h) for h in args.horizons_hours.split(",")]
    horizons_samples = {h: max(1, round(h / dt_hours)) for h in horizons_hours}
    trend_lag_samples = max(1, round(args.trend_lag_hours / dt_hours))
    vol_window_samples = max(1, round(args.vol_window_hours / dt_hours))
    origin_step = max(1, round(args.origin_step_hours / dt_hours))

    values = segment.values
    n = len(values)
    if n < args.train_samples + 2000:
        print(f"Segment too short for train+backtest split - stopping.")
        return

    wind_values = None
    if args.use_wind:
        nc_path = args.data_dir / f"{args.buoy}.nc"
        with xr.open_dataset(nc_path) as ds:
            lat, lon = get_scalar_latlon(ds)
        era5 = load_era5_for_buoy(lat, lon, args.era5_dir)
        wind_on_hs_grid = era5["wind_speed"].reindex(
            era5["wind_speed"].index.union(segment.index)).interpolate(method="time").reindex(segment.index)
        wind_values = wind_on_hs_grid.values
        coverage = (~np.isnan(wind_values)).mean()
        print(f"Wind feature enabled: {coverage*100:.1f}% coverage in segment")

    out_dir = default_paths("21_forecast_exceedance")
    all_results = []

    for pct, threshold in thresholds.items():
        for horizon_h, horizon_samples in horizons_samples.items():
            # --- Build training data ---
            train_X, train_y = [], []
            for origin in range(vol_window_samples, args.train_samples, max(1, origin_step // 2)):
                feat = build_features(values, origin, trend_lag_samples, vol_window_samples, wind_values)
                label = build_label(values, origin, horizon_samples, threshold)
                if feat is not None and label is not None:
                    train_X.append(feat)
                    train_y.append(label)
            train_X, train_y = np.array(train_X), np.array(train_y)

            if len(np.unique(train_y)) < 2:
                print(f"  p{pct:.0f}/{horizon_h}h: training labels are all one class - "
                      f"threshold/horizon combination not usable, skipping.")
                continue

            base_rate = train_y.mean()
            model = LogisticRegression().fit(train_X, train_y)

            # --- Backtest ---
            test_X, test_y = [], []
            for origin in range(args.train_samples, n, origin_step):
                feat = build_features(values, origin, trend_lag_samples, vol_window_samples, wind_values)
                label = build_label(values, origin, horizon_samples, threshold)
                if feat is not None and label is not None:
                    test_X.append(feat)
                    test_y.append(label)

            if len(test_y) < 30 or len(np.unique(test_y)) < 2:
                print(f"  p{pct:.0f}/{horizon_h}h: insufficient test data or single-class "
                      f"test labels - skipping.")
                continue

            test_X, test_y = np.array(test_X), np.array(test_y)
            probs = model.predict_proba(test_X)[:, 1]

            brier = brier_score_loss(test_y, probs)
            brier_baseline = brier_score_loss(test_y, np.full_like(probs, base_rate))
            brier_skill = 1 - brier / brier_baseline if brier_baseline > 0 else np.nan
            auc = roc_auc_score(test_y, probs)

            print(f"  p{pct:.0f} ({threshold:.2f}m) / {horizon_h}h: "
                  f"base_rate={base_rate:.3f}, n_test={len(test_y)}, "
                  f"Brier={brier:.4f} (baseline={brier_baseline:.4f}, "
                  f"skill={brier_skill:+.3f}), ROC-AUC={auc:.3f}")

            result_row = {
                "threshold_percentile": pct, "threshold_m": threshold,
                "horizon_hours": horizon_h, "base_rate": float(base_rate),
                "n_train": len(train_y), "n_test": len(test_y),
                "brier_score": float(brier), "brier_baseline": float(brier_baseline),
                "brier_skill_score": float(brier_skill), "roc_auc": float(auc),
                "coef_level": float(model.coef_[0][0]), "coef_trend": float(model.coef_[0][1]),
                "coef_volatility": float(model.coef_[0][2]),
            }
            if args.use_wind:
                result_row["coef_wind"] = float(model.coef_[0][3])
            all_results.append(result_row)

    if not all_results:
        print("No usable threshold/horizon combinations - stopping.")
        return

    results_df = pd.DataFrame(all_results)
    suffix = "_with_wind" if args.use_wind else ""
    results_df.to_csv(out_dir / f"{args.buoy}_{args.var}_exceedance_summary{suffix}.csv", index=False)

    # Compare against the Hs-only result, if it exists
    if args.use_wind:
        hs_only_path = out_dir / f"{args.buoy}_{args.var}_exceedance_summary.csv"
        if hs_only_path.exists():
            hs_only = pd.read_csv(hs_only_path)
            merged = results_df.merge(hs_only, on=["threshold_percentile", "horizon_hours"],
                                       suffixes=("_wind", "_hsonly"))
            merged["brier_skill_delta"] = merged["brier_skill_score_wind"] - merged["brier_skill_score_hsonly"]
            merged["auc_delta"] = merged["roc_auc_wind"] - merged["roc_auc_hsonly"]
            print("\nWind-augmented vs. Hs-only (positive delta = wind helps):")
            print(merged[["threshold_percentile", "horizon_hours", "brier_skill_delta",
                           "auc_delta"]].to_string(index=False))
        else:
            print(f"\nNo Hs-only result found at {hs_only_path} - run without --use-wind "
                  f"first for a comparison.")

    print(f"\n(brier_skill_score > 0 means better than the naive base-rate baseline; "
          f"roc_auc > 0.5 means better than random ranking. coef_* signs show "
          f"which features the model actually relies on.)")

    # --- Plot: ROC-AUC and Brier skill by horizon, one line per threshold ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for pct in thresholds:
        sub = results_df[results_df["threshold_percentile"] == pct].sort_values("horizon_hours")
        if len(sub) == 0:
            continue
        axes[0].plot(sub["horizon_hours"], sub["roc_auc"], marker="o", label=f"p{pct:.0f}")
        axes[1].plot(sub["horizon_hours"], sub["brier_skill_score"], marker="s", label=f"p{pct:.0f}")
    axes[0].axhline(0.5, color="gray", ls="--", alpha=0.5, label="random")
    axes[0].set_xlabel("exceedance window (hours)")
    axes[0].set_ylabel("ROC-AUC")
    axes[0].set_title(f"{args.buoy} — exceedance ROC-AUC")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    axes[1].axhline(0, color="gray", ls="--", alpha=0.5, label="tied with baseline")
    axes[1].set_xlabel("exceedance window (hours)")
    axes[1].set_ylabel("Brier skill score")
    axes[1].set_title(f"{args.buoy} — exceedance Brier skill vs. base rate")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_exceedance_skill.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
