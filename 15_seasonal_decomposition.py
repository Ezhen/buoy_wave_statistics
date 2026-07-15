"""
Seasonal decomposition (STL) - separate real winter-storm-season vs.
summer-calm seasonality from the M2 tide, which is the only
"seasonality" this pipeline has been able to detect on shorter records.

Runs on MONTHLY aggregates (mean and p95 Hs separately, same
typical-climate vs. storm-intensity split as Stage 14's Mann-Kendall),
not raw high-frequency data - a 12-month seasonal cycle needs monthly
resolution to estimate cleanly; raw 30-min data would be dominated by
the tidal/diurnal signal Stage 03b already handles, not the annual one
this stage is after.

Only meaningful with several full annual cycles - default minimum 3
years to let statsmodels' STL run at all (needs >= 2 periods), but
treat anything under ~10 years with real caution; this network's 6
long-record buoys (30+ years) are the ones this was really built for.

CAVEAT worth internalizing before reading output: STL always fits SOME
period-12-shaped component, even to pure noise with no real seasonality
(confirmed on a synthetic no-seasonal-signal control: seasonal
"variance fraction" came out 0.30-0.40 despite zero true seasonality,
vs. 0.99+ on a control with a real injected seasonal signal). The
peak-to-trough AMPLITUDE is a much more reliable signal-vs-noise
discriminator than the variance-fraction number - a small amplitude
(a few cm) with an inconsistent peak/trough month across variables is
likely fitting noise; a large, consistent amplitude with a physically
sensible peak month (winter) is the real thing.

Usage:
    python 15_seasonal_decomposition.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.seasonal import STL

from utils import default_paths


def build_monthly_series(hs: pd.Series, agg: str):
    if agg == "mean":
        monthly = hs.resample("MS").mean()
    else:
        monthly = hs.resample("MS").quantile(0.95)

    full_index = pd.date_range(monthly.index.min(), monthly.index.max(), freq="MS")
    monthly = monthly.reindex(full_index)

    n_missing = int(monthly.isna().sum())
    if n_missing > 0:
        monthly = monthly.interpolate(limit_direction="both")

    return monthly, n_missing


def run_stl(monthly: pd.Series, robust: bool = True):
    result = STL(monthly, period=12, robust=robust).fit()
    return result


def seasonal_by_calendar_month(seasonal: pd.Series):
    return seasonal.groupby(seasonal.index.month).mean()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--min-years", default=3.0, type=float,
                         help="refuse to run below this - STL needs multiple "
                              "full annual cycles to estimate seasonality at all")
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

    print(f"--- {args.buoy} / {args.var} seasonal decomposition (STL) ---")
    print(f"Record length: {record_years:.1f} years")
    if record_years < args.min_years:
        print(f"Record is under --min-years={args.min_years} - refusing to run. "
              f"STL needs multiple full annual cycles to estimate seasonality "
              f"at all, let alone reliably.")
        return
    if record_years < 10:
        print(f"NOTE: {record_years:.1f} years clears --min-years but is short "
              f"for a stable seasonal estimate (10+ years preferred). Treat the "
              f"seasonal component with real caution.")

    out_dir = default_paths("15_seasonal_decomposition")
    summary = {"record_years": record_years}

    for agg, label in [("mean", "monthly mean Hs (typical climate)"),
                        ("p95", "monthly p95 Hs (storm intensity)")]:
        print(f"\n[{label}]")
        monthly, n_missing = build_monthly_series(hs, agg)
        n_total = len(monthly)
        if n_missing > 0:
            frac = n_missing / n_total
            extra = " NOTE: this is a substantial fraction - treat results with extra caution." if frac > 0.1 else ""
            print(f"  {n_missing}/{n_total} months had no data - linearly "
                  f"interpolated to give STL a complete series (STL cannot "
                  f"handle NaN).{extra}")

        if n_total < 24:
            print(f"  Only {n_total} months - too few for STL (needs >= 24, "
                  f"2 full periods). Skipping.")
            summary[agg] = {"error": f"only {n_total} months, need >= 24"}
            continue

        result = run_stl(monthly)
        trend, seasonal, resid = result.trend, result.seasonal, result.resid

        seasonal_amplitude = float(seasonal.max() - seasonal.min())
        by_month = seasonal_by_calendar_month(seasonal)
        peak_month = int(by_month.idxmax())
        trough_month = int(by_month.idxmin())

        # variance explained by each component, as a rough diagnostic
        total_var = monthly.var()
        seasonal_var_frac = float(seasonal.var() / total_var) if total_var > 0 else float("nan")
        trend_var_frac = float(trend.var() / total_var) if total_var > 0 else float("nan")

        print(f"  Seasonal amplitude (peak-to-trough): {seasonal_amplitude:.3f} m")
        print(f"  Peak month: {peak_month} ({pd.Timestamp(2000, peak_month, 1).strftime('%B')}), "
              f"trough month: {trough_month} ({pd.Timestamp(2000, trough_month, 1).strftime('%B')})")
        print(f"  Variance fraction - seasonal: {seasonal_var_frac:.3f}, trend: {trend_var_frac:.3f}")

        summary[agg] = {
            "n_missing_months": n_missing,
            "n_total_months": n_total,
            "seasonal_amplitude_m": seasonal_amplitude,
            "peak_month": peak_month,
            "trough_month": trough_month,
            "seasonal_variance_fraction": seasonal_var_frac,
            "trend_variance_fraction": trend_var_frac,
            "seasonal_by_month": {int(m): float(v) for m, v in by_month.items()},
        }

        fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
        axes[0].plot(monthly.index, monthly.values, lw=0.8, color="steelblue")
        axes[0].set_ylabel("observed")
        axes[1].plot(trend.index, trend.values, lw=1, color="darkorange")
        axes[1].set_ylabel("trend")
        axes[2].plot(seasonal.index, seasonal.values, lw=0.8, color="firebrick")
        axes[2].set_ylabel("seasonal")
        axes[3].plot(resid.index, resid.values, lw=0.5, color="gray")
        axes[3].axhline(0, color="black", lw=0.5)
        axes[3].set_ylabel("residual")
        fig.suptitle(f"{args.buoy} — STL decomposition, {label}")
        fig.tight_layout()
        fig.savefig(out_dir / f"{args.buoy}_{args.var}_{agg}_stl.png", dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5))
        month_names = [pd.Timestamp(2000, m, 1).strftime("%b") for m in by_month.index]
        ax.bar(month_names, by_month.values, color="steelblue", alpha=0.7)
        ax.axhline(0, color="black", lw=0.5)
        ax.set_ylabel("mean seasonal effect (m)")
        ax.set_title(f"{args.buoy} — seasonal cycle by calendar month, {label}")
        fig.tight_layout()
        fig.savefig(out_dir / f"{args.buoy}_{args.var}_{agg}_seasonal_by_month.png", dpi=150)
        plt.close(fig)

    with open(out_dir / f"{args.buoy}_{args.var}_seasonal_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
