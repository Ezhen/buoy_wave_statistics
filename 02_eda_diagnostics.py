"""
Stage 1 - Exploratory diagnostics.

- Rolling mean/variance (volatility clustering check)
- ACF/PACF of the raw series
- Periodogram with the M2 tidal frequency (~12.42 h) marked

Usage:
    python 02_eda_diagnostics.py --buoy WesthinderBuoy --var VHM0
(reads the cleaned CSV from Stage 0)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from scipy.signal import periodogram

from utils import default_paths

M2_PERIOD_HOURS = 12.4206  # principal lunar semi-diurnal tide


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--rolling-window", default="3D",
                         help="pandas offset string, e.g. 3D for 3 days")
    args = parser.parse_args()

    in_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    s = pd.read_csv(in_path, index_col=0, parse_dates=True)[args.var]

    out_dir = default_paths("02_eda_diagnostics")

    # --- Rolling mean/variance ---
    roll_mean = s.rolling(args.rolling_window).mean()
    roll_var = s.rolling(args.rolling_window).var()

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    axes[0].plot(s.index, s.values, lw=0.6, color="steelblue")
    axes[0].set_ylabel(f"{args.var} (raw)")
    axes[1].plot(roll_mean.index, roll_mean.values, color="darkorange")
    axes[1].set_ylabel(f"rolling mean\n({args.rolling_window})")
    axes[2].plot(roll_var.index, roll_var.values, color="firebrick")
    axes[2].set_ylabel(f"rolling variance\n({args.rolling_window})")
    axes[2].set_xlabel("Time")
    fig.suptitle(f"{args.buoy} — {args.var}: level, rolling mean, rolling variance")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_rolling.png", dpi=150)
    plt.close(fig)

    # --- ACF / PACF ---
    s_dropna = s.dropna()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    plot_acf(s_dropna, lags=min(200, len(s_dropna) // 2 - 1), ax=axes[0])
    plot_pacf(s_dropna, lags=min(50, len(s_dropna) // 2 - 1), ax=axes[1], method="ywm")
    axes[0].set_title(f"ACF — {args.var}")
    axes[1].set_title(f"PACF — {args.var}")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_acf_pacf.png", dpi=150)
    plt.close(fig)

    # --- Periodogram with M2 tidal marker ---
    # infer sampling interval in hours
    dt_hours = (s.index[1] - s.index[0]).total_seconds() / 3600.0
    fs = 1.0 / dt_hours  # cycles per hour

    freqs, power = periodogram(s_dropna.values, fs=fs, detrend="linear")
    periods_hours = np.divide(1.0, freqs, out=np.full_like(freqs, np.inf), where=freqs > 0)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.semilogy(periods_hours, power)
    ax.axvline(M2_PERIOD_HOURS, color="red", ls="--",
               label=f"M2 tidal period ({M2_PERIOD_HOURS:.2f} h)")
    ax.set_xlim(0, 48)  # focus on sub-2-day periods where tidal signal would show
    ax.set_xlabel("Period (hours)")
    ax.set_ylabel("Power (log scale)")
    ax.set_title(f"{args.buoy} — {args.var} periodogram (check for M2 tidal peak)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_periodogram.png", dpi=150)
    plt.close(fig)

    # Report power near the M2 peak vs. surrounding baseline
    near_m2 = np.abs(periods_hours - M2_PERIOD_HOURS) < 0.5
    baseline = (periods_hours > 6) & (periods_hours < 24) & ~near_m2
    m2_ratio = None
    if near_m2.any() and baseline.any():
        m2_power = power[near_m2].max()
        baseline_power = np.median(power[baseline])
        m2_ratio = float(m2_power / baseline_power) if baseline_power > 0 else None
        if m2_ratio is not None:
            print(f"M2 peak power / local baseline power ratio: {m2_ratio:.2f}")
            print("(ratio >> 1 suggests a real tidal signature worth notch-filtering later)")

    import json
    with open(out_dir / f"{args.buoy}_{args.var}_eda_summary.json", "w") as f:
        json.dump({"m2_ratio_raw": m2_ratio, "sampling_interval_hours": dt_hours}, f, indent=2)

    print(f"\nSaved diagnostics to {out_dir}")


if __name__ == "__main__":
    main()
