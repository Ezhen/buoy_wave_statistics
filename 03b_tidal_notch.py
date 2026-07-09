"""
Stage 2b - Tidal notch (harmonic regression).

Slots between Stage 2 (stationarity) and Stage 3 (differencing). The
periodogram in Stage 02 showed a dominant M2 peak (~12.42h), and Stage 04/05
showed lag-1 differencing does nothing to it - a ~24.84-sample period isn't
an integer lag, so seasonal differencing can't target it cleanly either.

Instead: fit and subtract M2 (+ its first harmonic, since coastal tidal
signatures are rarely pure sinusoids) via least-squares harmonic regression,
done on the Box-Cox-transformed LEVEL series (variance-stabilize first,
same ordering rule as before - then detide, then Stage 04 only needs to
difference what's left).

Usage:
    python 03b_tidal_notch.py --buoy WesthinderBuoy --var VHM0 --harmonics 2
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import boxcox
from scipy.signal import periodogram

from utils import default_paths

M2_PERIOD_HOURS = 12.4206


def build_harmonic_design(t_hours: np.ndarray, period_hours: float, n_harmonics: int):
    """Columns: [1, cos(w t), sin(w t), cos(2w t), sin(2w t), ...]"""
    cols = [np.ones_like(t_hours)]
    for k in range(1, n_harmonics + 1):
        w = 2 * np.pi * k / period_hours
        cols.append(np.cos(w * t_hours))
        cols.append(np.sin(w * t_hours))
    return np.column_stack(cols)


def periodogram_m2_ratio(series: np.ndarray, dt_hours: float):
    fs = 1.0 / dt_hours
    freqs, power = periodogram(series, fs=fs, detrend="linear")
    periods = np.divide(1.0, freqs, out=np.full_like(freqs, np.inf), where=freqs > 0)
    near_m2 = np.abs(periods - M2_PERIOD_HOURS) < 0.5
    baseline = (periods > 6) & (periods < 24) & ~near_m2
    if near_m2.any() and baseline.any() and np.median(power[baseline]) > 0:
        return power[near_m2].max() / np.median(power[baseline])
    return np.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--harmonics", default=2, type=int,
                         help="number of M2 harmonics to fit (1=pure sinusoid, "
                              "2=+first overtone, etc.)")
    args = parser.parse_args()

    in_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    s = pd.read_csv(in_path, index_col=0, parse_dates=True)[args.var].dropna()

    # --- Variance-stabilize first (same rule as Stage 3/04) ---
    shifted = s + 1e-6 if (s <= 0).any() else s
    boxcox_vals, lam = boxcox(shifted.values)
    s_bc = pd.Series(boxcox_vals, index=s.index, name=f"{args.var}_boxcox")
    print(f"Box-Cox lambda = {lam:.4f}")

    # --- Harmonic regression fit ---
    t0 = s_bc.index[0]
    t_hours = (s_bc.index - t0).total_seconds().values / 3600.0
    X = build_harmonic_design(t_hours, M2_PERIOD_HOURS, args.harmonics)
    coefs, *_ = np.linalg.lstsq(X, s_bc.values, rcond=None)
    tidal_fit = X @ coefs
    s_detided = pd.Series(s_bc.values - tidal_fit, index=s_bc.index, name=f"{args.var}_detided")

    # --- Report M2 power before/after ---
    dt_hours = (s.index[1] - s.index[0]).total_seconds() / 3600.0
    ratio_before = periodogram_m2_ratio(s_bc.values, dt_hours)
    ratio_after = periodogram_m2_ratio(s_detided.values, dt_hours)
    print(f"M2 peak/baseline power ratio - before: {ratio_before:.2f}, after: {ratio_after:.2f}")
    if ratio_after > 3:
        print("WARNING: M2 power still elevated after notch - "
              "try --harmonics with a higher value, or check the sampling grid for gaps "
              "that would smear the harmonic fit.")
    else:
        print("M2 signature substantially reduced -> OK to proceed to Stage 3 differencing "
              "on the detided series.")

    out_dir = default_paths("03b_tidal_notch")
    s_detided.to_csv(out_dir / f"{args.buoy}_{args.var}_detided_boxcox.csv",
                      header=[s_detided.name])

    import json
    with open(out_dir / f"{args.buoy}_{args.var}_notch_summary.json", "w") as f:
        json.dump({
            "boxcox_lambda": float(lam),
            "m2_ratio_before": None if np.isnan(ratio_before) else float(ratio_before),
            "m2_ratio_after": None if np.isnan(ratio_after) else float(ratio_after),
            "harmonics": args.harmonics,
        }, f, indent=2)

    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    axes[0].plot(s_bc.index, s_bc.values, lw=0.6)
    axes[0].set_title(f"{args.var} — Box-Cox level (pre-detide)")
    axes[1].plot(s_bc.index, tidal_fit, lw=0.8, color="darkorange")
    axes[1].set_title(f"Fitted M2 + harmonics (n_harmonics={args.harmonics})")
    axes[2].plot(s_detided.index, s_detided.values, lw=0.6, color="firebrick")
    axes[2].set_title(f"{args.var} — detided (ready for Stage 3 differencing)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_tidal_notch.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved detided series + plot to {out_dir}")


if __name__ == "__main__":
    main()
