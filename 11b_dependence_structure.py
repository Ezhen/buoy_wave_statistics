"""
Dependence structure: estimate each buoy's integral (decorrelation)
timescale from the ACF, rather than assuming independence anywhere
downstream.

This exists because Ljung-Box fails at all 19 buoys - Hs is confirmed
autocorrelated everywhere - which means:
  - a naive/IID bootstrap on Hs quantiles will understate uncertainty
    (Stage 12, not yet built) unless block length is tied to a real
    persistence timescale instead of guessed
  - Fisher z CIs on Stage 11's pairwise correlations need an effective-N
    correction, not the raw sample size
  - Stage 08's EVA declustering window (currently a round 24h default)
    should be justified by an actual persistence estimate, not just
    "whatever got enough storm peaks"

Method: integral timescale via the standard formula
    tau = dt * (1 + 2 * sum_{k=1}^{M} rho(k))
where M is the first lag at which the ACF crosses zero (avoids summing
noise past the point where the signal has decorrelated). This is a
simplified version of the idea behind Politis-White optimal block length,
not the full algorithm - stated explicitly so it isn't over-trusted.

Usage:
    python 11b_dependence_structure.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf

from utils import default_paths


def integral_timescale(series: np.ndarray, dt_hours: float, max_lag: int):
    rho = acf(series, nlags=max_lag, fft=True)
    # first lag where ACF crosses zero (first index >0 with rho <= 0)
    zero_crossing = None
    for k in range(1, len(rho)):
        if rho[k] <= 0:
            zero_crossing = k
            break
    if zero_crossing is None:
        zero_crossing = len(rho) - 1  # never crossed within max_lag - use what we have

    tau_hours = dt_hours * (1 + 2 * np.sum(rho[1:zero_crossing]))
    return tau_hours, zero_crossing, rho


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--max-lag", default=300, type=int,
                         help="max lag (samples) to search for the ACF zero-crossing")
    args = parser.parse_args()

    in_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    s = pd.read_csv(in_path, index_col=0, parse_dates=True)[args.var].dropna()

    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"
    import json
    dt_hours = 0.5
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            dt_hours = json.load(f).get("sampling_interval_hours", 0.5)

    max_lag = min(args.max_lag, len(s) // 4)
    tau_hours, zero_crossing, rho = integral_timescale(s.values, dt_hours, max_lag)

    lag1_acf = float(rho[1])

    print(f"--- {args.buoy} / {args.var} dependence structure ---")
    print(f"Sampling interval: {dt_hours}h")
    print(f"Lag-1 ACF: {lag1_acf:.4f}")
    print(f"ACF first crosses zero at lag {zero_crossing} ({zero_crossing * dt_hours:.1f}h)")
    print(f"Integral timescale (persistence time): {tau_hours:.2f}h")

    suggested_block_samples = max(1, int(np.ceil(tau_hours / dt_hours)))
    suggested_decluster_hours = round(2 * tau_hours, 1)

    print(f"\nSuggested block bootstrap block length: {suggested_block_samples} samples "
          f"({suggested_block_samples * dt_hours:.1f}h) - use this instead of a default "
          f"guess when Stage 12's block bootstrap gets built.")
    print(f"Suggested EVA declustering window: {suggested_decluster_hours}h "
          f"(2x integral timescale) - compare against Stage 08's current default "
          f"before trusting its xi estimates at face value.")
    if suggested_decluster_hours > 48:
        print(f"NOTE: this exceeds Stage 08's current 48h default entirely - "
              f"if so, some 'independent' storm peaks there may be the same storm "
              f"counted twice.")
    elif suggested_decluster_hours > 24:
        print(f"NOTE: this exceeds the 24h option some Stage 08 runs have used - "
              f"same caution applies.")

    out_dir = default_paths("11b_dependence_structure")

    with open(out_dir / f"{args.buoy}_{args.var}_dependence_summary.json", "w") as f:
        json.dump({
            "sampling_interval_hours": dt_hours,
            "lag1_acf": lag1_acf,
            "acf_zero_crossing_lag": int(zero_crossing),
            "integral_timescale_hours": float(tau_hours),
            "suggested_block_length_samples": suggested_block_samples,
            "suggested_decluster_hours": suggested_decluster_hours,
        }, f, indent=2)

    fig, ax = plt.subplots(figsize=(10, 5))
    lags_hours = np.arange(len(rho)) * dt_hours
    ax.plot(lags_hours, rho, lw=1)
    ax.axhline(0, color="gray", ls="-", lw=0.8)
    ax.axvline(zero_crossing * dt_hours, color="firebrick", ls="--",
               label=f"first zero-crossing ({zero_crossing * dt_hours:.1f}h)")
    ax.axvline(tau_hours, color="darkorange", ls="--",
               label=f"integral timescale ({tau_hours:.1f}h)")
    ax.set_xlabel("lag (hours)")
    ax.set_ylabel("ACF")
    ax.set_title(f"{args.buoy} — {args.var} dependence structure")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_acf_timescale.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
