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
  - Stage 08's EVA declustering window (currently a round 24-48h default)
    should be justified by an actual persistence estimate, not just
    "whatever got enough storm peaks"

FIXED (v2) from the first version, after the first 19-buoy run exposed
two problems:
  1. Runs the ACF on the DETIDED series (Stage 03b output), not Stage 0's
     raw one. A strongly periodic (tidal) component makes the ACF
     oscillate and cross zero at fractions of the tidal period - that's
     the tide's own rhythm, not evidence of losing memory. Zeebrugge (the
     most tidally contaminated buoy in the network) showed the SHORTEST
     apparent persistence in the raw-series version, backwards from what
     genuine storm memory should look like. Falls back to raw with a
     loud warning if Stage 03b hasn't been run for this buoy yet.
  2. Cutoff is now a significance-band criterion (|rho(k)| under the
     1.96/sqrt(n) band, held for several consecutive lags) instead of a
     single zero-crossing - the old criterion let several buoys hit the
     --max-lag search ceiling and report a lower bound as if it were a
     measurement, and is generally noise-sensitive.

Still a simplified version of the idea behind Politis-White optimal block
length, not the full algorithm - stated explicitly so it isn't
over-trusted.

Usage:
    python 11b_dependence_structure.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf

from utils import default_paths


def integral_timescale(series: np.ndarray, dt_hours: float, max_lag: int, consecutive: int):
    n = len(series)
    rho = acf(series, nlags=max_lag, fft=True)
    band = 1.96 / np.sqrt(n)

    criterion_lag = None
    for k in range(1, len(rho) - consecutive + 1):
        if np.all(np.abs(rho[k:k + consecutive]) < band):
            criterion_lag = k
            break
    hit_ceiling = criterion_lag is None
    if hit_ceiling:
        criterion_lag = len(rho) - 1  # never stayed inside the band - use what we have as a lower bound

    tau_hours = dt_hours * (1 + 2 * np.sum(rho[1:criterion_lag]))
    return tau_hours, criterion_lag, rho, band, hit_ceiling


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--max-lag", default=500, type=int,
                         help="max lag (samples) to search for the ACF to settle inside the significance band")
    parser.add_argument("--consecutive", default=5, type=int,
                         help="number of consecutive lags the ACF must stay inside "
                              "the significance band to count as decorrelated")
    args = parser.parse_args()

    detided_path = Path("pipeline_out/03b_tidal_notch") / f"{args.buoy}_{args.var}_detided_boxcox.csv"
    raw_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"

    used_detided = detided_path.exists()
    if used_detided:
        s = pd.read_csv(detided_path, index_col=0, parse_dates=True).iloc[:, 0].dropna()
        print(f"Using detided series from Stage 2b: {detided_path}")
    else:
        print(f"WARNING: no Stage 2b (03b_tidal_notch) output found for {args.buoy} - "
              f"falling back to the RAW series. If this buoy has meaningful tidal "
              f"contamination (check its Stage 02 M2 ratio), the resulting timescale "
              f"will likely reflect tidal oscillation, not storm memory. Run "
              f"03b_tidal_notch.py first if that's a concern.")
        s = pd.read_csv(raw_path, index_col=0, parse_dates=True)[args.var].dropna()

    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"
    dt_hours = 0.5
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            dt_hours = json.load(f).get("sampling_interval_hours", 0.5)

    max_lag = min(args.max_lag, len(s) // 3)
    tau_hours, criterion_lag, rho, band, hit_ceiling = integral_timescale(
        s.values, dt_hours, max_lag, args.consecutive)

    lag1_acf = float(rho[1])

    print(f"\n--- {args.buoy} / {args.var} dependence structure ---")
    print(f"Input series: {'detided (Stage 2b)' if used_detided else 'RAW (Stage 0)'}")
    print(f"Sampling interval: {dt_hours}h")
    print(f"Significance band: +/-{band:.4f}")
    print(f"Lag-1 ACF: {lag1_acf:.4f}")
    if hit_ceiling:
        print(f"WARNING: ACF never stayed inside the significance band for "
              f"{args.consecutive} consecutive lags within --max-lag={max_lag} "
              f"({max_lag * dt_hours:.1f}h). Reported timescale is a LOWER BOUND, "
              f"not a measurement - rerun with a higher --max-lag.")
    else:
        print(f"ACF settles inside the significance band at lag {criterion_lag} "
              f"({criterion_lag * dt_hours:.1f}h)")
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
            "used_detided_input": used_detided,
            "sampling_interval_hours": dt_hours,
            "lag1_acf": lag1_acf,
            "acf_criterion_lag": int(criterion_lag),
            "hit_max_lag_ceiling": bool(hit_ceiling),
            "integral_timescale_hours": float(tau_hours),
            "suggested_block_length_samples": suggested_block_samples,
            "suggested_decluster_hours": suggested_decluster_hours,
        }, f, indent=2)

    fig, ax = plt.subplots(figsize=(10, 5))
    lags_hours = np.arange(len(rho)) * dt_hours
    ax.plot(lags_hours, rho, lw=1)
    ax.axhline(0, color="gray", ls="-", lw=0.8)
    ax.axhline(band, color="gray", ls=":", lw=0.8, label=f"significance band (+/-{band:.3f})")
    ax.axhline(-band, color="gray", ls=":", lw=0.8)
    ax.axvline(criterion_lag * dt_hours, color="firebrick", ls="--",
               label=f"{'lower bound (hit ceiling)' if hit_ceiling else 'settles in band'} "
                     f"({criterion_lag * dt_hours:.1f}h)")
    ax.axvline(tau_hours, color="darkorange", ls="--",
               label=f"integral timescale ({tau_hours:.1f}h)")
    ax.set_xlabel("lag (hours)")
    ax.set_ylabel("ACF")
    ax.set_title(f"{args.buoy} — {args.var} dependence structure "
                 f"({'detided' if used_detided else 'RAW - not detided'})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_acf_timescale.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
