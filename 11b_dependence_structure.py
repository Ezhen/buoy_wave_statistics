"""
Dependence structure: estimate each buoy's integral (decorrelation)
timescale from the ACF, rather than assuming independence anywhere
downstream.

This exists because Ljung-Box fails at all 19 buoys - Hs is confirmed
autocorrelated everywhere - which means:
  - a naive/IID bootstrap on Hs quantiles will understate uncertainty
    (Stage 12) unless block length is tied to a real persistence
    timescale instead of guessed
  - Fisher z CIs on Stage 11's pairwise correlations need an effective-N
    correction, not the raw sample size
  - Stage 08's EVA declustering window (currently a round 24-48h default)
    should be justified by an actual persistence estimate, not just
    "whatever got enough storm peaks"

v2 fixed: (1) runs on the DETIDED series, not raw - a periodic tidal
component makes the ACF cross zero at fractions of the tidal period,
not evidence of losing memory; (2) significance-band cutoff instead of
a single zero-crossing, to avoid both spurious early crossings and
silently reporting a search-ceiling lower bound as a real measurement.

v3 (this version) fixed: on a heavily fragmented multi-year record
(Westhinder: 1905 contiguous segments after gaps), "use the longest
segment" from v2 discarded ~93% of valid data - the longest single
stretch was only 6.9% of the record. Now aggregates the integral
timescale across EVERY qualifying segment (length-weighted mean),
instead of computing on one segment and throwing the rest away. The
diagnostic plot still shows the single longest segment's ACF (most
useful to look at), but the reported number is the aggregate.

Still a simplified version of the idea behind Politis-White optimal
block length, not the full algorithm - stated explicitly so it isn't
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

from utils import default_paths, all_contiguous_segments


def integral_timescale(series: np.ndarray, dt_hours: float, max_lag: int, consecutive: int):
    n = len(series)
    max_lag = min(max_lag, n // 3) if n >= 9 else max(1, n - 1)
    rho = acf(series, nlags=max_lag, fft=True)
    band = 1.96 / np.sqrt(n)

    criterion_lag = None
    for k in range(1, len(rho) - consecutive + 1):
        if np.all(np.abs(rho[k:k + consecutive]) < band):
            criterion_lag = k
            break
    hit_ceiling = criterion_lag is None
    if hit_ceiling:
        criterion_lag = len(rho) - 1

    tau_hours = dt_hours * (1 + 2 * np.sum(rho[1:criterion_lag]))
    return tau_hours, criterion_lag, rho, band, hit_ceiling


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--max-lag", default=2000, type=int,
                         help="max lag (samples) to search for the ACF to settle "
                              "inside the significance band, per segment")
    parser.add_argument("--consecutive", default=5, type=int,
                         help="number of consecutive lags the ACF must stay inside "
                              "the significance band to count as decorrelated")
    parser.add_argument("--min-segment-length", default=200, type=int,
                         help="minimum samples for a contiguous segment to be used "
                              "at all - too short to say anything meaningful about "
                              "a 50-100h+ persistence scale below this")
    parser.add_argument("--max-segments", default=50, type=int,
                         help="cap on how many (longest-first) segments to actually "
                              "process, in case a record is fragmented into hundreds")
    args = parser.parse_args()

    detided_path = Path("pipeline_out/03b_tidal_notch") / f"{args.buoy}_{args.var}_detided_boxcox.csv"
    raw_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"

    used_detided = detided_path.exists()
    if used_detided:
        s_full = pd.read_csv(detided_path, index_col=0, parse_dates=True).iloc[:, 0]
        print(f"Using detided series from Stage 2b: {detided_path}")
    else:
        print(f"WARNING: no Stage 2b (03b_tidal_notch) output found for {args.buoy} - "
              f"falling back to the RAW series. If this buoy has meaningful tidal "
              f"contamination (check its Stage 02 M2 ratio), the resulting timescale "
              f"will likely reflect tidal oscillation, not storm memory. Run "
              f"03b_tidal_notch.py first if that's a concern.")
        s_full = pd.read_csv(raw_path, index_col=0, parse_dates=True)[args.var]

    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"
    dt_hours = 0.5
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            dt_hours = json.load(f).get("sampling_interval_hours", 0.5)

    all_segments = all_contiguous_segments(s_full, min_length=args.min_segment_length)
    n_found = len(all_segments)
    total_valid = int(s_full.notna().sum())

    if n_found == 0:
        print(f"No contiguous segment reaches --min-segment-length={args.min_segment_length} "
              f"samples - cannot estimate a persistence timescale. Lower the threshold "
              f"or check Stage 0's gap report for this buoy.")
        return

    segments_used = all_segments[:args.max_segments]
    n_used = len(segments_used)
    coverage = sum(len(s) for s in segments_used)

    print(f"\n--- {args.buoy} / {args.var} dependence structure ---")
    print(f"Input series: {'detided (Stage 2b)' if used_detided else 'RAW (Stage 0)'}")
    print(f"Sampling interval: {dt_hours}h")
    if n_found > 1:
        print(f"Record fragmented into {n_found} contiguous segments >= "
              f"{args.min_segment_length} samples; longest covers "
              f"{100*len(all_segments[0])/total_valid:.1f}% of valid data alone.")
        print(f"Using {n_used} segment(s) (capped at --max-segments={args.max_segments}), "
              f"covering {coverage} samples ({100*coverage/total_valid:.1f}% of valid data) "
              f"- aggregating across all of them instead of discarding the rest.")

    # --- Per-segment integral timescale ---
    per_segment = []
    for seg in segments_used:
        tau, crit_lag, rho, band, hit_ceiling = integral_timescale(
            seg.values, dt_hours, args.max_lag, args.consecutive)
        per_segment.append({
            "n_samples": len(seg), "tau_hours": tau, "lag1_acf": float(rho[1]),
            "hit_ceiling": hit_ceiling, "start": str(seg.index[0]), "end": str(seg.index[-1]),
        })

    weights = np.array([p["n_samples"] for p in per_segment], dtype=float)
    taus = np.array([p["tau_hours"] for p in per_segment])
    lag1s = np.array([p["lag1_acf"] for p in per_segment])

    tau_hours = float(np.average(taus, weights=weights))
    lag1_acf = float(np.average(lag1s, weights=weights))
    n_ceiling_hits = sum(p["hit_ceiling"] for p in per_segment)

    print(f"\nLength-weighted mean lag-1 ACF: {lag1_acf:.4f}")
    print(f"Length-weighted mean integral timescale: {tau_hours:.2f}h "
          f"(per-segment range: {taus.min():.1f}h - {taus.max():.1f}h across "
          f"{n_used} segment(s))")
    if n_ceiling_hits > 0:
        print(f"WARNING: {n_ceiling_hits}/{n_used} segment(s) hit the --max-lag="
              f"{args.max_lag} search ceiling individually - their per-segment tau is "
              f"a lower bound, which pulls the weighted mean down somewhat. Raise "
              f"--max-lag further if this matters, though each segment is still capped "
              f"at its own length//3.")

    suggested_block_samples = max(1, int(np.ceil(tau_hours / dt_hours)))
    suggested_decluster_hours = round(2 * tau_hours, 1)

    print(f"\nSuggested block bootstrap block length: {suggested_block_samples} samples "
          f"({suggested_block_samples * dt_hours:.1f}h)")
    print(f"Suggested EVA declustering window: {suggested_decluster_hours}h "
          f"(2x integral timescale)")
    if suggested_decluster_hours > 48:
        print("NOTE: exceeds Stage 08's current 48h default - some 'independent' "
              "storm peaks there may be the same storm counted twice.")
    elif suggested_decluster_hours > 24:
        print("NOTE: exceeds the 24h option some Stage 08 runs have used - same caution.")

    out_dir = default_paths("11b_dependence_structure")

    with open(out_dir / f"{args.buoy}_{args.var}_dependence_summary.json", "w") as f:
        json.dump({
            "used_detided_input": used_detided,
            "sampling_interval_hours": dt_hours,
            "lag1_acf": lag1_acf,
            "integral_timescale_hours": tau_hours,
            "suggested_block_length_samples": suggested_block_samples,
            "suggested_decluster_hours": suggested_decluster_hours,
            "n_segments_found": n_found,
            "n_segments_used": n_used,
            "pct_valid_data_used": round(100 * coverage / total_valid, 1),
            "n_segments_hit_ceiling": int(n_ceiling_hits),
            "per_segment_tau_min_hours": float(taus.min()),
            "per_segment_tau_max_hours": float(taus.max()),
            "hit_max_lag_ceiling": bool(n_ceiling_hits > 0),  # kept for backward compat with consumers
        }, f, indent=2)

    # --- Diagnostic plot: the single longest segment's ACF (most informative one to look at) ---
    longest = segments_used[0]
    tau_l, crit_lag_l, rho_l, band_l, hit_ceiling_l = integral_timescale(
        longest.values, dt_hours, args.max_lag, args.consecutive)

    fig, ax = plt.subplots(figsize=(10, 5))
    lags_hours = np.arange(len(rho_l)) * dt_hours
    ax.plot(lags_hours, rho_l, lw=1)
    ax.axhline(0, color="gray", ls="-", lw=0.8)
    ax.axhline(band_l, color="gray", ls=":", lw=0.8, label=f"significance band (+/-{band_l:.3f})")
    ax.axhline(-band_l, color="gray", ls=":", lw=0.8)
    ax.axvline(tau_hours, color="darkorange", ls="--",
               label=f"weighted-mean integral timescale ({tau_hours:.1f}h, "
                     f"across {n_used} segments)")
    ax.set_xlabel("lag (hours)")
    ax.set_ylabel("ACF")
    ax.set_title(f"{args.buoy} — {args.var} ACF on longest segment "
                 f"({len(longest)} samples, {100*len(longest)/total_valid:.1f}% of valid data)\n"
                 f"({'detided' if used_detided else 'RAW - not detided'})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_acf_timescale.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
