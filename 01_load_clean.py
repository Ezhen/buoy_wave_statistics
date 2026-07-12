"""
Stage 0 - Data prep.

- Confirms/regularizes the sampling grid
- Interpolates short gaps (linear), flags long gaps instead of bridging them
- Sanity-bounds VHM0 (>= 0, drops absurd spikes)
- Saves a cleaned CSV + a gap report

Usage:
    python 01_load_clean.py --data-dir data --buoy WesthinderBuoy --var VHM0
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from utils import load_buoy_series, default_paths, count_raw_duplicate_timestamps, detect_available_variables

MAX_REASONABLE_VHM0 = 15.0  # meters; North Sea storm Hs rarely exceeds this
SHORT_GAP_LIMIT = 3         # interpolate gaps up to this many samples; longer -> segment


def infer_grid_freq(index: pd.DatetimeIndex) -> pd.Timedelta:
    diffs = pd.Series(index).diff().dropna()
    return diffs.mode().iloc[0]


def regularize_and_clean(s: pd.Series, var: str):
    freq = infer_grid_freq(s.index)
    full_index = pd.date_range(s.index.min(), s.index.max(), freq=freq)
    s_reg = s.reindex(full_index)

    n_missing_before = s_reg.isna().sum()

    # Sanity bounds (VHM0-specific; skip for other vars)
    if var == "VHM0":
        bad = (s_reg < 0) | (s_reg > MAX_REASONABLE_VHM0)
        n_bad = bad.sum()
        s_reg[bad] = np.nan
    else:
        n_bad = 0

    # Identify gap runs
    is_na = s_reg.isna()
    gap_id = (is_na != is_na.shift()).cumsum()
    gap_lengths = is_na.groupby(gap_id).transform("sum")
    short_gap_mask = is_na & (gap_lengths <= SHORT_GAP_LIMIT)
    long_gap_mask = is_na & (gap_lengths > SHORT_GAP_LIMIT)

    s_clean = s_reg.copy()
    s_clean[short_gap_mask] = s_reg.interpolate(method="linear")[short_gap_mask]

    # Longest gap, in samples and hours
    if is_na.any():
        longest_gap_samples = int(gap_lengths[is_na].max())
    else:
        longest_gap_samples = 0

    report = {
        "inferred_freq": str(freq),
        "sampling_interval_hours": freq.total_seconds() / 3600.0,
        "n_samples_regularized": len(s_reg),
        "n_missing_before_clean": int(n_missing_before),
        "n_sanity_flagged_bad": int(n_bad),
        "n_short_gap_interpolated": int(short_gap_mask.sum()),
        "n_long_gap_left_as_nan": int(long_gap_mask.sum()),
        "longest_gap_samples": longest_gap_samples,
        "longest_gap_hours": round(longest_gap_samples * freq.total_seconds() / 3600.0, 2),
        "pct_missing_after_clean": float(s_clean.isna().mean() * 100),
    }
    return s_clean, report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data", type=Path)
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    args = parser.parse_args()

    nc_path = args.data_dir / f"{args.buoy}.nc"
    n_duplicates_raw = count_raw_duplicate_timestamps(nc_path)
    available_vars = detect_available_variables(nc_path)
    s = load_buoy_series(nc_path, args.var)

    s_clean, report = regularize_and_clean(s, args.var)
    report["n_duplicate_timestamps_raw"] = n_duplicates_raw
    report["available_variables"] = available_vars
    report["record_years"] = round(
        (s_clean.index[-1] - s_clean.index[0]).total_seconds() / (3600 * 24 * 365.25), 4)

    out_dir = default_paths("01_load_clean")
    s_clean.to_csv(out_dir / f"{args.buoy}_{args.var}_clean.csv", header=[args.var])

    import json
    with open(out_dir / f"{args.buoy}_{args.var}_load_summary.json", "w") as f:
        json.dump(report, f, indent=2)

    print(f"--- {args.buoy} / {args.var} ---")
    for k, v in report.items():
        print(f"{k:32s}: {v}")

    if report["n_long_gap_left_as_nan"] > 0:
        print("\nNOTE: long gaps were left as NaN (not bridged). "
              "Segment the series around these before Stage 2 if they're substantial, "
              "rather than treating it as one continuous record.")


if __name__ == "__main__":
    main()
