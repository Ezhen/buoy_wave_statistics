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

from utils import (load_buoy_series, default_paths, count_raw_duplicate_timestamps,
                    detect_available_variables, get_provenance)

MAX_REASONABLE_VHM0 = 15.0  # meters; North Sea storm Hs rarely exceeds this
SHORT_GAP_LIMIT = 3         # interpolate gaps up to this many samples; longer -> segment
MIN_ERA_SAMPLES = 500       # a native-interval run must persist this long to count as a
                             # real rate-change era, not a one-off missed sample or blip


def infer_grid_freq(index: pd.DatetimeIndex) -> pd.Timedelta:
    """Kept for backward compatibility (existing callers, tests) - this
    is the OLD single-global-frequency behavior, now only used as the
    fallback path inside detect_sampling_eras when no sustained rate
    change is found."""
    diffs = pd.Series(index).diff().dropna()
    return diffs.mode().iloc[0]


def detect_sampling_eras(index: pd.DatetimeIndex, min_era_samples: int = MIN_ERA_SAMPLES):
    """Partition a raw (irregular) time index into contiguous eras of
    consistent native sampling interval, instead of assuming ONE
    frequency for the whole record.

    Found on real data, not a hypothetical: Zeebrugge sampled at 15min
    from 2009-2017, then 30min from 2018 onward. The old single-global-
    frequency approach picked 15min (more total raw samples in that
    era) and used it for the WHOLE record. For 2018+, every real
    30-minute reading left an empty adjacent 15-minute slot - a
    length-1 gap, under SHORT_GAP_LIMIT, silently linearly interpolated
    and reported as 0% missing. Confirmed: 144,018 of 605,091 total
    regularized samples (23.8% of the full record, concentrated almost
    entirely in the "best covered" recent decade) were fabricated
    interpolation, not real readings - invisible in every existing
    diagnostic because interpolated points don't count as missing. The
    reverse direction (coarse-to-fine) is a different failure mode -
    silent data DISCARD instead of fabrication, since extra samples
    that don't land on the coarser grid are dropped by reindex with no
    trace either.

    Method: round each raw diff to the nearest minute (tolerant of a
    few seconds of instrument jitter; fine enough to distinguish
    genuinely different intervals like 15 vs 20 vs 30 min), then find
    contiguous runs of that rounded value. A run only becomes a
    confirmed era if it persists for >= min_era_samples - a single
    missed sample (which doubles the local diff for one step) or a
    short return-to-normal blip must NOT be mistaken for a rate change;
    only a sustained shift should. Any short run that doesn't clear the
    threshold is absorbed into the PRECEDING confirmed era, since a
    short gap within an era is a missing-data event for the existing
    short/long gap logic to handle, not a new rate regime.

    Returns a list of (start_sample_idx, end_sample_idx_exclusive,
    interval_hours) tuples partitioning the full index with no gaps or
    overlaps.
    """
    n = len(index)
    if n < 2:
        return [(0, n, np.nan)]

    diffs_minutes = (index.to_series().diff().dt.total_seconds() / 60.0).round()
    vals = diffs_minutes.values[1:]  # drop leading NaN; vals[i] = diff between sample i, i+1
    n_diffs = len(vals)

    change_points = np.where(vals[1:] != vals[:-1])[0] + 1
    run_starts = np.concatenate(([0], change_points))
    run_ends = np.concatenate((change_points, [n_diffs]))  # exclusive, in diff-index space

    confirmed = [(s, e, vals[s]) for s, e in zip(run_starts, run_ends)
                 if (e - s) >= min_era_samples]

    if not confirmed:
        # No sustained rate change anywhere - single era at the global
        # mode, identical to the old behavior (regression-safe for
        # every buoy that doesn't have this problem).
        mode_val = pd.Series(vals).mode().iloc[0]
        return [(0, n, mode_val / 60.0)]

    # MERGE adjacent confirmed runs that share the same native interval.
    # A genuine long GAP (missing samples, buoy offline) creates its own
    # tiny, uniquely-valued diff-run - which would otherwise split two
    # runs of the SAME true rate into two separate "eras" that happen to
    # coincide. That's wrong: the rate never changed, there's just a gap
    # inside it, which the existing short/long-gap NaN logic is exactly
    # built to represent - PROVIDED the gap period stays inside a single
    # era's date_range span rather than falling in the dead zone between
    # two artificially-separated eras and being silently dropped. Found
    # on real data: Zeebrugge's real ~2.93-year 2010-2012 outage, sitting
    # between two 15-min-rate stretches, was originally mis-split into
    # two same-rate "eras" this way - erasing the gap period entirely
    # (fewer output rows, record_years shrunk by exactly the gap length)
    # instead of representing it as NaN. Only a genuine VALUE change
    # should start a new era.
    merged = [list(confirmed[0])]
    for s, e, mode_min in confirmed[1:]:
        if abs(mode_min - merged[-1][2]) < 1e-6:
            merged[-1][1] = e  # extend the current era's own run-end
        else:
            merged.append([s, e, mode_min])

    eras = []
    prev_boundary = 0
    for i, (s, e, mode_min) in enumerate(merged):
        era_end = n if i == len(merged) - 1 else e + 1
        eras.append((prev_boundary, era_end, mode_min / 60.0))
        prev_boundary = era_end

    return eras


def regularize_and_clean_one_era(s_era: pd.Series, var: str, freq_hours: float):
    """The original single-frequency regularize/clean logic, applied to
    ONE era's slice of the raw series at ITS OWN native frequency -
    everything downstream of era detection is otherwise unchanged,
    EXCEPT the timestamp-snap step below, which was missing from the
    original code too (this bug predates today's era-aware rewrite).

    Real buoy telemetry timestamps often carry a few seconds of jitter
    around the nominal grid point - detect_sampling_eras tolerates this
    fine (it rounds diffs to the nearest MINUTE before classifying a
    run), but a naive reindex against an exact-second grid does not:
    even 1 second of jitter makes a real sample miss its grid slot
    entirely, and reindex reports that slot as NaN with no indication
    the underlying reading was actually present. Found on real
    Zeebrugge data - one era showed 152,526 of 152,619 slots "missing"
    despite independent coverage checks confirming that era has
    substantial real data; reproduced exactly on synthetic data with
    only +/-3s of random jitter (85.9% apparent missingness from
    jitter alone, zero real gap). Fix: snap each raw timestamp to the
    nearest grid point at this era's own frequency BEFORE reindexing.
    A snap collision (two raw samples rounding to the same grid slot)
    is reported, not silently dropped without a trace."""
    freq = pd.Timedelta(hours=freq_hours)

    s_era = s_era.copy()
    s_era.index = s_era.index.round(freq)
    dup = s_era.index.duplicated(keep="first")
    n_snap_collisions = int(dup.sum())
    if n_snap_collisions > 0:
        s_era = s_era[~dup]

    full_index = pd.date_range(s_era.index.min(), s_era.index.max(), freq=freq)
    s_reg = s_era.reindex(full_index)

    n_missing_before = s_reg.isna().sum()

    if var == "VHM0":
        bad = (s_reg < 0) | (s_reg > MAX_REASONABLE_VHM0)
        n_bad = bad.sum()
        s_reg[bad] = np.nan
    else:
        n_bad = 0

    is_na = s_reg.isna()
    gap_id = (is_na != is_na.shift()).cumsum()
    gap_lengths = is_na.groupby(gap_id).transform("sum")
    short_gap_mask = is_na & (gap_lengths <= SHORT_GAP_LIMIT)
    long_gap_mask = is_na & (gap_lengths > SHORT_GAP_LIMIT)

    s_clean = s_reg.copy()
    s_clean[short_gap_mask] = s_reg.interpolate(method="linear")[short_gap_mask]

    longest_gap_samples = int(gap_lengths[is_na].max()) if is_na.any() else 0

    era_report = {
        "start": str(s_era.index.min()),
        "end": str(s_era.index.max()),
        "inferred_freq_hours": freq_hours,
        "n_samples_regularized": len(s_reg),
        "n_snap_collisions": n_snap_collisions,
        "n_missing_before_clean": int(n_missing_before),
        "n_sanity_flagged_bad": int(n_bad),
        "n_short_gap_interpolated": int(short_gap_mask.sum()),
        "n_long_gap_left_as_nan": int(long_gap_mask.sum()),
        "longest_gap_samples": longest_gap_samples,
        "longest_gap_hours": round(longest_gap_samples * freq_hours, 2),
    }
    return s_clean, era_report


def regularize_and_clean(s: pd.Series, var: str):
    """Era-aware regularization: detect any sustained native-sampling-
    rate change, then regularize each era at its OWN frequency instead
    of forcing one global frequency onto the whole record. For a buoy
    with no rate change (the common case), this produces output
    IDENTICAL to the old single-frequency function - verified as a
    regression test before this replaced it.
    """
    eras = detect_sampling_eras(s.index)

    era_reports = []
    era_series = []
    for start, end, freq_hours in eras:
        s_era = s.iloc[start:end]
        s_clean_era, era_report = regularize_and_clean_one_era(s_era, var, freq_hours)
        era_reports.append(era_report)
        era_series.append(s_clean_era)

    s_clean = pd.concat(era_series).sort_index()
    # Guard against any accidental overlap at an era boundary (shouldn't
    # happen given detect_sampling_eras partitions sample INDICES with
    # no overlap, but a duplicate TIMESTAMP straddling the boundary
    # would silently corrupt aggregate stats below if not caught).
    dup = s_clean.index.duplicated()
    if dup.any():
        s_clean = s_clean[~s_clean.index.duplicated(keep="first")]

    # Aggregate stats across eras, for backward compatibility with every
    # downstream stage that reads these top-level keys. The "dominant"
    # era (most raw samples) is reported as the top-level
    # inferred_freq/sampling_interval_hours, same field names as before -
    # but see multi_era_detected/eras below for the real picture; a
    # single number can no longer represent the whole record honestly
    # once more than one era exists.
    dominant = max(era_reports, key=lambda r: r["n_samples_regularized"])

    n_missing_before_total = sum(r["n_missing_before_clean"] for r in era_reports)
    n_bad_total = sum(r["n_sanity_flagged_bad"] for r in era_reports)
    n_short_total = sum(r["n_short_gap_interpolated"] for r in era_reports)
    n_long_total = sum(r["n_long_gap_left_as_nan"] for r in era_reports)
    longest_gap_hours_overall = max((r["longest_gap_hours"] for r in era_reports), default=0.0)
    longest_gap_samples_overall = max((r["longest_gap_samples"] for r in era_reports), default=0)

    report = {
        "inferred_freq": f"{dominant['inferred_freq_hours']} hours",
        "sampling_interval_hours": dominant["inferred_freq_hours"],
        "n_samples_regularized": len(s_clean),
        "n_missing_before_clean": n_missing_before_total,
        "n_sanity_flagged_bad": n_bad_total,
        "n_short_gap_interpolated": n_short_total,
        "n_long_gap_left_as_nan": n_long_total,
        "longest_gap_samples": longest_gap_samples_overall,
        "longest_gap_hours": longest_gap_hours_overall,
        "pct_missing_after_clean": float(s_clean.isna().mean() * 100),
        "multi_era_detected": len(eras) > 1,
        "eras": era_reports,
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
    report["provenance"] = get_provenance(input_path=nc_path)

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

    if report["multi_era_detected"]:
        print(f"\nNOTE: {len(report['eras'])} distinct native-sampling-rate eras detected - "
              f"this buoy's sampling interval changed mid-record (see 'eras' in the saved "
              f"JSON for exact boundaries/frequencies). The top-level "
              f"'sampling_interval_hours' above ({report['sampling_interval_hours']}h) is "
              f"only the DOMINANT era's value, for backward compatibility with stages that "
              f"read that single field as a lag/dt parameter (11b, 13, 24) - those stages "
              f"are NOT yet era-aware and may use the wrong dt for the non-dominant era(s) "
              f"of this specific buoy's record.")


if __name__ == "__main__":
    main()
