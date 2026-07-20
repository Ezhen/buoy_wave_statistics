"""
Checks whether a buoy's RAW native sampling interval changes over the
record, by computing the modal time-gap per calendar year directly from
the .nc file's TIME array - bypassing Stage 01's single global
regularization grid entirely, since that grid uses ONE frequency
(infer_grid_freq's mode of ALL time diffs) for the whole record. If the
true native interval changes partway through (e.g. an instrument or
telemetry upgrade), a single global frequency is structurally wrong on
one side of the transition - and would explain both an oddly-small
"longest contiguous segment" AND why that segment lands specifically in
the post-change period, without needing any other explanation.

Also reports what Stage 01 actually chose (from its saved
load_summary.json), so a mismatch between "true per-year native
interval" and "what Stage 01 assumed for the whole record" is visible
directly, not inferred indirectly the way Stage 26's segment arithmetic
first surfaced it.

Usage:
    python check_sampling_interval_history.py --data-dir data_multiyear --buoy ZeebruggeZandopvangkadeBuoy --var VHM0
"""

import argparse
import json
import sys
from pathlib import Path

# tools/ scripts sit one level below repo root, where utils.py and all
# pipeline stages live; add that root to sys.path so `from utils import`
# resolves regardless of where the script is invoked from.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import xarray as xr

from utils import resolve_coord_name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data_multiyear", type=Path)
    parser.add_argument("--buoy", default="ZeebruggeZandopvangkadeBuoy")
    parser.add_argument("--var", default="VHM0")
    args = parser.parse_args()

    nc_path = args.data_dir / f"{args.buoy}.nc"
    if not nc_path.exists():
        alt = args.data_dir / f"NO_TS_MO_{args.buoy}.nc"
        nc_path = alt if alt.exists() else nc_path
    if not nc_path.exists():
        raise FileNotFoundError(f"Could not find {args.buoy}.nc under {args.data_dir}")

    with xr.open_dataset(nc_path) as ds:
        time_name = resolve_coord_name(ds, "TIME")
        time = pd.to_datetime(ds[time_name].values)

    time = pd.Series(time).sort_values().reset_index(drop=True)
    diffs_hours = time.diff().dt.total_seconds().dropna() / 3600.0
    years = time.iloc[1:].dt.year.values  # diff i corresponds to the later timestamp

    print(f"=== {nc_path} - raw TIME native interval by year ===\n")
    print(f"Overall modal interval (this is what Stage 01's infer_grid_freq would "
          f"pick for the WHOLE record): {diffs_hours.mode().iloc[0]:.4f}h\n")

    df = pd.DataFrame({"year": years, "diff_hours": diffs_hours.values})
    by_year = df.groupby("year")["diff_hours"].agg(
        modal_interval_hours=lambda x: x.mode().iloc[0],
        n_samples="count",
        pct_at_modal=lambda x: round(100 * (x == x.mode().iloc[0]).mean(), 1),
    )
    print(by_year.to_string())

    # Flag any year whose modal interval differs from the overall mode -
    # a real per-year regime change, not just a couple of odd samples.
    overall_mode = diffs_hours.mode().iloc[0]
    changed_years = by_year.index[
        (by_year["modal_interval_hours"] - overall_mode).abs() > 1e-6
    ].tolist()
    print(f"\nYears whose MODAL interval differs from the overall/global mode "
          f"({overall_mode:.4f}h): {changed_years}")

    # Cross-reference against what Stage 01 actually assumed, if it's been run.
    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            summary = json.load(f)
        stage01_freq_hours = summary.get("sampling_interval_hours")
        print(f"\nStage 01 (01_load_clean.py) actually used: "
              f"{summary.get('inferred_freq')} = {stage01_freq_hours}h for the WHOLE record")
        if changed_years:
            print(f"-> Every sample in a year NOT at this interval would be treated as "
                  f"'missing' by Stage 01's regularization (reindexed against a grid "
                  f"built at {stage01_freq_hours}h), regardless of whether it was "
                  f"actually recorded - this would silently distort coverage %, gap "
                  f"detection, and every downstream stage that reads Stage 01's output.")
    else:
        print(f"\n{load_summary_path} not found - run Stage 01 for this buoy/var to "
              f"compare against what it actually assumed.")


if __name__ == "__main__":
    main()
