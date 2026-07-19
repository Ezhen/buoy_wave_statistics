"""
Batch-runs Stage 25 (change-point detection) across every buoy that
qualifies for the long-record scope (same min-record-years gate Stage
14/25 already use via their own --min-years=10.0 default).

Stages 14/15/16/17-25 are deliberately NOT wired into run_all_buoys.py's
STAGES list - they're run manually per buoy per the README. That also
means there's no existing "the 6 long-record buoys" list anywhere in the
repo; it discovers this from pipeline_out/01_load_clean/*_load_summary.json
(Stage 01 must already have been run for every buoy you want considered).

For each qualifying buoy, runs Stage 25 twice:
  1. Full record.
  2. Trimmed to skip the first --edge-trim-years of THAT BUOY'S OWN
     record start - not a fixed calendar year. Westhinder's manual
     --start-year 1996 trim was specific to its 1990 start; the plan
     notes the other long-record buoys start anywhere 1990-1997, so a
     flat 1996 cutoff would over-trim some and under-trim others.
     Skipping N years from each buoy's own start generalizes the same
     edge-artifact check (a short pre-record segment gives PELT's cost
     model an unstable mean estimate, which can register as a spurious
     change point near the record start) without assuming every buoy
     needs it - it's a check, not a correction applied blindly.

Usage:
    python run_changepoint_batch.py --var VHM0 --min-years 10 --edge-trim-years 5
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


def discover_long_record_buoys(var: str, min_years: float):
    load_dir = Path("pipeline_out/01_load_clean")
    suffix = f"_{var}_load_summary.json"
    buoys = []
    for p in sorted(load_dir.glob(f"*{suffix}")):
        buoy = p.name[: -len(suffix)]
        with open(p) as f:
            info = json.load(f)
        years = info.get("record_years", 0.0)
        if years >= min_years:
            buoys.append((buoy, years))
    return buoys


def record_start_year(buoy: str, var: str) -> int:
    clean_path = Path("pipeline_out/01_load_clean") / f"{buoy}_{var}_clean.csv"
    first_ts = pd.read_csv(clean_path, index_col=0, parse_dates=True, nrows=1).index[0]
    return int(first_ts.year)


def run_stage25(buoy: str, var: str, start_year=None):
    cmd = [sys.executable, "25_changepoint_detection.py", "--buoy", buoy, "--var", var]
    if start_year is not None:
        cmd += ["--start-year", str(start_year)]
    return subprocess.run(cmd, capture_output=True, text=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--min-years", default=10.0, type=float,
                         help="matches Stage 25's own --min-years default - keep these in "
                              "sync or a buoy could be selected here and then refuse to run")
    parser.add_argument("--edge-trim-years", default=5, type=int)
    parser.add_argument("--log-file", default="pipeline_out/changepoint_batch.log", type=Path)
    args = parser.parse_args()

    buoys = discover_long_record_buoys(args.var, args.min_years)
    if not buoys:
        print(f"No buoys found with >= {args.min_years} years of {args.var} in "
              f"pipeline_out/01_load_clean/ - run Stage 01 for all buoys first "
              f"(e.g. via run_all_buoys.py, which always runs Stage 01 regardless "
              f"of a buoy's Core/Advanced tier eligibility for the rest of its list).")
        return

    print(f"Found {len(buoys)} long-record buoy(s) (>= {args.min_years}yr): "
          f"{[b for b, _ in buoys]}")

    args.log_file.parent.mkdir(parents=True, exist_ok=True)
    results = {}
    with open(args.log_file, "a") as log:
        for buoy, years in buoys:
            print(f"\n=== {buoy} ({years:.1f}yr) ===")
            log.write(f"\n{'=' * 70}\n{buoy} ({years:.1f}yr) - full record\n{'=' * 70}\n")
            proc = run_stage25(buoy, args.var)
            log.write(proc.stdout)
            if proc.returncode != 0:
                print(f"  FAILED (full record) - see {args.log_file}")
                log.write(f"[stderr]\n{proc.stderr}\n")
                results[buoy] = "FAILED (full record)"
                continue
            print(proc.stdout)

            start = record_start_year(buoy, args.var) + args.edge_trim_years
            print(f"  -- trimmed re-run, --start-year {start} (skips first "
                  f"{args.edge_trim_years}yr of {buoy}'s own record start) --")
            log.write(f"\n{'-' * 70}\n{buoy} - trimmed to start {start}\n{'-' * 70}\n")
            proc2 = run_stage25(buoy, args.var, start_year=start)
            log.write(proc2.stdout)
            if proc2.returncode != 0:
                print(f"  FAILED (trimmed) - see {args.log_file}")
                log.write(f"[stderr]\n{proc2.stderr}\n")
                results[buoy] = "FAILED (trimmed)"
                continue
            print(proc2.stdout)
            results[buoy] = "ok"

    print("\n=== Batch summary ===")
    for buoy, status in results.items():
        print(f"  {buoy}: {status}")
    print(f"\nFull log: {args.log_file}")


if __name__ == "__main__":
    main()
