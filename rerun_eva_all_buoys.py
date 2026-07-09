"""
Re-run ONLY Stage 08 (extreme value analysis) across all buoys, at a lower
threshold percentile - without re-running stages 01-07. Stage 08 only
depends on Stage 0's cleaned series (pipeline_out/01_load_clean/), so
this is safe to re-run standalone once that already exists.

Use this after a batch run showed most buoys had too few storm peaks
(<10) at the default 95th-percentile threshold.

Usage:
    python rerun_eva_all_buoys.py --threshold-percentile 85
    python rerun_eva_all_buoys.py --threshold-percentile 90 --min-separation-hours 24
"""

import argparse
import subprocess
import sys
from pathlib import Path

STAGE0_DIR = Path("pipeline_out/01_load_clean")


def discover_buoys(var: str):
    suffix = f"_{var}_clean.csv"
    return sorted(p.name[: -len(suffix)] for p in STAGE0_DIR.glob(f"*{suffix}"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--threshold-percentile", default=85.0, type=float)
    parser.add_argument("--min-separation-hours", default=48.0, type=float)
    parser.add_argument("--return-periods-years", default="1,5,10,25,50")
    parser.add_argument("--buoys", default=None,
                         help="comma-separated buoy names; default = every buoy with Stage 0 output")
    args = parser.parse_args()

    buoys = [b.strip() for b in args.buoys.split(",")] if args.buoys else discover_buoys(args.var)
    if not buoys:
        print(f"No Stage 0 output found for var={args.var} - run the main pipeline first.")
        return

    print(f"Re-running Stage 08 for {len(buoys)} buoy(s) at "
          f"threshold-percentile={args.threshold_percentile}, "
          f"min-separation-hours={args.min_separation_hours}")

    reliable_count = 0
    for buoy in buoys:
        cmd = [
            sys.executable, "08_extreme_value_analysis.py",
            "--buoy", buoy, "--var", args.var,
            "--threshold-percentile", str(args.threshold_percentile),
            "--min-separation-hours", str(args.min_separation_hours),
            "--return-periods-years", args.return_periods_years,
        ]
        print(f"\n=== {buoy} ===")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        print(proc.stdout)
        if proc.stderr:
            print("[stderr]", proc.stderr)
        if proc.returncode != 0:
            print(f"  FAILED (exit code {proc.returncode})")
            continue

    print("\nDone. Run summarize_results.py again to refresh the comparison table "
          "with the updated EVA columns.")


if __name__ == "__main__":
    main()
