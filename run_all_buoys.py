"""
Batch-run the validated pipeline (stages 01 -> 08) across every buoy in
data/, logging everything to one file. Continues past a failed buoy
(e.g. a station with no wave sensor) instead of aborting the whole batch.

Usage:
    python run_all_buoys.py --data-dir data --var VHM0
    python run_all_buoys.py --data-dir data --var VHM0 --buoys WesthinderBuoy,A2Buoy

Then:
    python summarize_results.py --var VHM0
to build the cross-buoy comparison table.
"""

import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path

STAGES = [
    ("01_load_clean.py", ["--data-dir", "{data_dir}", "--buoy", "{buoy}", "--var", "{var}"]),
    ("02_eda_diagnostics.py", ["--buoy", "{buoy}", "--var", "{var}"]),
    ("11b_dependence_structure.py", ["--buoy", "{buoy}", "--var", "{var}"]),
    ("03_stationarity_tests.py", ["--buoy", "{buoy}", "--var", "{var}"]),
    ("03b_tidal_notch.py", ["--buoy", "{buoy}", "--var", "{var}", "--harmonics", "2"]),
    ("04_transform_detrend.py", ["--buoy", "{buoy}", "--var", "{var}", "--diff-order", "1"]),
    ("05_whiteness_check.py", ["--buoy", "{buoy}", "--var", "{var}"]),
    ("06_distribution_fit.py", ["--buoy", "{buoy}", "--var", "{var}"]),
    ("07_arch_lm_test.py", ["--buoy", "{buoy}", "--var", "{var}"]),
    ("08_extreme_value_analysis.py", ["--buoy", "{buoy}", "--var", "{var}"]),
    ("10_regime_identification.py", ["--data-dir", "{data_dir}", "--buoy", "{buoy}", "--var", "{var}"]),
]
# 11b runs right after Stage 0/02 (only needs the cleaned level series) so its
# persistence-time numbers are on disk before Stage 08 runs, even though
# Stage 08 doesn't auto-consume them yet - check 11b's output and pass
# --min-separation-hours to Stage 08 manually if its suggestion differs a lot
# from Stage 08's default.
#
# 09_cross_variable_analysis.py is intentionally NOT in the default loop:
# only CadzandBoei and Deurlo carry VTPK/VMDR (confirmed via ncdump across all
# 19 files) - the other 17 always exit with "not enough variables." Run it
# standalone against those two specifically:
#   python 09_cross_variable_analysis.py --data-dir data --buoy CadzandBoei
#   python 09_cross_variable_analysis.py --data-dir data --buoy Deurlo
# Also dropped --include-period from Stage 10's default args here for the
# same reason - it silently falls back to Hs-only for 17/19 buoys anyway;
# run it explicitly with --include-period for CadzandBoei/Deurlo if wanted.


def discover_buoys(data_dir: Path):
    return sorted(p.stem for p in data_dir.glob("*.nc"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data", type=Path)
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--buoys", default=None,
                         help="comma-separated buoy names; default = every .nc in --data-dir")
    parser.add_argument("--log-file", default="pipeline_out/batch_run.log", type=Path)
    args = parser.parse_args()

    args.log_file.parent.mkdir(parents=True, exist_ok=True)

    buoys = [b.strip() for b in args.buoys.split(",")] if args.buoys else discover_buoys(args.data_dir)
    if not buoys:
        print(f"No .nc files found in {args.data_dir}")
        return

    print(f"Running pipeline for {len(buoys)} buoy(s): {buoys}")
    results = {}

    with open(args.log_file, "a") as log:
        log.write(f"\n{'=' * 70}\nBATCH RUN started {datetime.now().isoformat()}\n{'=' * 70}\n")

        for buoy in buoys:
            log.write(f"\n--- {buoy} ---\n")
            print(f"\n=== {buoy} ===")
            failed_at = None

            for script, arg_template in STAGES:
                filled = [a.format(data_dir=args.data_dir, buoy=buoy, var=args.var)
                          for a in arg_template]
                cmd = [sys.executable, script] + filled
                log.write(f"\n$ {' '.join(cmd)}\n")

                proc = subprocess.run(cmd, capture_output=True, text=True)
                log.write(proc.stdout)
                if proc.stderr:
                    log.write("\n[stderr]\n" + proc.stderr)

                if proc.returncode != 0:
                    print(f"  FAILED at {script} (see {args.log_file})")
                    log.write(f"\n*** {script} exited with code {proc.returncode} - "
                              f"stopping remaining stages for {buoy} ***\n")
                    failed_at = script
                    break
                else:
                    print(f"  ok: {script}")

            results[buoy] = failed_at if failed_at else "ok"

        log.write(f"\n{'=' * 70}\nBATCH RUN finished {datetime.now().isoformat()}\n{'=' * 70}\n")
        log.write("\nSummary:\n")
        for buoy, status in results.items():
            log.write(f"  {buoy}: {status}\n")

    print("\n=== Batch summary ===")
    for buoy, status in results.items():
        print(f"  {buoy}: {status}")
    print(f"\nFull log: {args.log_file}")

    n_ok = sum(1 for s in results.values() if s == "ok")
    if n_ok >= 3:
        print(f"\nRunning Stage 11 (spatial statistics) once across all {n_ok} successful buoys...")
        cmd = [sys.executable, "11_spatial_statistics.py", "--data-dir", str(args.data_dir),
               "--var", args.var]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        with open(args.log_file, "a") as log:
            log.write(f"\n--- Stage 11 (network-wide) ---\n$ {' '.join(cmd)}\n{proc.stdout}\n")
            if proc.stderr:
                log.write(f"[stderr]\n{proc.stderr}\n")
        print(proc.stdout)
        if proc.returncode != 0:
            print(f"Stage 11 failed (exit {proc.returncode}) - see {args.log_file}")
    else:
        print(f"\nOnly {n_ok} buoy(s) succeeded - skipping Stage 11 (needs >=3).")

    print("Run summarize_results.py next to build the cross-buoy comparison table.")


if __name__ == "__main__":
    main()
