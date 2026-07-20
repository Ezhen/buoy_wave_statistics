"""
Batch-run the validated pipeline across every buoy in data/, gating each
stage by an explicit, declared requirement instead of the three different
ad hoc patterns that had accumulated (Stage 10 silently falling back,
Stage 09 manually excluded from the loop, a script crashing on a missing
variable). Stage 0 always runs first per buoy - it's what discovers which
variables and how much record length that buoy actually has, which every
other stage's eligibility check reads.

Tiers:
  - Core: no requirements, always runs.
  - Advanced: requires specific variables (e.g. VTPK for period-based
    analysis) - most buoys in this network only carry VHM0, so most
    Advanced-gated stages will be skipped for most buoys. That's expected,
    not a failure.
  - Multi-year: requires a minimum record length. Nothing uses this yet
    (no multi-year data downloaded) - the mechanism is ready for when
    Mann-Kendall/seasonal-STL/robust-EVA stages get built.

A SKIPPED stage (requirement not met) is logged distinctly from a FAILED
one (the script actually errored) - skips are expected/normal, failures
are not, and the batch summary keeps them separate.

Usage:
    python run_all_buoys.py --data-dir data --var VHM0
    python run_all_buoys.py --data-dir data --var VHM0 --buoys WesthinderBuoy,A2Buoy

Then:
    python summarize_results.py --var VHM0
to build the cross-buoy comparison table.
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from utils import stage_eligible

# Each entry: (script, arg_template, requirements)
# requirements is one of:
#   {}                                   -> Core, always eligible
#   {"variables_any": ["VTPK", "VMDR"]}  -> eligible if buoy has at least one
#   {"variables_all": ["VTPK"]}          -> eligible if buoy has all listed
#   {"min_record_years": 2.0}            -> eligible once record is long enough
STAGES = [
    ("02_eda_diagnostics.py", ["--buoy", "{buoy}", "--var", "{var}"], {}),
    ("03_stationarity_tests.py", ["--buoy", "{buoy}", "--var", "{var}"], {}),
    ("03b_tidal_notch.py", ["--buoy", "{buoy}", "--var", "{var}", "--harmonics", "2"], {}),
    ("11b_dependence_structure.py", ["--buoy", "{buoy}", "--var", "{var}"], {}),
    ("04_transform_detrend.py", ["--buoy", "{buoy}", "--var", "{var}", "--diff-order", "1"], {}),
    ("05_whiteness_check.py", ["--buoy", "{buoy}", "--var", "{var}"], {}),
    ("06_distribution_fit.py", ["--buoy", "{buoy}", "--var", "{var}"], {}),
    ("07_arch_lm_test.py", ["--buoy", "{buoy}", "--var", "{var}"], {}),
    ("08_extreme_value_analysis.py", ["--buoy", "{buoy}", "--var", "{var}"], {}),
    # 08's --min-separation-hours is decided dynamically per buoy below
    # from Stage 11b's own persistence estimate, not this static default -
    # see build_stage08_args().
    ("12_confidence_intervals.py", ["--buoy", "{buoy}", "--var", "{var}"], {}),
    # 10's --include-period flag is likewise decided dynamically - see
    # build_stage10_args().
    ("10_regime_identification.py", ["--data-dir", "{data_dir}", "--buoy", "{buoy}", "--var", "{var}"], {}),
    ("13_stability_analysis.py", ["--buoy", "{buoy}", "--var", "{var}"], {}),
    ("09_cross_variable_analysis.py", ["--data-dir", "{data_dir}", "--buoy", "{buoy}"],
     {"variables_any": ["VTPK", "VMDR"]}),
]


def discover_buoys(data_dir: Path):
    return sorted(p.stem for p in data_dir.glob("*.nc"))


def build_stage10_args(base_args: list, buoy_info: dict):
    eligible, _ = stage_eligible({"variables_all": ["VTPK"]}, buoy_info)
    if eligible:
        return base_args + ["--include-period"]
    return base_args


def build_stage08_args(base_args: list, buoy: str, var: str):
    """Inject --min-separation-hours from Stage 11b's own per-buoy
    persistence estimate, instead of leaving Stage 08 on its round
    24-48h default network-wide. Found missing after the first
    multi-year batch run: 11b was already computing a real
    per-buoy-justified window (e.g. 184.8h for A2Buoy) but Stage 08
    never consumed it - every buoy's EVA ran on the same generic
    default regardless. Concrete symptom that motivated the fix:
    A2Buoy's GPD xi CI at the default window was [-0.198, 0.065],
    straddling zero - uninformative about tail boundedness, the same
    failure mode too-few-peaks caused earlier, just from a
    too-short-window cause this time."""
    dep_path = Path("pipeline_out/11b_dependence_structure") / f"{buoy}_{var}_dependence_summary.json"
    if dep_path.exists():
        with open(dep_path) as f:
            dep = json.load(f)
        decluster_hours = dep.get("suggested_decluster_hours")
        if decluster_hours:
            return base_args + ["--min-separation-hours", str(decluster_hours)]
    return base_args  # no 11b output yet - falls back to Stage 08's own default


def run_one(script, args_filled, log):
    cmd = [sys.executable, script] + args_filled
    log.write(f"\n$ {' '.join(cmd)}\n")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log.write(proc.stdout)
    if proc.stderr:
        log.write("\n[stderr]\n" + proc.stderr)
    return proc


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

            # Stage 0 always runs first, unconditionally - it's what
            # discovers available_variables/record_years for the gate.
            stage0_args = ["--data-dir", str(args.data_dir), "--buoy", buoy, "--var", args.var]
            proc = run_one("01_load_clean.py", stage0_args, log)
            if proc.returncode != 0:
                print(f"  FAILED at 01_load_clean.py (see {args.log_file})")
                log.write(f"\n*** 01_load_clean.py exited with code {proc.returncode} - "
                          f"stopping remaining stages for {buoy} ***\n")
                results[buoy] = "01_load_clean.py"
                continue
            print("  ok: 01_load_clean.py")

            load_summary_path = (Path("pipeline_out/01_load_clean")
                                  / f"{buoy}_{args.var}_load_summary.json")
            buoy_info = {}
            if load_summary_path.exists():
                with open(load_summary_path) as f:
                    buoy_info = json.load(f)

            failed_at = None
            for script, arg_template, requirements in STAGES:
                if script == "10_regime_identification.py":
                    filled = build_stage10_args(
                        [a.format(data_dir=args.data_dir, buoy=buoy, var=args.var)
                         for a in arg_template],
                        buoy_info)
                elif script == "08_extreme_value_analysis.py":
                    filled = build_stage08_args(
                        [a.format(data_dir=args.data_dir, buoy=buoy, var=args.var)
                         for a in arg_template],
                        buoy, args.var)
                else:
                    filled = [a.format(data_dir=args.data_dir, buoy=buoy, var=args.var)
                              for a in arg_template]

                eligible, reason = stage_eligible(requirements, buoy_info)
                if not eligible:
                    print(f"  SKIPPED: {script} ({reason})")
                    log.write(f"\n$ (skipped) {script} - {reason}\n")
                    continue

                proc = run_one(script, filled, log)
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

        print("\nRunning Stage 12b (correlation confidence intervals)...")
        cmd = [sys.executable, "12b_correlation_confidence.py", "--var", args.var]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        with open(args.log_file, "a") as log:
            log.write(f"\n--- Stage 12b (network-wide) ---\n$ {' '.join(cmd)}\n{proc.stdout}\n")
            if proc.stderr:
                log.write(f"[stderr]\n{proc.stderr}\n")
        print(proc.stdout)
        if proc.returncode != 0:
            print(f"Stage 12b failed (exit {proc.returncode}) - see {args.log_file}")
    else:
        print(f"\nOnly {n_ok} buoy(s) succeeded - skipping Stage 11/12b (need >=3).")

    print("Run summarize_results.py next to build the cross-buoy comparison table.")


if __name__ == "__main__":
    main()
