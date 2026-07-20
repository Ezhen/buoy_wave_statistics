"""
Generic single-stage batch runner: given a stage script and a variable,
discovers which buoys are eligible (via buoy_registry.json's raw facts
+ stage_registry.json's requirements, using the SAME utils.stage_eligible
logic run_all_buoys.py uses for its own default Stage 02-13 batch) and
runs the stage for each eligible buoy.

Exists to replace the pattern of writing a new bespoke discovery
script every time a stage needs batch invocation - run_changepoint_batch.py
did this once, by hand, for Stage 25 specifically; this generalizes
that so the next stage needing the same treatment doesn't need its own
copy of the same discovery logic.

Does NOT change run_all_buoys.py's own default behavior. Stages 14+
are deliberately excluded from that default batch (see its own
docstring - "standalone... run manually," an intentional choice, not
an oversight) - this tool is the explicit, on-demand alternative for
running one of those stages across its eligible buoys, not a way to
fold them into the automatic default.

Usage:
    python tools/build_buoy_registry.py --var VHM0     # if not already built
    python tools/build_stage_registry.py                # if not already built
    python tools/run_stage.py --stage 25_changepoint_detection.py --var VHM0
    python tools/run_stage.py --stage 14_mann_kendall_trend.py --var VHM0 --buoys A2Buoy,WesthinderBuoy
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from utils import stage_eligible  # noqa: E402


def load_registry(path: Path, builder_script: str, build_args: list):
    """Load a registry JSON, building it fresh first if missing - both
    registries are cheap to regenerate and can go stale (e.g. after a
    Stage 01 rerun changes a buoy's record_years), so default to
    freshness rather than trusting a possibly-stale file silently."""
    if not path.exists():
        print(f"{path.name} not found - building it first...")
        subprocess.run([sys.executable, builder_script] + build_args, check=True, cwd=REPO_ROOT)
    with open(path) as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True,
                         help="stage script filename, e.g. 25_changepoint_detection.py")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--buoys", default=None,
                         help="comma-separated buoy names to consider; default = every buoy "
                              "in buoy_registry.json")
    parser.add_argument("--extra-args", default="",
                         help="additional CLI args passed through to the stage script verbatim, "
                              "e.g. --extra-args '--n-components 12'")
    parser.add_argument("--rebuild-registries", action="store_true",
                         help="force-rebuild both registries even if already present")
    parser.add_argument("--log-file", default=Path("pipeline_out/run_stage.log"), type=Path)
    args = parser.parse_args()

    buoy_registry_path = REPO_ROOT / "buoy_registry.json"
    stage_registry_path = REPO_ROOT / "stage_registry.json"

    if args.rebuild_registries:
        buoy_registry_path.unlink(missing_ok=True)
        stage_registry_path.unlink(missing_ok=True)

    buoy_registry = load_registry(
        buoy_registry_path, str(REPO_ROOT / "tools" / "build_buoy_registry.py"),
        ["--var", args.var],
    )
    stage_registry = load_registry(
        stage_registry_path, str(REPO_ROOT / "tools" / "build_stage_registry.py"), [],
    )

    if args.stage not in stage_registry:
        print(f"'{args.stage}' not found in stage_registry.json - check the exact filename "
              f"(available: {sorted(stage_registry.keys())})")
        return
    requirements = stage_registry[args.stage]["requirements"]

    candidate_buoys = (
        [b.strip() for b in args.buoys.split(",")] if args.buoys
        else sorted(buoy_registry.keys())
    )

    eligible_buoys = []
    for buoy in candidate_buoys:
        buoy_info = buoy_registry.get(buoy, {}).get(args.var)
        if buoy_info is None:
            print(f"  SKIPPED {buoy}: no {args.var} entry in buoy_registry.json "
                  f"(run Stage 01 for this buoy/var first)")
            continue
        ok, reason = stage_eligible(requirements, buoy_info)
        if ok:
            eligible_buoys.append(buoy)
        else:
            print(f"  SKIPPED {buoy}: {reason}")

    print(f"\n{len(eligible_buoys)}/{len(candidate_buoys)} buoy(s) eligible for {args.stage} "
          f"(requirements: {requirements or 'none'}): {eligible_buoys}")

    if not eligible_buoys:
        return

    extra = args.extra_args.split() if args.extra_args else []
    args.log_file.parent.mkdir(parents=True, exist_ok=True)
    results = {}

    with open(args.log_file, "a") as log:
        for buoy in eligible_buoys:
            cmd = [sys.executable, args.stage, "--buoy", buoy, "--var", args.var] + extra
            print(f"\n=== {buoy} ===")
            log.write(f"\n=== {buoy} ===\n$ {' '.join(cmd)}\n")
            proc = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT)
            log.write(proc.stdout)
            if proc.stderr:
                log.write(f"\n[stderr]\n{proc.stderr}\n")
            if proc.returncode != 0:
                print(f"  FAILED (see {args.log_file})")
                results[buoy] = "FAILED"
            else:
                print(proc.stdout)
                results[buoy] = "ok"

    print("\n=== Summary ===")
    for buoy, status in results.items():
        print(f"  {buoy}: {status}")
    print(f"\nFull log: {args.log_file}")


if __name__ == "__main__":
    main()
