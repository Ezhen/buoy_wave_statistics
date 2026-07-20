"""
Aggregates Stage 01's per-buoy `*_load_summary.json` files into one
`buoy_registry.json` - a single place to check "does this buoy qualify
for stage X" instead of the current situation, where that logic is
independently reimplemented in at least three places
(`run_changepoint_batch.py`'s own discovery, this session's ad hoc
Phase-A shell loop, and README prose listing "6 long-record buoys" with
no buoy names actually written down anywhere in the repo).

Per-stage eligibility thresholds below are read directly from each
stage script's own argparse default, NOT assumed to be a single shared
"long-record" cutoff - confirmed they genuinely differ:
  - Stage 14 (Mann-Kendall trend):      --min-years default 10.0
  - Stage 15 (seasonal decomposition):  --min-years default 3.0
  - Stage 25 (change-point detection):  --min-years default 10.0
If any of these stage defaults change, update STAGE_MIN_YEARS below to
match - this file does not read the other scripts' source at runtime,
so it can silently drift out of sync if a threshold changes elsewhere
and this isn't updated too. Grep the other scripts' `--min-years`
argparse lines periodically to check for drift.

Stage 09 (cross-variable analysis) eligibility isn't a year threshold -
it's presence of VTPK and/or VMDR in `available_variables`, which is
reported directly per buoy below rather than as a separate flag.

The "known_quirks" section is manually curated, NOT derived from
load_summary.json - it exists for real per-buoy facts that don't live
in Stage 01's output at all (e.g. Zeebrugge's edge-artifact-prone
record start, motivating the --start-year trims used throughout the
change-point investigation). Extend it by hand as new quirks are found;
don't expect this script to discover them automatically.

Usage:
    python tools/build_buoy_registry.py --var VHM0
    python tools/build_buoy_registry.py --var VHM0 --min-years-override 15
"""

import argparse
import json
from pathlib import Path

STAGE_MIN_YEARS = {
    "stage14_mann_kendall": 10.0,
    "stage15_seasonal_decomposition": 3.0,
    "stage25_changepoint": 10.0,
}

# Manually curated - see module docstring. Keys are buoy names.
KNOWN_QUIRKS = {
    "ZeebruggeZandopvangkadeBuoy": [
        "Record start (1990s-era portion for some vars) prone to PELT "
        "edge artifacts in Stage 25 - every change-point investigation "
        "this session used a --start-year trim relative to this buoy's "
        "own record start, not a shared calendar year.",
        "Structurally distinct from the rest of the network on 5 "
        "independent lines of evidence (see METHODS.md Section 6) - "
        "singleton spatial cluster, unstable distribution-fit, "
        "collapsed wind-coupling, persistent tidal-notch failure, "
        "harbor-interior siting.",
    ],
    "AkkaertSouthwestBuoy": [
        "Anomalously high snap-collision rate in Stage 01 (~2.9% of "
        "its record, ~10-100x every other buoy) - flagged for "
        "follow-up, not yet root-caused. Possibly irregular native "
        "sampling cadence rather than ordinary telemetry jitter.",
    ],
    "CadzandBoei": [
        "Distinct instrument class - finer native sampling (VTPK/VMDR "
        "sensors), reconfirmed via fingerprint clustering.",
    ],
    "Deurlo": [
        "Distinct instrument class - finer native sampling (VTPK/VMDR "
        "sensors), same class as CadzandBoei.",
    ],
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", default="VHM0",
                         help="which variable's load_summary.json files to aggregate")
    parser.add_argument("--load-clean-dir", default=Path("pipeline_out/01_load_clean"),
                         type=Path)
    parser.add_argument("--out", default=Path("buoy_registry.json"), type=Path)
    args = parser.parse_args()

    suffix = f"_{args.var}_load_summary.json"
    files = sorted(args.load_clean_dir.glob(f"*{suffix}"))
    if not files:
        print(f"No {suffix} files found under {args.load_clean_dir} - "
              f"run Stage 01 for the buoys you want in the registry first.")
        return

    registry = {}
    for f in files:
        buoy = f.name[: -len(suffix)]
        with open(f) as fh:
            summary = json.load(fh)

        record_years = summary.get("record_years", 0.0)
        entry = {
            "record_years": record_years,
            "available_variables": summary.get("available_variables", [args.var]),
            "sampling_interval_hours": summary.get("sampling_interval_hours"),
            "multi_era_detected": summary.get("multi_era_detected", False),
            "n_eras": len(summary.get("eras", [1])),
            "pct_missing_after_clean": summary.get("pct_missing_after_clean"),
            "eligibility": {
                stage: record_years >= min_years
                for stage, min_years in STAGE_MIN_YEARS.items()
            },
            "has_vtpk_or_vmdr": any(
                v in summary.get("available_variables", [])
                for v in ("VTPK", "VMDR")
            ),
            "known_quirks": KNOWN_QUIRKS.get(buoy, []),
        }
        registry.setdefault(buoy, {})[args.var] = entry

    with open(args.out, "w") as fh:
        json.dump(registry, fh, indent=2)

    n_eligible_14 = sum(
        1 for b in registry.values()
        for v in b.values() if v["eligibility"]["stage14_mann_kendall"]
    )
    n_multi_era = sum(
        1 for b in registry.values()
        for v in b.values() if v["multi_era_detected"]
    )
    print(f"Wrote {args.out}: {len(registry)} buoy(s), {len(files)} {args.var} entries.")
    print(f"  Clears Stage 14/25 (>={STAGE_MIN_YEARS['stage14_mann_kendall']}yr): {n_eligible_14}")
    print(f"  Multi-era detected: {n_multi_era}")
    print(f"  Known quirks on file for: {sorted(k for k in KNOWN_QUIRKS if k in registry)}")


if __name__ == "__main__":
    main()
