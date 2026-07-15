"""
Run Stage 16 (wind-wave coupling) across every buoy with Stage 0 output,
then build one cross-buoy comparison table - same pattern as
rerun_eva_all_buoys.py / summarize_results.py.

Usage:
    python run_all_wind_coupling.py --data-dir data_multiyear --var VHM0
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


def discover_buoys(var: str):
    load_dir = Path("pipeline_out/01_load_clean")
    suffix = f"_{var}_clean.csv"
    if not load_dir.exists():
        return []
    return sorted(p.name[: -len(suffix)] for p in load_dir.glob(f"*{suffix}"))


def run_one(buoy: str, data_dir: Path, var: str, era5_dir: str):
    cmd = [sys.executable, "16_wind_wave_coupling.py", "--data-dir", str(data_dir),
           "--buoy", buoy, "--var", var, "--era5-dir", era5_dir]
    print(f"\n=== {buoy} ===")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        print(f"  FAILED (exit {proc.returncode})")
        if proc.stderr:
            print(proc.stderr[-2000:])
        return False
    return True


def build_summary_table(buoys, var: str, out_path: Path):
    rows = []
    for buoy in buoys:
        summary_path = (Path("pipeline_out/16_wind_wave_coupling")
                         / f"{buoy}_{var}_wind_coupling_summary.json")
        if not summary_path.exists():
            continue
        with open(summary_path) as f:
            s = json.load(f)
        rows.append({
            "buoy": buoy,
            "lat": s.get("buoy_lat"),
            "lon": s.get("buoy_lon"),
            "n_overlap_samples": s.get("n_overlap_samples"),
            "wind_hs_best_lag_hours": s.get("wind_hs_best_lag_hours"),
            "wind_hs_best_corr": s.get("wind_hs_best_corr"),
            "mslp_hs_best_lag_hours": s.get("mslp_hs_best_lag_hours"),
            "mslp_hs_best_corr": s.get("mslp_hs_best_corr"),
            "r2_wind_linear": s.get("r2_wind_linear"),
            "r2_wind_quadratic": s.get("r2_wind_quadratic"),
            "mean_wind_wave_dir_diff_deg": s.get("mean_wind_wave_dir_diff_deg"),
        })

    if not rows:
        print("No wind-coupling summaries found - nothing to aggregate.")
        return

    df = pd.DataFrame(rows).sort_values("buoy")
    df.to_csv(out_path, index=False)

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 200)
    print(f"\n{'=' * 70}\nCross-buoy wind-wave coupling summary\n{'=' * 70}")
    print(df.to_string(index=False))
    print(f"\nSaved: {out_path}")

    if df["wind_hs_best_lag_hours"].notna().any():
        print(f"\nWind lag range across network: "
              f"{df['wind_hs_best_lag_hours'].min():.0f}h to "
              f"{df['wind_hs_best_lag_hours'].max():.0f}h "
              f"(median {df['wind_hs_best_lag_hours'].median():.0f}h)")
    if df["r2_wind_linear"].notna().any():
        print(f"R^2 (linear) range: {df['r2_wind_linear'].min():.3f} to "
              f"{df['r2_wind_linear'].max():.3f} "
              f"(median {df['r2_wind_linear'].median():.3f})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data_multiyear", type=Path)
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--era5-dir", default="meteo_era5")
    parser.add_argument("--buoys", default=None,
                         help="comma-separated buoy names; default = every buoy "
                              "with Stage 0 output")
    args = parser.parse_args()

    buoys = [b.strip() for b in args.buoys.split(",")] if args.buoys else discover_buoys(args.var)
    if not buoys:
        print(f"No Stage 0 output found for var={args.var} - run the main pipeline first.")
        return

    print(f"Running Stage 16 for {len(buoys)} buoy(s): {buoys}")
    succeeded = []
    for buoy in buoys:
        ok = run_one(buoy, args.data_dir, args.var, args.era5_dir)
        if ok:
            succeeded.append(buoy)

    print(f"\n{len(succeeded)}/{len(buoys)} succeeded.")
    build_summary_table(succeeded, args.var,
                         Path("pipeline_out/16_wind_wave_coupling") / f"{args.var}_wind_coupling_comparison.csv")


if __name__ == "__main__":
    main()
