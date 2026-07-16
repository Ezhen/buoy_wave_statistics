"""
Build one comparison table across every buoy that's been through the
pipeline, by reading the small summary files each stage writes alongside
its plots (rather than re-parsing console output).

Usage:
    python summarize_results.py --var VHM0 --out pipeline_out/bcz_comparison_summary.csv
"""

import argparse
import json
from pathlib import Path

import pandas as pd

STAGE_DIR = Path("pipeline_out")


def discover_buoys(var: str):
    load_dir = STAGE_DIR / "01_load_clean"
    suffix = f"_{var}_load_summary.json"
    if not load_dir.exists():
        return []
    return sorted(p.name[: -len(suffix)] for p in load_dir.glob(f"*{suffix}"))


def safe_load_json(path: Path):
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def safe_load_csv(path: Path):
    if not path.exists():
        return None
    return pd.read_csv(path)


def summarize_one(buoy: str, var: str) -> dict:
    row = {"buoy": buoy}

    load = safe_load_json(STAGE_DIR / "01_load_clean" / f"{buoy}_{var}_load_summary.json")
    row["n_samples"] = load.get("n_samples_regularized")
    row["sampling_interval_hours"] = load.get("sampling_interval_hours")
    row["pct_missing_after_clean"] = load.get("pct_missing_after_clean")

    eda = safe_load_json(STAGE_DIR / "02_eda_diagnostics" / f"{buoy}_{var}_eda_summary.json")
    row["m2_ratio_raw"] = eda.get("m2_ratio_raw")

    dep = safe_load_json(STAGE_DIR / "11b_dependence_structure" / f"{buoy}_{var}_dependence_summary.json")
    row["lag1_acf"] = dep.get("lag1_acf")
    row["integral_timescale_hours"] = dep.get("integral_timescale_hours")
    row["suggested_decluster_hours"] = dep.get("suggested_decluster_hours")

    stat = safe_load_json(STAGE_DIR / "03_stationarity_tests" / f"{buoy}_{var}_stationarity.json")
    row["adf_pvalue_raw"] = stat.get("adf_pvalue")
    row["kpss_pvalue_raw"] = stat.get("kpss_pvalue")
    row["adf_kpss_agree"] = stat.get("tests_agree")

    notch = safe_load_json(STAGE_DIR / "03b_tidal_notch" / f"{buoy}_{var}_notch_summary.json")
    row["boxcox_lambda"] = notch.get("boxcox_lambda")
    row["m2_ratio_after_notch"] = notch.get("m2_ratio_after")

    detrend = safe_load_json(STAGE_DIR / "04_transform_detrend" / f"{buoy}_{var}_detrend_summary.json")
    row["adf_pvalue_after_detrend"] = detrend.get("adf_pvalue_after")

    lb = safe_load_csv(STAGE_DIR / "05_whiteness_check" / f"{buoy}_{var}_ljungbox.csv")
    if lb is not None and "lb_pvalue" in lb.columns:
        row["ljungbox_min_pvalue"] = lb["lb_pvalue"].min()
        row["ljungbox_all_white"] = bool((lb["lb_pvalue"] > 0.05).all())
        interval = row.get("sampling_interval_hours")
        if interval and lb.iloc[:, 0].notna().all():
            lags_hours = (lb.iloc[:, 0].astype(float) * interval).round(1).tolist()
            row["ljungbox_lags_hours"] = ",".join(f"{h}h" for h in lags_hours)
    else:
        row["ljungbox_min_pvalue"] = None
        row["ljungbox_all_white"] = None

    fit = safe_load_csv(STAGE_DIR / "06_distribution_fit" / f"{buoy}_{var}_fit_summary.csv")
    if fit is not None and len(fit):
        best = fit.loc[fit["ks_stat"].idxmin()]
        row["best_distribution"] = best["distribution"]
        row["best_dist_ks_stat"] = best["ks_stat"]
    else:
        row["best_distribution"] = None
        row["best_dist_ks_stat"] = None

    arch = safe_load_csv(STAGE_DIR / "07_arch_lm_test" / f"{buoy}_{var}_arch_lm.csv")
    row["arch_effects_detected"] = bool((arch["lm_pvalue"] < 0.05).any()) if arch is not None else None

    eva = safe_load_json(STAGE_DIR / "08_extreme_value_analysis" / f"{buoy}_{var}_eva_summary.json")
    row["record_years"] = eva.get("record_years")
    row["n_storm_peaks"] = eva.get("n_peaks")
    row["gpd_shape_xi"] = eva.get("gpd_shape")
    row["eva_fit_reliable"] = eva.get("fit_reliable")

    # Added for the GPD xi external-literature cross-check (README/PLAN):
    # is the network's most extreme xi values (down to -1.31) coming from
    # buoys with reliable (narrow-CI, many-peak) estimates, or from the
    # same small-peak-count buoys already known to have unreliable CIs?
    conf = safe_load_json(STAGE_DIR / "12_confidence_intervals" / f"{buoy}_{var}_confidence_summary.json")
    xi_ci = conf.get("gpd_xi_ci")
    if xi_ci:
        row["gpd_xi_ci_low"] = xi_ci.get("ci_low")
        row["gpd_xi_ci_high"] = xi_ci.get("ci_high")
        if xi_ci.get("ci_low") is not None and xi_ci.get("ci_high") is not None:
            row["gpd_xi_ci_width"] = xi_ci["ci_high"] - xi_ci["ci_low"]
            row["gpd_xi_ci_crosses_zero"] = bool(xi_ci["ci_low"] < 0 < xi_ci["ci_high"])
        else:
            row["gpd_xi_ci_width"] = None
            row["gpd_xi_ci_crosses_zero"] = None
    else:
        row["gpd_xi_ci_low"] = None
        row["gpd_xi_ci_high"] = None
        row["gpd_xi_ci_width"] = None
        row["gpd_xi_ci_crosses_zero"] = None

    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--out", default="pipeline_out/bcz_comparison_summary.csv", type=Path)
    args = parser.parse_args()

    buoys = discover_buoys(args.var)
    if not buoys:
        print(f"No Stage 0 output found for var={args.var} - run the pipeline first.")
        return

    rows = [summarize_one(buoy, args.var) for buoy in buoys]
    df = pd.DataFrame(rows)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 200)
    print(df.to_string(index=False))

    print(f"\nSaved: {args.out}")

    if "eva_fit_reliable" in df.columns:
        unreliable = df[df["eva_fit_reliable"] == False]  # noqa: E712
        if len(unreliable):
            print(f"\nNote: {len(unreliable)} buoy(s) have <10 storm peaks - "
                  f"treat their gpd_shape_xi as noisy, not comparable across the network yet:")
            print(", ".join(unreliable["buoy"].tolist()))

    if "sampling_interval_hours" in df.columns:
        distinct_intervals = df["sampling_interval_hours"].dropna().unique()
        if len(distinct_intervals) > 1:
            print(f"\nWARNING: sampling intervals differ across buoys: {sorted(distinct_intervals)}. "
                  f"ljungbox_min_pvalue and ljungbox_lags_hours are NOT directly comparable "
                  f"across buoys with different intervals - check ljungbox_lags_hours per "
                  f"buoy before comparing whiteness results.")
            mismatched = df.loc[df["sampling_interval_hours"] != df["sampling_interval_hours"].mode()[0],
                                 ["buoy", "sampling_interval_hours"]]
            if len(mismatched):
                print(mismatched.to_string(index=False))

    # GPD xi external-literature cross-check: this network's range (-1.31
    # to -0.04) extends far more negative than a published shallow-water
    # North Sea reference (Caires 2011, xi=-0.12 to -0.13). Is the most
    # extreme end genuine site-specific depth-limiting, or a small-peak-
    # count reliability artifact? Sorted by |xi| so the most extreme
    # values sit next to the exact numbers needed to judge that.
    if "gpd_shape_xi" in df.columns and df["gpd_shape_xi"].notna().any():
        xi_check = df[["buoy", "gpd_shape_xi", "n_storm_peaks", "gpd_xi_ci_width",
                        "gpd_xi_ci_crosses_zero", "eva_fit_reliable"]].dropna(subset=["gpd_shape_xi"])
        xi_check = xi_check.reindex(xi_check["gpd_shape_xi"].abs().sort_values(ascending=False).index)
        print(f"\n{'=' * 70}\nGPD xi cross-check: most extreme first (vs. published shallow-water "
              f"North Sea reference of xi~-0.12 to -0.13, Caires 2011)\n{'=' * 70}")
        print(xi_check.to_string(index=False))
        extreme = xi_check[xi_check["gpd_shape_xi"] < -0.5]
        if len(extreme):
            wide_ci = extreme[(extreme["gpd_xi_ci_width"] > 0.5) | (extreme["gpd_xi_ci_crosses_zero"] == True)]  # noqa: E712
            if len(wide_ci):
                print(f"\nNOTE: {len(wide_ci)} of the most extreme xi value(s) also have a wide "
                      f"or zero-crossing CI ({', '.join(wide_ci['buoy'].tolist())}) - treat these "
                      f"specific extreme values as likely small-sample artifacts, not established "
                      f"physical findings, until checked further.")
            narrow_ci = extreme[~extreme.index.isin(wide_ci.index)] if len(wide_ci) else extreme
            if len(narrow_ci):
                print(f"\nNOTE: {len(narrow_ci)} extreme xi value(s) have a NARROW, non-zero-"
                      f"crossing CI ({', '.join(narrow_ci['buoy'].tolist())}) - these look like "
                      f"genuine, reliably-estimated site-specific depth-limiting, not noise. "
                      f"Worth reporting with real confidence.")


if __name__ == "__main__":
    main()
