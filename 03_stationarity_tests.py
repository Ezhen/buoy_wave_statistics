"""
Stage 2 - Stationarity.

Runs ADF (H0: unit root / non-stationary) and KPSS (H0: stationary) side by
side. Agreement between the two is much stronger evidence than either alone.

Usage:
    python 03_stationarity_tests.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
from pathlib import Path
import json

import pandas as pd
from statsmodels.tsa.stattools import adfuller, kpss

from utils import default_paths, longest_contiguous_segment


def run_tests(series: pd.Series, label: str):
    # ADF/KPSS's own internal lag construction is invalid across a
    # spliced gap - restrict to the longest contiguous segment rather
    # than the full gap-spliced record.
    s, seg_meta = longest_contiguous_segment(series)
    if seg_meta["n_segments"] > 1:
        print(f"Record has {seg_meta['n_segments']} contiguous segments - ADF/KPSS "
              f"run on the longest one only ({seg_meta['pct_of_valid_used']}% of "
              f"valid samples), not the full gap-spliced record.")

    adf_stat, adf_p, adf_lags, adf_nobs, adf_crit, _ = adfuller(s, autolag="AIC")
    kpss_stat, kpss_p, kpss_lags, kpss_crit = kpss(s, regression="c", nlags="auto")

    adf_says = "stationary" if adf_p < 0.05 else "non-stationary (unit root not rejected)"
    kpss_says = "non-stationary" if kpss_p < 0.05 else "stationary (fails to reject)"

    agree = (adf_p < 0.05) == (kpss_p >= 0.05)

    result = {
        "label": label,
        "n": int(len(s)),
        "adf_stat": float(adf_stat),
        "adf_pvalue": float(adf_p),
        "adf_verdict": adf_says,
        "kpss_stat": float(kpss_stat),
        "kpss_pvalue": float(kpss_p),
        "kpss_verdict": kpss_says,
        "tests_agree": bool(agree),
        "n_gap_segments": seg_meta["n_segments"],
        "longest_segment_pct_of_valid": seg_meta["pct_of_valid_used"],
    }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    args = parser.parse_args()

    in_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    s = pd.read_csv(in_path, index_col=0, parse_dates=True)[args.var]

    result = run_tests(s, label="raw")

    out_dir = default_paths("03_stationarity_tests")
    out_path = out_dir / f"{args.buoy}_{args.var}_stationarity.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"--- {args.buoy} / {args.var} (raw series, n={result['n']}) ---")
    print(f"ADF:  stat={result['adf_stat']:.3f}  p={result['adf_pvalue']:.4f}  -> {result['adf_verdict']}")
    print(f"KPSS: stat={result['kpss_stat']:.3f}  p={result['kpss_pvalue']:.4f}  -> {result['kpss_verdict']}")
    print(f"Tests agree: {result['tests_agree']}")
    if not result["tests_agree"]:
        print("Disagreement -> ambiguous case, common with regime-switching series. "
              "Proceed to Stage 3 (transform + differencing) regardless.")

    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
