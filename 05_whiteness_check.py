"""
Stage 4 - Residual whiteness check.

Confirms the residual from Stage 3 has no leftover autocorrelation before
moving on to distribution testing. If this fails, the transform/detrend
step wasn't sufficient - go back to Stage 3 (e.g. check the periodogram
again for a missed tidal peak) rather than proceeding.

Usage:
    python 05_whiteness_check.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
from pathlib import Path

import pandas as pd
from statsmodels.stats.diagnostic import acorr_ljungbox

from utils import default_paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--lags", default="6,12,24,48", help="comma-separated lags to test")
    args = parser.parse_args()

    in_path = Path("pipeline_out/04_transform_detrend") / f"{args.buoy}_{args.var}_residual.csv"
    resid = pd.read_csv(in_path, index_col=0, parse_dates=True)["residual"].dropna()

    lags = [int(x) for x in args.lags.split(",")]
    lb = acorr_ljungbox(resid, lags=lags, return_df=True)

    out_dir = default_paths("05_whiteness_check")
    out_path = out_dir / f"{args.buoy}_{args.var}_ljungbox.csv"
    lb.to_csv(out_path)

    print(f"--- Ljung-Box on residual ({args.buoy} / {args.var}) ---")
    print(lb)

    all_white = (lb["lb_pvalue"] > 0.05).all()
    if all_white:
        print("\nResidual looks white across tested lags -> OK to proceed to distribution testing.")
    else:
        failed = lb.index[lb["lb_pvalue"] <= 0.05].tolist()
        print(f"\nResidual still has significant autocorrelation at lag(s) {failed}.")
        print("Don't proceed to distribution testing yet - revisit Stage 3 "
              "(differencing order, or a missed periodic component like the tidal peak).")

    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
