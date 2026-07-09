"""
ARCH-LM test (Engle) - formal test for volatility clustering.

Tests whether squared residuals are themselves autocorrelated, i.e.
whether local variance is time-varying (calm vs. storm regimes) rather
than constant. Confirms/refutes what the Stage 02 rolling-variance plot
suggested visually, independent of the forecasting question.

Runs on the Stage 3 (04_transform_detrend) residual - the detided,
differenced series - since ARCH-LM is a test on the residual of a fitted
mean model, not on the raw level series.

Usage:
    python 07_arch_lm_test.py --buoy WesthinderBuoy --var VHM0 --lags 6,12,24,48
"""

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.diagnostic import het_arch
from statsmodels.graphics.tsaplots import plot_acf

from utils import default_paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--lags", default="6,12,24,48",
                         help="comma-separated lag orders to test")
    args = parser.parse_args()

    in_path = Path("pipeline_out/04_transform_detrend") / f"{args.buoy}_{args.var}_residual.csv"
    resid = pd.read_csv(in_path, index_col=0, parse_dates=True)["residual"].dropna()

    out_dir = default_paths("07_arch_lm_test")
    lags = [int(x) for x in args.lags.split(",")]

    rows = []
    print(f"--- ARCH-LM test ({args.buoy} / {args.var} residual, n={len(resid)}) ---")
    for lag in lags:
        lm_stat, lm_p, f_stat, f_p = het_arch(resid, nlags=lag)
        rows.append({"lag": lag, "lm_stat": lm_stat, "lm_pvalue": lm_p,
                     "f_stat": f_stat, "f_pvalue": f_p})
        print(f"lag={lag:3d}: LM stat={lm_stat:8.2f}  p={lm_p:.4g}   "
              f"F stat={f_stat:8.2f}  p={f_p:.4g}")

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / f"{args.buoy}_{args.var}_arch_lm.csv", index=False)

    any_sig = (df["lm_pvalue"] < 0.05).any()
    if any_sig:
        print("\nARCH effects detected at one or more lags -> volatility clustering "
              "confirmed formally. Matches the storm-driven rolling-variance spikes "
              "seen in Stage 02. If you later build a forecasting stage, this is the "
              "evidence for preferring ARMA-GARCH over plain ARMA.")
    else:
        print("\nNo significant ARCH effects at tested lags -> variance looks locally "
              "stable at this sampling resolution; plain ARMA may be sufficient later.")

    # --- Diagnostic plots: squared residual + its ACF ---
    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    axes[0].plot(resid.index, resid.values ** 2, lw=0.6, color="firebrick")
    axes[0].set_ylabel("squared residual")
    axes[0].set_title(f"{args.buoy} — squared residual (local-variance proxy)")
    plot_acf(resid.values ** 2, lags=min(96, len(resid) // 2 - 1), ax=axes[1])
    axes[1].set_title("ACF of squared residual (persistent decay = ARCH signature)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_arch_diagnostic.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
