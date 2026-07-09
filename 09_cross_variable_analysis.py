"""
Cross-variable relationships: Hs, peak period, mean period, and direction
together instead of Hs in isolation.

- Correlation matrix (Pearson). VMDR is a DIRECTION (circular quantity) -
  correlating it linearly against Hs/Tp is not meaningful, so it's
  decomposed into sin/cos components first, which is the standard way to
  bring a circular variable into a linear-correlation framework.
- Lagged cross-correlation between Hs and Tp: does peak period shift
  before Hs rises (swell arriving ahead of a storm) or together with it?
- PCA on the continuous wave variables, to see how many independent
  "directions of variation" the sea state actually has.

Usage:
    python 09_cross_variable_analysis.py --buoy WesthinderBuoy --data-dir data
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import load_buoy_dataframe, default_paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data", type=Path)
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--max-lag-samples", default=48, type=int,
                         help="max lag (in samples) to test for Hs/Tp cross-correlation")
    args = parser.parse_args()

    nc_path = args.data_dir / f"{args.buoy}.nc"
    df = load_buoy_dataframe(nc_path)
    print(f"Variables available for {args.buoy}: {list(df.columns)}")

    out_dir = default_paths("09_cross_variable_analysis")

    if "VMDR" in df.columns:
        theta = np.deg2rad(df["VMDR"])
        df["VMDR_sin"] = np.sin(theta)
        df["VMDR_cos"] = np.cos(theta)

    corr_cols = [c for c in df.columns if c != "VMDR"]  # drop raw angle, keep sin/cos
    df_corr_input = df[corr_cols].dropna()

    if len(df_corr_input) < 10 or df_corr_input.shape[1] < 2:
        print("Not enough overlapping variables/samples for correlation analysis. Stopping.")
        return

    corr = df_corr_input.corr(method="pearson")
    corr.to_csv(out_dir / f"{args.buoy}_correlation_matrix.csv")
    print("\nCorrelation matrix:")
    print(corr.round(3))

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(corr.values, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(corr.columns)))
    ax.set_yticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticklabels(corr.columns)
    for i in range(len(corr)):
        for j in range(len(corr)):
            ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, label="Pearson r")
    ax.set_title(f"{args.buoy} — cross-variable correlation")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_correlation_heatmap.png", dpi=150)
    plt.close(fig)

    # --- Lagged cross-correlation: does Tp lead Hs? ---
    if "VHM0" in df.columns and "VTPK" in df.columns:
        pair = df[["VHM0", "VTPK"]].dropna()
        hs = (pair["VHM0"] - pair["VHM0"].mean()).values
        tp = (pair["VTPK"] - pair["VTPK"].mean()).values

        lags = range(-args.max_lag_samples, args.max_lag_samples + 1)
        ccf_vals = []
        for lag in lags:
            if lag < 0:
                a, b = hs[-lag:], tp[:lag] if lag != 0 else tp
            elif lag > 0:
                a, b = hs[:-lag], tp[lag:]
            else:
                a, b = hs, tp
            if len(a) < 10:
                ccf_vals.append(np.nan)
                continue
            ccf_vals.append(np.corrcoef(a, b)[0, 1])

        ccf_vals = np.array(ccf_vals)
        best_lag = list(lags)[int(np.nanargmax(np.abs(ccf_vals)))]
        best_corr = ccf_vals[int(np.nanargmax(np.abs(ccf_vals)))]

        print(f"\nHs/Tp cross-correlation: strongest |r|={abs(best_corr):.3f} at lag={best_lag} samples")
        if best_lag < 0:
            print("  (negative lag: Tp changes BEFORE Hs — consistent with swell/period "
                  "shifting ahead of a storm's height rise)")
        elif best_lag > 0:
            print("  (positive lag: Tp changes AFTER Hs)")
        else:
            print("  (zero lag: Hs and Tp move together, no lead/lag)")

        fig, ax = plt.subplots(figsize=(9, 4))
        ax.stem(list(lags), ccf_vals)
        ax.axvline(0, color="gray", ls="--")
        ax.set_xlabel("lag (samples, Tp relative to Hs)")
        ax.set_ylabel("cross-correlation")
        ax.set_title(f"{args.buoy} — Hs/Tp cross-correlation function")
        fig.tight_layout()
        fig.savefig(out_dir / f"{args.buoy}_hs_tp_ccf.png", dpi=150)
        plt.close(fig)

        pd.DataFrame({"lag": list(lags), "ccf": ccf_vals}).to_csv(
            out_dir / f"{args.buoy}_hs_tp_ccf.csv", index=False)
    else:
        best_lag, best_corr = None, None

    # --- PCA on the continuous variables ---
    pca_cols = [c for c in ["VHM0", "VTPK", "VTM02", "VMDR_sin", "VMDR_cos"] if c in df.columns]
    pca_input = df[pca_cols].dropna()
    pca_summary = {}
    if len(pca_cols) >= 2 and len(pca_input) > 20:
        X = (pca_input - pca_input.mean()) / pca_input.std()
        cov = np.cov(X.values, rowvar=False)
        eigvals, eigvecs = np.linalg.eigh(cov)
        order = np.argsort(eigvals)[::-1]
        eigvals, eigvecs = eigvals[order], eigvecs[:, order]
        explained = eigvals / eigvals.sum()

        print("\nPCA explained variance ratio:", np.round(explained, 3).tolist())
        loadings = pd.DataFrame(eigvecs, index=pca_cols,
                                 columns=[f"PC{i+1}" for i in range(len(pca_cols))])
        loadings.to_csv(out_dir / f"{args.buoy}_pca_loadings.csv")
        print("PC1 loadings (dominant direction of joint variation):")
        print(loadings["PC1"].round(3))

        pca_summary = {
            "explained_variance_ratio": explained.tolist(),
            "pc1_loadings": loadings["PC1"].to_dict(),
        }
    else:
        print("\nNot enough variables/samples for PCA - skipping.")

    import json
    with open(out_dir / f"{args.buoy}_cross_var_summary.json", "w") as f:
        json.dump({
            "variables_available": list(df.columns),
            "hs_tp_best_lag_samples": best_lag,
            "hs_tp_best_corr": None if best_corr is None else float(best_corr),
            "pca": pca_summary,
        }, f, indent=2, default=str)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
