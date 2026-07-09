"""
Regime identification: classify each timestep into a sea-state regime
(calm / moderate / energetic / storm) via Gaussian Mixture Model
clustering on Hs (optionally + Tp), rather than eyeballing it from a
rolling-variance plot.

This is NOT forecasting - it labels what already happened, using only
the observed record. Answers: "what fraction of time does this buoy
spend in each regime?"

Usage:
    python 10_regime_identification.py --buoy WesthinderBuoy --var VHM0 \
        --n-regimes 4 --include-period
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

from utils import default_paths, load_buoy_dataframe

REGIME_LABELS_4 = ["calm", "moderate", "energetic", "storm"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data", type=Path)
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--n-regimes", default=4, type=int)
    parser.add_argument("--include-period", action="store_true",
                         help="cluster on [Hs, Tp] jointly instead of Hs alone")
    parser.add_argument("--random-state", default=0, type=int)
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    s = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var]

    out_dir = default_paths("10_regime_identification")

    use_period = args.include_period
    if use_period:
        nc_path = args.data_dir / f"{args.buoy}.nc"
        df_raw = load_buoy_dataframe(nc_path, varnames=("VTPK",))
        if "VTPK" not in df_raw.columns:
            print(f"WARNING: --include-period was requested but {args.buoy} has no VTPK "
                  f"variable (this buoy likely only reports Hs) - falling back to Hs-only "
                  f"clustering instead of crashing.")
            use_period = False

    if use_period:
        features = pd.concat([s.rename("VHM0"), df_raw["VTPK"]], axis=1).dropna()
        feature_cols = ["VHM0", "VTPK"]
    else:
        features = s.dropna().to_frame(name=args.var)
        feature_cols = [args.var]

    X = features[feature_cols].values
    X_std = (X - X.mean(axis=0)) / X.std(axis=0)

    gmm = GaussianMixture(n_components=args.n_regimes, random_state=args.random_state,
                           n_init=5)
    labels_raw = gmm.fit_predict(X_std)

    # Relabel clusters by mean Hs, ascending, so "0" is always calmest
    mean_hs_per_cluster = pd.Series(X[:, 0], index=features.index).groupby(labels_raw).mean()
    order = mean_hs_per_cluster.sort_values().index.tolist()
    relabel_map = {old: new for new, old in enumerate(order)}
    labels = np.array([relabel_map[l] for l in labels_raw])

    if args.n_regimes == 4:
        names = REGIME_LABELS_4
    else:
        names = [f"regime_{i}" for i in range(args.n_regimes)]

    regime_series = pd.Series(labels, index=features.index, name="regime")
    regime_series.to_csv(out_dir / f"{args.buoy}_{args.var}_regime_labels.csv", header=["regime"])

    fractions = regime_series.value_counts(normalize=True).sort_index()
    print(f"--- {args.buoy} regime identification ({'Hs+Tp' if use_period else 'Hs only'}, "
          f"k={args.n_regimes}) ---")
    summary_rows = []
    for i in range(args.n_regimes):
        name = names[i] if i < len(names) else f"regime_{i}"
        frac = fractions.get(i, 0.0)
        mean_hs = features.loc[regime_series == i, feature_cols[0]].mean()
        print(f"  {name:12s}: {frac*100:5.1f}% of time, mean {feature_cols[0]}={mean_hs:.2f}")
        summary_rows.append({"regime_id": i, "regime_name": name,
                              "fraction_of_time": frac, "mean_hs": mean_hs})

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(out_dir / f"{args.buoy}_{args.var}_regime_summary.csv", index=False)

    # --- Plot: series colored by regime ---
    fig, ax = plt.subplots(figsize=(12, 5))
    cmap = plt.cm.viridis
    colors = [cmap(i / max(1, args.n_regimes - 1)) for i in range(args.n_regimes)]
    for i in range(args.n_regimes):
        mask = regime_series == i
        name = names[i] if i < len(names) else f"regime_{i}"
        ax.scatter(features.index[mask], features.loc[mask, feature_cols[0]],
                   s=4, color=colors[i], label=name)
    ax.set_ylabel(feature_cols[0])
    ax.set_title(f"{args.buoy} — sea-state regime timeline")
    ax.legend(markerscale=3, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_regime_timeline.png", dpi=150)
    plt.close(fig)

    if use_period:
        fig, ax = plt.subplots(figsize=(7, 6))
        for i in range(args.n_regimes):
            mask = regime_series == i
            name = names[i] if i < len(names) else f"regime_{i}"
            ax.scatter(features.loc[mask, "VHM0"], features.loc[mask, "VTPK"],
                       s=6, color=colors[i], label=name, alpha=0.6)
        ax.set_xlabel("Hs (m)")
        ax.set_ylabel("Tp (s)")
        ax.set_title(f"{args.buoy} — regimes in Hs/Tp space")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / f"{args.buoy}_{args.var}_regime_scatter.png", dpi=150)
        plt.close(fig)

    import json
    with open(out_dir / f"{args.buoy}_{args.var}_regime_meta.json", "w") as f:
        json.dump({
            "n_regimes": args.n_regimes,
            "included_period": use_period,
            "fractions": {names[i] if i < len(names) else f"regime_{i}": float(fractions.get(i, 0.0))
                          for i in range(args.n_regimes)},
        }, f, indent=2)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
