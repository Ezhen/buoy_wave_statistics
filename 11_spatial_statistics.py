"""
Spatial statistics: the payoff of having 19 buoys instead of one.

- Pairwise correlation of Hs across all buoys, on their overlapping time
  window.
- Correlation vs. inter-buoy distance (do nearby buoys move together more
  than distant ones? - the basic spatial-coherence question).
- Hierarchical clustering of buoys by (1 - correlation) as a distance
  metric, to see which stations group together.

This runs AFTER Stage 0 (and needs the raw .nc files for lat/lon).

Usage:
    python 11_spatial_statistics.py --data-dir data --var VHM0
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import xarray as xr

from utils import default_paths, haversine_km, get_scalar_latlon


def discover_buoys(var: str):
    load_dir = Path("pipeline_out/01_load_clean")
    suffix = f"_{var}_clean.csv"
    return sorted(p.name[: -len(suffix)] for p in load_dir.glob(f"*{suffix}"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data", type=Path)
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--n-clusters", default=4, type=int)
    args = parser.parse_args()

    buoys = discover_buoys(args.var)
    if len(buoys) < 3:
        print(f"Only {len(buoys)} buoy(s) with Stage 0 output - need at least 3 "
              f"for spatial statistics. Run Stage 0 for more buoys first.")
        return

    out_dir = default_paths("11_spatial_statistics")

    # --- Load all series, aligned ---
    series = {}
    latlon = {}
    for buoy in buoys:
        s = pd.read_csv(Path("pipeline_out/01_load_clean") / f"{buoy}_{args.var}_clean.csv",
                         index_col=0, parse_dates=True)[args.var]
        series[buoy] = s

        nc_path = args.data_dir / f"{buoy}.nc"
        if nc_path.exists():
            with xr.open_dataset(nc_path) as ds:
                latlon[buoy] = get_scalar_latlon(ds)

    df = pd.DataFrame(series)  # outer join on time index
    print(f"Loaded {len(buoys)} buoys, {df.shape[0]} timesteps (union of all time indices)")
    overlap = df.dropna()
    print(f"{len(overlap)} timesteps with ALL buoys present simultaneously")

    if len(overlap) < 20:
        print("WARNING: very little full overlap across all buoys - correlations below "
              "are computed pairwise, using each pair's own overlap, not the all-buoy overlap.")

    # --- Pairwise correlation (each pair uses its own overlapping timestamps) ---
    corr = df.corr(method="pearson", min_periods=20)
    corr.to_csv(out_dir / f"{args.var}_pairwise_correlation.csv")

    fig, ax = plt.subplots(figsize=(10, 9))
    im = ax.imshow(corr.values, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(corr.columns)))
    ax.set_yticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=90, fontsize=7)
    ax.set_yticklabels(corr.columns, fontsize=7)
    fig.colorbar(im, label="Pearson r")
    ax.set_title(f"{args.var} — pairwise correlation across BCZ buoys")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.var}_correlation_heatmap.png", dpi=150)
    plt.close(fig)

    # --- Correlation vs distance ---
    rows = []
    for i, b1 in enumerate(buoys):
        for b2 in buoys[i + 1:]:
            if b1 not in latlon or b2 not in latlon:
                continue
            r = corr.loc[b1, b2]
            if pd.isna(r):
                continue
            dist = haversine_km(*latlon[b1], *latlon[b2])
            rows.append({"buoy1": b1, "buoy2": b2, "distance_km": dist, "correlation": r})

    dist_corr_df = pd.DataFrame(rows)
    dist_corr_df.to_csv(out_dir / f"{args.var}_correlation_vs_distance.csv", index=False)

    if len(dist_corr_df):
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(dist_corr_df["distance_km"], dist_corr_df["correlation"], alpha=0.6)
        if len(dist_corr_df) > 3:
            fit = np.polyfit(dist_corr_df["distance_km"], dist_corr_df["correlation"], 1)
            xs = np.linspace(0, dist_corr_df["distance_km"].max(), 50)
            ax.plot(xs, np.polyval(fit, xs), color="firebrick", lw=2,
                    label=f"linear fit (slope={fit[0]:.4f}/km)")
            ax.legend()
        ax.set_xlabel("Inter-buoy distance (km)")
        ax.set_ylabel(f"{args.var} correlation")
        ax.set_title("Spatial coherence: correlation vs. distance")
        fig.tight_layout()
        fig.savefig(out_dir / f"{args.var}_correlation_vs_distance.png", dpi=150)
        plt.close(fig)

        spearman_r = dist_corr_df[["distance_km", "correlation"]].corr(method="spearman").iloc[0, 1]
        print(f"\nCorrelation vs. distance: Spearman r={spearman_r:.3f} "
              f"({'decreasing' if spearman_r < 0 else 'not clearly decreasing'} with distance)")

    # --- Hierarchical clustering of buoys by (1 - correlation) ---
    valid_buoys = [b for b in buoys if corr[b].notna().sum() > 1]
    corr_valid = corr.loc[valid_buoys, valid_buoys].fillna(0)
    dist_matrix = 1 - corr_valid.values
    np.fill_diagonal(dist_matrix, 0)
    dist_matrix = (dist_matrix + dist_matrix.T) / 2  # enforce symmetry against float noise
    condensed = squareform(dist_matrix, checks=False)
    Z = linkage(condensed, method="average")

    fig, ax = plt.subplots(figsize=(10, 6))
    dendrogram(Z, labels=valid_buoys, ax=ax, leaf_rotation=90)
    ax.set_title(f"{args.var} — buoy clustering by (1 - correlation)")
    ax.set_ylabel("distance (1 - r)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.var}_dendrogram.png", dpi=150)
    plt.close(fig)

    cluster_labels = fcluster(Z, t=args.n_clusters, criterion="maxclust")
    cluster_df = pd.DataFrame({"buoy": valid_buoys, "cluster": cluster_labels})
    cluster_df.to_csv(out_dir / f"{args.var}_buoy_clusters.csv", index=False)
    print(f"\nBuoy clusters (k={args.n_clusters}):")
    for c in sorted(set(cluster_labels)):
        members = cluster_df.loc[cluster_df["cluster"] == c, "buoy"].tolist()
        print(f"  cluster {c}: {members}")

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
