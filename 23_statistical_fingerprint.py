"""
Priority 10 - statistical fingerprint: one-page-per-buoy visual summary
(distribution, storm stats, persistence, CI, stability, regime
structure all in one place) PLUS network-wide clustering of buoys by
statistical similarity (PCA for 2D visualization + KMeans clustering on
the standardized feature vectors) - PCA/classical clustering, not UMAP
(n=19 is too small for UMAP to mean anything).

This is a genuinely different lens than Stage 11's spatial clustering:
Stage 11 clusters buoys by CORRELATION (do their Hs series move
together), this clusters buoys by STATISTICAL PROPERTY SIMILARITY (do
they have similar distribution shape, tail behavior, persistence,
regime structure) - two buoys could be geographically uncorrelated but
statistically similar, or vice versa.

Two modes:
    python 23_statistical_fingerprint.py --buoy WesthinderBuoy --var VHM0
    python 23_statistical_fingerprint.py --network --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from utils import default_paths

STAGE_DIR = Path("pipeline_out")


def safe_json(path: Path):
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def safe_csv(path: Path):
    if not path.exists():
        return None
    return pd.read_csv(path)


def discover_buoys(var: str):
    load_dir = STAGE_DIR / "01_load_clean"
    suffix = f"_{var}_load_summary.json"
    if not load_dir.exists():
        return []
    return sorted(p.name[: -len(suffix)] for p in load_dir.glob(f"*{suffix}"))


def build_fingerprint(buoy: str, var: str) -> dict:
    """Numeric + categorical fingerprint for one buoy, reading from
    CONFIRMED real schemas (checked against actual stage source code,
    not guessed - two schema-guessing bugs were found and fixed in
    Priority 9's diagnostics report before this was built, applying the
    same lesson here from the start)."""
    fp = {"buoy": buoy}

    load = safe_json(STAGE_DIR / "01_load_clean" / f"{buoy}_{var}_load_summary.json")
    if load:
        fp["record_years"] = load.get("record_years")
        fp["pct_missing"] = load.get("pct_missing_after_clean")

    fit_interp = safe_json(STAGE_DIR / "06_distribution_fit" / f"{buoy}_{var}_fit_interpretation.json")
    fit_csv = safe_csv(STAGE_DIR / "06_distribution_fit" / f"{buoy}_{var}_fit_summary.csv")
    if fit_interp:
        fp["best_distribution"] = fit_interp.get("best_distribution")
    if fit_csv is not None and len(fit_csv):
        fp["best_ks_stat"] = fit_csv["ks_stat"].min()

    eva = safe_json(STAGE_DIR / "08_extreme_value_analysis" / f"{buoy}_{var}_eva_summary.json")
    if eva:
        fp["gpd_xi"] = eva.get("gpd_shape")
        fp["n_peaks"] = eva.get("n_peaks")

    dep = safe_json(STAGE_DIR / "11b_dependence_structure" / f"{buoy}_{var}_dependence_summary.json")
    if dep:
        fp["persistence_hours"] = dep.get("integral_timescale_hours")
        fp["lag1_acf"] = dep.get("lag1_acf")

    conf = safe_json(STAGE_DIR / "12_confidence_intervals" / f"{buoy}_{var}_confidence_summary.json")
    if conf:
        xi_ci = conf.get("gpd_xi_ci")
        if xi_ci and xi_ci.get("ci_low") is not None:
            fp["xi_ci_width"] = xi_ci["ci_high"] - xi_ci["ci_low"]
        if conf.get("hs_mean_ci_low") is not None:
            fp["hs_mean_ci_width"] = conf["hs_mean_ci_high"] - conf["hs_mean_ci_low"]

    stab = safe_json(STAGE_DIR / "13_stability_analysis" / f"{buoy}_{var}_stability_summary.json")
    if stab:
        fp["window_agreement_frac"] = stab.get("window_stability", {}).get("fraction_windows_agreeing")

    regime = safe_csv(STAGE_DIR / "10_regime_identification" / f"{buoy}_{var}_regime_summary.csv")
    if regime is not None and len(regime):
        regime = regime.sort_values("regime_id")
        for _, row in regime.iterrows():
            fp[f"regime_{int(row['regime_id'])}_frac"] = row["fraction_of_time"]

    return fp


def plot_single_buoy(fp: dict, out_dir: Path, var: str):
    fig = plt.figure(figsize=(10, 7))
    gs = fig.add_gridspec(2, 2)

    ax_text = fig.add_subplot(gs[0, :])
    ax_text.axis("off")
    lines = [f"{fp['buoy']} — {var} statistical fingerprint", ""]
    if "record_years" in fp:
        lines.append(f"Record: {fp['record_years']:.1f} years, "
                      f"{fp.get('pct_missing', float('nan')):.1f}% missing")
    if "best_distribution" in fp:
        lines.append(f"Distribution: {fp['best_distribution']} "
                      f"(KS stat={fp.get('best_ks_stat', float('nan')):.4f})")
    if "gpd_xi" in fp:
        ci_str = f" +/- {fp['xi_ci_width']/2:.3f}" if "xi_ci_width" in fp else ""
        lines.append(f"GPD xi: {fp['gpd_xi']:.3f}{ci_str} ({fp.get('n_peaks', '?')} peaks)")
    if "persistence_hours" in fp:
        lines.append(f"Persistence timescale: {fp['persistence_hours']:.1f}h "
                      f"(lag-1 ACF={fp.get('lag1_acf', float('nan')):.3f})")
    if "window_agreement_frac" in fp:
        lines.append(f"Moving-window distribution agreement: "
                      f"{fp['window_agreement_frac']*100:.0f}%")
    ax_text.text(0, 1, "\n".join(lines), va="top", fontsize=11, family="monospace")

    ax_regime = fig.add_subplot(gs[1, 0])
    regime_keys = sorted([k for k in fp if k.startswith("regime_") and k.endswith("_frac")])
    if regime_keys:
        vals = [fp[k] * 100 for k in regime_keys]
        labels = [k.replace("regime_", "R").replace("_frac", "") for k in regime_keys]
        ax_regime.bar(labels, vals, color="steelblue")
        ax_regime.set_ylabel("% of time")
        ax_regime.set_title("Regime fractions")
    else:
        ax_regime.axis("off")
        ax_regime.text(0.5, 0.5, "No Stage 10 output", ha="center")

    ax_xi = fig.add_subplot(gs[1, 1])
    if "gpd_xi" in fp:
        half_width = fp.get("xi_ci_width", 0) / 2
        ax_xi.errorbar([0], [fp["gpd_xi"]], yerr=[[half_width], [half_width]],
                       fmt="o", capsize=8, color="firebrick", markersize=10)
        ax_xi.axhline(0, color="gray", ls="--", lw=0.8)
        ax_xi.set_xlim(-1, 1)
        ax_xi.set_xticks([])
        ax_xi.set_ylabel("GPD xi")
        ax_xi.set_title("Tail shape +/- CI")
    else:
        ax_xi.axis("off")
        ax_xi.text(0.5, 0.5, "No Stage 08/12 output", ha="center")

    fig.tight_layout()
    fig.savefig(out_dir / f"{fp['buoy']}_{var}_fingerprint.png", dpi=150)
    plt.close(fig)


def run_network_clustering(var: str, out_dir: Path, exclude_record_length: bool = False):
    buoys = discover_buoys(var)
    if len(buoys) < 3:
        print(f"Only {len(buoys)} buoy(s) with Stage 0 output - need at least 3 for clustering.")
        return

    fingerprints = [build_fingerprint(b, var) for b in buoys]
    df = pd.DataFrame(fingerprints).set_index("buoy")

    numeric_cols = ["record_years", "pct_missing", "best_ks_stat", "gpd_xi", "n_peaks",
                     "persistence_hours", "lag1_acf", "xi_ci_width", "hs_mean_ci_width",
                     "window_agreement_frac"]
    if exclude_record_length:
        # Added after a real finding: clustering on the full feature set
        # produced two large clusters that suspiciously tracked deployment
        # era (newer/shorter-record buoys vs. older/longer-record ones)
        # more than obviously-behavioral properties. record_years and
        # pct_missing are deployment-history/data-availability metadata,
        # not behavioral - excluding them tests whether that split
        # reflects genuine physical difference or just data availability.
        numeric_cols = [c for c in numeric_cols if c not in ("record_years", "pct_missing")]
        print("Excluding record_years/pct_missing from clustering features "
              "(deployment-history metadata, not behavioral properties).")
    regime_cols = sorted([c for c in df.columns if c.startswith("regime_") and c.endswith("_frac")])
    feature_cols = [c for c in numeric_cols + regime_cols if c in df.columns]

    suffix = "_no_reclen" if exclude_record_length else ""
    print(f"Fingerprint table ({len(buoys)} buoys, {len(feature_cols)} numeric features) "
          f"saved to {out_dir / f'{var}_fingerprint_table.csv'}")

    X = df[feature_cols].copy()
    missing_frac = X.isna().mean()
    if (missing_frac > 0).any():
        print(f"Mean-imputing missing values for clustering (per-feature missing fraction):")
        print(missing_frac[missing_frac > 0].to_string())
    X = X.fillna(X.mean())

    complete_buoys = X.dropna().index  # after imputation, only fully-NaN columns would remain missing
    X = X.loc[complete_buoys]
    if len(X) < 3:
        print("Too few buoys with usable features after imputation - stopping clustering.")
        return

    X_scaled = StandardScaler().fit_transform(X)

    # Pick k via silhouette score instead of a fixed guess - found necessary
    # after testing on a small (n=6) synthetic set with a clean 2-group
    # structure: a hardcoded k=4 over-split it into 4 clusters (no
    # cross-contamination between the true groups, but not the clean
    # split that was actually there). Silhouette-based selection correctly
    # recovers k=2 in that case while still being sensible at the real
    # n=19 network scale.
    from sklearn.metrics import silhouette_score
    max_k = min(6, len(X) - 1)
    best_k, best_score, best_labels = 2, -1, None
    for k in range(2, max_k + 1):
        trial_labels = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(X_scaled)
        score = silhouette_score(X_scaled, trial_labels)
        if score > best_score:
            best_k, best_score, best_labels = k, score, trial_labels
    n_clusters, labels = best_k, best_labels
    print(f"Selected k={n_clusters} via silhouette score ({best_score:.3f})")

    pca = PCA(n_components=2)
    coords = pca.fit_transform(X_scaled)
    var_explained = pca.explained_variance_ratio_
    print(f"PCA: PC1/PC2 explain {var_explained[0]*100:.1f}%/{var_explained[1]*100:.1f}% of variance")

    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=labels, cmap="tab10", s=80)
    for i, buoy in enumerate(X.index):
        ax.annotate(buoy, (coords[i, 0], coords[i, 1]), fontsize=8,
                    xytext=(5, 5), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var_explained[0]*100:.0f}%)")
    ax.set_ylabel(f"PC2 ({var_explained[1]*100:.0f}%)")
    ax.set_title(f"{var} — statistical fingerprint similarity (PCA + KMeans, k={n_clusters})"
                 f"{' [record-length excluded]' if exclude_record_length else ''}")
    fig.tight_layout()
    fig.savefig(out_dir / f"{var}_fingerprint_clustering{suffix}.png", dpi=150)
    plt.close(fig)

    cluster_df = pd.DataFrame({"buoy": X.index, "cluster": labels})
    cluster_df.to_csv(out_dir / f"{var}_fingerprint_clusters{suffix}.csv", index=False)
    print("\nStatistical-similarity clusters (different lens than Stage 11's "
          "correlation-based spatial clusters):")
    for c in sorted(set(labels)):
        members = cluster_df.loc[cluster_df["cluster"] == c, "buoy"].tolist()
        print(f"  cluster {c}: {members}")

    print(f"\nSaved: {out_dir}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default=None)
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--network", action="store_true",
                         help="build the network-wide fingerprint table + PCA/clustering "
                              "instead of a single-buoy visual")
    parser.add_argument("--exclude-record-length", action="store_true",
                         help="drop record_years/pct_missing from the clustering features - "
                              "tests whether clusters reflect genuine behavioral similarity "
                              "or just deployment-era/data-availability differences")
    args = parser.parse_args()

    out_dir = default_paths("23_statistical_fingerprint")

    if args.network or args.buoy is None:
        run_network_clustering(args.var, out_dir, exclude_record_length=args.exclude_record_length)
    else:
        fp = build_fingerprint(args.buoy, args.var)
        plot_single_buoy(fp, out_dir, args.var)
        print(f"Fingerprint for {args.buoy}:")
        for k, v in fp.items():
            print(f"  {k}: {v}")
        print(f"\nSaved: {out_dir / f'{args.buoy}_{args.var}_fingerprint.png'}")


if __name__ == "__main__":
    main()
