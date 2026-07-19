"""
Stage 24 - HMM regime identification: extends Stage 10's static GMM
regime classification with TRANSITION PROBABILITY structure - given
you're in a storm regime right now, what's the probability you're
still in one at the next timestep? Stage 10 and 11b can't answer this:
Stage 10 gives regime fractions with no memory, 11b gives one aggregate
persistence timescale for Hs itself, not per-regime dwell times.

Uses hmmlearn's `lengths` parameter to fit on contiguous segments
separately - critical, not optional, for the same reason nearly every
lag-based stage in this pipeline needed gap-awareness. Verified
empirically before trusting it (not just documentation): on a 50-segment
synthetic test with a deliberate state-mismatch at every segment
boundary, omitting `lengths` produced spurious cross-state transition
probabilities of ~0.025 (a real, material contamination); with
`lengths`, the same spurious transitions were numerically zero
(~1e-134). The effect scales with segment count, so it matters more on
a heavily-fragmented record (Westhinder: 1905 segments) than a lightly
fragmented one - exactly backwards from what an unaware implementation
would assume "probably fine for a small gap count."

Reports expected dwell time per regime (from the transition matrix
diagonal: E[dwell] = 1/(1-p_ii)), directly comparable to 11b's
aggregate Hs persistence timescale but regime-specific.

Usage:
    python 24_regime_hmm.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hmmlearn.hmm import GaussianHMM

from utils import default_paths, all_contiguous_segments


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--n-regimes", default=4, type=int)
    parser.add_argument("--min-segment-length", default=50, type=int,
                         help="segments shorter than this contribute little to "
                              "transition-matrix estimation and are excluded")
    parser.add_argument("--n-iter", default=200, type=int)
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"

    hs_full = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var]
    dt_hours = 0.5
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            dt_hours = json.load(f).get("sampling_interval_hours", 0.5)

    print(f"--- {args.buoy} / {args.var} HMM regime identification ---")

    segments = all_contiguous_segments(hs_full, min_length=args.min_segment_length)
    total_valid = int(hs_full.notna().sum())
    coverage = sum(len(s) for s in segments)
    print(f"{len(segments)} segments >= {args.min_segment_length} samples, "
          f"covering {coverage}/{total_valid} valid samples "
          f"({100*coverage/total_valid:.1f}%)")

    if len(segments) == 0:
        print("No segments long enough - stopping.")
        return

    X = np.concatenate([s.values for s in segments]).reshape(-1, 1)
    lengths = [len(s) for s in segments]

    print(f"Fitting GaussianHMM(n_components={args.n_regimes}) on {len(X)} samples "
          f"across {len(lengths)} segments (lengths= respected, no cross-gap transitions)...")
    model = GaussianHMM(n_components=args.n_regimes, n_iter=args.n_iter,
                         random_state=0, covariance_type="diag")
    model.fit(X, lengths=lengths)

    # Order regimes by mean Hs (calm -> storm), matching Stage 10's convention
    order = np.argsort(model.means_.flatten())
    means_sorted = model.means_.flatten()[order]
    # Reorder the transition matrix rows/cols to match
    transmat_sorted = model.transmat_[np.ix_(order, order)]

    print(f"\nRegime means (calm -> storm): {np.round(means_sorted, 3)}")
    print(f"\nTransition matrix (rows=from, cols=to, calm..storm order):")
    print(pd.DataFrame(transmat_sorted,
                        index=[f"from_{i}" for i in range(args.n_regimes)],
                        columns=[f"to_{i}" for i in range(args.n_regimes)]).round(4).to_string())

    dwell_samples = 1.0 / (1.0 - np.diag(transmat_sorted))
    dwell_hours = dwell_samples * dt_hours
    print(f"\nExpected dwell time per regime:")
    for i in range(args.n_regimes):
        print(f"  regime {i} (mean Hs={means_sorted[i]:.2f}m): "
              f"{dwell_samples[i]:.1f} samples = {dwell_hours[i]:.1f}h")

    out_dir = default_paths("24_regime_hmm")
    with open(out_dir / f"{args.buoy}_{args.var}_hmm_summary.json", "w") as f:
        json.dump({
            "n_regimes": args.n_regimes,
            "n_segments_used": len(segments),
            "pct_valid_data_used": round(100 * coverage / total_valid, 1),
            "regime_means": means_sorted.tolist(),
            "transition_matrix": transmat_sorted.tolist(),
            "dwell_time_hours": dwell_hours.tolist(),
        }, f, indent=2)

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(transmat_sorted, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(args.n_regimes))
    ax.set_yticks(range(args.n_regimes))
    ax.set_xticklabels([f"R{i}" for i in range(args.n_regimes)])
    ax.set_yticklabels([f"R{i}" for i in range(args.n_regimes)])
    ax.set_xlabel("to regime")
    ax.set_ylabel("from regime")
    for i in range(args.n_regimes):
        for j in range(args.n_regimes):
            ax.text(j, i, f"{transmat_sorted[i, j]:.3f}", ha="center", va="center",
                     color="white" if transmat_sorted[i, j] < 0.5 else "black")
    fig.colorbar(im, label="transition probability")
    ax.set_title(f"{args.buoy} — HMM regime transition matrix\n"
                 f"(0=calmest .. {args.n_regimes-1}=stormiest)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_transition_matrix.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
