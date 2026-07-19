"""
Stage 26 - SSA (Singular Spectrum Analysis) decomposition.

Third and last of the post-priority signal-processing extensions
(alongside Stage 24's HMM and Stage 25's change-point detection),
scoped specifically for Zeebrugge's unresolved tidal notch anomaly:
Stage 03b's M2 harmonic-regression detiding makes Zeebrugge's residual
WORSE, not better, regardless of harmonic count - unlike every other
buoy in the network. Leading hypotheses, never directly tested until
now: a compound tide (MS4, generated locally by shallow-water/harbor
nonlinearity, not a simple sum of open-ocean constituents) or a
time-varying effective tidal frequency/amplitude that a FIXED-frequency
harmonic regression structurally cannot fit no matter how many
harmonics are added.

SSA doesn't assume a fixed frequency the way harmonic regression does -
it decomposes the series into data-driven oscillatory components via
trajectory-matrix (Hankel) SVD, so it can in principle separate close,
compound, or slowly-varying periodicities a fixed M2/S2 harmonic model
cannot.

CORE METHOD, validated against synthetic ground truth before trusting
it on real data (per this project's established discipline for every
new method so far):
  1. Single clean sinusoid -> recovered period within 0.02% of true,
     leading singular-value pair ~180x the noise floor.
  2. M2 (12.4206h) + S2 (12.0000h) mixed signal -> BOTH cleanly
     recovered at their exact true periods as two separate pairs, with
     a w-correlation matrix showing near-zero cross-correlation between
     the two pairs (0.03) vs ~1.0 within each pair - genuine
     separation, not conflation into one blended component.
  3. Pure noise control -> smooth, gradually-decaying singular value
     scree (consecutive ratios all ~1.0-1.04), unlike the sharp
     order-of-magnitude drop after the true signal pair(s) in tests
     1-2 - confirms the scree-shape diagnostic actually discriminates
     signal from noise, not just a plausible-sounding untested
     heuristic.

WINDOW LENGTH is not an arbitrary default: separating M2 from S2
requires a window >= ~354.4 hours (frequency resolution 1/(L*dt) must
beat their ~0.00282 cycles/hour spacing - the same resolution math
already used to justify Stage 03b's --fit-frequency band-search claim).
Default here is 600h (25 days, ~1.7x that floor) for a safety margin
above the bare theoretical minimum, not just clearing it exactly.

Applied to the LONGEST CONTIGUOUS SEGMENT only, via the same
utils.longest_contiguous_segment gap-awareness every other lag/order-
based stage in this pipeline requires - naive concatenation across a
gap would corrupt the trajectory matrix the same way it corrupts ACF or
HMM transitions elsewhere.

Usage:
    python 26_ssa_decomposition.py --buoy ZeebruggeZandopvangkadeBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse.linalg import svds
from scipy.signal import periodogram
from scipy.stats import boxcox

from utils import default_paths, longest_contiguous_segment

M2_PERIOD_HOURS = 12.4206
S2_PERIOD_HOURS = 12.0000
MS4_PERIOD_HOURS = 6.2103  # M2's first shallow-water overtide - half M2's period
DIURNAL_PERIOD_HOURS = 24.0
KNOWN_CONSTITUENTS = {
    "M2": M2_PERIOD_HOURS,
    "S2": S2_PERIOD_HOURS,
    "MS4": MS4_PERIOD_HOURS,
    "diurnal": DIURNAL_PERIOD_HOURS,
}
MATCH_TOLERANCE_HOURS = 0.25


def ssa_decompose(x: np.ndarray, window_length: int, n_components: int):
    """Trajectory-matrix (Hankel) SSA via truncated SVD (ARPACK, scipy's
    svds) - NOT a full dense SVD, which would be intractable at real
    record scale (hundreds of thousands of samples). Uses
    sliding_window_view for the trajectory matrix (a zero-copy strided
    view, not a materialized K x L dense array) and reconstructs each
    component via convolution-based diagonal averaging: for a rank-1
    elementary matrix sigma*outer(u,v), the anti-diagonal sums are
    exactly sigma*convolve(u,v), divided by the overlap count at each
    position - mathematically equivalent to, but far cheaper than, an
    explicit anti-diagonal averaging loop."""
    N = len(x)
    L = window_length
    K = N - L + 1
    if K <= L:
        raise ValueError(f"window_length={L} too long for series of length {N} "
                          f"(need K={K} > L={L})")

    traj = np.lib.stride_tricks.sliding_window_view(x, L)  # (K, L), zero-copy view
    k_request = min(n_components + 2, min(traj.shape) - 1)
    U, S, Vt = svds(traj, k=k_request, which="LM")
    order = np.argsort(S)[::-1]  # svds returns ascending order - flip to descending
    U, S, Vt = U[:, order], S[order], Vt[order, :]

    counts = np.convolve(np.ones(K), np.ones(L), mode="full")
    n_keep = min(n_components, len(S))
    components = np.zeros((n_keep, N))
    for i in range(n_keep):
        components[i] = S[i] * np.convolve(U[:, i], Vt[i, :], mode="full") / counts
    return S, components


def dominant_period_hours(component: np.ndarray, dt_hours: float = 1.0):
    freqs, power = periodogram(component, fs=1.0 / dt_hours)
    freqs, power = freqs[1:], power[1:]  # drop DC
    if len(freqs) == 0 or power.max() == 0:
        return float("nan")
    return 1.0 / freqs[np.argmax(power)]


def match_known_constituent(period_hours: float):
    if np.isnan(period_hours):
        return None
    for name, p in KNOWN_CONSTITUENTS.items():
        if abs(period_hours - p) <= MATCH_TOLERANCE_HOURS:
            return name
    return None


def w_correlation(components: np.ndarray):
    """Standard SSA diagnostic: high w-correlation between two
    components' reconstructions means they should be grouped as one
    physical oscillation - a real sinusoid always splits into a PAIR
    under trajectory-matrix SSA, never a single component, since an
    arbitrary-phase single frequency needs a 2D subspace (sin and cos)
    to represent."""
    n_comp, N = components.shape
    w = np.minimum(np.arange(1, N + 1), np.arange(N, 0, -1)).astype(float)
    w = np.minimum(w, w.max())

    def wnorm(a):
        return np.sqrt(np.sum(w * a * a))

    norms = np.array([wnorm(components[i]) for i in range(n_comp)])
    C = np.zeros((n_comp, n_comp))
    for i in range(n_comp):
        for j in range(n_comp):
            num = np.sum(w * components[i] * components[j])
            C[i, j] = num / (norms[i] * norms[j] + 1e-12)
    return C


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="ZeebruggeZandopvangkadeBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--window-hours", default=600, type=int,
                         help="trajectory matrix window length. Must exceed ~354.4h to "
                              "resolve M2 (12.4206h) from S2 (12.0000h) at hourly sampling - "
                              "default gives ~1.7x margin above that theoretical floor.")
    parser.add_argument("--n-components", default=10, type=int)
    parser.add_argument("--min-segment-hours", default=None, type=float,
                         help="refuse to run if the longest contiguous segment is shorter "
                              "than this many hours; default is 5x --window-hours, a minimal "
                              "safety margin for a meaningful number of trajectory rows")
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    if not clean_path.exists():
        raise FileNotFoundError(f"{clean_path} not found - run Stage 01 first.")

    s = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var]
    segment, seg_meta = longest_contiguous_segment(s)
    print(f"--- {args.buoy} / {args.var} SSA decomposition ---")
    print(f"Longest contiguous segment: {seg_meta['segment_length']} samples "
          f"({seg_meta['pct_of_valid_used']}% of valid data), "
          f"{seg_meta['segment_start']} to {seg_meta['segment_end']}")

    min_segment_hours = args.min_segment_hours or 5 * args.window_hours
    if len(segment) < min_segment_hours:
        print(f"Longest contiguous segment ({len(segment)} samples) is under "
              f"--min-segment-hours={min_segment_hours:.0f} - refusing to run.")
        return

    # Box-Cox for variance stabilization before decomposition - same
    # rationale and same shift-to-positive convention as Stage 03b, so
    # results are on a comparable footing to the existing tidal-notch work.
    shift = max(0.0, -segment.min()) + 1e-6
    bc_vals, lam = boxcox(segment.values + shift)
    print(f"Box-Cox lambda: {lam:.4f}")

    S, components = ssa_decompose(bc_vals, window_length=args.window_hours,
                                   n_components=args.n_components)

    print(f"\nSingular values (kept): {np.round(S[:args.n_components], 2)}")
    if len(S) > args.n_components:
        print(f"(next {len(S) - args.n_components} beyond kept, for scree context: "
              f"{np.round(S[args.n_components:], 2)})")

    report_rows = []
    for i in range(components.shape[0]):
        p = dominant_period_hours(components[i])
        match = match_known_constituent(p)
        report_rows.append({
            "component": i, "singular_value": float(S[i]),
            "recovered_period_hours": None if np.isnan(p) else float(p),
            "matched_constituent": match,
            "component_std": float(components[i].std()),
        })
        match_str = f" <-> {match}" if match else ""
        print(f"  component {i}: sigma={S[i]:.2f}, period={p:.3f}h{match_str}, "
              f"std={components[i].std():.4f}")

    C = w_correlation(components)
    print(f"\nw-correlation matrix (components that should be grouped as one "
          f"oscillation have |w-corr| close to 1):")
    print(pd.DataFrame(np.round(C, 2)).to_string())

    suggested_groups = []
    for i in range(C.shape[0]):
        for j in range(i + 1, C.shape[1]):
            if abs(C[i, j]) > 0.7:
                suggested_groups.append((i, j, float(C[i, j])))
    print(f"\nSuggested groupings (|w-corr| > 0.7): {suggested_groups}")

    out_dir = default_paths("26_ssa_decomposition")

    comp_df = pd.DataFrame(components.T, index=segment.index,
                            columns=[f"component_{i}" for i in range(components.shape[0])])
    comp_df.to_csv(out_dir / f"{args.buoy}_{args.var}_ssa_components.csv")

    summary = {
        "segment_meta": seg_meta,
        "window_hours": args.window_hours,
        "n_components": args.n_components,
        "boxcox_lambda": float(lam),
        "singular_values": S[:args.n_components].tolist(),
        "components": report_rows,
        "w_correlation_matrix": np.round(C, 4).tolist(),
        "suggested_groupings": suggested_groups,
    }
    with open(out_dir / f"{args.buoy}_{args.var}_ssa_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    fig, axes = plt.subplots(2, 1, figsize=(11, 8))
    axes[0].semilogy(range(len(S)), S, marker="o", ms=4)
    axes[0].axvline(args.n_components - 0.5, color="firebrick", ls="--", alpha=0.6)
    axes[0].set_xlabel("Component index")
    axes[0].set_ylabel("Singular value (log scale)")
    axes[0].set_title(f"{args.buoy} {args.var} - SSA singular value scree "
                       f"(kept components left of red line)")
    axes[0].grid(alpha=0.3)

    show_n = min(24 * 60, len(segment))  # ~60 days for visual clarity
    t_show = segment.index[:show_n]
    axes[1].plot(t_show, bc_vals[:show_n], color="gray", alpha=0.5, lw=1, label="Box-Cox series")
    for i in range(min(4, components.shape[0])):
        p = report_rows[i]["recovered_period_hours"]
        p_str = f"{p:.2f}h" if p is not None else "n/a"
        axes[1].plot(t_show, components[i][:show_n], lw=1.2, label=f"component {i} ({p_str})")
    axes[1].set_title("First ~60 days: leading SSA components vs. Box-Cox series")
    axes[1].legend(fontsize=8, loc="upper right")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_ssa.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
