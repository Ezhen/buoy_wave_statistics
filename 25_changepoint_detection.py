"""
Stage 25 - change-point detection: directly addresses a gap Mann-Kendall
(Stage 14) structurally cannot fill. Mann-Kendall tests for a smooth
monotonic trend; if the real story at a buoy is a discrete step-change
(mooring relocation, sensor upgrade, processing-method change at a
specific era boundary) rather than gradual drift, Mann-Kendall can
report "no trend" correctly while still missing a genuine regime shift
entirely - these are different hypotheses, not the same question at
different sensitivity.

Uses PELT (Pruned Exact Linear Time, via `ruptures`) with an L2
(mean-shift) cost model, not full MCMC-based Bayesian change-point
inference - a deliberate scope decision, not a shortcut taken without
thinking about it: PELT answers "does a change point exist, and
roughly where" (the actual question here) without the added complexity
of prior specification and convergence diagnostics a fully Bayesian
version would need for the same practical answer.

Penalty selection is a real, acknowledged design choice (higher
penalty = fewer detected change points) - reported across a sweep of
multipliers rather than presented as a single definitive value, so the
sensitivity of the result to this choice is visible rather than hidden
behind one arbitrary number.

Runs on annual mean AND annual p95 Hs, recomputed directly from Stage
0's output (Stage 14 does not save its own annual series, only the
trend-test result) - same 6-long-record-buoy scope as Stage 14/15.

Usage:
    python 25_changepoint_detection.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ruptures as rpt

from utils import default_paths


def annual_coverage(hs_full: pd.Series) -> pd.DataFrame:
    """Per-year sample coverage, computed from the FULL regularized grid
    (including NaN gap rows left by Stage 01) - must run on this BEFORE
    any .dropna(), or the gap information is gone (a dropna'd series has
    no record of how many samples were missing, only which ones survived).

    Coverage is measured against the record's own regularized grid (how
    many samples for THIS year exist in the full index at the inferred
    sampling interval), not a theoretical 365.25-day calendar count - so
    a real gap and a legitimately partial first/last year (buoy deployed
    or decommissioned mid-year) aren't conflated. The latter is flagged
    separately via `is_boundary_year` rather than scored as a "gap."

    Also reports STORM-SEASON (Oct-Mar) coverage separately from the
    annual figure. Total annual coverage alone can't distinguish a gap
    that falls in the calm season (no effect on annual mean/p95) from
    one that falls in the storm season (biases both low) - two years
    can have near-identical annual coverage and completely different
    effects on the statistic, so the annual number by itself is not
    sufficient evidence either way.
    """
    years = hs_full.index.year
    is_storm_season = hs_full.index.month.isin([10, 11, 12, 1, 2, 3])

    n_expected = hs_full.groupby(years).size()
    n_valid = hs_full.groupby(years).apply(lambda s: int(s.notna().sum()))

    storm = hs_full[is_storm_season]
    storm_years = storm.index.year
    n_storm_expected = storm.groupby(storm_years).size()
    n_storm_valid = storm.groupby(storm_years).apply(lambda s: int(s.notna().sum()))

    all_years = sorted(years.unique())
    df = pd.DataFrame({
        "year": n_expected.index.values,
        "n_expected": n_expected.values,
        "n_valid": n_valid.values,
    })
    df["coverage_pct"] = (100 * df["n_valid"] / df["n_expected"]).round(1)
    df = df.set_index("year")
    df["n_storm_expected"] = n_storm_expected
    df["n_storm_valid"] = n_storm_valid
    df["storm_coverage_pct"] = (100 * df["n_storm_valid"] / df["n_storm_expected"]).round(1)
    df["is_boundary_year"] = df.index.isin([all_years[0], all_years[-1]])
    return df


def detect_changepoints(series: np.ndarray, penalty_multiplier: float = 1.0):
    """PELT with an L2 cost model (mean-shift detection), penalty scaled
    by a BIC-style heuristic (variance * log(n)) - standard for L2 cost,
    not an arbitrary constant."""
    n = len(series)
    sigma2 = np.var(series)
    penalty = penalty_multiplier * sigma2 * np.log(n)
    algo = rpt.Pelt(model="l2").fit(series)
    breakpoints = algo.predict(pen=penalty)
    return breakpoints[:-1]  # ruptures includes n as the final "breakpoint" by convention


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--min-years", default=10.0, type=float)
    parser.add_argument("--start-year", default=None, type=int,
                         help="trim the record to start at this year - added to test "
                              "whether an early-record change point is a genuine signal "
                              "or an edge artifact from a short 'before' segment (PELT's "
                              "mean estimate for a short segment is less stable, which can "
                              "look like a real shift even without one)")
    parser.add_argument("--min-coverage-pct", default=80.0, type=float,
                         help="years with sample coverage below this threshold are flagged "
                              "as possible data-completeness artifacts rather than real "
                              "signal - a year with a coverage gap in the stormy season "
                              "biases annual mean/p95 low without any real physical change")
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"

    hs_full = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var]
    if args.start_year is not None:
        n_before = len(hs_full)
        hs_full = hs_full[hs_full.index.year >= args.start_year]
        print(f"Trimmed to start at {args.start_year}: {n_before} -> {len(hs_full)} samples")

    coverage = annual_coverage(hs_full)
    hs = hs_full.dropna()

    record_years = None
    if args.start_year is None and load_summary_path.exists():
        with open(load_summary_path) as f:
            record_years = json.load(f).get("record_years")
    if record_years is None:
        record_years = (hs.index[-1] - hs.index[0]).days / 365.25

    print(f"--- {args.buoy} / {args.var} change-point detection ---")
    print(f"Record length: {record_years:.1f} years")
    if record_years < args.min_years:
        print(f"Record is under --min-years={args.min_years} - refusing to run.")
        return

    annual_mean = hs.groupby(hs.index.year).mean()
    annual_p95 = hs.groupby(hs.index.year).quantile(0.95)
    years = annual_mean.index.values

    out_dir = default_paths("25_changepoint_detection")
    summary = {"record_years": record_years, "years": years.tolist()}

    penalty_multipliers = [0.5, 1.0, 1.5, 2.0, 4.0]

    fig, axes = plt.subplots(3, 1, figsize=(11, 11))
    for ax, series, label, key in [
        (axes[0], annual_mean, "Annual mean Hs", "mean"),
        (axes[1], annual_p95, "Annual p95 Hs (storm intensity)", "p95"),
    ]:
        print(f"\n[{label}]")
        sweep = {}
        for pm in penalty_multipliers:
            bkps_idx = detect_changepoints(series.values, penalty_multiplier=pm)
            bkps_years = [int(years[i]) for i in bkps_idx]
            sweep[pm] = bkps_years
            print(f"  penalty x{pm}: change point(s) at {bkps_years if bkps_years else '(none detected)'}")

        summary[f"{key}_changepoints_by_penalty"] = {str(k): v for k, v in sweep.items()}

        ax.plot(years, series.values, marker="o", ms=4, lw=1, color="steelblue")
        # Mark the penalty=1.0 (default) result specifically
        default_bkps = sweep[1.0]
        for cp_year in default_bkps:
            ax.axvline(cp_year, color="firebrick", ls="--", alpha=0.7)
        ax.set_ylabel(label)
        ax.set_xlabel("Year")
        ax.set_title(f"{label} - change points at default penalty: "
                     f"{default_bkps if default_bkps else 'none'}")
        ax.grid(alpha=0.3)

    # Coverage subplot - storm-season (Oct-Mar) coverage is the bar since
    # that's the variable that actually explains an anomalous mean/p95;
    # annual coverage shown as a marker for context, since the two can
    # diverge (a year can have fine annual coverage but a gap concentrated
    # in the storm months, or vice versa).
    ax_cov = axes[2]
    cov_years = coverage.index.values
    colors = ["firebrick" if (c < args.min_coverage_pct and not b) else "steelblue"
              for c, b in zip(coverage["storm_coverage_pct"], coverage["is_boundary_year"])]
    ax_cov.bar(cov_years, coverage["storm_coverage_pct"].values, color=colors, width=0.8,
               label="Storm-season (Oct-Mar) coverage")
    ax_cov.plot(cov_years, coverage["coverage_pct"].values, marker="d", ms=4, lw=0,
                color="black", alpha=0.7, label="Annual coverage (context)")
    ax_cov.axhline(args.min_coverage_pct, color="black", ls=":", lw=1, alpha=0.6)
    ax_cov.set_ylabel("Coverage %")
    ax_cov.set_xlabel("Year")
    ax_cov.set_ylim(0, 105)
    ax_cov.set_title(f"Storm-season coverage (bars, red below {args.min_coverage_pct:.0f}%) "
                      f"vs annual coverage (diamonds)")
    ax_cov.legend(loc="lower left", fontsize=8)
    ax_cov.grid(alpha=0.3)

    fig.suptitle(f"{args.buoy} — {args.var} change-point detection (PELT, L2 cost)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_changepoints.png", dpi=150)
    plt.close(fig)

    summary["annual_coverage"] = {
        str(y): {"n_expected": int(r.n_expected), "n_valid": int(r.n_valid),
                  "coverage_pct": float(r.coverage_pct),
                  "n_storm_expected": int(r.n_storm_expected), "n_storm_valid": int(r.n_storm_valid),
                  "storm_coverage_pct": float(r.storm_coverage_pct),
                  "is_boundary_year": bool(r.is_boundary_year)}
        for y, r in coverage.iterrows()
    }

    with open(out_dir / f"{args.buoy}_{args.var}_changepoint_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Cross-reference: scan every year INSIDE each detected segment (not
    # just years adjacent to a change point boundary) against storm-season
    # coverage specifically. Annual coverage alone isn't sufficient here -
    # real Westhinder data showed two segments with near-identical average
    # annual coverage (~80%) where only one had an anomalous mean, so total
    # volume of missing data doesn't discriminate; a gap concentrated in
    # the storm season (Oct-Mar) biases mean/p95 low, a gap in the calm
    # season doesn't. PELT marks segment BOUNDARIES, not the interior point
    # actually driving an anomalous segment mean, so a low-coverage year
    # can sit anywhere inside a segment, not just at its edges - checking
    # only cp_year +/- 1 misses that by construction.
    low_storm_cov_years = set(coverage.index[(coverage["storm_coverage_pct"] < args.min_coverage_pct)
                                              & (~coverage["is_boundary_year"])])

    default_mean_bkps = sorted(summary["mean_changepoints_by_penalty"]["1.0"])
    default_p95_bkps = sorted(summary["p95_changepoints_by_penalty"]["1.0"])
    all_bkps = sorted(set(default_mean_bkps) | set(default_p95_bkps))

    print(f"\n[Storm-season (Oct-Mar) coverage check, --min-coverage-pct={args.min_coverage_pct:.0f}]")
    if low_storm_cov_years:
        print(f"  Low storm-season coverage years: {sorted(low_storm_cov_years)}")
    else:
        print(f"  No low storm-season coverage years found (excluding boundary years).")

    if all_bkps:
        segment_bounds = [int(years[0])] + all_bkps + [int(years[-1]) + 1]
        print(f"\n  Segments at default (1.0x) penalty, bounded by: {all_bkps}")
        for lo, hi in zip(segment_bounds[:-1], segment_bounds[1:]):
            seg_years = [y for y in coverage.index if lo <= y < hi]
            seg_low = sorted(set(seg_years) & low_storm_cov_years)
            avg_annual = coverage.loc[seg_years, "coverage_pct"].mean()
            avg_storm = coverage.loc[seg_years, "storm_coverage_pct"].mean()
            print(f"    [{lo}-{hi-1}] avg annual coverage {avg_annual:.1f}%, "
                  f"avg storm-season coverage {avg_storm:.1f}%"
                  + (f" -- LOW STORM-SEASON COVERAGE YEARS INSIDE SEGMENT: {seg_low}"
                     if seg_low else ""))
    else:
        print(f"  No change points at default penalty - nothing to segment-check.")

    print(f"\n(Penalty is a real design choice - higher penalty finds fewer, more "
          f"conservative change points. The sweep above shows how sensitive the "
          f"result is; a change point that persists across most/all penalty levels "
          f"is a more robust finding than one that only appears at the most "
          f"permissive setting.)")
    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
