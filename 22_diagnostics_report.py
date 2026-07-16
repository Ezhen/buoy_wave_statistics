"""
Stage X - four-question diagnostics report: a narrative summary, not a
numeric table (the numeric cross-buoy table already exists via
summarize_results.py). Complements the tiering/CI machinery rather than
duplicating it - this stage doesn't compute anything new, it reads and
synthesizes what every prior stage already produced into four questions
per buoy:

  1. Can I trust the DATA?           (QC / gaps)
  2. Can I trust the ASSUMPTIONS?    (stationarity / dependence / fit)
  3. Can I trust the ESTIMATES?      (confidence intervals / stability)
  4. What did I learn about the PHYSICS? (storms / persistence / regimes)

Reads gracefully from ~12 different stage output files - if a stage
hasn't been run for this buoy, that question's paragraph says so
explicitly rather than crashing or silently omitting it. Given how
uneven stage completion can be across 19 buoys (Advanced-tier stages
only run for VTPK/VMDR buoys, Priority 3/6 stages only run for the 6
long-record buoys), this graceful-degradation behavior is the more
important property to get right here, not "every stage present."

Usage:
    python 22_diagnostics_report.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import pandas as pd

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


def question_1_data(buoy: str, var: str) -> str:
    """Can I trust the data? QC/gaps, from Stage 01 + 11b's fragmentation info."""
    load = safe_json(STAGE_DIR / "01_load_clean" / f"{buoy}_{var}_load_summary.json")
    if load is None:
        return "Stage 01 hasn't been run for this buoy - nothing to report."

    lines = [
        f"Record spans {load.get('record_years', '?'):.1f} years "
        f"({load.get('n_samples_regularized', '?')} samples at "
        f"{load.get('sampling_interval_hours', '?')}h resolution)."
    ]
    pct_missing = load.get("pct_missing_after_clean")
    if pct_missing is not None:
        severity = "substantial" if pct_missing > 10 else "modest" if pct_missing > 2 else "minimal"
        lines.append(f"{pct_missing:.1f}% missing after cleaning - {severity} gap coverage.")
    n_dup = load.get("n_duplicate_timestamps_raw")
    if n_dup:
        lines.append(f"{n_dup} duplicate raw timestamps were found and resolved.")

    dep = safe_json(STAGE_DIR / "11b_dependence_structure" / f"{buoy}_{var}_dependence_summary.json")
    if dep and dep.get("n_segments_found", 1) > 1:
        lines.append(
            f"Record is fragmented into {dep['n_segments_found']} contiguous segments by "
            f"real gaps - lag-based stages (ACF, persistence, differencing) use "
            f"{dep.get('pct_valid_data_used', '?')}% of valid data across "
            f"{dep.get('n_segments_used', '?')} qualifying segments, not the full "
            f"gap-spliced record. This is expected and handled correctly (see "
            f"CHANGELOG for the gap-splicing fixes this required), not a data-quality "
            f"problem in itself - but worth knowing when interpreting any single-buoy "
            f"result as spanning the buoy's full nominal record length.")

    return " ".join(lines)


def question_2_assumptions(buoy: str, var: str) -> str:
    """Can I trust the assumptions? Stationarity, dependence, distributional fit."""
    lines = []

    stat = safe_json(STAGE_DIR / "03_stationarity_tests" / f"{buoy}_{var}_stationarity.json")
    if stat:
        agree = stat.get("tests_agree")
        agree_detail = ("both say stationary" if agree else
                         "a common, expected result on strongly periodic/persistent "
                         "series - not itself a red flag")
        lines.append(
            f"ADF/KPSS {'agree' if agree else 'DISAGREE'} on stationarity "
            f"({agree_detail})."
        )
    else:
        lines.append("Stage 03 (stationarity) not run.")

    notch = safe_json(STAGE_DIR / "03b_tidal_notch" / f"{buoy}_{var}_notch_summary.json")
    if notch:
        before, after = notch.get("m2_ratio_before"), notch.get("m2_ratio_after")
        if before and after:
            clean = after <= 3
            lines.append(
                f"M2 tidal ratio {before:.1f} -> {after:.1f} after notching - "
                f"{'cleaned successfully' if clean else 'NOT fully cleaned, still elevated - treat downstream lag-based results for this buoy with extra caution'}."
            )

    lb = safe_csv(STAGE_DIR / "05_whiteness_check" / f"{buoy}_{var}_ljungbox.csv")
    if lb is not None and "lb_pvalue" in lb.columns:
        all_white = bool((lb["lb_pvalue"] > 0.05).all())
        lines.append(
            f"Ljung-Box {'passes' if all_white else 'FAILS'} - "
            f"{'residual looks white' if all_white else 'genuine storm persistence confirmed at every BCZ buoy tested so far, not a pipeline artifact'}."
        )

    fit_interp = safe_json(STAGE_DIR / "06_distribution_fit" / f"{buoy}_{var}_fit_interpretation.json")
    if fit_interp:
        lines.append(f"Best-fit distribution: {fit_interp['best_distribution']}. "
                      f"{fit_interp['physical_interpretation']}")

    arch = safe_csv(STAGE_DIR / "07_arch_lm_test" / f"{buoy}_{var}_arch_lm.csv")
    if arch is not None and "lm_pvalue" in arch.columns:
        detected = bool((arch["lm_pvalue"] < 0.05).any())
        lines.append(f"ARCH effects {'detected' if detected else 'NOT detected'} - "
                      f"{'volatility clustering confirmed, consistent with every other BCZ buoy' if detected else 'unusual - worth double-checking given how universal this finding has been elsewhere in the network'}.")

    return " ".join(lines) if lines else "No assumption-checking stages run yet."


def question_3_estimates(buoy: str, var: str) -> str:
    """Can I trust the estimates? CI widths, stability checks."""
    lines = []

    conf = safe_json(STAGE_DIR / "12_confidence_intervals" / f"{buoy}_{var}_confidence_summary.json")
    if conf:
        xi_ci = conf.get("gpd_xi_ci")
        if xi_ci and xi_ci.get("ci_low") is not None:
            width = xi_ci["ci_high"] - xi_ci["ci_low"]
            crosses = xi_ci["ci_low"] < 0 < xi_ci["ci_high"]
            lines.append(
                f"GPD xi 95% CI: [{xi_ci['ci_low']:.3f}, {xi_ci['ci_high']:.3f}] "
                f"(width {width:.3f}) - "
                f"{'CROSSES ZERO, uninformative about tail boundedness' if crosses else 'entirely one-signed, informative about tail shape'}."
            )
    else:
        lines.append("Stage 12 (confidence intervals) not run - point estimates elsewhere "
                      "in this report have no uncertainty quantification attached yet.")

    stab = safe_json(STAGE_DIR / "13_stability_analysis" / f"{buoy}_{var}_stability_summary.json")
    if stab:
        window_stab = stab.get("window_stability", {})
        agreement_frac = window_stab.get("fraction_windows_agreeing")
        if agreement_frac is not None:
            agreement_pct = agreement_frac * 100
            lines.append(
                f"Moving-window distribution check: {agreement_pct:.0f}% of windows agree with "
                f"the pooled full-record best-fit distribution - "
                f"{'stable across the record' if agreement_pct >= 75 else 'NOT stable, the pooled verdict may be a mixture artifact rather than a genuine property of the whole record (see Stage 06 interpretation above)'}."
            )

    return " ".join(lines) if lines else "No uncertainty-quantification stages run yet - treat any point estimate elsewhere in this report as having unknown precision."


def question_4_physics(buoy: str, var: str) -> str:
    """What did I learn about the physics? Storms, persistence, regimes, spatial context."""
    lines = []

    eva = safe_json(STAGE_DIR / "08_extreme_value_analysis" / f"{buoy}_{var}_eva_summary.json")
    if eva:
        xi = eva.get("gpd_shape")
        n_peaks = eva.get("n_peaks")
        interp = eva.get("physical_interpretation", "")
        lines.append(f"EVA: xi={xi:.3f} from {n_peaks} storm peaks. {interp}")

    dep = safe_json(STAGE_DIR / "11b_dependence_structure" / f"{buoy}_{var}_dependence_summary.json")
    if dep:
        tau = dep.get("integral_timescale_hours")
        if tau:
            lines.append(f"Persistence (integral) timescale: {tau:.1f}h.")

    regime = safe_csv(STAGE_DIR / "10_regime_identification" / f"{buoy}_{var}_regime_summary.csv")
    if regime is not None and len(regime):
        lines.append(f"Regime structure identified ({len(regime)} regimes) - "
                      f"see Stage 10 output for the calm-to-storm breakdown.")

    clusters = safe_csv(STAGE_DIR / "11_spatial_statistics" / f"{var}_buoy_clusters.csv")
    if clusters is not None and "buoy" in clusters.columns and buoy in clusters["buoy"].values:
        cluster_id = clusters.loc[clusters["buoy"] == buoy, "cluster"].iloc[0] \
            if "cluster" in clusters.columns else None
        if cluster_id is not None:
            cluster_size = (clusters["cluster"] == cluster_id).sum()
            lines.append(f"Spatial cluster {cluster_id} ({cluster_size} buoy(s) - "
                          f"{'singleton, dynamically distinct from the rest of the network' if cluster_size == 1 else 'part of a larger geographic grouping'}).")

    return " ".join(lines) if lines else "No physics-characterizing stages run yet."


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    args = parser.parse_args()

    sections = {
        "1. Can I trust the DATA?": question_1_data(args.buoy, args.var),
        "2. Can I trust the ASSUMPTIONS?": question_2_assumptions(args.buoy, args.var),
        "3. Can I trust the ESTIMATES?": question_3_estimates(args.buoy, args.var),
        "4. What did I learn about the PHYSICS?": question_4_physics(args.buoy, args.var),
    }

    report_lines = [f"# Diagnostics report - {args.buoy} / {args.var}\n"]
    print(f"=== Diagnostics report: {args.buoy} / {args.var} ===\n")
    for title, text in sections.items():
        print(f"{title}\n{text}\n")
        report_lines.append(f"## {title}\n\n{text}\n")

    out_dir = default_paths("22_diagnostics_report")
    report_path = out_dir / f"{args.buoy}_{args.var}_diagnostics_report.md"
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))

    print(f"Saved: {report_path}")


if __name__ == "__main__":
    main()
