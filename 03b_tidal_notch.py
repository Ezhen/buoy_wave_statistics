"""
Stage 2b - Tidal notch (harmonic regression).

Slots between Stage 2 (stationarity) and Stage 3 (differencing). The
periodogram in Stage 02 showed a dominant M2 peak (~12.42h), and Stage 04/05
showed lag-1 differencing does nothing to it - a ~24.84-sample period isn't
an integer lag, so seasonal differencing can't target it cleanly either.

Instead: fit and subtract M2 (+ its first harmonic, since coastal tidal
signatures are rarely pure sinusoids) via least-squares harmonic regression,
done on the Box-Cox-transformed LEVEL series (variance-stabilize first,
same ordering rule as before - then detide, then Stage 04 only needs to
difference what's left).

--fit-frequency: instead of assuming the dominant tidal period is exactly
M2 (12.4206h), find the actual spectral peak in a band around it and use
THAT period for the harmonic regression. Motivated by Zeebrugge: its
notch never fully cleaned regardless of harmonic count (35-36x baseline
residual power at both 2 and 3 harmonics, barely different) - consistent
with the fundamental frequency assumption being wrong, not insufficient
harmonic count (a shallow/harbor site can have a nearby S2 constituent
at 12.0000h, or a shoaling-shifted effective frequency). Only reliable
with enough record length to resolve M2 from S2 - frequency resolution
scales as 1/T; confirmed at Zeebrugge's 14.3-year record this is ~355x
finer than the M2/S2 spacing needs, comfortably sufficient (was not
true on the original 2-month record, where the two constituents
couldn't be told apart at all).

Usage:
    python 03b_tidal_notch.py --buoy WesthinderBuoy --var VHM0 --harmonics 2
    python 03b_tidal_notch.py --buoy ZeebruggeZandopvangkadeBuoy --var VHM0 --fit-frequency
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import boxcox
from scipy.signal import periodogram

from utils import default_paths, longest_contiguous_segment

M2_PERIOD_HOURS = 12.4206
S2_PERIOD_HOURS = 12.0000  # principal solar semi-diurnal - the usual nearby confuser


DIURNAL_PERIOD_HOURS = 24.0


def build_harmonic_design(t_hours: np.ndarray, period_hours: float, n_harmonics: int):
    """Columns: [1, cos(w t), sin(w t), cos(2w t), sin(2w t), ...]"""
    cols = [np.ones_like(t_hours)]
    for k in range(1, n_harmonics + 1):
        w = 2 * np.pi * k / period_hours
        cols.append(np.cos(w * t_hours))
        cols.append(np.sin(w * t_hours))
    return np.column_stack(cols)


def build_multi_harmonic_design(t_hours: np.ndarray, components: list):
    """Generalization of build_harmonic_design to multiple simultaneous
    fundamental frequencies, fit jointly (one shared intercept, then each
    component's own harmonic columns appended). `components` is a list
    of (period_hours, n_harmonics) tuples.

    Added to test the land-sea-breeze hypothesis: Stage 02 found a
    substantial ~24h diurnal periodogram peak in Hs at every buoy tested
    (85-289x baseline), and Stage 16 confirmed ERA5 wind itself shows the
    same daily cycle at the three "normal" open-water buoys (75-289x) -
    but NOT proportionally at Zeebrugge, where the diurnal signal is far
    stronger (1986x) despite the weakest wind/Hs coupling of the network,
    suggesting a different/additional local driver there. This joint fit
    is offered as an option (--diurnal-harmonics, default 0/off) rather
    than the new default, specifically to keep it separate from
    Zeebrugge's already-established M2-only treatment - land-sea breeze
    is a real physical driver at the other sites but not necessarily
    appropriate to notch the same way at a site where it isn't the
    dominant explanation.

    A fixed-frequency, fixed-phase sinusoid is a reasonable first-order
    model for land-sea breeze specifically because solar heating (its
    physical driver) has a very stable phase relative to local solar
    time even though its AMPLITUDE varies seasonally (stronger land-sea
    temperature contrast in summer) - this is a real simplification
    (one fixed amplitude across the whole multi-year record won't
    capture that seasonal modulation), stated explicitly rather than
    implied to be a complete treatment.
    """
    cols = [np.ones_like(t_hours)]
    for period_hours, n_harmonics in components:
        for k in range(1, n_harmonics + 1):
            w = 2 * np.pi * k / period_hours
            cols.append(np.cos(w * t_hours))
            cols.append(np.sin(w * t_hours))
    return np.column_stack(cols)


def find_dominant_period_in_band(values: np.ndarray, dt_hours: float,
                                   period_min: float, period_max: float):
    """Find the actual spectral peak period within [period_min, period_max]
    hours, rather than assuming it's exactly M2. Frequency resolution is
    1/T (T = series duration) - the caller is responsible for confirming
    that's fine enough to matter before trusting the result (see the
    module docstring's frequency-resolution note)."""
    fs = 1.0 / dt_hours
    freqs, power = periodogram(values, fs=fs, detrend="linear")
    periods = np.divide(1.0, freqs, out=np.full_like(freqs, np.inf), where=freqs > 0)
    mask = (periods >= period_min) & (periods <= period_max)
    if not mask.any():
        return None, None
    idx = np.argmax(power[mask])
    return periods[mask][idx], power[mask][idx]


def periodogram_m2_ratio(series: np.ndarray, dt_hours: float, center_period: float = M2_PERIOD_HOURS):
    fs = 1.0 / dt_hours
    freqs, power = periodogram(series, fs=fs, detrend="linear")
    periods = np.divide(1.0, freqs, out=np.full_like(freqs, np.inf), where=freqs > 0)
    near_center = np.abs(periods - center_period) < 0.5
    baseline = (periods > 6) & (periods < 24) & ~near_center
    if near_center.any() and baseline.any() and np.median(power[baseline]) > 0:
        return power[near_center].max() / np.median(power[baseline])
    return np.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--harmonics", default=2, type=int,
                         help="number of harmonics to fit (1=pure sinusoid, "
                              "2=+first overtone, etc.)")
    parser.add_argument("--fit-frequency", action="store_true",
                         help="find the actual dominant period near M2 instead of "
                              "assuming exactly 12.4206h - see module docstring")
    parser.add_argument("--search-period-min", default=11.5, type=float)
    parser.add_argument("--search-period-max", default=13.5, type=float)
    parser.add_argument("--diurnal-harmonics", default=0, type=int,
                         help="if > 0, jointly fit a 24h diurnal component with "
                              "this many harmonics, alongside M2 - see "
                              "build_multi_harmonic_design docstring for why this "
                              "is opt-in, not the new default")
    args = parser.parse_args()

    in_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    s_full = pd.read_csv(in_path, index_col=0, parse_dates=True)[args.var]
    valid_mask = s_full.notna()
    s_valid = s_full[valid_mask]

    if valid_mask.sum() < len(s_full):
        print(f"NOTE: {len(s_full) - int(valid_mask.sum())} gap sample(s) present "
              f"({100 * (1 - valid_mask.mean()):.1f}% of the grid). Fitting on valid "
              f"data only, but preserving gap positions as NaN in the output so "
              f"downstream lag-based stages (differencing, ACF) don't silently splice "
              f"pre-gap and post-gap periods together as if they were adjacent.")

    # --- Variance-stabilize first (same rule as Stage 3/04), fit on valid data only ---
    shifted = s_valid + 1e-6 if (s_valid <= 0).any() else s_valid
    boxcox_vals, lam = boxcox(shifted.values)
    print(f"Box-Cox lambda = {lam:.4f}")

    # Reindex onto the FULL grid - gap positions become NaN here, not dropped
    s_bc_valid = pd.Series(boxcox_vals, index=s_valid.index)
    s_bc_full = s_bc_valid.reindex(s_full.index)
    s_bc_full.name = f"{args.var}_boxcox"

    dt_hours = (s_full.index[1] - s_full.index[0]).total_seconds() / 3600.0

    # --- Optionally fit the actual dominant tidal period instead of
    #     assuming M2 exactly - on the longest contiguous segment, since
    #     a periodogram is also invalid across a spliced gap. ---
    fitted_period = None
    if args.fit_frequency:
        seg_for_fit, seg_fit_meta = longest_contiguous_segment(s_bc_full)
        record_hours = (seg_for_fit.index[-1] - seg_for_fit.index[0]).total_seconds() / 3600.0
        freq_resolution = 1.0 / record_hours
        freq_needed = abs(1 / S2_PERIOD_HOURS - 1 / M2_PERIOD_HOURS)
        resolution_ratio = freq_needed / freq_resolution
        print(f"Frequency-fitting on longest segment ({len(seg_for_fit)} samples, "
              f"{record_hours:.0f}h): resolution is {resolution_ratio:.1f}x finer "
              f"than the M2/S2 spacing needs "
              f"({'sufficient' if resolution_ratio > 5 else 'MARGINAL - interpret with caution'}).")

        fitted_period, fitted_power = find_dominant_period_in_band(
            seg_for_fit.values, dt_hours, args.search_period_min, args.search_period_max)
        if fitted_period is None:
            print(f"No spectral peak found in [{args.search_period_min}, "
                  f"{args.search_period_max}]h - falling back to nominal M2.")
        else:
            delta_from_m2 = fitted_period - M2_PERIOD_HOURS
            delta_from_s2 = fitted_period - S2_PERIOD_HOURS
            closer_to = "M2" if abs(delta_from_m2) < abs(delta_from_s2) else "S2"
            print(f"Fitted dominant period: {fitted_period:.4f}h (nominal M2={M2_PERIOD_HOURS}h, "
                  f"delta={delta_from_m2:+.4f}h; nominal S2={S2_PERIOD_HOURS}h, "
                  f"delta={delta_from_s2:+.4f}h) - closer to {closer_to}.")

    period_used = fitted_period if fitted_period is not None else M2_PERIOD_HOURS

    # --- Harmonic regression fit ---
    # Design matrix built for the FULL grid (purely a function of elapsed
    # time, so this is fine even at gap positions), but the fit itself
    # only uses rows where data actually exists.
    t0 = s_full.index[0]
    t_hours_full = (s_full.index - t0).total_seconds().values / 3600.0
    components = [(period_used, args.harmonics)]
    if args.diurnal_harmonics > 0:
        components.append((DIURNAL_PERIOD_HOURS, args.diurnal_harmonics))
        print(f"Jointly fitting M2 (period={period_used:.4f}h, {args.harmonics} harmonics) "
              f"+ diurnal (24h, {args.diurnal_harmonics} harmonics)")
    X_full = build_multi_harmonic_design(t_hours_full, components)
    X_fit = X_full[valid_mask.values]
    y_fit = s_bc_valid.values
    coefs, *_ = np.linalg.lstsq(X_fit, y_fit, rcond=None)
    tidal_fit_full = X_full @ coefs
    s_detided = s_bc_full - pd.Series(tidal_fit_full, index=s_full.index)
    s_detided.name = f"{args.var}_detided"

    # --- Report power before/after at the period actually used, on the
    #     longest contiguous segment only ---
    seg_bc, seg_meta = longest_contiguous_segment(s_bc_full)
    seg_detided, _ = longest_contiguous_segment(s_detided)
    if seg_meta["n_segments"] > 1:
        print(f"Record has {seg_meta['n_segments']} contiguous segments after gaps; "
              f"ratio computed on the longest one only "
              f"({seg_meta['pct_of_valid_used']}% of valid samples, "
              f"{seg_meta['segment_start']} to {seg_meta['segment_end']}).")

    ratio_before = periodogram_m2_ratio(seg_bc.values, dt_hours, period_used)
    ratio_after = periodogram_m2_ratio(seg_detided.values, dt_hours, period_used)
    print(f"Peak/baseline power ratio at {period_used:.4f}h - before: {ratio_before:.2f}, "
          f"after: {ratio_after:.2f}")

    diurnal_ratio_before, diurnal_ratio_after = None, None
    if args.diurnal_harmonics > 0:
        diurnal_ratio_before = periodogram_m2_ratio(seg_bc.values, dt_hours, DIURNAL_PERIOD_HOURS)
        diurnal_ratio_after = periodogram_m2_ratio(seg_detided.values, dt_hours, DIURNAL_PERIOD_HOURS)
        print(f"Peak/baseline power ratio at 24.0h (diurnal) - before: "
              f"{diurnal_ratio_before:.2f}, after: {diurnal_ratio_after:.2f}")
    m2_clean = ratio_after <= 3
    if m2_clean:
        print("M2: substantially reduced -> OK.")
    else:
        suggestion = ("try --harmonics with a higher value, or check the sampling "
                       "grid for gaps that would smear the harmonic fit."
                       if args.fit_frequency else
                       "try --fit-frequency (the fundamental may not be exactly M2 "
                       "at this site), or --harmonics with a higher value.")
        print(f"WARNING: M2 power still elevated after notch - {suggestion}")

    if args.diurnal_harmonics > 0:
        diurnal_clean = diurnal_ratio_after <= 3
        if diurnal_clean:
            print("Diurnal: substantially reduced -> OK.")
        else:
            reduction_pct = 100 * (1 - diurnal_ratio_after / diurnal_ratio_before) if diurnal_ratio_before else None
            print(f"WARNING: diurnal power still elevated after notch "
                  f"({'%.0f%% reduced' % reduction_pct if reduction_pct is not None else 'no reduction'}, "
                  f"still {diurnal_ratio_after:.1f}x baseline, threshold is 3x) - "
                  f"increasing --diurnal-harmonics is unlikely to help much further "
                  f"(confirmed on real data: partial, inconsistent reduction of 5-80% "
                  f"across 3 buoys tested, none reached the clean threshold). More "
                  f"likely cause: land-sea breeze phase shifts seasonally (sunrise/"
                  f"sunset timing, storm interference) - a single fixed-phase "
                  f"sinusoid across a multi-year record can't track that. Would need "
                  f"a genuinely different approach (time-varying envelope, seasonal "
                  f"deseasonalizing) to do meaningfully better, not more harmonics "
                  f"at the same fixed phase.")

    if m2_clean and (args.diurnal_harmonics == 0 or diurnal_clean):
        print("-> OK to proceed to Stage 3 differencing on the detided series.")
    else:
        print("-> proceed with caution; see warnings above before trusting "
              "downstream stages built on this detided series.")

    out_dir = default_paths("03b_tidal_notch")
    s_detided.to_csv(out_dir / f"{args.buoy}_{args.var}_detided_boxcox.csv",
                      header=[s_detided.name])

    import json
    with open(out_dir / f"{args.buoy}_{args.var}_notch_summary.json", "w") as f:
        json.dump({
            "boxcox_lambda": float(lam),
            "period_used_hours": float(period_used),
            "fitted_frequency": args.fit_frequency,
            "fitted_period_hours": None if fitted_period is None else float(fitted_period),
            "m2_ratio_before": None if np.isnan(ratio_before) else float(ratio_before),
            "m2_ratio_after": None if np.isnan(ratio_after) else float(ratio_after),
            "harmonics": args.harmonics,
            "diurnal_harmonics": args.diurnal_harmonics,
            "diurnal_ratio_before": (None if diurnal_ratio_before is None or np.isnan(diurnal_ratio_before)
                                      else float(diurnal_ratio_before)),
            "diurnal_ratio_after": (None if diurnal_ratio_after is None or np.isnan(diurnal_ratio_after)
                                     else float(diurnal_ratio_after)),
            "n_gap_segments": seg_meta["n_segments"],
            "longest_segment_pct_of_valid": seg_meta["pct_of_valid_used"],
        }, f, indent=2)

    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    axes[0].plot(s_bc_full.index, s_bc_full.values, lw=0.6)
    axes[0].set_title(f"{args.var} — Box-Cox level (pre-detide)")
    axes[1].plot(s_full.index, tidal_fit_full, lw=0.8, color="darkorange")
    axes[1].set_title(f"Fitted {period_used:.4f}h period + harmonics "
                       f"(n_harmonics={args.harmonics}"
                       f"{', fitted' if args.fit_frequency else ', nominal M2'})")
    axes[2].plot(s_detided.index, s_detided.values, lw=0.6, color="firebrick")
    axes[2].set_title(f"{args.var} — detided (ready for Stage 3 differencing, "
                       f"gaps preserved as NaN)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_tidal_notch.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved detided series + plot to {out_dir}")


if __name__ == "__main__":
    main()
