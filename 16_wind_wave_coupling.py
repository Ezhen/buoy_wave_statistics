"""
Wind-wave coupling: first real cross-variable analysis with actual wind
data, pairing ERA5 meteo with buoy Hs. Extends Stage 09's cross-variable
framework - genuinely useful across the whole network now, not just
CadzandBoei/Deurlo, since wind is external data rather than something
only 2 buoys happen to measure.

Data alignment handled by utils.load_era5_for_buoy(): ERA5 is gridded
(0.25deg), buoys are points - nearest-grid-cell extraction, concatenated
across every monthly ERA5 file. Buoy Hs (10-30 min) is downsampled to
3-hourly bins (mean AND max, matching ERA5's native resolution) rather
than upsampling ERA5, and only the genuinely overlapping time period is
used (buoy records may start before ERA5's downloaded range).

Analyses:
  - Wind-speed/Hs cross-correlation at multiple lags - does wind lead Hs
    by a predictable amount? Physical expectation: near-zero or small
    positive lag, unlike a slower-responding quantity.
  - MSLP/Hs cross-correlation - pressure drops are a known leading
    indicator of wind speed increases; report the best lag.
  - Directional alignment (only for buoys with VMDR - CadzandBoei/Deurlo):
    circular difference between wind direction and wave direction, split
    by Stage 10's regime label if available - tighter alignment during
    storm regimes would indicate wind-sea; decoupling would indicate
    swell dominance.
  - Simple regression (Hs ~ wind_speed, Hs ~ wind_speed^2 - wave energy
    scales roughly with wind^2 in fetch-limited growth) as a
    "how much does wind alone explain" R^2 summary.
  - Diurnal (24h) power ratio in ERA5's own wind_speed series, checked
    on wind's full downloaded range (not the buoy-overlap period, for
    more statistical power on the wind signal itself) - added after
    Stage 02 found a substantial ~24h periodogram peak in Hs at every
    buoy tested (85-179x baseline), too consistent across sites to be
    tidal/site-specific. Tests the land-sea-breeze hypothesis directly:
    if wind itself shows the same daily cycle, that's the likely
    physical driver, not an astronomical tidal constituent - which
    would call for a different fix than the fixed-phase harmonic notch
    used for M2 (a stochastic, weather-conditional pattern shouldn't be
    notched the same way a precise tidal constituent is).

Usage:
    python 16_wind_wave_coupling.py --data-dir data_multiyear --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr

from utils import (default_paths, load_era5_for_buoy, get_scalar_latlon,
                    load_buoy_dataframe, integral_timescale)


def cross_correlation(x: np.ndarray, y: np.ndarray, max_lag: int):
    """CCF of y relative to x, lags in samples. Verified empirically
    (not just derived by eye - lag-sign conventions are easy to get
    backwards): POSITIVE lag = x leads y (x's current value shows up in
    y at a later time)."""
    x = (x - x.mean()) / x.std()
    y = (y - y.mean()) / y.std()
    lags = range(-max_lag, max_lag + 1)
    ccf = []
    for lag in lags:
        if lag < 0:
            a, b = x[-lag:], y[:lag] if lag != 0 else y
        elif lag > 0:
            a, b = x[:-lag], y[lag:]
        else:
            a, b = x, y
        ccf.append(np.corrcoef(a, b)[0, 1] if len(a) > 10 else np.nan)
    return np.array(list(lags)), np.array(ccf)


def circular_diff_deg(a: np.ndarray, b: np.ndarray):
    """Smallest angular difference between two sets of compass bearings, 0-180."""
    d = np.abs(a - b) % 360
    return np.minimum(d, 360 - d)


def diurnal_power_ratio(values: np.ndarray, dt_hours: float, period_hours: float = 24.0):
    """Peak periodogram power near `period_hours` vs. local baseline -
    same logic as Stage 02's M2/diurnal checks, reimplemented here
    (not refactored into a shared utility, to avoid touching
    already-tested code) for checking ERA5 wind_speed's own periodicity."""
    from scipy.signal import periodogram
    fs = 1.0 / dt_hours
    freqs, power = periodogram(values, fs=fs, detrend="linear")
    periods = np.divide(1.0, freqs, out=np.full_like(freqs, np.inf), where=freqs > 0)
    near = np.abs(periods - period_hours) < 1.0
    baseline = (periods > 6) & (periods < 48) & ~near
    if near.any() and baseline.any() and np.median(power[baseline]) > 0:
        return float(power[near].max() / np.median(power[baseline]))
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data_multiyear", type=Path)
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--era5-dir", default="meteo_era5")
    parser.add_argument("--max-lag-hours", default=36, type=int,
                         help="max lag to test in the CCF, in hours (converted to "
                              "3-hourly steps internally)")
    args = parser.parse_args()

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    hs = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var].dropna()

    nc_path = args.data_dir / f"{args.buoy}.nc"
    with xr.open_dataset(nc_path) as ds:
        lat, lon = get_scalar_latlon(ds)
    print(f"--- {args.buoy} / {args.var} wind-wave coupling ---")
    print(f"Buoy location: lat={lat:.3f}, lon={lon:.3f}")

    era5 = load_era5_for_buoy(lat, lon, args.era5_dir)
    print(f"ERA5: {len(era5)} 3-hourly samples, {era5.index.min()} -> {era5.index.max()}")

    # --- Align: downsample Hs to 3-hourly (mean and max), inner join with ERA5 ---
    hs_3h_mean = hs.resample("3h").mean()
    hs_3h_max = hs.resample("3h").max()

    paired = pd.DataFrame({"hs_mean": hs_3h_mean, "hs_max": hs_3h_max}).join(era5, how="inner")
    paired = paired.dropna(subset=["hs_mean", "wind_speed"])
    print(f"Overlapping period: {len(paired)} samples, "
          f"{paired.index.min()} -> {paired.index.max()} "
          f"({(paired.index.max() - paired.index.min()).days / 365.25:.1f} years)")

    if len(paired) < 100:
        print("Too little overlap to analyze meaningfully - stopping.")
        return

    out_dir = default_paths("16_wind_wave_coupling")
    summary = {"buoy_lat": lat, "buoy_lon": lon, "n_overlap_samples": len(paired),
               "overlap_start": str(paired.index.min()), "overlap_end": str(paired.index.max())}

    # --- Wind speed / Hs cross-correlation ---
    max_lag_samples = args.max_lag_hours // 3
    lags, ccf_wind = cross_correlation(paired["wind_speed"].values, paired["hs_mean"].values, max_lag_samples)
    best_idx = np.nanargmax(np.abs(ccf_wind))
    best_lag_wind = int(lags[best_idx]) * 3  # back to hours
    print(f"\nWind speed / Hs: strongest |r|={abs(ccf_wind[best_idx]):.3f} at lag={best_lag_wind}h "
          f"({'wind leads' if best_lag_wind > 0 else 'wind lags' if best_lag_wind < 0 else 'contemporaneous'})")
    summary["wind_hs_best_lag_hours"] = best_lag_wind
    summary["wind_hs_best_corr"] = float(ccf_wind[best_idx])

    # --- MSLP / Hs cross-correlation ---
    if "msl" in paired.columns:
        lags_p, ccf_msl = cross_correlation(paired["msl"].values, paired["hs_mean"].values, max_lag_samples)
        best_idx_p = np.nanargmax(np.abs(ccf_msl))
        best_lag_msl = int(lags_p[best_idx_p]) * 3
        print(f"MSLP / Hs: strongest |r|={abs(ccf_msl[best_idx_p]):.3f} at lag={best_lag_msl}h "
              f"({'MSLP leads' if best_lag_msl > 0 else 'MSLP lags' if best_lag_msl < 0 else 'contemporaneous'})")
        if best_lag_msl > 0 and ccf_msl[best_idx_p] < 0:
            print("  (positive lag + negative correlation: pressure DROP precedes Hs RISE - "
                  "matches the expected leading-indicator mechanism)")
        summary["mslp_hs_best_lag_hours"] = best_lag_msl
        summary["mslp_hs_best_corr"] = float(ccf_msl[best_idx_p])

    # --- Simple regression: how much does wind alone explain ---
    ws = paired["wind_speed"].values
    hsv = paired["hs_mean"].values
    r2_linear = np.corrcoef(ws, hsv)[0, 1] ** 2
    r2_quad = np.corrcoef(ws ** 2, hsv)[0, 1] ** 2
    print(f"\nHs ~ wind_speed:  R^2={r2_linear:.3f}")
    print(f"Hs ~ wind_speed^2: R^2={r2_quad:.3f}")
    summary["r2_wind_linear"] = float(r2_linear)
    summary["r2_wind_quadratic"] = float(r2_quad)

    # --- Diurnal (24h) check on wind's OWN full range (land-sea-breeze hypothesis) ---
    wind_diurnal_ratio = diurnal_power_ratio(era5["wind_speed"].dropna().values, dt_hours=3.0)
    if wind_diurnal_ratio is not None:
        print(f"\nERA5 wind_speed diurnal (24h) power ratio (full ERA5 range, "
              f"n={era5['wind_speed'].notna().sum()}): {wind_diurnal_ratio:.2f}")
        if wind_diurnal_ratio > 3:
            print("  (ratio > 3: wind itself shows a real daily cycle - land-sea-breeze "
                  "is a plausible driver of the diurnal signal Stage 02 found in Hs at "
                  "every buoy, distinct from an astronomical tidal constituent)")
        else:
            print("  (ratio <= 3: wind does not show a strong daily cycle here - the "
                  "diurnal Hs signal likely has a different driver, e.g. a weak "
                  "astronomical S1 tidal constituent)")
        summary["wind_diurnal_ratio"] = wind_diurnal_ratio

    # --- Wind's own persistence timescale - added to explain a real finding:
    #     ARMAX (Stage E) showed skill staying positive out to 24h at
    #     Westhinder, not decaying/reversing beyond the 3h wind-Hs lag the
    #     way a synthetic test (with fast, ~20h-period oscillating wind)
    #     predicted. Same ACF-based method as 11b, just applied to wind
    #     instead of Hs.
    #
    #     MUST deseasonalize first - raw wind_speed gave tau=343h (~14
    #     days), far beyond any plausible synoptic timescale - the SAME
    #     category of mistake 11b itself made originally (before the
    #     raw-vs-detided fix): an ACF-based integral timescale can't
    #     distinguish genuine short-term memory decay from real-but-slow
    #     seasonal cycling (wind is windier in winter, calmer in summer).
    #
    #     FIRST deseasonalizing attempt (30-day rolling mean subtraction)
    #     was ALSO wrong, in a different way - produced a NEGATIVE
    #     timescale (-11.0h) on real data, a mathematical impossibility.
    #     Root cause: subtracting a rolling mean is a known source of
    #     spurious NEGATIVE autocorrelation in the residual (a mechanical
    #     artifact of local-average filtering, not a real property of
    #     wind - nearby points tend to deviate from their rolling average
    #     in opposite directions). Fixed by reusing a method already
    #     TRUSTED in this exact pipeline instead of introducing a new,
    #     less-tested one: harmonic regression at the annual period -
    #     literally the same approach Stage 03b already uses for the M2
    #     tidal notch, just fit at 365.25 days instead of 12.4206 hours.
    #     A fitted deterministic sinusoid doesn't mechanically induce
    #     this artifact the way a local rolling filter does. ---
    wind_valid = era5["wind_speed"].dropna()
    if len(wind_valid) > 30:
        t_hours_wind = (wind_valid.index - wind_valid.index[0]).total_seconds().values / 3600.0
        annual_period_hours = 365.25 * 24
        n_harmonics = 3
        design_cols = [np.ones_like(t_hours_wind)]
        for k in range(1, n_harmonics + 1):
            w = 2 * np.pi * k / annual_period_hours
            design_cols.append(np.cos(w * t_hours_wind))
            design_cols.append(np.sin(w * t_hours_wind))
        X_seasonal = np.column_stack(design_cols)
        coefs, *_ = np.linalg.lstsq(X_seasonal, wind_valid.values, rcond=None)
        seasonal_fit = X_seasonal @ coefs
        wind_deseasonalized = wind_valid.values - seasonal_fit

        tau_raw, _, _, _, hit_ceiling_raw = integral_timescale(
            wind_valid.values, dt_hours=3.0, max_lag=2000, consecutive=5)
        tau_hours, crit_lag, rho, band, hit_ceiling = integral_timescale(
            wind_deseasonalized, dt_hours=3.0, max_lag=2000, consecutive=5)

        print(f"\nERA5 wind_speed integral (persistence) timescale:")
        print(f"  RAW (seasonally contaminated, not the right number to use): "
              f"{tau_raw:.1f}h")
        print(f"  DESEASONALIZED (annual harmonic regression, {n_harmonics} "
              f"harmonics, removed): {tau_hours:.1f}h "
              f"{'(hit search ceiling - LOWER BOUND)' if hit_ceiling else ''}")
        if tau_hours < 0:
            print(f"  WARNING: negative timescale is mathematically impossible - "
                  f"the deseasonalizing method is still introducing an artifact. "
                  f"Do not use this number.")
        print(f"  Compare the DESEASONALIZED value to where ARMAX skill vs. plain "
              f"ARMA actually peaks/decays (Stage 20's output) - if wind's own "
              f"short-term persistence extends out to a similar timescale, that "
              f"directly explains why 'persist last known wind' stays a reasonable "
              f"approximation well beyond the raw wind-Hs lag, rather than "
              f"degrading sharply the way the synthetic validation test (fast, "
              f"~20h-period oscillating wind, deliberately unlike real weather) "
              f"predicted.")
        summary["wind_persistence_timescale_hours_raw"] = float(tau_raw)
        summary["wind_persistence_timescale_hours_deseasonalized"] = float(tau_hours)
        summary["wind_persistence_hit_ceiling"] = bool(hit_ceiling)

    # --- Directional alignment, only if this buoy has VMDR ---
    buoy_df = load_buoy_dataframe(nc_path, varnames=("VMDR",))
    if "VMDR" in buoy_df.columns:
        vmdr_3h = buoy_df["VMDR"].resample("3h").mean()
        dir_paired = pd.DataFrame({"wave_dir": vmdr_3h}).join(
            era5[["wind_dir_from_deg"]], how="inner").dropna()
        if len(dir_paired) > 50:
            diffs = circular_diff_deg(dir_paired["wave_dir"].values,
                                       dir_paired["wind_dir_from_deg"].values)
            print(f"\nWind/wave direction alignment ({len(dir_paired)} samples): "
                  f"mean angular difference={diffs.mean():.1f} deg "
                  f"({'well-aligned (wind-sea)' if diffs.mean() < 45 else 'decoupled (swell-influenced)'})")
            summary["mean_wind_wave_dir_diff_deg"] = float(diffs.mean())

            regime_path = Path("pipeline_out/10_regime_identification") / f"{args.buoy}_{args.var}_regime_labels.csv"
            if regime_path.exists():
                regimes = pd.read_csv(regime_path, index_col=0, parse_dates=True)["regime"]
                regimes_3h = regimes.reindex(dir_paired.index, method="nearest", tolerance=pd.Timedelta("3h"))
                by_regime = pd.Series(diffs, index=dir_paired.index).groupby(regimes_3h).mean()
                print("  Mean angular difference by regime (0=calmest):")
                for r, v in by_regime.items():
                    print(f"    regime {int(r)}: {v:.1f} deg")
                summary["dir_diff_by_regime"] = {int(k): float(v) for k, v in by_regime.items()}
    else:
        print(f"\n{args.buoy} has no VMDR - skipping directional alignment "
              f"(only CadzandBoei/Deurlo carry this variable).")

    with open(out_dir / f"{args.buoy}_{args.var}_wind_coupling_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # --- Plots ---
    fig, axes = plt.subplots(1, 2 if "msl" in paired.columns else 1, figsize=(14, 5))
    axes = np.atleast_1d(axes)
    axes[0].stem(lags * 3, ccf_wind)
    axes[0].axvline(0, color="gray", ls="--")
    axes[0].set_xlabel("lag (hours; positive = wind leads Hs)")
    axes[0].set_ylabel("cross-correlation")
    axes[0].set_title("Wind speed / Hs CCF")
    if "msl" in paired.columns:
        axes[1].stem(lags_p * 3, ccf_msl)
        axes[1].axvline(0, color="gray", ls="--")
        axes[1].set_xlabel("lag (hours; positive = MSLP leads Hs)")
        axes[1].set_title("MSLP / Hs CCF")
    fig.suptitle(f"{args.buoy} — wind-wave coupling CCFs")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_ccf.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(ws, hsv, s=4, alpha=0.3)
    ax.set_xlabel("ERA5 wind speed (m/s)")
    ax.set_ylabel(f"{args.var} (m)")
    ax.set_title(f"{args.buoy} — Hs vs. wind speed (R^2={r2_linear:.3f})")
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_scatter.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
