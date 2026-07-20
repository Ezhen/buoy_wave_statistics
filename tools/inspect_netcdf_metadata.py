"""
Exploratory - dump everything the raw NetCDF file itself carries, not
just TIME + the value column every other stage reads.

CMEMS in-situ files routinely carry global attributes (history, source,
platform_code, wmo_platform_code, instrument, comment, references,
qc_manual) and per-observation QC flag variables (<VAR>_QC) - none of
this has been inspected anywhere in this pipeline so far. A sensor swap,
reprocessing pass, or QC-methodology change is exactly the kind of thing
that could show up in a `history` string or a shift in QC flag
composition at a specific year, even with no public documentation of it
anywhere else.

Also checks whether LATITUDE/LONGITUDE are stored per-timestep (some
CMEMS products do, for moored platforms that get redeployed) - if so, a
real mooring relocation would show up directly as a position shift.

Usage:
    python inspect_netcdf_metadata.py --data-dir data_multiyear --buoy WesthinderBuoy --var VHM0
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import xarray as xr

from utils import resolve_coord_name, pick_var, VAR_CANDIDATES


def describe_attrs(attrs: dict, indent: str = ""):
    if not attrs:
        print(f"{indent}(none)")
        return
    for k, v in attrs.items():
        v_str = str(v)
        if len(v_str) > 300:
            v_str = v_str[:300] + " ... [truncated]"
        print(f"{indent}{k}: {v_str}")


def qc_flag_report(ds, var_name: str, time: pd.DatetimeIndex):
    """CMEMS convention names the QC companion variable '<var>_QC' with
    integer flag codes - but the exact code meanings live in that
    variable's OWN flag_values/flag_meanings attributes, not assumed
    from the general CMEMS convention, since a provider silently
    deviating from it is exactly the kind of thing worth catching here
    rather than papering over with an assumed mapping."""
    qc_name = f"{var_name}_QC"
    if qc_name not in ds.variables:
        print(f"  No {qc_name} variable found in this file.")
        return

    qc = np.asarray(ds[qc_name].values)
    if qc.ndim > 1:
        qc = qc.reshape(qc.shape[0], -1)[:, 0]

    val = np.asarray(ds[var_name].values)
    if val.ndim > 1:
        val = val.reshape(val.shape[0], -1)[:, 0]

    print(f"\n  --- {qc_name} ---")
    print(f"  Attributes (flag_values/flag_meanings define what the codes below mean):")
    describe_attrs(ds[qc_name].attrs, indent="    ")

    s = pd.Series(qc, index=time[: len(qc)])
    print(f"\n  Overall flag value counts:")
    print(s.value_counts().sort_index().to_string())

    print(f"\n  Flag value counts by year (rows=year, cols=flag code):")
    by_year = s.groupby(s.index.year).value_counts().unstack(fill_value=0)
    print(by_year.to_string())

    # A QC-methodology change often shows up as a step in the FRACTION of
    # non-good flags per year, not necessarily total sample count (which
    # Stage 25's coverage check already covers) - flag any year where the
    # "not the modal/most-common flag code" share jumps.
    modal_flag = s.mode().iloc[0]
    non_modal_frac = s.groupby(s.index.year).apply(lambda x: (x != modal_flag).mean())
    print(f"\n  Fraction NOT flag={modal_flag} (the most common code), by year:")
    print(non_modal_frac.round(3).to_string())

    # CRITICAL CHECK: nothing in this pipeline has ever filtered on QC -
    # Stage 01's cleaning only reacts to actual NaN in the value array.
    # A row flagged e.g. "missing_value" or "bad_data" by the provider
    # but NOT stored as NaN in {var_name} itself would silently survive
    # as "valid" data all the way through the pipeline. Check every
    # non-modal flag code, not just the missing-value one - a "bad_data"
    # or "value_changed" flag surviving as a finite number would be an
    # even more direct correctness problem than a missing-value one.
    print(f"\n  --- Does a non-{modal_flag} QC flag actually correspond to NaN in {var_name}? ---")
    n = min(len(qc), len(val))
    qc_n, val_n = qc[:n], val[:n]
    for flag_code in sorted(np.unique(qc_n)):
        if flag_code == modal_flag:
            continue
        mask = qc_n == flag_code
        n_flagged = int(mask.sum())
        n_nan = int(np.isnan(val_n[mask]).sum())
        n_finite = n_flagged - n_nan
        status = "OK - all NaN" if n_finite == 0 else "FLAG: finite values present despite non-good QC code"
        print(f"    flag={flag_code}: {n_flagged} rows, {n_nan} are NaN in {var_name}, "
              f"{n_finite} are FINITE (present as 'valid' data downstream) - {status}")
        if n_finite > 0:
            finite_vals = val_n[mask][~np.isnan(val_n[mask])]
            print(f"      finite value stats at flag={flag_code}: "
                  f"min={finite_vals.min():.3f}, max={finite_vals.max():.3f}, "
                  f"mean={finite_vals.mean():.3f}")


def position_report(ds, time: pd.DatetimeIndex):
    try:
        lat_name = resolve_coord_name(ds, "LATITUDE")
        lon_name = resolve_coord_name(ds, "LONGITUDE")
    except KeyError:
        print("\n--- No LATITUDE/LONGITUDE found. ---")
        return

    lat = np.asarray(ds[lat_name].values).reshape(-1)
    lon = np.asarray(ds[lon_name].values).reshape(-1)

    if lat.size <= 1 and lon.size <= 1:
        print(f"\n--- Position: single scalar LATITUDE/LONGITUDE "
              f"({float(lat[0]):.5f}, {float(lon[0]):.5f}) - file does not store "
              f"per-timestep position, so a mooring relocation wouldn't be "
              f"visible here even if one happened. ---")
        return

    n = min(len(lat), len(lon), len(time))
    lat_s = pd.Series(lat[:n], index=time[:n])
    lon_s = pd.Series(lon[:n], index=time[:n])
    print(f"\n--- Position over time (mooring relocation check) ---")
    print(f"  LATITUDE range: {lat_s.min():.5f} to {lat_s.max():.5f}")
    print(f"  LONGITUDE range: {lon_s.min():.5f} to {lon_s.max():.5f}")

    lat_mode = lat_s.round(3).mode().iloc[0]
    lon_mode = lon_s.round(3).mode().iloc[0]
    shifted = (lat_s.round(3) != lat_mode) | (lon_s.round(3) != lon_mode)
    if shifted.any():
        shifted_years = sorted(set(shifted.index[shifted].year))
        print(f"  POSITION SHIFT detected in year(s): {shifted_years} "
              f"(relative to modal position {lat_mode}, {lon_mode})")
    else:
        print(f"  No position shift detected - single fixed location for the whole record.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data_multiyear", type=Path)
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    args = parser.parse_args()

    nc_path = args.data_dir / f"{args.buoy}.nc"
    if not nc_path.exists():
        # Multi-year history download uses NO_TS_MO_<name>.nc per the
        # session handoff's naming-convention fix - try that before
        # crashing on a guessed filename.
        alt = args.data_dir / f"NO_TS_MO_{args.buoy}.nc"
        if alt.exists():
            nc_path = alt
        else:
            raise FileNotFoundError(f"Neither {nc_path} nor {alt} exists.")

    with xr.open_dataset(nc_path) as ds:
        print(f"=== {nc_path} ===\n")

        print("--- Global attributes ---")
        describe_attrs(ds.attrs)

        print(f"\n--- Dimensions ---")
        for d, size in ds.sizes.items():
            print(f"  {d}: {size}")

        print(f"\n--- All variables in file ---")
        for v in ds.variables:
            print(f"  {v}  dims={ds[v].dims}  dtype={ds[v].dtype}")

        time_name = resolve_coord_name(ds, "TIME")
        time = pd.to_datetime(ds[time_name].values)

        actual_var = pick_var(ds, VAR_CANDIDATES.get(args.var, [args.var]))
        if actual_var is None:
            print(f"\n{args.var} not found in this file. Available: {list(ds.variables)}")
            return

        print(f"\n--- {actual_var} variable attributes ---")
        describe_attrs(ds[actual_var].attrs)

        qc_flag_report(ds, actual_var, time)

        other_qc = [v for v in ds.variables if v.endswith("_QC") and v != f"{actual_var}_QC"]
        if other_qc:
            print(f"\n--- Other QC variables present in file (not inspected in "
                  f"detail here): {other_qc} ---")

        position_report(ds, time)


if __name__ == "__main__":
    main()
