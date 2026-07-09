"""
Plot each Belgian coastal buoy's wave time series, plus a combined map of
buoy locations. Run this on your HPC environment where the .nc files live.

Usage:
    python plot_belgian_buoys.py --data-dir data --out-dir plots

Requires: xarray, netCDF4, matplotlib, pandas, numpy
    pip install xarray netCDF4 matplotlib pandas numpy
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Variables of interest, with fallback names (CMEMS in-situ files vary slightly)
VAR_CANDIDATES = {
    "VHM0": ["VHM0"],                  # significant wave height
    "VTPK": ["VTPK"],                  # peak period
    "VTM02": ["VTM02"],                # mean period
    "VMDR": ["VMDR"],                  # mean wave direction
}


def pick_var(ds, candidates):
    for name in candidates:
        if name in ds.variables:
            return name
    return None


def get_scalar_latlon(ds):
    lat = float(np.asarray(ds["latitude"].values).flat[0])
    lon = float(np.asarray(ds["longitude"].values).flat[0])
    return lat, lon


def plot_one_buoy(nc_path: Path, out_dir: Path):
    name = nc_path.stem
    with xr.open_dataset(nc_path) as ds:
        try:
            lat, lon = get_scalar_latlon(ds)
        except Exception:
            lat, lon = np.nan, np.nan

        # Find a time coordinate
        time_dim = "time" if "time" in ds.variables else None
        if time_dim is None:
            print(f"[skip] {name}: no TIME variable found")
            return None

        time = pd.to_datetime(ds["time"].values)

        present_vars = {k: pick_var(ds, v) for k, v in VAR_CANDIDATES.items()}
        present_vars = {k: v for k, v in present_vars.items() if v is not None}

        if not present_vars:
            print(f"[skip] {name}: none of {list(VAR_CANDIDATES)} found")
            return None

        n = len(present_vars)
        fig, axes = plt.subplots(n, 1, figsize=(11, 2.5 * n), sharex=True)
        if n == 1:
            axes = [axes]

        for ax, (label, varname) in zip(axes, present_vars.items()):
            data = ds[varname].values
            # collapse extra dims (e.g. DEPTH) if present, take first level
            if data.ndim > 1:
                data = data.reshape(data.shape[0], -1)[:, 0]
            ax.plot(time, data, lw=0.8)
            units = ds[varname].attrs.get("units", "")
            ax.set_ylabel(f"{label}\n[{units}]")
            ax.grid(alpha=0.3)

        axes[-1].set_xlabel("Time")
        fig.suptitle(f"{name}  (lat={lat:.3f}, lon={lon:.3f})")
        fig.tight_layout()
        out_path = out_dir / f"{name}_timeseries.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"[ok] {name}: saved {out_path}")

        return {"name": name, "lat": lat, "lon": lon}


def plot_map(locations, out_dir: Path):
    locations = [loc for loc in locations if loc is not None and not np.isnan(loc["lat"])]
    if not locations:
        print("[map] no valid locations to plot")
        return

    fig, ax = plt.subplots(figsize=(8, 9))
    lats = [loc["lat"] for loc in locations]
    lons = [loc["lon"] for loc in locations]
    ax.scatter(lons, lats, c="crimson", s=40, zorder=3)

    for loc in locations:
        ax.annotate(
            loc["name"].replace("Buoy", "").replace(".nc", ""),
            (loc["lon"], loc["lat"]),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=8,
        )

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Belgian coastal zone buoy locations")
    ax.grid(alpha=0.3)
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    out_path = out_dir / "buoy_locations_map.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[ok] map: saved {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data", type=Path)
    parser.add_argument("--out-dir", default="plots", type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(exist_ok=True, parents=True)

    nc_files = sorted(args.data_dir.glob("*.nc"))
    if not nc_files:
        print(f"No .nc files found in {args.data_dir}")
        return

    locations = []
    for nc_path in nc_files:
        loc = plot_one_buoy(nc_path, args.out_dir)
        locations.append(loc)

    plot_map(locations, args.out_dir)


if __name__ == "__main__":
    main()
