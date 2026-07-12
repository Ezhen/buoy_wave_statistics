"""
Download ERA5 reanalysis meteo (wind, pressure, temperature) over the
Belgian coastal zone, to pair with the wave buoy data.

Variables pulled (see discussion of why each matters):
  - 10m u/v wind components -> speed + direction, the dominant physical
    driver of wave generation
  - mean sea-level pressure -> leading indicator, pressure drops precede
    wind speed increases by hours; useful if a forecasting stage wants an
    early predictor rather than a contemporaneous one
  - 2m air temperature + sea surface temperature -> their DIFFERENCE
    affects atmospheric stability / drag coefficient, a second-order
    correction on wind-driven wave growth, not a standalone driver

REVISED after hitting CDS's request-size limit on a single full year of
hourly data: now chunks by MONTH (not year) and defaults to 3-HOURLY
sampling instead of hourly - Stage 11b already established storm
persistence on the order of 50-100h, so 3-hourly (8 samples/day)
resolves that scale comfortably without requesting resolution this
pipeline has no current use for. Together this is roughly a 96x
reduction in request size vs. the first attempt. Use --hour-step 1 if
you later need hourly for a specific short window.

RECOMMENDATION: don't request the full historical range blindly. Run
--start-year --end-year for 1-2 recent years first to validate
credentials/variable names/queue behavior, THEN launch the full range
(matched to whatever the multi-year buoy download actually covers) as a
background job once you know real per-buoy coverage.

Setup (one-time):
    pip install cdsapi
    # ~/.cdsapirc with your CDS (Climate Data Store) credentials -
    # cds.climate.copernicus.eu, a THIRD separate Copernicus portal/account,
    # distinct from both CDSE and CMEMS.

Usage:
    python download_era5_meteo.py --start-year 2025 --end-year 2025
    python download_era5_meteo.py --start-year 1990 --end-year 2026 --hour-step 3
"""

import argparse
from pathlib import Path

import cdsapi
import xarray as xr

OUTPUT_DIR = Path("meteo_era5")
OUTPUT_DIR.mkdir(exist_ok=True)

# [North, West, South, East] - ERA5's area format, NOT the same order as
# the CMEMS bbox dict used in the buoy download scripts. Slightly padded
# beyond the buoy bounding box so coastal buoys aren't sitting right on
# a grid edge.
AREA = [51.7, 2.0, 51.0, 3.6]

VARIABLES = [
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "mean_sea_level_pressure",
    "2m_temperature",
    "sea_surface_temperature",
]

DAYS = [f"{d:02d}" for d in range(1, 32)]
MONTHS = [f"{m:02d}" for m in range(1, 13)]


def hours_list(step: int):
    return [f"{h:02d}:00" for h in range(0, 24, step)]


def download_month(client: cdsapi.Client, year: int, month: str, area, variables, hours):
    out_path = OUTPUT_DIR / f"era5_belgium_{year}_{month}.nc"
    if out_path.exists():
        print(f"{out_path} already exists - skipping (delete it to re-download)")
        return out_path

    print(f"Requesting ERA5 {year}-{month}...")
    client.retrieve(
        "reanalysis-era5-single-levels",
        {
            "product_type": "reanalysis",
            "format": "netcdf",
            "variable": variables,
            "year": str(year),
            "month": month,
            "day": DAYS,
            "time": hours,
            "area": area,
        },
        str(out_path),
    )
    print(f"Saved: {out_path}")
    return out_path


def report_coverage(paths):
    print(f"\n{'=' * 70}\nDownloaded chunks\n{'=' * 70}")
    for p in sorted(paths):
        if not p.exists():
            continue
        try:
            with xr.open_dataset(p) as ds:
                time_dim = "time" if "time" in ds.dims else "valid_time"
                n = ds.dims.get(time_dim, "?")
                print(f"  {p.name}: {n} timesteps")
        except Exception as e:
            print(f"  {p.name}: could not read ({e})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start-year", required=True, type=int)
    parser.add_argument("--end-year", required=True, type=int)
    parser.add_argument("--hour-step", default=3, type=int,
                         help="hours between samples, e.g. 3 = 3-hourly (default). "
                              "Use 1 for full hourly if a specific window needs it - "
                              "will likely need --start-year==--end-year and a single "
                              "month at a time to stay under the request-size limit.")
    args = parser.parse_args()

    n_years = args.end_year - args.start_year + 1
    if n_years > 10:
        print(f"WARNING: requesting {n_years} years ({n_years * 12} monthly chunks). "
              f"Each is its own CDS queue request - consider running this as a "
              f"background/detached job (nohup, screen, or your HPC scheduler) "
              f"rather than interactively.")

    client = cdsapi.Client()
    hours = hours_list(args.hour_step)

    downloaded = []
    for year in range(args.start_year, args.end_year + 1):
        for month in MONTHS:
            path = download_month(client, year, month, AREA, VARIABLES, hours)
            downloaded.append(path)

    report_coverage(downloaded)


if __name__ == "__main__":
    main()
