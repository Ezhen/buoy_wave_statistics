"""
Download the FULL multi-year history for the Belgian coastal wave buoys.

Two prior corrections got us here:
  1. subset() only serves the "latest" 30-day window for in-situ products.
     The full archive is a different dataset PART - "history" - served
     through a different API: the Files service (copernicusmarine.get()).
  2. get()'s actual signature (confirmed via `inspect.signature` and
     `copernicusmarine get --help` on the real installed version, not
     guessed from documentation for a possibly-different release):
     show_outputnames doesn't exist; force_download is deprecated with
     no direct replacement. What DOES exist and is exactly what we need:
     `create_file_list` - "only create a file containing the names of
     the targeted files instead of downloading them." That's the preview
     step, confirmed from --help text, not inferred.

Two-step workflow, run as two separate commands (deliberately - don't
skip straight to downloading 21GB of the whole North West Shelf region
when we only need ~19 Belgian platforms out of it):

STEP 1 - list only, no download:
    python download_belgian_wave_buoys_history.py --list

Writes every matching filename in the "history" part to a .txt file.
Inspect it for the Belgian platforms (likely named by WMO/platform code,
not the friendly names like "WesthinderBuoy" - those were assigned by
subset()'s own metadata layer, not the raw file naming convention).

STEP 2 - real download, filtered to just the Belgian platforms:
    python download_belgian_wave_buoys_history.py --download --regex "PATTERN"

where PATTERN is built from what step 1's file list actually shows -
not guessed in advance. If you can't find an obvious pattern (e.g. no
shared prefix/country code), fall back to --file-list pointing at a
manually curated .txt of just the ~19 file paths you want, copied out
of the full listing - --file-list downloads exactly the paths given,
no pattern-matching needed.

Setup (one-time, same account as the other scripts):
    pip install copernicusmarine
    copernicusmarine login

Usage:
    python download_belgian_wave_buoys_history.py --list
    python download_belgian_wave_buoys_history.py --download --regex "PATTERN"
    python download_belgian_wave_buoys_history.py --download --file-list belgian_files.txt
"""

import argparse
import shutil
from pathlib import Path

import copernicusmarine
import xarray as xr

from utils import resolve_coord_name

RAW_DIR = Path("data_multiyear_raw")       # everything get() downloads, unfiltered
FILTERED_DIR = Path("data_multiyear")      # only the Belgian-coast subset, after lat/lon check
FILE_LIST_PATH = Path("history_file_list.txt")
RAW_DIR.mkdir(exist_ok=True)
FILTERED_DIR.mkdir(exist_ok=True)

DATASET_ID = "cmems_obs-ins_nws_phybgcwav_mynrt_na_irr"  # confirmed via an earlier run's own output
DATASET_PART = "history"

BBOX = dict(
    minimum_longitude=2.2,
    maximum_longitude=3.4,
    minimum_latitude=51.1,
    maximum_latitude=51.6,
)


def list_files():
    """Step 1: write matching filenames to a .txt, download nothing."""
    print(f"Listing files under dataset_part='{DATASET_PART}' for {DATASET_ID} "
          f"-> {FILE_LIST_PATH} (no data downloaded)")
    result = copernicusmarine.get(
        dataset_id=DATASET_ID,
        dataset_part=DATASET_PART,
        create_file_list=str(FILE_LIST_PATH),
    )
    print("get() result:", result)
    if FILE_LIST_PATH.exists():
        lines = FILE_LIST_PATH.read_text().splitlines()
        print(f"\n{len(lines)} file(s) listed. First 20:")
        for line in lines[:20]:
            print(f"  {line}")
        print(f"\nFull list in {FILE_LIST_PATH} - inspect it for the Belgian platform "
              f"naming pattern before building --regex or a --file-list subset.")


def download_filtered(regex, file_list):
    """Step 2: real download, scoped down via regex or an explicit file list."""
    if not regex and not file_list:
        print("Refusing to download the full 'history' part with no filter - "
              "that's the 21GB whole-region pull. Pass --regex or --file-list "
              "(built from step 1's output) to scope this down first.")
        return

    kwargs = dict(
        dataset_id=DATASET_ID,
        dataset_part=DATASET_PART,
        output_directory=str(RAW_DIR),
        no_directories=True,
    )
    if regex:
        kwargs["regex"] = regex
        print(f"Downloading files matching regex: {regex}")
    if file_list:
        kwargs["file_list"] = file_list
        print(f"Downloading exactly the paths listed in: {file_list}")

    result = copernicusmarine.get(**kwargs)
    print("get() result:", result)


def filter_to_belgian_coast():
    """Extra safety net: even after a scoped download, keep only files whose
    LATITUDE/LONGITUDE actually fall in the Belgian coastal bbox - copy
    rather than move, so a mistake here doesn't lose anything raw."""
    nc_files = sorted(RAW_DIR.glob("*.nc"))
    if not nc_files:
        print("No .nc files in the raw download dir yet - run --download first.")
        return
    print(f"\n{len(nc_files)} file(s) downloaded - checking coordinates...")

    kept = []
    for nc_path in nc_files:
        try:
            with xr.open_dataset(nc_path) as ds:
                lat = float(ds[resolve_coord_name(ds, "LATITUDE")].values.flat[0])
                lon = float(ds[resolve_coord_name(ds, "LONGITUDE")].values.flat[0])
            if (BBOX["minimum_latitude"] <= lat <= BBOX["maximum_latitude"]
                    and BBOX["minimum_longitude"] <= lon <= BBOX["maximum_longitude"]):
                # Native history filenames look like NO_TS_MO_WesthinderBuoy.nc -
                # strip the network/platform-type prefix so the rest of the
                # pipeline (which expects "<buoy>.nc") can find it without
                # a manual rename step.
                clean_name = nc_path.name
                if "_" in clean_name:
                    clean_name = clean_name.split("_")[-1]
                dest = FILTERED_DIR / clean_name
                shutil.copy2(nc_path, dest)
                kept.append((clean_name, lat, lon))
        except Exception as e:
            print(f"  Could not read {nc_path.name}: {e}")

    print(f"\n{len(kept)} file(s) fall inside the Belgian coastal bbox, copied to {FILTERED_DIR}/:")
    for name, lat, lon in kept:
        print(f"  {name}  (lat={lat:.3f}, lon={lon:.3f})")

    if not kept:
        print("\nWARNING: nothing matched the bbox - check the regex/file-list actually "
              "targeted Belgian platforms, not some other region.")


def report_coverage():
    nc_files = sorted(FILTERED_DIR.glob("*.nc"))
    if not nc_files:
        return
    print(f"\n{'=' * 70}\nActual per-buoy coverage (filtered set)\n{'=' * 70}")
    for nc_path in nc_files:
        try:
            with xr.open_dataset(nc_path) as ds:
                time = ds[resolve_coord_name(ds, "TIME")].values
                print(f"  {nc_path.name:45s} {str(time.min())[:10]} -> "
                      f"{str(time.max())[:10]}  (n={len(time)})")
        except Exception as e:
            print(f"  {nc_path.name}: could not read ({e})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--list", action="store_true",
                         help="Step 1: write matching filenames to a .txt, download nothing.")
    parser.add_argument("--download", action="store_true",
                         help="Step 2: actually download, requires --regex or --file-list.")
    parser.add_argument("--regex", default=None)
    parser.add_argument("--file-list", default=None)
    args = parser.parse_args()

    if args.list:
        list_files()
        return

    if args.download:
        download_filtered(args.regex, args.file_list)
        filter_to_belgian_coast()
        report_coverage()
        return

    print("Specify --list (step 1, preview filenames) or --download (step 2, "
          "requires --regex or --file-list). See the script's docstring.")


if __name__ == "__main__":
    main()
