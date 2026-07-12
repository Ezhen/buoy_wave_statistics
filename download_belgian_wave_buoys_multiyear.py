"""
Pull the FULL available history for the Belgian coastal wave buoys.

CORRECTION from the first version of this script: there is no separate
MY (delayed-mode-only) product for the North West Shelf region. Unlike
the global product, which splits into INSITU_GLO_PHYBGCWAV_DISCRETE_
MYNRT_013_030 (NRT) and a distinct INSITU_GLO_WAV_DISCRETE_MY_013_045
(delayed-mode), NWS has a single COMBINED product -
INSITU_NWS_PHYBGCWAV_DISCRETE_MYNRT_013_036 - the same one the original
download_belgian_wave_buoys.py already used. There's no different
product to switch to.

What was actually different: the original script's subset() call didn't
specify a date range, so it likely defaulted to a recent window (~2
months). This version explicitly requests a wide historical range
instead, on the SAME product, to get whatever multi-year history it
actually holds per buoy - which may still turn out to be uneven or
short for some buoys; that's a real finding to check via
report_actual_coverage() below, not something to assume either way
beforehand.

Saves to data_multiyear/ (NOT data/), so the existing NRT-window pipeline
run and its outputs stay untouched and comparable.

Setup (one-time, same account as the original script):
    pip install copernicusmarine
    copernicusmarine login

Usage:
    python download_belgian_wave_buoys_multiyear.py
"""

import copernicusmarine
import xarray as xr
from pathlib import Path

OUTPUT_DIR = Path("data_multiyear")
OUTPUT_DIR.mkdir(exist_ok=True)

# Same product as the original NRT-window script - NWS doesn't have a
# separate MY product to switch to (see correction note above).
PRODUCT_ID = "INSITU_NWS_PHYBGCWAV_DISCRETE_MYNRT_013_036"

BBOX = dict(
    minimum_longitude=2.2,
    maximum_longitude=3.4,
    minimum_latitude=51.1,
    maximum_latitude=51.6,
)

VARIABLES = ["VHM0", "VTPK", "VTM02", "VMDR"]

# Deliberately far back - if a buoy's history doesn't go this far, the
# product simply won't return data before its actual deployment start.
# No harm in asking too early; there IS harm in asking too late and
# silently truncating real history, which is what likely happened
# without any start_datetime at all.
START_DATETIME = "1990-01-01T00:00:00"


def find_wave_dataset_id(product_id: str) -> str:
    """Look up the actual dataset id instead of hardcoding a guess -
    same reasoning as the original NRT script, after the earlier
    dataset-id mixup."""
    catalogue = copernicusmarine.describe(product_id=product_id, contains=[product_id])
    candidates = []
    for product in catalogue.products:
        for dataset in product.datasets:
            candidates.append(dataset.dataset_id)

    if not candidates:
        raise RuntimeError(
            f"No datasets found under {product_id}. Browse "
            f"https://data.marine.copernicus.eu and search 'INSITU_NWS_PHYBGCWAV' "
            f"to confirm the current product id before assuming the script is wrong."
        )

    preferred = [c for c in candidates if "history" in c.lower()]
    chosen = preferred[0] if preferred else candidates[0]
    print(f"Found {len(candidates)} dataset(s) under {product_id}:")
    for c in candidates:
        marker = "  <-- using this one" if c == chosen else ""
        print(f"  - {c}{marker}")
    return chosen


def report_actual_coverage(output_dir: Path):
    """After downloading, check what date range each buoy's file actually
    contains - this is the real answer to 'how many years did we get',
    more reliable than trusting catalogue metadata per platform."""
    nc_files = sorted(output_dir.glob("*.nc"))
    if not nc_files:
        print("No .nc files found after download - something went wrong upstream.")
        return

    print(f"\n{'='*70}\nActual per-buoy coverage after download\n{'='*70}")
    rows = []
    for nc_path in nc_files:
        try:
            with xr.open_dataset(nc_path) as ds:
                time = ds["TIME"].values
                start, end = time.min(), time.max()
                n = len(time)
                rows.append((nc_path.stem, str(start)[:10], str(end)[:10], n))
        except Exception as e:
            rows.append((nc_path.stem, "ERROR", str(e), 0))

    for name, start, end, n in sorted(rows):
        print(f"  {name:35s} {start} -> {end}  (n={n})")

    print("\nCheck the spread above before assuming uniform history - if some "
          "buoys are much shorter than others, that's real deployment history, "
          "not a download bug.")


def main():
    dataset_id = find_wave_dataset_id(PRODUCT_ID)

    print(f"\nSubsetting {dataset_id} over the Belgian coastal bounding box, "
          f"from {START_DATETIME} onward - explicitly wide, since the original "
          f"script's unspecified date range likely defaulted to a recent window.")

    result = copernicusmarine.subset(
        dataset_id=dataset_id,
        variables=VARIABLES,
        minimum_longitude=BBOX["minimum_longitude"],
        maximum_longitude=BBOX["maximum_longitude"],
        minimum_latitude=BBOX["minimum_latitude"],
        maximum_latitude=BBOX["maximum_latitude"],
        start_datetime=START_DATETIME,
        output_directory=str(OUTPUT_DIR),
    )
    print("Saved:", result)

    report_actual_coverage(OUTPUT_DIR)


if __name__ == "__main__":
    main()
