"""
Download in-situ wave buoy time series for the Belgian Coastal Zone from CMEMS.

Requires: pip install copernicusmarine
Auth: copernicusmarine login   (uses your data.marine.copernicus.eu account,
      separate from any CDSE credentials)
"""

import copernicusmarine
from pathlib import Path

OUTPUT_DIR = Path("data")
OUTPUT_DIR.mkdir(exist_ok=True)

# North West Shelf regional in-situ product - covers the Belgian coast.
# (The global product's regional coverage there is thinner; NWS is the right one.)
PRODUCT_ID = "INSITU_NWS_PHYBGCWAV_DISCRETE_MYNRT_013_036"

BBOX = dict(
    minimum_longitude=2.2,
    maximum_longitude=3.4,
    minimum_latitude=51.1,
    maximum_latitude=51.6,
)

VARIABLES = ["VHM0", "VTPK", "VTM02", "VMDR"]  # Hs, peak period, mean period, mean direction


def find_hourly_wave_dataset_id(product_id: str) -> str:
    """Look up the actual dataset id instead of hardcoding a guess."""
    catalogue = copernicusmarine.describe(product_id=product_id, contains=[product_id])
    candidates = []
    for product in catalogue.products:
        for dataset in product.datasets:
            candidates.append(dataset.dataset_id)

    if not candidates:
        raise RuntimeError(
            f"No datasets found under {product_id}. "
            "Browse https://data.marine.copernicus.eu and search the product id "
            "manually to confirm it's still current."
        )

    # Prefer the hourly ("PT1H" / "history") time series over monthly-file layers
    preferred = [c for c in candidates if "PT1H" in c or "history" in c.lower()]
    chosen = preferred[0] if preferred else candidates[0]
    print(f"Found {len(candidates)} dataset(s) under {product_id}:")
    for c in candidates:
        marker = "  <-- using this one" if c == chosen else ""
        print(f"  - {c}{marker}")
    return chosen


def main():
    dataset_id = find_hourly_wave_dataset_id(PRODUCT_ID)

    result = copernicusmarine.subset(
        dataset_id=dataset_id,
        variables=VARIABLES,
        minimum_longitude=BBOX["minimum_longitude"],
        maximum_longitude=BBOX["maximum_longitude"],
        minimum_latitude=BBOX["minimum_latitude"],
        maximum_latitude=BBOX["maximum_latitude"],
        output_directory=str(OUTPUT_DIR),
        output_filename="belgian_coastal_wave_buoys.nc",
    )
    print("Saved:", result)


if __name__ == "__main__":
    main()
