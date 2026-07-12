"""Shared helpers for the Belgian buoy characterization pipeline."""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr

VAR_CANDIDATES = {
    "VHM0": ["VHM0"],    # significant wave height (m)
    "VTPK": ["VTPK"],    # peak period (s)
    "VTM02": ["VTM02"],  # mean period (s)
    "VMDR": ["VMDR"],    # mean wave direction (deg)
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


def load_buoy_series(nc_path: Path, varname: str = "VHM0") -> pd.Series:
    """Load a single variable from a buoy .nc file as a pandas Series indexed by time."""
    with xr.open_dataset(nc_path) as ds:
        actual_var = pick_var(ds, VAR_CANDIDATES.get(varname, [varname]))
        if actual_var is None:
            raise ValueError(f"{varname} not found in {nc_path.name}. "
                              f"Available: {list(ds.variables)}")
        time = pd.to_datetime(ds["time"].values)
        data = ds[actual_var].values
        if data.ndim > 1:
            data = data.reshape(data.shape[0], -1)[:, 0]
        s = pd.Series(data, index=time, name=actual_var)
        s = s[~s.index.duplicated(keep="first")].sort_index()
        return s


def load_buoy_dataframe(nc_path: Path, varnames=("VHM0", "VTPK", "VTM02", "VMDR")) -> pd.DataFrame:
    """Load multiple variables from one buoy file, aligned by time. Silently
    skips any variable not present in this particular buoy's file."""
    with xr.open_dataset(nc_path) as ds:
        time = pd.to_datetime(ds["time"].values)
        cols = {}
        for v in varnames:
            actual = pick_var(ds, VAR_CANDIDATES.get(v, [v]))
            if actual is None:
                continue
            data = ds[actual].values
            if data.ndim > 1:
                data = data.reshape(data.shape[0], -1)[:, 0]
            cols[v] = data
        df = pd.DataFrame(cols, index=time)
        df = df[~df.index.duplicated(keep="first")].sort_index()
        return df


def buoy_name(nc_path: Path) -> str:
    return nc_path.stem


def haversine_km(lat1, lon1, lat2, lon2):
    import math
    R = 6371.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    return 2 * R * math.asin(min(1.0, np.sqrt(a)))


def count_raw_duplicate_timestamps(nc_path: Path) -> int:
    """Count duplicate TIME entries before load_buoy_series silently drops them."""
    with xr.open_dataset(nc_path) as ds:
        time = pd.to_datetime(ds["time"].values)
    return int(pd.Index(time).duplicated().sum())


def detect_available_variables(nc_path: Path, candidates=("VHM0", "VTPK", "VTM02", "VMDR")):
    """Which of the standard wave variables does this buoy's file actually have?
    Used by the tiering gate to decide which Advanced stages a buoy qualifies for."""
    with xr.open_dataset(nc_path) as ds:
        return [v for v in candidates if pick_var(ds, VAR_CANDIDATES.get(v, [v])) is not None]


def resolve_block_length(buoy: str, var: str, n_samples: int):
    """Same resolution Stage 12 uses: prefer Stage 11b's persistence-based
    block length, fall back to a sqrt(n) heuristic with a loud caveat."""
    dep_path = Path("pipeline_out/11b_dependence_structure") / f"{buoy}_{var}_dependence_summary.json"
    if dep_path.exists():
        with open(dep_path) as f:
            dep = json.load(f)
        return dep["suggested_block_length_samples"], "Stage 11b integral timescale", dep.get("hit_max_lag_ceiling", False)
    block_length = max(1, int(np.sqrt(n_samples)))
    return block_length, "sqrt(n) fallback - Stage 11b not found, this is a WEAK default", False


def default_paths(stage_out: str):
    """Standard output dir for a pipeline stage."""
    out = Path("pipeline_out") / stage_out
    out.mkdir(parents=True, exist_ok=True)
    return out
