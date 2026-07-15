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


def resolve_coord_name(ds, logical_name: str) -> str:
    """Find the actual variable/coord name in ds for a logical name like
    'TIME', 'LATITUDE', 'LONGITUDE' - case-insensitively. Exists because
    subset() (NRT) and get() (multi-year history) apparently return
    different capitalization for the same logical coordinates - rather
    than rename files by hand every download, resolve it at read time."""
    candidates = list(ds.variables) + list(ds.coords)
    for c in candidates:
        if c.upper() == logical_name.upper():
            return c
    raise KeyError(f"No variable matching '{logical_name}' (any case) found. "
                    f"Available: {list(ds.variables)}")


def get_scalar_latlon(ds):
    lat_name = resolve_coord_name(ds, "LATITUDE")
    lon_name = resolve_coord_name(ds, "LONGITUDE")
    lat = float(np.asarray(ds[lat_name].values).flat[0])
    lon = float(np.asarray(ds[lon_name].values).flat[0])
    return lat, lon


def load_buoy_series(nc_path: Path, varname: str = "VHM0") -> pd.Series:
    """Load a single variable from a buoy .nc file as a pandas Series indexed by time."""
    with xr.open_dataset(nc_path) as ds:
        actual_var = pick_var(ds, VAR_CANDIDATES.get(varname, [varname]))
        if actual_var is None:
            raise ValueError(f"{varname} not found in {nc_path.name}. "
                              f"Available: {list(ds.variables)}")
        time = pd.to_datetime(ds[resolve_coord_name(ds, "TIME")].values)
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
        time = pd.to_datetime(ds[resolve_coord_name(ds, "TIME")].values)
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
        time = pd.to_datetime(ds[resolve_coord_name(ds, "TIME")].values)
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


def longest_contiguous_segment(series: pd.Series):
    """Given a series that may contain NaN gaps (Stage 0 leaves long gaps
    unfilled rather than bridging them), return the longest contiguous
    non-NaN run, plus metadata about how much of the record it represents.

    Exists because lag/order-based methods (ACF, differencing, anything
    that assumes consecutive samples are actually temporally adjacent)
    are invalid across a spliced gap - a naive .dropna() silently
    concatenates the pre-gap and post-gap periods as if they were
    adjacent in time, which corrupts cumulative statistics like ACF far
    more than it corrupts a single point estimate. Value-only stages
    (distribution fits, EVA peak-finding) don't need this - order doesn't
    matter for those.
    """
    is_valid = series.notna()
    if not is_valid.any():
        empty = series.iloc[0:0]
        return empty, {"n_segments": 0, "segment_length": 0,
                        "total_valid_samples": 0, "pct_of_valid_used": 0.0}

    seg_id = (is_valid != is_valid.shift()).cumsum()
    segments = [group for _, group in series[is_valid].groupby(seg_id[is_valid])]

    longest = max(segments, key=len)
    total_valid = int(is_valid.sum())
    meta = {
        "n_segments": len(segments),
        "segment_length": len(longest),
        "total_valid_samples": total_valid,
        "pct_of_valid_used": round(100 * len(longest) / total_valid, 1),
        "segment_start": str(longest.index[0]),
        "segment_end": str(longest.index[-1]),
    }
    return longest, meta


def all_contiguous_segments(series: pd.Series, min_length: int = 1):
    """Like longest_contiguous_segment but returns every contiguous
    non-NaN run at least min_length samples long, sorted longest-first.

    Use when a statistic should aggregate across all usable stretches of
    a fragmented record rather than discarding everything except the
    single longest segment - e.g. a 36-year buoy record with hundreds of
    short outages scattered throughout can have its longest single
    segment cover under 10% of total valid data, which throws away real
    information a persistence estimate could otherwise use.
    """
    is_valid = series.notna()
    if not is_valid.any():
        return []
    seg_id = (is_valid != is_valid.shift()).cumsum()
    segments = [group for _, group in series[is_valid].groupby(seg_id[is_valid])]
    segments = [s for s in segments if len(s) >= min_length]
    segments.sort(key=len, reverse=True)
    return segments


def segments_by_time_gap(series: pd.Series, dt_hours: float, gap_multiplier: float = 1.5):
    """Split a time-indexed series into contiguous segments based on
    actual elapsed time between consecutive index entries - for series
    where gap positions are entirely ABSENT (missing rows) rather than
    NaN-marked on a regular grid. Stage 10's regime labels are like this:
    written only for valid timestamps, so the row-to-row index gap
    itself is the only signal a real gap happened.

    A break is detected wherever consecutive timestamps are more than
    gap_multiplier * dt_hours apart (default 1.5x the nominal sampling
    interval - generous enough to not trip on normal jitter, tight
    enough to catch a real gap).
    """
    if len(series) < 2:
        return [series] if len(series) else []
    idx = series.index
    deltas_hours = (idx[1:] - idx[:-1]).total_seconds() / 3600.0
    threshold = dt_hours * gap_multiplier
    break_points = np.where(deltas_hours > threshold)[0] + 1
    segments = []
    start = 0
    for bp in break_points:
        segments.append(series.iloc[start:bp])
        start = bp
    segments.append(series.iloc[start:])
    return segments


def load_era5_for_buoy(lat: float, lon: float, era5_dir: str = "meteo_era5"):
    """Load ERA5 meteo for the grid cell nearest a buoy's location,
    concatenated across every monthly file in era5_dir. Returns a
    DataFrame indexed by time with u10, v10, msl, t2m, sst (raw ERA5
    names) plus derived wind_speed, wind_dir_from_deg, and
    air_sea_temp_diff_c.

    Wind direction convention (easy to get backwards, so stated
    explicitly): meteorological "FROM" convention - the compass bearing
    the wind is blowing FROM, not the direction the wind vector points
    TOWARD. u10/v10 give the vector the wind blows TOWARD (eastward/
    northward components), so:
        bearing_toward = atan2(u10, v10)   (bearing from north)
        wind_dir_from = (bearing_toward + 180) mod 360
    """
    era5_files = sorted(Path(era5_dir).glob("era5_belgium_*.nc"))
    if not era5_files:
        raise FileNotFoundError(f"No era5_belgium_*.nc files found in {era5_dir}")

    frames = []
    for fp in era5_files:
        with xr.open_dataset(fp) as ds:
            time_name = "time" if "time" in ds.variables or "time" in ds.dims else "valid_time"
            lat_name = resolve_coord_name(ds, "latitude")
            lon_name = resolve_coord_name(ds, "longitude")

            point = ds.sel({lat_name: lat, lon_name: lon}, method="nearest")

            data = {}
            for var in ["u10", "v10", "msl", "t2m", "sst"]:
                if var in point.variables:
                    data[var] = point[var].values

            time_vals = pd.to_datetime(point[time_name].values)
            df = pd.DataFrame(data, index=time_vals)
            frames.append(df)

    era5 = pd.concat(frames).sort_index()
    era5 = era5[~era5.index.duplicated(keep="first")]

    if "u10" in era5.columns and "v10" in era5.columns:
        era5["wind_speed"] = np.sqrt(era5["u10"] ** 2 + era5["v10"] ** 2)
        bearing_toward = np.degrees(np.arctan2(era5["u10"], era5["v10"]))
        era5["wind_dir_from_deg"] = (bearing_toward + 180) % 360

    if "t2m" in era5.columns and "sst" in era5.columns:
        era5["air_sea_temp_diff_c"] = era5["t2m"] - era5["sst"]  # Kelvin difference == Celsius difference

    return era5


def integral_timescale(series, dt_hours: float, max_lag: int, consecutive: int):
    """ACF-based integral (decorrelation) timescale, with a significance-
    band cutoff instead of a single zero-crossing (avoids both spurious
    early crossings and silently reporting a search-ceiling lower bound
    as a real measurement). Moved here from 11b_dependence_structure.py
    so other analyses can reuse it (e.g. checking ERA5 wind's own
    persistence, not just Hs) without duplicating the logic - Python
    module names can't start with a digit, so direct import from
    "11b_dependence_structure" isn't possible."""
    from statsmodels.tsa.stattools import acf
    n = len(series)
    max_lag = min(max_lag, n // 3) if n >= 9 else max(1, n - 1)
    rho = acf(series, nlags=max_lag, fft=True)
    band = 1.96 / np.sqrt(n)

    criterion_lag = None
    for k in range(1, len(rho) - consecutive + 1):
        if np.all(np.abs(rho[k:k + consecutive]) < band):
            criterion_lag = k
            break
    hit_ceiling = criterion_lag is None
    if hit_ceiling:
        criterion_lag = len(rho) - 1

    tau_hours = dt_hours * (1 + 2 * np.sum(rho[1:criterion_lag]))
    return tau_hours, criterion_lag, rho, band, hit_ceiling


def default_paths(stage_out: str):
    """Standard output dir for a pipeline stage."""
    out = Path("pipeline_out") / stage_out
    out.mkdir(parents=True, exist_ok=True)
    return out
