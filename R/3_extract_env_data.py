#!/usr/bin/env python3
"""
Script 3: Environmental data extraction — T. teniotis fall migration
Via Open-Meteo Historical Weather API (ERA5 backend)

WHY OPEN-METEO OVER CDS/cdsapi:
  - No API key, no account, no queue
  - No local NetCDF storage (terabytes for ERA5 GRIB files)
  - Point queries only: exactly the coordinates and dates you need
  - ERA5-seamless: 9 km resolution from 2017+, consistent ERA5 back to 1940
  - Free for non-commercial research use

VARIABLES EXTRACTED (matching noctule pipeline script 3/8):
  Hourly surface:
    wind_u_component (u10), wind_v_component (v10) at 10 m and 100 m
    surface_pressure, mean_sea_level_pressure
    temperature_2m, dewpoint_2m, relative_humidity_2m
    precipitation, boundary_layer_height (PBL)
    is_day (for dusk/dawn context)

  Derived (computed here):
    wind_speed_10m, wind_speed_100m
    wind_direction_10m, wind_direction_100m
    tailwind (component along the bat's bearing)
    crosswind (perpendicular component)
    pressure_tendency (1-hr delta MSLP — migration trigger proxy)

TIME WINDOW STRATEGY:
  T. teniotis migrates nocturnally. We want atmospheric conditions at dusk
  (departure window) and during the night. For each daily bat position we
  extract the 3-hour block around sunset ± 1 hr and the full night
  (sunset to sunrise next day). Summary statistics (mean, max, range) are
  computed per window.

INPUT:
  CSV exported from R bats_daily with columns:
    individual_local_identifier, timestamp (UTC), lon, lat,
    behavior, daily_displacement_km, bearing (optional)

OUTPUT:
  CSV: bats_daily_env.csv  — one row per individual-day with env columns appended
"""

import math
import time
import warnings
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

# ── Configuration ─────────────────────────────────────────────────────────────

BASE_URL   = "https://archive-api.open-meteo.com/v1/archive"
MODEL      = "era5_seamless"   # ERA5 + ERA5-Land seamless; use "era5" for strict consistency
TIMEZONE   = "UTC"
MAX_RETRIES = 4
RETRY_DELAY = 5   # seconds between retries (rate-limit friendly)
BATCH_DELAY = 0.2 # seconds between API calls (Open-Meteo: ~10 k free requests/day)

# Variables to request from Open-Meteo
HOURLY_VARS = [
    "temperature_2m",
    "relative_humidity_2m",
    "dewpoint_2m",
    "surface_pressure",
    "pressure_msl",           # mean sea level pressure
    "wind_speed_10m",
    "wind_speed_100m",
    "wind_direction_10m",
    "wind_direction_100m",
    "precipitation",
    "boundary_layer_height",  # PBL — key for assessing soaring/gliding conditions
    "is_day",                 # 0/1 flag — helps isolate nocturnal flight hours
    "cloud_cover",
]

# File paths — adjust to your Dropbox structure
INPUT_CSV   = Path("../../../Dropbox/MPI/Tadarida/Data/bats_daily_export.csv")
OUTPUT_CSV  = Path("../../../Dropbox/MPI/Tadarida/Data/bats_daily_env.csv")
CACHE_DIR   = Path("../../../Dropbox/MPI/Tadarida/Data/env_cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# ── Helpers ───────────────────────────────────────────────────────────────────

def cache_path(lat: float, lon: float, date_str: str) -> Path:
    """One JSON cache file per (lat, lon, date) rounded to 2 decimal places."""
    return CACHE_DIR / f"{lat:.2f}_{lon:.2f}_{date_str}.json"


def fetch_hourly(lat: float, lon: float, date_str: str) -> pd.DataFrame | None:
    """
    Fetch one day of hourly ERA5 data centred on `date_str` (YYYY-MM-DD).
    Returns a DataFrame with columns: time + all HOURLY_VARS.
    Returns None on persistent failure.
    """
    cp = cache_path(round(lat, 2), round(lon, 2), date_str)
    if cp.exists():
        return pd.read_json(cp)

    # Fetch date ± 1 day to allow full nocturnal window spanning midnight
    d      = datetime.strptime(date_str, "%Y-%m-%d")
    d_prev = (d - timedelta(days=1)).strftime("%Y-%m-%d")
    d_next = (d + timedelta(days=1)).strftime("%Y-%m-%d")

    params = {
        "latitude":   round(lat, 4),
        "longitude":  round(lon, 4),
        "start_date": d_prev,
        "end_date":   d_next,
        "hourly":     ",".join(HOURLY_VARS),
        "models":     MODEL,
        "timezone":   TIMEZONE,
        "wind_speed_unit": "ms",   # m/s matches ERA5 native; convert later if needed
    }

    for attempt in range(MAX_RETRIES):
        try:
            r = requests.get(BASE_URL, params=params, timeout=30)
            r.raise_for_status()
            data = r.json()
            df = pd.DataFrame({"time": data["hourly"]["time"]})
            for v in HOURLY_VARS:
                df[v] = data["hourly"].get(v, np.nan)
            df["time"] = pd.to_datetime(df["time"], utc=True)
            df.to_json(cp)
            return df
        except (requests.RequestException, KeyError) as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                warnings.warn(f"Failed to fetch {date_str} @ ({lat:.3f},{lon:.3f}): {e}")
                return None


def wind_components(speed: pd.Series, direction: pd.Series):
    """
    Convert meteorological wind (speed m/s, direction degrees from N)
    to U (eastward) and V (northward) components.
    Meteorological convention: direction = direction *from* which wind blows.
    """
    dir_rad = np.deg2rad(direction)
    u = -speed * np.sin(dir_rad)
    v = -speed * np.cos(dir_rad)
    return u, v


def tailwind_crosswind(wind_speed: pd.Series,
                       wind_dir: pd.Series,
                       bat_bearing: float) -> tuple[pd.Series, pd.Series]:
    """
    Decompose wind into tailwind (positive = helpful) and crosswind components
    relative to the bat's movement bearing.

    bat_bearing: direction of travel in degrees from N (0–360).
    """
    u_wind, v_wind = wind_components(wind_speed, wind_dir)
    # Unit vector of bat travel direction
    bat_rad  = np.deg2rad(bat_bearing)
    u_bat    =  np.sin(bat_rad)   # eastward component of travel
    v_bat    =  np.cos(bat_rad)   # northward component of travel
    # Tailwind: dot product of wind vector with travel unit vector
    tailwind  = u_wind * u_bat + v_wind * v_bat
    # Crosswind: magnitude of cross product (scalar in 2D)
    crosswind = u_wind * v_bat - v_wind * u_bat
    return tailwind, crosswind


def pressure_tendency(df_hourly: pd.DataFrame,
                      target_hours: list[int],
                      delta_h: int = 3) -> float:
    """
    3-hour pressure tendency at dusk: MSLP(t) - MSLP(t - delta_h).
    Falling pressure → approaching fronts → often precedes migration.
    """
    sub = df_hourly[df_hourly["time"].dt.hour.isin(target_hours)].copy()
    if sub.empty or len(sub) < 2:
        return np.nan
    sub = sub.sort_values("time").reset_index(drop=True)
    # Use the last available dusk hour pair
    return float(sub["pressure_msl"].iloc[-1] - sub["pressure_msl"].iloc[0])


def summarise_window(df_h: pd.DataFrame,
                     window_mask: pd.Series,
                     bat_bearing: float | None,
                     prefix: str) -> dict:
    """
    Compute summary statistics for a subset of hourly data (defined by window_mask).
    Returns a flat dict with keys like `{prefix}_mean_wind_speed_10m` etc.
    """
    sub = df_h[window_mask].copy()
    out = {}

    if sub.empty:
        for v in HOURLY_VARS + ["tailwind_10m", "crosswind_10m",
                                 "tailwind_100m", "crosswind_100m",
                                 "pressure_tendency_3h"]:
            for stat in ["mean", "max", "min", "range"]:
                out[f"{prefix}_{stat}_{v}"] = np.nan
        return out

    numeric_vars = [v for v in HOURLY_VARS if v != "is_day"]

    for v in numeric_vars:
        col = sub[v].dropna()
        out[f"{prefix}_mean_{v}"]  = col.mean()
        out[f"{prefix}_max_{v}"]   = col.max()
        out[f"{prefix}_min_{v}"]   = col.min()
        out[f"{prefix}_range_{v}"] = col.max() - col.min()

    # Tailwind / crosswind (only if bearing known)
    for height, spd_col, dir_col in [
        ("10m",  "wind_speed_10m",  "wind_direction_10m"),
        ("100m", "wind_speed_100m", "wind_direction_100m"),
    ]:
        if bat_bearing is not None and not sub[spd_col].isna().all():
            tw, cw = tailwind_crosswind(sub[spd_col], sub[dir_col], bat_bearing)
            out[f"{prefix}_mean_tailwind_{height}"]  = float(tw.mean())
            out[f"{prefix}_mean_crosswind_{height}"] = float(cw.mean())
        else:
            out[f"{prefix}_mean_tailwind_{height}"]  = np.nan
            out[f"{prefix}_mean_crosswind_{height}"] = np.nan

    # Pressure tendency across the window
    if len(sub) >= 2:
        mslp_sorted = sub.sort_values("time")["pressure_msl"].dropna()
        out[f"{prefix}_pressure_tendency"] = (
            float(mslp_sorted.iloc[-1] - mslp_sorted.iloc[0])
            if len(mslp_sorted) >= 2 else np.nan
        )
    else:
        out[f"{prefix}_pressure_tendency"] = np.nan

    return out


# ── Main extraction loop ──────────────────────────────────────────────────────

def main():
    print(f"Reading input: {INPUT_CSV}")
    bats = pd.read_csv(INPUT_CSV, parse_dates=["timestamp"])

    # Ensure timestamp is UTC-aware
    if bats["timestamp"].dt.tz is None:
        bats["timestamp"] = bats["timestamp"].dt.tz_localize("UTC")

    # Extract coordinates — handle sf WKT geometry column if present
    if "geometry" in bats.columns and "lon" not in bats.columns:
        from shapely import wkt
        bats["geometry"] = bats["geometry"].apply(wkt.loads)
        bats["lon"] = bats["geometry"].apply(lambda g: g.x)
        bats["lat"] = bats["geometry"].apply(lambda g: g.y)

    required = {"individual_local_identifier", "timestamp", "lon", "lat"}
    missing  = required - set(bats.columns)
    if missing:
        raise ValueError(f"Input CSV missing columns: {missing}")

    print(f"  {len(bats)} individual-days across "
          f"{bats['individual_local_identifier'].nunique()} individuals")

    results = []

    for _, row in tqdm(bats.iterrows(), total=len(bats), desc="Extracting env data"):
        ts      = row["timestamp"]
        lat     = row["lat"]
        lon     = row["lon"]
        bearing = row.get("bearing", None)   # movement bearing if computed in Script 2
        if pd.isna(bearing):
            bearing = None

        date_str = ts.strftime("%Y-%m-%d")

        # Fetch 3-day hourly window (prev / date / next)
        df_h = fetch_hourly(lat, lon, date_str)
        time.sleep(BATCH_DELAY)

        env_row = {
            "individual_local_identifier": row["individual_local_identifier"],
            "timestamp": ts,
        }

        if df_h is None:
            results.append(env_row)
            continue

        # ── Window 1: Dusk ± 1 hour (hours 18–22 UTC as broad proxy for dusk in
        #    southern Europe/N Africa during Aug–Nov; adjust if needed)
        #    A proper solar calculation is preferable — see note below. ──────────
        dusk_hours = [18, 19, 20, 21, 22]
        same_day   = df_h["time"].dt.date == ts.date()
        dusk_mask  = same_day & df_h["time"].dt.hour.isin(dusk_hours)

        env_row.update(summarise_window(df_h, dusk_mask, bearing, prefix="dusk"))

        # ── Window 2: Nocturnal flight window (21:00–05:00 UTC next day) ─────────
        ts_date = ts.replace(hour=0, minute=0, second=0, microsecond=0)
        night_start = ts_date.replace(hour=21)
        night_end   = ts_date.replace(hour=5) + timedelta(days=1)
        night_mask  = (df_h["time"] >= night_start) & (df_h["time"] <= night_end)

        env_row.update(summarise_window(df_h, night_mask, bearing, prefix="night"))

        # ── Scalar: 3-hour dusk pressure tendency ───────────────────────────────
        env_row["dusk_pressure_tendency_3h"] = pressure_tendency(
            df_h, target_hours=dusk_hours, delta_h=3
        )

        results.append(env_row)

    env_df = pd.DataFrame(results)

    # Merge back onto bats_daily
    out = bats.merge(
        env_df,
        on=["individual_local_identifier", "timestamp"],
        how="left"
    )

    out.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved: {OUTPUT_CSV}  ({len(out)} rows, {len(out.columns)} columns)")
    print(f"Cache: {len(list(CACHE_DIR.glob('*.json')))} files in {CACHE_DIR}")


if __name__ == "__main__":
    main()
