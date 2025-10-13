#!/usr/bin/env python3
"""
Generate a YAML config using default stripe parameters, room dimensions,
and UE positions from 'ue_locations_5681.nc'.

Usage:
  python generate_config.py [output.yaml]

If no argument is given, defaults to 'generated_config.yaml'.
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import xarray as xr
import yaml


# ---- Load parameters from office_rt_config.yml ----
rt_config_path = Path("configurations/office_rt_config.yml")
with open(rt_config_path, "r") as f:
    rt_params = yaml.safe_load(f)

N_RUs = rt_params.get("N_RUs", 42)
N_stripes = rt_params.get("N_stripes", 13)
space_between_RUs = rt_params.get("space_between_RUs", 0.5)
space_between_stripes = rt_params.get("space_between_stripes", 0.5)
stripe_start_pos = tuple(rt_params.get("stripe_start_pos", [2.0, 2.5, 3.5]))
ROOM = rt_params.get("room", {"x": 10, "y": 25, "z": 5.45})

# Load additional params for sub_thz and sub10GHz from YAML, following example.yml structure
sub_thz_params = rt_params.get("subTHz_config", {})
sub10GHz_params = rt_params.get("sub10GHz_config", {})
antenna_params = rt_params.get("antenna_config", {})

output_file = "configurations/office_config.yml"
UE_NC_PATH = Path("configurations/ue_locations_5681.nc")


def _select_1d_numeric_var(ds: xr.Dataset, names: list[str]) -> xr.DataArray | None:
    """Find a numeric 1D variable by candidate names; reduce by selecting first index on extra dims."""
    for name in names:
        if name in ds.variables:
            da = ds[name]
            if da.ndim > 1:
                sel = {d: 0 for d in da.dims if da.sizes[d] > 1}
                da = da.isel(**sel).squeeze(drop=True)
            if da.ndim == 1 and np.issubdtype(da.dtype, np.number):
                return da
    return None


def load_ue_positions(nc_path: Path) -> list[dict[str, float]]:
    if not nc_path.exists():
        raise FileNotFoundError(f"UE NetCDF file not found: {nc_path}")

    ds = xr.open_dataset(nc_path)

    # Common candidate names; extend if your file differs
    x_candidates = ["ue_x", "x_ue", "x", "X", "pos_x", "ueX"]
    y_candidates = ["ue_y", "y_ue", "y", "Y", "pos_y", "ueY"]
    z_candidates = ["ue_z", "z_ue", "z", "Z", "pos_z", "ueZ", "height"]

    x_da = _select_1d_numeric_var(ds, x_candidates)
    y_da = _select_1d_numeric_var(ds, y_candidates)
    z_da = _select_1d_numeric_var(ds, z_candidates)

    if x_da is None or y_da is None:
        raise ValueError("Could not find UE x/y arrays in NetCDF.")

    x_arr = np.asarray(x_da.values, dtype=float)
    y_arr = np.asarray(y_da.values, dtype=float)
    z_arr = np.asarray(
        (np.zeros_like(x_arr) if z_da is None else z_da.values), dtype=float
    )

    # Load invalid_point if it exists
    mask = None
    if "invalid_point" in ds.variables:
        invalid_da = ds["invalid_point"]
        if invalid_da.ndim > 1:
            sel = {d: 0 for d in invalid_da.dims if invalid_da.sizes[d] > 1}
            invalid_da = invalid_da.isel(**sel).squeeze(drop=True)
        mask = ~np.asarray(invalid_da.values, dtype=bool)
        # Only keep where invalid_point == True
        if mask.shape != x_arr.shape:
            mask = np.broadcast_to(mask, x_arr.shape)
    else:
        # Keep all if not present
        mask = np.ones_like(x_arr, dtype=bool)

    # Apply filter
    x_arr = x_arr[mask]
    y_arr = y_arr[mask]
    z_arr = z_arr[mask]

    return [
        {"x": float(x), "y": float(y), "z": float(z)}
        for x, y, z in zip(x_arr.ravel(), y_arr.ravel(), z_arr.ravel())
    ]


def build_radio_stripes() -> list[list[dict[str, dict[str, float]]]]:
    """Build radio_stripes with entries like {'radio_unit': {'x':..., 'y':..., 'z':...}}."""
    stripes: list[list[dict[str, dict[str, float]]]] = []
    x0, y0, z0 = stripe_start_pos
    for stripe_idx in range(N_stripes):
        stripe_list: list[dict[str, dict[str, float]]] = []
        stripe_list.append(
            {
                "central_unit": {
                    "x": x0 + stripe_idx * space_between_stripes,
                    "y": y0 - space_between_RUs,
                    "z": z0,
                }
            }
        )
        for ru_idx in range(N_RUs):
            x = x0 + stripe_idx * space_between_stripes
            y = y0 + ru_idx * space_between_RUs
            z = z0
            stripe_list.append({"radio_unit": {"x": x, "y": y, "z": z}})
        stripes.append(stripe_list)
    return stripes


def main():
    # ---- Load parameters from office_co_config.yml ----
    co_config_path = Path("configurations/office_co_config.yml")
    with open(co_config_path, "r") as f:
        co_params = yaml.safe_load(f)
    wf_config_path = Path("configurations/office_wf_config.yml")
    with open(wf_config_path, "r") as f:
        wf_params = yaml.safe_load(f)

    output_yaml = Path(output_file)

    ue_positions = load_ue_positions(UE_NC_PATH)
    radio_stripes = build_radio_stripes()

    config = {
        "stripe_config": {
            "N_RUs": N_RUs,
            "N_stripes": N_stripes,
            "space_between_RUs": space_between_RUs,
            "space_between_stripes": space_between_stripes,
            "stripe_start_pos": list(stripe_start_pos),
        },
        "component_config":co_params,
        "waveform_config": wf_params,
        "room": ROOM,
        "radio_stripes": radio_stripes,
        "ue_positions": ue_positions,
        "sub_thz": sub_thz_params,
        "sub10GHz": sub10GHz_params,
        "antenna": antenna_params,
        "central_unit_fiber_length": 0.5,
    }

    output_yaml.parent.mkdir(parents=True, exist_ok=True)
    with open(output_yaml, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False, default_flow_style=False)

    print(f"Wrote YAML config to: {output_yaml.resolve()}")


if __name__ == "__main__":
    main()
