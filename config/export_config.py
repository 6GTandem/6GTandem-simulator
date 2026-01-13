#!/usr/bin/env python3
"""
Generate a YAML config using default stripe parameters, room dimensions,
and UE positions from 'ue_locations_5681.nc'.

Usage:
  python generate_config.py [output.yaml]

If no argument is given, defaults to 'generated_config.yaml'.
"""

from __future__ import annotations
from pathlib import Path, PurePath

import argparse
import numpy as np
import xarray as xr
import yaml


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


def load_ue_positions(ds: xr.Dataset) -> list[dict[str, float]]:
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
    z_arr = np.asarray((np.zeros_like(x_arr) if z_da is None else z_da.values), dtype=float)

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
        {"x": float(x), "y": float(y), "z": float(z)} for x, y, z in zip(x_arr.ravel(), y_arr.ravel(), z_arr.ravel())
    ]


def build_radio_stripes(
    stripe_start_pos, stripe_end_pos, N_stripes, space_between_stripes, space_between_RUs, N_RUs, stripe_direction
) -> list[list[dict[str, dict[str, float]]]]:
    """Build radio_stripes with entries like {'radio_unit': {'x':..., 'y':..., 'z':...}}."""
    stripes: list[list[dict[str, dict[str, float]]]] = []

    # Compute the RU positions.
    x_pos = np.linspace(stripe_start_pos[0], stripe_end_pos[0], N_RUs).tolist()
    y_pos = np.linspace(stripe_start_pos[1], stripe_end_pos[1], N_RUs).tolist()
    z_pos = np.linspace(stripe_start_pos[2], stripe_end_pos[2], N_RUs).tolist()

    for stripe_idx in range(N_stripes):
        stripe_list: list[dict[str, dict[str, float]]] = []
        rux = x_pos[0]
        ruy = y_pos[0]

        if stripe_direction == "y":
            rux -= space_between_RUs
            ruy += stripe_idx * space_between_stripes
        elif stripe_direction == "x":
            ruy -= space_between_RUs
            rux += stripe_idx * space_between_stripes

        stripe_list.append(
            {
                "central_unit": {
                    "x": rux,
                    "y": ruy,
                    "z": z_pos[0],
                }
            }
        )
        for ru_idx in range(N_RUs):
            rux = x_pos[ru_idx]
            ruy = y_pos[ru_idx]

            if stripe_direction == "y":
                rux += stripe_idx * space_between_stripes
            elif stripe_direction == "x":
                ruy += stripe_idx * space_between_stripes

            stripe_list.append({"radio_unit": {"x": rux, "y": ruy, "z": z_pos[ru_idx]}})
        stripes.append(stripe_list)

    return stripes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a YAML config file for the simulation.")
    parser.add_argument(
        "simulation_environment", type=str, help="Name of the simulation environment present in the sionna datasets."
    )
    parser.add_argument("output_file", type=str, help="Output YAML file name.")
    args = parser.parse_args()

    simulation_environment = args.simulation_environment
    output_file = args.output_file

    # ---- Load parameters from office_rt_config.yml ----
    dataset_path = Path("wireless_channel/sionna_dataset")
    rt_config_path = PurePath(dataset_path, simulation_environment, "config.yaml")
    with open(rt_config_path, "r") as f:
        rt_params = yaml.safe_load(f)

    stripe_params = rt_params.get("stripe_config")
    N_RUs = stripe_params.get("N_RUs")
    N_stripes = stripe_params.get("N_stripes")
    space_between_RUs = stripe_params.get("space_between_RUs")
    stripe_direction = stripe_params.get("stripe_direction")
    space_between_stripes = stripe_params.get("space_between_stripes")
    stripe_start_pos = tuple(stripe_params.get("stripe_start_pos"))
    stripe_end_pos = tuple(stripe_params.get("stripe_end_pos"))

    # Make the room area the same as the area in which the UEs are plotted.
    ue_locations_config = rt_params.get("ue_locations_config")
    ue_area = ue_locations_config.get("ue_area")
    # If ue_area are defined in zones we have to extract it differently.
    if ue_area is None:
        zones = ue_locations_config.get("zones")
        x, y = 0, 0

        # Collect all the x and y values of the UE zones.
        uex = []
        uey = []
        uez = []
        for zone_name, zone in zones.items():
            xstart, xstop, ystart, ystop, zstart, zstop = zone["bounds"]

            uex.extend([xstart, xstop])
            uey.extend([ystart, ystop])
            uez.extend([zstart, zstop])
    else:
        xstart, xstop, ystart, ystop = ue_area
        uex = [xstart, xstop]
        uey = [ystart, ystop]
        uez = [ue_locations_config["z_height"]]
        
    # Collect all the bounds of the stripe.
    x_start, y_start, z_start = stripe_start_pos
    x_stop, y_stop, z_stop = stripe_end_pos

    # Look for the largest x, y and z value and use this for the room bounds.
    x = float(np.max([x_start, x_stop] + uex))
    y = float(np.max([y_start, y_stop] + uey))
    z = float(np.max([z_start, z_stop] + uez))

    ue_area = (x, y, z)

    # Load additional params for sub_thz and sub10GHz from YAML, following example.yml structure
    sub_thz_params = rt_params.get("subTHz_config")
    sub10GHz_params = rt_params.get("sub10GHz_config")
    antenna_params = rt_params.get("antenna_config")

    nc_path = PurePath(dataset_path, simulation_environment, "ue_locations/ue_locations.nc")
    ds = xr.open_dataset(nc_path)
    ue_positions = load_ue_positions(ds)
    radio_stripes = build_radio_stripes(
        stripe_start_pos, stripe_end_pos, N_stripes, space_between_stripes, space_between_RUs, N_RUs, stripe_direction
    )

    config = {
        "stripe_config": stripe_params,
        "room": ue_area,
        "radio_stripes": radio_stripes,
        "ue_positions": ue_positions,
        "sub_thz": sub_thz_params,
        "sub10GHz": sub10GHz_params,
        "antenna": antenna_params,
        "central_unit_fiber_length": 1.0,
    }

    output_yaml = Path(output_file)
    output_yaml.parent.mkdir(parents=True, exist_ok=True)
    with open(output_yaml, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False, default_flow_style=False)

    print(f"Wrote YAML config to: {output_yaml.resolve()}")
