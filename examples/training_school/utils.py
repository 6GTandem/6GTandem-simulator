"""
Utility helpers for the 6GTandem training-school scripts.

The main entry point for newcomers is :func:`build_radio_stripe_config`, which
generates a complete environment configuration dictionary from a small set of
high-level parameters (room size, number of stripes, number of RUs per stripe).

The returned dict matches the YAML structure that the simulator reads from the
``environments/`` folder, so you can either pass it directly in Python or dump
it to a YAML file and load it the normal way.
"""

from __future__ import annotations

import random

import numpy as np


def select_active_radio_unit(
    ue_position: dict,
    radio_stripes: list[list[dict]],
    mode: str = "distance",
) -> tuple[int, int, dict, float]:
    """Select the active RU for a UE.

    Parameters
    ----------
    ue_position : dict
        UE coordinates with keys ``x``, ``y``, ``z``.
    radio_stripes : list[list[dict]]
        Stripe configuration where each stripe contains one ``central_unit``
        entry followed by one or more ``radio_unit`` entries.
    mode : str
        Selection mode. Supported values:
        - ``"distance"``: choose the RU with minimum Euclidean distance
          to the UE (default).
        - ``"random"``: choose a random RU uniformly across all stripes.

    Returns
    -------
    tuple[int, int, dict, float]
        ``(stripe_idx, ru_idx, ru_position, distance_m)`` for the selected RU.
    """
    normalized_mode = str(mode).strip().lower()
    if normalized_mode not in ("distance", "random"):
        raise ValueError(f"Unsupported RU selection mode '{mode}'. Supported: distance, random")

    if not radio_stripes:
        raise ValueError("radio_stripes is empty.")

    # Collect all RU candidates across all stripes.
    all_rus: list[tuple[int, int, dict]] = []
    for stripe_idx, stripe in enumerate(radio_stripes):
        ru_idx = 0
        for entry in stripe:
            ru_pos = entry.get("radio_unit")
            if ru_pos is None:
                continue
            all_rus.append((stripe_idx, ru_idx, ru_pos))
            ru_idx += 1

    if not all_rus:
        raise ValueError("No radio_unit entries found in radio_stripes.")

    if normalized_mode == "random":
        stripe_idx, ru_idx, ru_position = random.choice(all_rus)
        ue_x = float(ue_position["x"])
        ue_y = float(ue_position["y"])
        ue_z = float(ue_position["z"])
        dx = float(ru_position["x"]) - ue_x
        dy = float(ru_position["y"]) - ue_y
        dz = float(ru_position["z"]) - ue_z
        distance_m = (dx * dx + dy * dy + dz * dz) ** 0.5
        return stripe_idx, ru_idx, ru_position, float(distance_m)

    # mode == "distance"
    ue_x = float(ue_position["x"])
    ue_y = float(ue_position["y"])
    ue_z = float(ue_position["z"])

    best = None
    for stripe_idx, ru_idx, ru_pos in all_rus:
        dx = float(ru_pos["x"]) - ue_x
        dy = float(ru_pos["y"]) - ue_y
        dz = float(ru_pos["z"]) - ue_z
        distance = (dx * dx + dy * dy + dz * dz) ** 0.5

        candidate = (distance, stripe_idx, ru_idx, ru_pos)
        if best is None or candidate[0] < best[0]:
            best = candidate

    distance_m, stripe_idx, ru_idx, ru_position = best
    return stripe_idx, ru_idx, ru_position, float(distance_m)


def build_radio_stripe_config(
    room_x: float,
    room_y: float,
    room_z: float,
    n_stripes: int,
    n_rus_per_stripe: int,
    *,
    ru_spacing_m: float = 1.0,
    stripe_spacing_m: float = 1.0,
    y_first_stripe: float = 1.0,
    x_cu: float = 1.0,
    x_first_stripe: float = 2.0,
    n_antennas: int = 4,
    fc: float = 157.75e9,
    bw: float = 3e9,
    num_subcarriers: int = 4096,
) -> dict:
    """
    Generate a full environment configuration dictionary for a RadioStripe
    deployment.

    Layout convention (same as the reference environments in environments/)
    -------------------------------------------------------------------
    - Stripes run along the *Y-axis* (lengthwise in the room).
        - Stripes are spaced ``stripe_spacing_m`` apart in *X*, starting from
            ``x_first_stripe``.
        - Radio Units along each stripe are spaced ``ru_spacing_m`` apart in *Y*,
            starting from ``y_first_stripe`` (distance from the front wall at y=0).
    - All RUs (and the Central Unit) are ceiling-mounted at z = ``room_z``.
    - The Central Unit for each stripe sits at x = ``x_cu`` (against the wall),
      at the same Y as the first RU of that stripe.
    - UE positions are *not* set inside this function — the caller adds them
      separately under the ``"ue_positions"`` key.

    Visual layout (top view, example: 3 stripes, 4 RUs each)::

        Y-axis →

        [CU]──[RU]──[RU]──[RU]──[RU]   stripe 0   (x=2)
        [CU]──[RU]──[RU]──[RU]──[RU]   stripe 1   (x=3)
        [CU]──[RU]──[RU]──[RU]──[RU]   stripe 2   (x=4)

        x=1  x=2  x=3  x=4  x=5  …
              ↑ stripe X positions

    Parameters
    ----------
    room_x, room_y, room_z : float
        Room dimensions in metres.  ``room_z`` is the ceiling height.
    n_stripes : int
        Number of RadioStripes.  Each stripe runs the full Y-length of the room
        (bounded by ``y_first_stripe`` at the front and the last RU at the back).
    n_rus_per_stripe : int
        Number of Radio Units placed along each stripe.
    ru_spacing_m : float
        Distance between adjacent Radio Units along each stripe in metres.
        Default is 1.0 m.
    stripe_spacing_m : float
        Distance between adjacent stripes in metres. Default is 1.0 m.
    y_first_stripe : float
        Gap (in metres) between the front wall (y=0) and the first RU.
        Default is 1.0 m.
    x_cu : float
        X position of the Central Unit for every stripe.  The CU represents
        the fibre-connected head node and is conventionally placed on the wall.
        Default is 1.0 m.
    x_first_stripe : float
        X position of the *first* RadioStripe.  Must satisfy
        ``x_first_stripe > x_cu`` so the CU is on the wall side.
        Default is 2.0 m.
    n_antennas : int
        Number of antenna elements per Radio Unit.  Default is 4.
    fc : float
        Carrier frequency in Hz.  Default is 157.75 GHz (sub-THz D-band).
    bw : float
        Signal bandwidth in Hz.  Default is 3 GHz.
    num_subcarriers : int
        Number of OFDM subcarriers.  Default is 4096.

    Returns
    -------
    dict
        Configuration dictionary with the following top-level keys:

        ``stripe_config``
            Metadata about the deployment geometry (used by the simulator
            to reconstruct stripe/RU indices).
        ``room``
            Room bounding-box dimensions.
        ``radio_stripes``
            List-of-lists.  ``radio_stripes[s]`` is the ordered list of
            units on stripe *s*: first entry is the ``central_unit``, then
            one ``radio_unit`` entry per RU.
        ``antenna``
            Antenna model and polarisation settings.
        ``sub_thz``
            Physical-layer parameters (carrier frequency, bandwidth, …).
        ``ue_positions``
            Empty list — add ``{"x": …, "y": …, "z": …}`` dicts here.
    """

    # ------------------------------------------------------------------ #
    # Validate inputs early so errors are obvious.                        #
    # ------------------------------------------------------------------ #
    if n_stripes < 1:
        raise ValueError(f"n_stripes must be >= 1, got {n_stripes}")
    if n_rus_per_stripe < 1:
        raise ValueError(f"n_rus_per_stripe must be >= 1, got {n_rus_per_stripe}")
    # if x_first_stripe <= x_cu:
    #     raise ValueError(
    #         f"x_first_stripe ({x_first_stripe}) must be greater than x_cu ({x_cu}) "
    #         "so that the Central Unit sits on the wall side of each stripe."
    #     )
    if ru_spacing_m <= 0:
        raise ValueError(f"ru_spacing_m must be > 0, got {ru_spacing_m}")
    if stripe_spacing_m <= 0:
        raise ValueError(f"stripe_spacing_m must be > 0, got {stripe_spacing_m}")

    # ------------------------------------------------------------------ #
    # 1.  Compute RU and stripe positions.                                #
    # ------------------------------------------------------------------ #
    # The CU sits at y_first_stripe; RUs start one ru_spacing_m further along Y.
    cu_y = y_first_stripe
    ru_y_positions = [y_first_stripe + (i + 1) * ru_spacing_m for i in range(n_rus_per_stripe)]

    # X coordinates of the n_stripes stripes.
    stripe_x_positions = [x_first_stripe + s * stripe_spacing_m for s in range(n_stripes)]

    # ------------------------------------------------------------------ #
    # 2.  Build stripe_config metadata block.                             #
    # ------------------------------------------------------------------ #
    # The simulator uses this block to quickly look up topology parameters
    # without having to parse the full radio_stripes list.
    stripe_cfg = {
        "N_RUs": n_rus_per_stripe,
        "N_stripes": n_stripes,
        "space_between_RUs": float(ru_spacing_m),
        "space_between_stripes": float(stripe_spacing_m),
        "stripe_direction": "y",          # stripes are oriented along the Y axis
        # stripe_start_pos / stripe_end_pos describe the first stripe's extent
        "stripe_start_pos": [stripe_x_positions[0], ru_y_positions[0], room_z],
        "stripe_end_pos": [stripe_x_positions[0], ru_y_positions[-1], room_z],
        "array_direction": "x",           # antenna array elements are spaced along X
    }

    # ------------------------------------------------------------------ #
    # 3.  Build radio_stripes list.                                       #
    # ------------------------------------------------------------------ #
    # radio_stripes is a list-of-lists (one inner list per stripe).
    # Each inner list starts with the Central Unit followed by the Radio Units.
    #
    # The Central Unit is placed one fiber length before the first RU
    # along the stripe direction (Y), connected via fiber.
    # All units are at the ceiling height z=room_z.
    radio_stripes = []
    for s_idx, stripe_x in enumerate(stripe_x_positions):
        stripe_entries = []

        # Central Unit — one per stripe, at the start of the stripe.
        stripe_entries.append({
            "central_unit": {
                "x": float(stripe_x),
                "y": float(cu_y),
                "z": float(room_z),
            }
        })

        # Radio Units — ceiling-mounted, spaced 1 m apart along Y.
        for ru_y in ru_y_positions:
            stripe_entries.append({
                "radio_unit": {
                    "x": float(stripe_x),
                    "y": float(ru_y),
                    "z": float(room_z),
                }
            })

        radio_stripes.append(stripe_entries)

    # ------------------------------------------------------------------ #
    # 4.  Antenna configuration.                                          #
    # ------------------------------------------------------------------ #
    # "tr38901" → directional pattern from 3GPP TR 38.901:
    #   - Main-lobe gain: 20 dBi, HPBW ≈ 65°.
    #   - Default boresight: RU points downward (-Z), UE points upward (+Z).
    # "isotropic" → omnidirectional (useful for debugging path-loss only).
    antenna = {
        "N_antennas": n_antennas,
        "pattern": "tr38901",
        "polarization": "V",
    }

    # ------------------------------------------------------------------ #
    # 5.  Sub-THz physical-layer parameters.                              #
    # ------------------------------------------------------------------ #
    sub_thz = {
        "fc": float(fc),
        "bw": float(bw),
        "num_subcarriers": num_subcarriers,
        "doppler": False,
    }

    # ------------------------------------------------------------------ #
    # 6.  Assemble the complete config dict.                              #
    # ------------------------------------------------------------------ #
    config = {
        "stripe_config": stripe_cfg,
        "room": {"x": float(room_x), "y": float(room_y), "z": float(room_z)},
        "radio_stripes": radio_stripes,
        "antenna": antenna,
        "sub_thz": sub_thz,
        # No UEs are defined by this function.  Add them as:
        #   config["ue_positions"].append({"x": …, "y": …, "z": …})
        "ue_positions": [],
    }

    # ------------------------------------------------------------------ #
    # 7.  Assert all units are within room bounds.                       #
    # ------------------------------------------------------------------ #
    for stripe in radio_stripes:
        for entry in stripe:
            unit = entry.get("central_unit") or entry.get("radio_unit")
            if unit is not None:
                x, y, z = unit["x"], unit["y"], unit["z"]
                assert 0 <= x <= room_x, f"Unit x={x} out of bounds (0, {room_x})"
                assert 0 <= y <= room_y, f"Unit y={y} out of bounds (0, {room_y})"
                assert 0 <= z <= room_z, f"Unit z={z} out of bounds (0, {room_z})"

    return config


def generate_random_ue_positions(
    n_users: int,
    room_x: float,
    room_y: float,
    *,
    z_height: float = 1.5,
    wall_margin: float = 0.1,
    seed: int | None = None,
) -> list[dict]:
    """Generate random UE positions inside the room bounds.

    Parameters
    ----------
    n_users : int
        Number of UE positions to generate.
    room_x, room_y : float
        Room dimensions in metres.
    z_height : float
        Fixed UE height in metres.
    wall_margin : float
        Keep UEs at least this many metres away from each wall.
    seed : int | None
        Optional RNG seed for reproducibility.
    """
    if n_users < 0:
        raise ValueError(f"n_users must be >= 0, got {n_users}")
    if room_x <= 0 or room_y <= 0:
        raise ValueError("room_x and room_y must be > 0")
    if wall_margin < 0:
        raise ValueError(f"wall_margin must be >= 0, got {wall_margin}")
    if wall_margin * 2 >= room_x or wall_margin * 2 >= room_y:
        raise ValueError(
            "wall_margin is too large for the room dimensions; need "
            "2*wall_margin < room_x and room_y"
        )

    rng = random.Random(seed)
    x_min, x_max = wall_margin, room_x - wall_margin
    y_min, y_max = wall_margin, room_y - wall_margin

    positions = []
    for _ in range(n_users):
        positions.append(
            {
                "x": rng.uniform(x_min, x_max),
                "y": rng.uniform(y_min, y_max),
                "z": float(z_height),
            }
        )
    return positions


def generate_grid_ue_positions(
    n_x: int,
    n_y: int,
    room_x: float,
    room_y: float,
    *,
    z_height: float = 1.5,
    wall_margin: float = 0.1,
) -> list[dict]:
    """Generate a uniform grid of UE positions inside the room bounds.

    Parameters
    ----------
    n_x, n_y : int
        Number of grid points along x and y axes.
    room_x, room_y : float
        Room dimensions in metres.
    z_height : float
        Fixed UE height in metres.
    wall_margin : float
        Keep UEs at least this many metres away from each wall.
    """
    if n_x < 1 or n_y < 1:
        raise ValueError(f"n_x and n_y must be >= 1, got n_x={n_x}, n_y={n_y}")
    if room_x <= 0 or room_y <= 0:
        raise ValueError("room_x and room_y must be > 0")
    if wall_margin < 0:
        raise ValueError(f"wall_margin must be >= 0, got {wall_margin}")
    if wall_margin * 2 >= room_x or wall_margin * 2 >= room_y:
        raise ValueError(
            "wall_margin is too large for the room dimensions; need "
            "2*wall_margin < room_x and room_y"
        )

    x_min, x_max = wall_margin, room_x - wall_margin
    y_min, y_max = wall_margin, room_y - wall_margin

    xs = np.linspace(x_min, x_max, n_x) if n_x > 1 else [0.5 * (x_min + x_max)]
    ys = np.linspace(y_min, y_max, n_y) if n_y > 1 else [0.5 * (y_min + y_max)]

    positions = []
    for y in ys:
        for x in xs:
            positions.append(
                {
                    "x": float(x),
                    "y": float(y),
                    "z": float(z_height),
                }
            )
    return positions
