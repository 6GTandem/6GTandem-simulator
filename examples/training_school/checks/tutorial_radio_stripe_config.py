# %% [markdown]
# # Tutorial: Generating a RadioStripe Environment Configuration
#
# **Goal of this script:**
# Show a newcomer how to *programmatically* create the environment configuration
# that describes a RadioStripe deployment: room dimensions, stripe positions,
# and Radio Unit (RU) coordinates.
#
# You will:
#   1. Define your scenario as a few plain Python variables (no YAML editing).
#   2. Call a helper function that builds the full configuration dictionary.
#   3. Inspect the result — stripe count, RU coordinates, physical-layer params.
#   4. Save the generated config to a YAML file so other scripts can load it.
#   5. Reload and verify the YAML round-trip.
#
# > **How to run:**
# > Open in VS Code and click *Run Cell* (▷) above each `# %%` marker, **or**
# > run the whole file from a terminal:
# >   python tutorial_radio_stripe_config.py
#
# No prior knowledge of the codebase is required beyond basic Python.

# %%
# ===========================================================================
# STEP 0 — Boilerplate: resolve the repository root so imports work
# ===========================================================================
# The simulator code lives in several sub-packages (wireless_channel/,
# sub_THz_stripe/, …).  We need the *repository root* on sys.path so that
# Python can find them regardless of where you launched this script from.

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml  # PyYAML — used to dump/load YAML config files

# Walk up from this file's location until we find the "wireless_channel" folder,
# which is a reliable marker for the repository root.
root = Path(__file__).resolve().parent
while not (root / "wireless_channel").exists() and root != root.parent:
    root = root.parent

if not (root / "wireless_channel").exists():
    raise RuntimeError(
        "Could not locate the repository root.  "
        "Make sure you are running this script from inside the 6GTandem-simulator repo."
    )

os.chdir(root)
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

# Avoid local module shadowing (examples/training_school/utils.py vs repo utils.py)
# because plotter imports repo-level utils with `from utils import ...`.
script_dir = str(Path(__file__).resolve().parent)
sys.path = [p for p in sys.path if p not in ("", script_dir)]
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

print(f"Repository root: {root}")

# Our training-school utilities (build_radio_stripe_config lives here).
from examples.training_school.utils import build_radio_stripe_config, generate_random_ue_positions
from examples.training_school.utils import select_active_radio_unit
from plotter import plotter
from sub_THz_stripe.amplifier.amplifier import Amplifier
from sub_THz_stripe.central_unit.central_unit import CentralUnit
from sub_THz_stripe.combiner.combiner import Combiner
from sub_THz_stripe.coupler.coupler import Coupler
from sub_THz_stripe.phase_shifter.phase_shifter import PhaseShifter
from sub_THz_stripe.radio_unit.radio_unit import RadioUnit
from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.splitter.splitter import Splitter
from wireless_channel.subTHz_channel import build_channel
from wireless_channel.waveforms import Waveform

# %%
# ===========================================================================
# STEP 1 — Define your scenario
# ===========================================================================
# These are the ONLY values you need to change to describe a new deployment.
# Everything else (stripe positions, RU coordinates, …) is derived from them.
# ---------------------------------------------------------------------------

# --- Room geometry ----------------------------------------------------------
ROOM_X = 10.0   # metres — room width  (X-axis, cross-stripe direction)
ROOM_Y = 20.0   # metres — room length (Y-axis, along-stripe direction)
ROOM_Z = 3.5    # metres — ceiling height (Z-axis)
#                           RUs are mounted at the ceiling: z = ROOM_Z
#                           UEs are assumed to be at floor level: z ≈ 1.5 m

# --- Stripe topology --------------------------------------------------------
N_STRIPES = 3           # number of RadioStripes running across the room
N_RUS_PER_STRIPE = 5   # number of Radio Units per stripe
RU_SPACING_M = 3     # distance between adjacent RUs along a stripe (Y direction)
STRIPE_SPACING_M = 2.0 # distance between adjacent stripes (X direction)

# With the defaults above the layout looks like this (top view):
#
#   Y →
#   ┌──────────────────────────────────────────────────┐
#   │                      (ROOM_X=10 m)                │
#   │  CU─RU─RU─RU─RU─RU   stripe 0   (x=2 m)        │
#   │  CU─RU─RU─RU─RU─RU   stripe 1   (x=3 m)        │
#   │  CU─RU─RU─RU─RU─RU   stripe 2   (x=4 m)        │
#   │                                                   │
#   │          ★  ← UE_X=5, UE_Y=10                   │
#   └──────────────────────────────────────────────────┘
#       ROOM_Y = 20 m

# --- User Equipment (UE) location ------------------------------------------
# Generate random UE positions inside the room.
N_RANDOM_USERS = 20
UE_Z = 1.5
UE_WALL_MARGIN_M = 0.1
UE_RANDOM_SEED = 42

# --- Plotting ---------------------------------------------------------------
plot_room = True

# --- Component-config coupling ----------------------------------------------
# If True, the tutorial component config will be updated so the fiber segment
# length between adjacent RUs matches RU_SPACING_M.
UPDATE_COMPONENT_FIBER_LENGTH = True

# %%
# ===========================================================================
# STEP 2 — Generate the configuration dictionary
# ===========================================================================
# build_radio_stripe_config() does all the geometry arithmetic:
#   • places stripes at x = 2, 2+STRIPE_SPACING_M, ...
#   • places RUs  at y = 1, 1+RU_SPACING_M, ...
#   • mounts everything at the ceiling (z = ROOM_Z)
#   • puts the Central Unit (fibre head-end) against the left wall at x = 1 m
#
# The function returns a plain Python dict — no files are touched yet.

env_config = build_radio_stripe_config(
    room_x=ROOM_X,
    room_y=ROOM_Y,
    room_z=ROOM_Z,
    n_stripes=N_STRIPES,
    n_rus_per_stripe=N_RUS_PER_STRIPE,
    ru_spacing_m=RU_SPACING_M,
    stripe_spacing_m=STRIPE_SPACING_M,
)

print("\n" + "=" * 60)
print("Generated environment configuration")
print("=" * 60)

# %%
# ===========================================================================
# STEP 2b — Add the UE to the configuration
# ===========================================================================
# build_radio_stripe_config() returns an empty ue_positions list on purpose —
# it only handles the infrastructure side (stripes, RUs, room).
# Here we generate N random users within room bounds.
env_config["ue_positions"] = generate_random_ue_positions(
    n_users=N_RANDOM_USERS,
    room_x=ROOM_X,
    room_y=ROOM_Y,
    z_height=UE_Z,
    wall_margin=UE_WALL_MARGIN_M,
    seed=UE_RANDOM_SEED,
)

print(
    f"Generated {len(env_config['ue_positions'])} random UE positions "
    f"(z={UE_Z:.2f} m, wall_margin={UE_WALL_MARGIN_M:.2f} m, seed={UE_RANDOM_SEED})"
)

# Optional top-view room plot (same helper used in other tutorials).
if plot_room:
    plotter.plot_room(env_config)
    plt.show()

# %%
# ===========================================================================
# STEP 3 — Inspect the configuration
# ===========================================================================
# Let's print a human-readable summary so you can verify the geometry without
# having to dig into the raw dict.

# --- Room summary -----------------------------------------------------------
room = env_config["room"]
print(f"\nRoom dimensions:  {room['x']} m (X) × {room['y']} m (Y) × {room['z']} m (Z)")

# --- Stripe / RU summary ----------------------------------------------------
sc = env_config["stripe_config"]
print(f"\nStripe topology:")
print(f"  Number of stripes       : {sc['N_stripes']}")
print(f"  RUs per stripe          : {sc['N_RUs']}")
print(f"  Spacing between stripes : {sc['space_between_stripes']} m  (X direction)")
print(f"  Spacing between RUs     : {sc['space_between_RUs']} m  (Y direction)")
print(f"  Stripe orientation      : {sc['stripe_direction']}-axis")
print(f"  First RU position       : {sc['stripe_start_pos']}")
print(f"  Last  RU position       : {sc['stripe_end_pos']}")

# --- Per-stripe RU coordinates ----------------------------------------------
print(f"\nRadio Unit coordinates (all at ceiling z={ROOM_Z} m):")

radio_stripes = env_config["radio_stripes"]
for s_idx, stripe in enumerate(radio_stripes):
    # The first entry of each stripe is the Central Unit; the rest are RUs.
    cu_entry = stripe[0]
    cu_xyz = cu_entry["central_unit"]

    ru_entries = stripe[1:]   # everything after the Central Unit
    ru_positions = [list(e["radio_unit"].values()) for e in ru_entries]

    # Summarise as x fixed, y ranging from first to last
    x_coord = ru_positions[0][0]
    y_coords = [pos[1] for pos in ru_positions]
    print(
        f"  Stripe {s_idx}:  CU at ({cu_xyz['x']:.1f}, {cu_xyz['y']:.1f}, {cu_xyz['z']:.1f})  |"
        f"  {len(ru_entries)} RUs at x={x_coord:.1f},  y = {y_coords[0]:.1f} … {y_coords[-1]:.1f},  z={ROOM_Z}"
    )

# --- Physical-layer summary -------------------------------------------------
st = env_config["sub_thz"]
print(f"\nSub-THz physical layer:")
print(f"  Carrier frequency : {st['fc']/1e9:.2f} GHz")
print(f"  Bandwidth         : {st['bw']/1e9:.1f} GHz")
print(f"  Subcarriers       : {st['num_subcarriers']}")
print(f"  Doppler           : {st['doppler']}")

# --- Antenna summary --------------------------------------------------------
ant = env_config["antenna"]
print(f"\nAntenna:")
print(f"  Elements per RU : {ant['N_antennas']}")
print(f"  Pattern model   : {ant['pattern']}")
print(f"  Polarisation    : {ant['polarization']}")

# --- UE positions -----------------------------------------------------------
ues = env_config["ue_positions"]
print(f"\nUE positions: {len(ues)} defined")
for i, ue in enumerate(ues):
    print(f"  UE {i}: x={ue['x']:.2f} m,  y={ue['y']:.2f} m,  z={ue['z']:.2f} m")

# %%
# ===========================================================================
# STEP 4 — Save the configuration to a YAML file
# ===========================================================================
# Saving to YAML makes the config reusable by any script that loads env configs
# the normal way (yaml.safe_load).

# Output path: configs/tutorial_env_empty.yaml (in the same folder as this script)
configs_dir = Path(__file__).parent / "configs"
configs_dir.mkdir(exist_ok=True)   # create the folder if it doesn't already exist
output_path = configs_dir / "tutorial_env_empty.yaml"

with open(output_path, "w") as f:
    yaml.dump(
        env_config,
        f,
        default_flow_style=False,   # human-readable block style
        sort_keys=False,            # keep original insertion order
        allow_unicode=True,
    )

print(f"\nConfiguration saved to: {output_path.resolve()}")

# Optionally keep component config aligned with the chosen RU spacing.
if UPDATE_COMPONENT_FIBER_LENGTH:
    component_config_path = configs_dir / "tutorial_component_config.yaml"
    with open(component_config_path, "r", encoding="utf8") as f:
        component_config = yaml.safe_load(f)

    component_config.setdefault("fiber", {})
    component_config["fiber"]["length"] = float(RU_SPACING_M)

    with open(component_config_path, "w", encoding="utf8") as f:
        yaml.dump(
            component_config,
            f,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
        )

    print(
        "Updated component config: "
        f"fiber.length={component_config['fiber']['length']} m at {component_config_path.resolve()}"
    )

# %%
# ===========================================================================
# STEP 5 — Reload and verify the YAML round-trip
# ===========================================================================
# This confirms that the saved YAML is valid and can be loaded back without
# any data loss.  Any script that needs the env config would start from here.

with open(output_path) as f:
    reloaded = yaml.safe_load(f)

# Quick sanity checks
assert reloaded["stripe_config"]["N_stripes"] == N_STRIPES, "Stripe count mismatch after reload!"
assert reloaded["stripe_config"]["N_RUs"] == N_RUS_PER_STRIPE, "RU count mismatch after reload!"
assert reloaded["room"]["x"] == ROOM_X
assert reloaded["room"]["y"] == ROOM_Y
assert reloaded["room"]["z"] == ROOM_Z
assert len(reloaded["radio_stripes"]) == N_STRIPES
assert all(
    len(stripe) == N_RUS_PER_STRIPE + 1   # +1 for the Central Unit
    for stripe in reloaded["radio_stripes"]
)

print("\nReload verification: PASSED")
print(f"  Stripes in file : {reloaded['stripe_config']['N_stripes']}")
print(f"  RUs per stripe  : {reloaded['stripe_config']['N_RUs']}")
print(f"  Entries/stripe  : {len(reloaded['radio_stripes'][0])}  (1 Central Unit + {N_RUS_PER_STRIPE} Radio Units)")

# %%
# ===========================================================================
# STEP 6 — Generate the signal (bits -> QAM -> OFDM)
# ===========================================================================
# This section mirrors the "generate the signal" part of basic_ul_tutorial.ipynb.
# We keep everything local to this script so newcomers can run one file end-to-end.

waveform_config_path = configs_dir / "tutorial_waveform_config.yaml"
with open(waveform_config_path, "r", encoding="utf8") as f:
    waveform_config = yaml.safe_load(f)

freq_band_config = env_config["sub_thz"]
wf = Waveform.from_config(waveform_config, freq_band_config)

# Generate random bit stream and map it to QAM symbols.
bits = wf.generate_bits()
qam = wf.qam_modulate()

# OFDM modulation to time-domain signal.
ofdm_time = wf.ofdm_modulate()  # shape: n_ofdm_symbols x (fft_size + cp_length)
x_combined_freq = wf.ofdm_time_to_freq(ofdm_time)
xsubc = wf.extract_subcarriers(x_combined_freq)

print("\nSignal generation summary")
print(f"  Bit stream length          : {bits.size}")
print(f"  Number of QAM symbols      : {qam.size}")
print(f"  OFDM time-domain shape     : {ofdm_time.shape}")
print(f"  Subcarrier grid shape      : {xsubc.shape}")

# Plot 1: generated bits (first 100 bits)
num_bits_to_plot = 100
plt.figure(figsize=(12, 3))
plt.step(np.arange(num_bits_to_plot), bits[:num_bits_to_plot], where="mid")
plt.ylim(-0.1, 1.1)
plt.xlim(0, num_bits_to_plot - 1)
plt.xlabel("Bit index")
plt.ylabel("Bit value")
plt.title("Generated bit sequence")
plt.grid(True, alpha=0.3)
plt.show()

# Plot 2: QAM constellation
plt.figure(figsize=(6, 6))
plt.scatter(qam.real, qam.imag, s=10, alpha=0.7)
plt.axhline(0, color="gray", linewidth=0.8)
plt.axvline(0, color="gray", linewidth=0.8)
plt.xlabel("In-phase (I)")
plt.ylabel("Quadrature (Q)")
plt.title("QAM symbols constellation")
plt.grid(True, alpha=0.3)
plt.axis("equal")
plt.show()

# Plot 3: OFDM time-domain waveform for first symbol
num_samples_to_plot = 600
start_idx = 0
end_idx = min(start_idx + num_samples_to_plot, ofdm_time.shape[-1])

ofdm_symbol = ofdm_time[0]  # first OFDM symbol
n = np.arange(start_idx, end_idx)

plt.figure(figsize=(12, 4))
plt.plot(n, ofdm_symbol[start_idx:end_idx].real, label="Real", linewidth=1.2)
plt.plot(n, ofdm_symbol[start_idx:end_idx].imag, label="Imag", linewidth=1.2, alpha=0.8)
plt.xlabel("Sample index")
plt.ylabel("Amplitude")
plt.title(f"OFDM time-domain signal (samples {start_idx} to {end_idx-1})")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# Plot 4: PSD of the OFDM waveform
ofdm_signal = ofdm_time.reshape(-1)  # flatten (n_symbols, n_samples) -> 1D
wf.plot_psd(ofdm_signal)
plt.show()

# %%
# ===========================================================================
# STEP 7 — Generate user and send uplink signal over LOS
# ===========================================================================
# This section mirrors the "generate user and send uplink signal" part of
# basic_ul_tutorial.ipynb using the topology/config generated above.

# Select which generated UE to use for this uplink example.
UE_INDEX_FOR_UPLINK = 5

# Beam and receive-selection knobs for this tutorial section.
UE_BEAM = 0
RU_BEAM = 0
USE_AUTO_ACTIVE_RU = True
ACTIVE_RU_SELECTION_MODE = "distance"  # default: closest RU to selected UE
MANUAL_ACTIVE_STRIPE_IDX = 0
MANUAL_ACTIVE_RU_IDX = 0

component_config_path = configs_dir / "tutorial_component_config.yaml"
with open(component_config_path, "r", encoding="utf8") as f:
    component_config = yaml.safe_load(f)

# 1) Build the UE transmit chain (same hardware model as a central unit).
ue = CentralUnit()
ofdm_time_after_ue = ue.run(ofdm_time)

amp = Amplifier(bw=env_config["sub_thz"]["bw"])
coup = Coupler(wf=wf)
split = Splitter(env_config["antenna"]["N_antennas"])
comb = Combiner()
ps = PhaseShifter(
    num_shifters=env_config["antenna"]["N_antennas"],
    resolution=int(component_config["phase_shifter"]["resolution"]),
)
ue_ru = RadioUnit(
    x=0,
    y=0,
    z=0,
    boost_amp=amp,
    antenna_amp=amp,
    coup_in=coup,
    coup_out=coup,
    splitter=split,
    combiner=comb,
    pshift=ps,
)

# 2) Build all stripes from the generated environment + component config.
stripes: list[RadioStripe] = []
for stripe_cfg in env_config["radio_stripes"]:
    stripes.append(
        RadioStripe.from_config_locations(
            stripe_cfg,
            component_config,
            env_config["antenna"]["N_antennas"],
            wf,
        )
    )

# 3) Generate UE transmit signal toward the stripe.
ue_shift = ue_ru.phase_shifter.get_phases(UE_BEAM)
iq_data_tx, _ = ue_ru.transmit(ofdm_time_after_ue, ue_shift)

ue_positions = env_config["ue_positions"]
if UE_INDEX_FOR_UPLINK < 0 or UE_INDEX_FOR_UPLINK >= len(ue_positions):
    raise ValueError(
        f"UE_INDEX_FOR_UPLINK={UE_INDEX_FOR_UPLINK} out of range [0, {len(ue_positions)-1}]"
    )
ue_pos = ue_positions[UE_INDEX_FOR_UPLINK]
print(f"\nSelected UE for uplink: index={UE_INDEX_FOR_UPLINK}, position={ue_pos}")

if USE_AUTO_ACTIVE_RU:
    ACTIVE_STRIPE_IDX, ACTIVE_RU_IDX, active_ru_pos, active_ru_distance = select_active_radio_unit(
        ue_position=ue_pos,
        radio_stripes=env_config["radio_stripes"],
        mode=ACTIVE_RU_SELECTION_MODE,
    )
    print(
        "Auto-selected active RU: "
        f"mode={ACTIVE_RU_SELECTION_MODE}, stripe={ACTIVE_STRIPE_IDX}, ru={ACTIVE_RU_IDX}, "
        f"distance={active_ru_distance:.3f} m, ru_pos={active_ru_pos}"
    )
else:
    ACTIVE_STRIPE_IDX = MANUAL_ACTIVE_STRIPE_IDX
    ACTIVE_RU_IDX = MANUAL_ACTIVE_RU_IDX
    active_ru_pos = env_config["radio_stripes"][ACTIVE_STRIPE_IDX][ACTIVE_RU_IDX + 1]["radio_unit"]
    print(f"Manual active RU: stripe={ACTIVE_STRIPE_IDX}, ru={ACTIVE_RU_IDX}")

# Plot room again with active selection highlighted.
if plot_room:
    plotter.plot_room_with_active_selection(
        env_config,
        active_ue_pos=ue_pos,
        active_ru_pos=active_ru_pos,
    )
    plt.show()

# 4) Build a LOS wireless channel and transmit to one selected RU.
channel = build_channel(
    channel_model="los",
    ue_coordinates=ue_pos,
    sim_env="office_space_perpendicular",
    component_config=component_config,
    stripe_positions=env_config["radio_stripes"],
    waveform=wf,
    Nr_ue_antennas=env_config["antenna"]["N_antennas"],
    Nr_ru_antennas=env_config["antenna"]["N_antennas"],
    debug=False,
    los_normalize_gain=False,
)

iq_data_rx = channel.transmit_ul_id(iq_data_tx, ACTIVE_STRIPE_IDX, ACTIVE_RU_IDX, wf)

# 5) Receive through the selected stripe RU and process back to subcarriers.
stripe = stripes[ACTIVE_STRIPE_IDX]
stripe.active_unit = ACTIVE_RU_IDX
ru_shift = stripe.radio_units[0].phase_shifter.get_phases(RU_BEAM)
y, _ = stripe.receive(iq_data_rx, ru_shift)

y_combined_freq = wf.ofdm_time_to_freq(y)
rx_subc = wf.extract_subcarriers(y_combined_freq)

# Plot received pilots first to visualise the raw received quality.
wf.plot_constellation(rx_subc[0, wf.pilot_indices], title="Received pilots (before equalization)")

# LS channel estimation + one-tap equalization + constellation plot.
H_est = wf.channel_estimate_ls(rx_subc, interp_mode="phase")
eq_subc = wf.equalize_one_tap(rx_subc, H_est)
wf.plot_constellation(eq_subc, symbols_tx=xsubc, title="Equalized received symbols")

print(
    "Uplink demo complete: "
    f"stripe={ACTIVE_STRIPE_IDX}, ru={ACTIVE_RU_IDX}, ue_index={UE_INDEX_FOR_UPLINK}, "
    f"ru_selection_mode={ACTIVE_RU_SELECTION_MODE if USE_AUTO_ACTIVE_RU else 'manual'}"
)

# %%
# ===========================================================================
# What's next?
# ===========================================================================
# This script only created and inspected the deployment geometry.
# The next tutorial steps will:
#
#   Step 2: Add UE positions and visualise the layout with a top-view plot.
#   Step 3: Load the waveform and component configs alongside the env config.
#   Step 4: Build the channel model and compute path gains + delays.
#   Step 5: Run a full uplink simulation and evaluate throughput.
#
# See tutorial.py for a complete end-to-end example.

print("\nDone — see tutorial.py for the next steps.")
