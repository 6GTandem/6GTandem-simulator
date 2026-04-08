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
#   6. Generate an OFDM uplink signal and transmit it over a LOS channel.
#   7. Receive, equalise, and visualise the result.
#
# > **How to run:**
# > Open in VS Code and click *Run Cell* (▷) above each `# %%` marker, **or**
# > run the whole file from a terminal:
# >   python anotated_tutorial.py
#
# No prior knowledge of the codebase is required beyond basic Python.

# %% [markdown]
# ## Background: The 6GTandem Simulator
#
# This simulator models a **sub-THz RadioStripe uplink** communication system —
# the kind of indoor wireless network envisioned for 6G. It lets you:
#
# - **Design a deployment**: place RadioStripes on the ceiling of a room and
#   scatter User Equipment (UEs) on the floor.
# - **Generate and transmit signals**: create OFDM waveforms, push them through
#   realistic (or idealised) hardware models, and propagate them over a wireless
#   channel.
# - **Evaluate performance**: measure Bit Error Rate (BER), Error Vector
#   Magnitude (EVM), and inspect constellation diagrams.
#
# ### Repository layout — the packages you will use
#
# | Package / folder | What it contains |
# |---|---|
# | `sub_THz_stripe/` | Hardware component models: amplifiers, couplers, phase shifters, splitters, combiners, fiber links, and the `RadioStripe` / `RadioUnit` / `CentralUnit` classes that chain them together. |
# | `wireless_channel/` | Channel models (`Channel` class) and the `Waveform` class for OFDM signal generation, modulation, and demodulation. |
# | `plotter/` | Visualisation helpers — room layout plots, constellation diagrams, PSD plots. |
# | `environments/` | Pre-built YAML configuration files for various indoor deployments (office, arena, industry hall, …). |
# | `examples/training_school/` | **This tutorial** and its supporting utilities and config files. |
#
# ### The three configuration files
#
# The simulator is **config-driven**. Three YAML files describe a simulation:
#
# 1. **Environment config** — room dimensions, stripe/RU positions, UE locations,
#    carrier frequency, bandwidth, antenna model.
# 2. **Waveform config** — OFDM parameters: QAM order, number of symbols,
#    cyclic-prefix length, pilot spacing, transmit power.
# 3. **Component config** — hardware impairment models: amplifier gain/noise
#    figure, fiber length, coupler damping, phase-shifter resolution.
#
# This tutorial focuses on **creating the environment config from scratch**.
# The waveform and component configs are loaded from pre-existing files in
# `examples/training_school/configs/`.
#
# ### How the `# %%` cell markers work
#
# This `.py` file is a **VS Code Interactive Python** script. Each `# %%` line
# marks the start of a new cell that you can run independently (like a Jupyter
# notebook cell). Markdown cells start with `# %% [markdown]` and are rendered
# as formatted text in the interactive window. You can also run the entire file
# as a normal Python script from the terminal.

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

# %% [markdown]
# ## What is a RadioStripe?
#
# A **RadioStripe** is a distributed antenna system designed for indoor
# coverage at very high frequencies (sub-THz / D-band, around 150 GHz).
# The key idea:
#
# - A thin **stripe** (like an LED strip) is mounted on the ceiling.
# - Along the stripe sit multiple **Radio Units (RUs)** — small antenna
#   elements spaced a few metres apart.
# - All RUs on one stripe are **daisy-chained** by optical fibre and
#   couplers back to a single **Central Unit (CU)** that performs
#   baseband processing.
# - On the floor, **User Equipment (UEs)** — phones, laptops, sensors —
#   transmit uplink signals that are received by the nearest RU(s) and
#   forwarded along the stripe to the CU.
#
# At sub-THz frequencies the available bandwidth is enormous (several GHz),
# enabling multi-Gbps data rates. However, signals are highly directional and
# attenuate quickly with distance, so dense RU placement is essential. The
# RadioStripe architecture achieves this without running a separate cable to
# every antenna.
#
# ### Geometry conventions used in this simulator
#
# - **X-axis** — across the room (stripes are spaced in X).
# - **Y-axis** — along the room (RUs within a stripe are spaced in Y).
# - **Z-axis** — height.  RUs sit at ceiling height; UEs at ~1.5 m.
# - The **Central Unit** for each stripe is placed against one wall
#   (small X), and the RUs extend away from it along Y.
#
# The diagram below (from the scenario parameters) shows a top-down view of a
# typical deployment:
#
# ```
#   Y →
#   ┌──────────────────────────────────────────────────┐
#   │                      (ROOM_X)                     │
#   │  CU─RU─RU─RU─RU─RU   stripe 0                   │
#   │  CU─RU─RU─RU─RU─RU   stripe 1                   │
#   │  CU─RU─RU─RU─RU─RU   stripe 2                   │
#   │                                                   │
#   │          ★  UE                                   │
#   └──────────────────────────────────────────────────┘
#       ROOM_Y
# ```

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

# %% [markdown]
# ## How the Configuration Dictionary Works
#
# The helper function `build_radio_stripe_config()` (defined in
# `examples/training_school/utils.py`) takes the room dimensions and stripe
# topology you defined above and returns a **plain Python dictionary**. This
# dict has the same structure as the YAML files stored in `environments/`, so
# you can either use it directly in Python or save it to disk for later.
#
# ### Key sections of the config dict
#
# | Key | Contents |
# |---|---|
# | `room` | `x`, `y`, `z` — room dimensions in metres. |
# | `radio_stripes` | A **list of lists**. Each inner list represents one stripe: the first entry is the Central Unit (`{"central_unit": {"x":…, "y":…, "z":…}}`), followed by one dict per Radio Unit (`{"radio_unit": {"x":…, "y":…, "z":…}}`). |
# | `stripe_config` | Metadata about the topology: number of stripes/RUs, spacing, direction. |
# | `sub_thz` | Physical-layer parameters: carrier frequency (`fc`), bandwidth (`bw`), number of subcarriers, Doppler setting. Default: **157.75 GHz** carrier, **3 GHz** bandwidth, **4096** subcarriers. |
# | `antenna` | Antenna model per RU: number of elements (4), pattern (`TR38901`), polarisation (`V`). |
# | `ue_positions` | List of `{"x":…, "y":…, "z":…}` dicts — one per UE. Initially empty; we add them in Step 2b. |
#
# The function only populates the infrastructure side (room + stripes).
# UE positions are added separately so you can swap in different user
# distributions without regenerating the entire config.

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

# %% [markdown]
# ## YAML Configs — The Simulator's Common Language
#
# Every simulation script in this repository loads its configuration from
# **YAML files** (`.yaml`). YAML is a human-readable data format — think of it
# as a cleaner version of JSON with indentation instead of braces.
#
# The config dictionary we just built in Python has the **exact same structure**
# as the YAML files in the `environments/` folder. So you can either:
#
# - Use the dict directly in Python (as we did above), **or**
# - Save it to a `.yaml` file and load it in any other script with
#   `yaml.safe_load()`.
#
# This makes configs **portable** and **version-controllable** — you can share
# a YAML file with a colleague and they can reproduce your exact deployment.
#
# ### Component config coupling
#
# Besides the environment config, the simulator also needs a **component config**
# that describes the hardware models (amplifier gain, fiber length, coupler
# damping, etc.). One important coupling: the **fiber length** between adjacent
# RUs should match the **RU spacing** you chose for the deployment. When
# `UPDATE_COMPONENT_FIBER_LENGTH = True`, this script automatically updates the
# component config so these stay in sync.

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

# %% [markdown]
# ## From Bits to Air: The OFDM Transmit Chain
#
# Now that we have a deployment, let's generate a signal and transmit it. The
# simulator uses **CP-OFDM** (Cyclic-Prefix Orthogonal Frequency Division
# Multiplexing), the same modulation scheme used in 4G/5G and proposed for 6G.
#
# ### OFDM in a nutshell
#
# Instead of sending one high-rate data stream on a single carrier, OFDM
# splits the available bandwidth into thousands of narrow **subcarriers**.
# Each subcarrier carries a low-rate QAM symbol in parallel. This makes the
# system robust against frequency-selective fading: if some subcarriers are
# in a deep fade, others can still carry data. A **cyclic prefix (CP)** is
# prepended to each OFDM symbol to absorb multipath delay spread and prevent
# inter-symbol interference.
#
# ### The transmit processing chain
#
# ```
#  Random bits  →  QAM mapping  →  OFDM subcarrier allocation  →  IFFT + CP  →  time-domain signal
# ```
#
# 1. **Bit generation** — a random bit stream (the "payload").
# 2. **QAM mapping** — groups of $\log_2(M)$ bits are mapped to complex
#    constellation points (e.g. QPSK has $M{=}4$, so 2 bits per symbol).
# 3. **Subcarrier allocation** — QAM symbols and known **pilot** symbols are
#    placed onto an OFDM resource grid.
# 4. **IFFT + cyclic prefix** — the frequency-domain grid is converted to a
#    time-domain waveform by an Inverse FFT, and the CP is prepended.
#
# ### The `Waveform` object
#
# In the simulator, all of this is encapsulated in the `Waveform` class
# (from `wireless_channel/waveforms.py`). You create one from a config dict:
#
# ```python
# wf = Waveform.from_config(waveform_config, freq_band_config)
# ```
#
# Key waveform parameters (from `tutorial_waveform_config.yaml`):
#
# | Parameter | Default | Meaning |
# |---|---|---|
# | `qam_order` | 4 (QPSK) | Constellation size $M$ |
# | `n_ofdm_symbols` | 2 | Number of OFDM symbols per transmission |
# | `cp_length` | 200 | Cyclic prefix length in samples |
# | `pilot_spacing` | 4 | Insert a pilot every N subcarriers |
# | `oversampling_factor` | 4 | Oversampling for analog signal representation |
# | `tx_power` | 30 dBm | Transmit power |
#
# ### What the plots below show
#
# - **Bit sequence** — the raw 0/1 payload (first 100 bits).
# - **QAM constellation** — the complex symbols; for QPSK you should see 4
#   clusters at (±1, ±1).
# - **OFDM time-domain waveform** — looks noise-like because many subcarriers
#   add up (central limit theorem).
# - **Power Spectral Density (PSD)** — should be flat within the occupied
#   bandwidth and roll off sharply outside.

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

# %% [markdown]
# ## Uplink Signal Flow: UE → Channel → RadioStripe
#
# With the OFDM signal ready, we now simulate the **complete uplink path**.
# Here is what happens, step by step:
#
# ### 1. UE transmit chain
#
# The UE side is modelled with two simulator objects:
#
# - **`CentralUnit`** — represents the UE's baseband processor. In this
#   simplified model it is essentially a pass-through, but it is the
#   place where digital pre-processing (e.g. DPD) would be added.
# - **`RadioUnit`** — models the UE's analog front-end: power amplifier,
#   splitter (to multiple antenna elements), phase shifters (for beam
#   steering), and coupler. The `transmit()` method takes the baseband
#   waveform and returns the IQ signal as it would appear at the antenna.
#
# The phase shifts are set to **zero** here (broadside beam — no steering).
# You can experiment with `UE_BEAM` to steer the beam.
#
# ### 2. Wireless channel
#
# The `build_channel()` factory creates a channel object. Two models are
# available:
#
# | Model | Description |
# |---|---|
# | `"los"` | Deterministic **Line-of-Sight** channel: free-space path loss + distance-dependent phase per subcarrier. Fast, predictable, good for debugging. |
# | `"sionna"` | **Ray-traced** channel from pre-computed Sionna datasets stored in `wireless_channel/sionna_dataset/`. Includes multipath reflections, realistic for performance evaluation. |
#
# `channel.transmit_ul_id()` takes the UE's IQ signal and returns the received
# signal at a specific (stripe, RU) pair, applying the channel response
# (amplitude and phase per subcarrier per antenna pair).
#
# ### 3. Active RU selection
#
# In a RadioStripe, the signal enters the daisy chain at the **active RU** —
# the one closest to the transmitting UE. From there it travels hop by hop
# (through fiber and couplers) to the Central Unit. The helper
# `select_active_radio_unit()` picks the RU with the shortest Euclidean
# distance to the UE.
#
# ### 4. Stripe reception
#
# `stripe.receive()` routes the received signal from the active RU through
# the daisy chain (couplers, fiber, amplifiers) to the Central Unit. The
# output is the time-domain signal at the CU output.
#
# ### 5. OFDM demodulation and equalisation
#
# The received time-domain signal is processed back to the frequency domain:
#
# ```
# time signal → FFT → remove CP → extract subcarriers → LS channel estimate → one-tap equalisation
# ```
#
# - **LS (Least Squares) channel estimation** — uses the known pilot symbols
#   to estimate the channel response $H[k]$ at pilot subcarriers, then
#   interpolates to all data subcarriers.
# - **One-tap equalisation** — divides each received subcarrier by the
#   estimated channel: $\hat{X}[k] = Y[k] / \hat{H}[k]$.
#
# After equalisation the received symbols should cluster back around the
# original QAM constellation points. The quality is measured by:
#
# - **BER** (Bit Error Rate) — fraction of bits that flipped.
# - **EVM** (Error Vector Magnitude) — RMS distance between received and ideal
#   constellation points, expressed as a percentage.

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
y, imdata_stripe = stripe.receive(iq_data_rx, ru_shift)

y_combined_freq = wf.ofdm_time_to_freq(y)
rx_subc = wf.extract_subcarriers(y_combined_freq)

# Plot received pilots first to visualise the raw received quality.
wf.plot_constellation(rx_subc[0, wf.pilot_indices], title="Received pilots (before equalization)")

# LS channel estimation + one-tap equalization + constellation plot.
H_est = wf.channel_estimate_ls(rx_subc, interp_mode="phase")
eq_subc = wf.equalize_one_tap(rx_subc, H_est)
wf.plot_constellation(eq_subc, symbols_tx=xsubc, title="Equalized received symbols")

# 6) Compute BER and EVM to quantify the link quality.
#    - Demap the equalized subcarriers back to QAM symbols, then to bits.
#    - Compare received bits against the original transmitted bits.
y_qam = wf.demap_data_from_grid(eq_subc).flatten()
y_bits = wf.qam_to_bits(y_qam)
ber = wf.compute_ber(bits, y_bits)

# EVM: root-mean-square distance between received and ideal constellation points,
# normalised by the average TX symbol power.  Lower is better.
evm_percent = wf.compute_evm(xsubc, eq_subc)

print(f"\n{'='*40}")
print(f"  BER : {ber:.3e}")
print(f"  EVM : {evm_percent:.2f} %")
print(f"{'='*40}")

print(
    "\nUplink demo complete: "
    f"stripe={ACTIVE_STRIPE_IDX}, ru={ACTIVE_RU_IDX}, ue_index={UE_INDEX_FOR_UPLINK}, "
    f"ru_selection_mode={ACTIVE_RU_SELECTION_MODE if USE_AUTO_ACTIVE_RU else 'manual'}"
)

# %% [markdown]
# ## Zooming In on the Hardware Impairments
#
# In the previous step the signal travelled from the UE, through the wireless
# channel, and into the RadioStripe — but we treated the stripe as a black box.
# Let's open it up and look at what happens to the signal **inside each
# component** of each Radio Unit it passes through.
#
# ### How a signal travels through a RadioStripe (uplink)
#
# When the active RU receives the signal from the air, it goes through the
# **receive chain**:
#
# ```
# Antenna → Antenna Amplifier → Phase Shifter → Combiner → Coupler (in) → Fiber
# ```
#
# The signal then enters the **booster chain** of each subsequent RU on its
# way to the Central Unit:
#
# ```
# Fiber → Coupler (in) → Boost Amplifier → Coupler (out) → Fiber → … → CU
# ```
#
# Each component can distort the signal:
# - **Amplifiers** introduce nonlinearity (gain compression at high input
#   power) and additive noise.
# - **Couplers** introduce insertion loss and frequency-dependent filtering.
# - **Fiber** adds chromatic dispersion and attenuation.
# - **Phase shifters** apply quantised phase shifts (finite resolution).
# - **Combiners** sum the antenna branches (ideally lossless).
#
# ### AM/AM plots
#
# An **AM/AM (amplitude-to-amplitude) plot** shows the relationship between
# input and output signal amplitude for a given component. For an ideal
# component the plot is a straight line through the origin. Deviations reveal:
#
# - **Gain** — the slope of the line.  A slope > 1 means amplification.
# - **Compression** — the line bends downward at high input amplitudes,
#   meaning the component clips or saturates.
# - **Attenuation** — the line has slope < 1, meaning signal power is lost.
#
# Below we plot AM/AM curves for every component the signal traverses, from
# the active RU's receive chain through all booster RUs until the signal
# reaches the Central Unit.  The `stripe.receive()` call captured
# **intermediate data** (`imdata_stripe`) — a snapshot of the signal after
# each processing step — which we now visualise.

# %%
# ===========================================================================
# STEP 8 — Visualise hardware impairments per component per RU
# ===========================================================================

# ---- Parse imdata_stripe into labelled (input, output, title) triples ----
# imdata_stripe layout (from radiostripe.py & radio_unit.py):
#   [0]   = raw antenna input
#   Active RU receive chain (4 items):
#     [1] = antenna_amp output
#     [2] = phase_shifter output
#     [3] = combiner output
#     [4] = coupler_in output
#   [5]   = fiber output (after active RU)
#   Per booster RU (3 items each):
#     coupler_in output, boost_amp output, coupler_out output
#     then fiber output

active_ru_label = f"Active RU {ACTIVE_RU_IDX} (stripe {ACTIVE_STRIPE_IDX})"

# Collect all transitions with descriptive labels.
transitions = []

# Active RU receive chain
transitions.append((imdata_stripe[0], imdata_stripe[1], f"{active_ru_label} — Antenna Amplifier"))
transitions.append((imdata_stripe[1], imdata_stripe[2], f"{active_ru_label} — Phase Shifter"))
transitions.append((imdata_stripe[2], imdata_stripe[3], f"{active_ru_label} — Combiner"))
transitions.append((imdata_stripe[3], imdata_stripe[4], f"{active_ru_label} — Coupler (in)"))
transitions.append((imdata_stripe[4], imdata_stripe[5], f"{active_ru_label} — Fiber"))

# Booster RUs (from active_unit-1 down to 0, towards the CU)
offset = 6  # first booster data starts at index 6
booster_rus = list(range(ACTIVE_RU_IDX - 1, -1, -1))
for booster_count, ru_idx in enumerate(booster_rus):
    base = offset + booster_count * 4  # 3 component outputs + 1 fiber output
    if base + 3 >= len(imdata_stripe):
        break
    ru_label = f"Booster RU {ru_idx} (stripe {ACTIVE_STRIPE_IDX})"
    transitions.append((imdata_stripe[base - 1], imdata_stripe[base],     f"{ru_label} — Coupler (in)"))
    transitions.append((imdata_stripe[base],     imdata_stripe[base + 1], f"{ru_label} — Boost Amplifier"))
    transitions.append((imdata_stripe[base + 1], imdata_stripe[base + 2], f"{ru_label} — Coupler (out)"))
    transitions.append((imdata_stripe[base + 2], imdata_stripe[base + 3], f"{ru_label} — Fiber"))

print(f"\nPlotting {len(transitions)} AM/AM transitions through the stripe …")

# ---- Plot in a grid of subplots ----
n_plots = len(transitions)
n_cols = min(4, n_plots)
n_rows = int(np.ceil(n_plots / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 3.5 * n_rows))
# Flatten axes for easy indexing (handle single-row edge case).
if n_plots == 1:
    axes = np.array([axes])
axes = np.atleast_2d(axes).reshape(-1)

for idx, (sig_in, sig_out, title) in enumerate(transitions):
    plotter.plot_amam_transition(sig_in, sig_out, title, ax=axes[idx])

# Hide unused subplot slots.
for idx in range(n_plots, len(axes)):
    axes[idx].set_visible(False)

fig.suptitle("AM/AM per Component — Signal Path Through the RadioStripe", fontsize=12, y=1.02)
plt.tight_layout()
plt.show()

print(f"  Total components traversed: {len(transitions)}")
print(f"  Active RU: {ACTIVE_RU_IDX}  →  Booster RUs: {booster_rus}  →  Central Unit")

# %% [markdown]
# ## Power Spectral Density Through the Stripe
#
# The AM/AM plots above reveal amplitude distortion per component. Another
# useful diagnostic is the **Power Spectral Density (PSD)**: it shows how the
# signal's frequency content changes as it passes through each hardware stage.
#
# For each component we plot:
# - **Input PSD** (blue) — spectral shape entering the component.
# - **Output PSD** (orange) — spectral shape leaving the component.
#
# Things to look for:
# - **Flat in-band response** — ideally both curves overlap within the
#   occupied bandwidth.
# - **Spectral regrowth** — nonlinear amplifiers can widen the spectrum,
#   raising the out-of-band floor.
# - **Frequency-dependent loss** — couplers and fiber may attenuate some
#   frequencies more than others, tilting the PSD.
# - **Noise floor rise** — each active component adds noise, which shows up
#   as a higher baseline outside the signal bandwidth.

# %%
# ===========================================================================
# STEP 9 — PSD before and after each component
# ===========================================================================

print(f"\nPlotting {len(transitions)} PSD transitions through the stripe …")

n_plots_psd = len(transitions)
n_cols_psd = min(4, n_plots_psd)
n_rows_psd = int(np.ceil(n_plots_psd / n_cols_psd))

fig_psd, axes_psd = plt.subplots(n_rows_psd, n_cols_psd, figsize=(4.5 * n_cols_psd, 3.5 * n_rows_psd))
if n_plots_psd == 1:
    axes_psd = np.array([axes_psd])
axes_psd = np.atleast_2d(axes_psd).reshape(-1)

for idx, (sig_in, sig_out, title) in enumerate(transitions):
    plotter.plot_psd_transition(sig_in, sig_out, title, fs=wf.fs, ax=axes_psd[idx])

for idx in range(n_plots_psd, len(axes_psd)):
    axes_psd[idx].set_visible(False)

fig_psd.suptitle("PSD per Component — Signal Path Through the RadioStripe", fontsize=12, y=1.02)
plt.tight_layout()
plt.show()

print(f"  Total components traversed: {len(transitions)}")

# %% [markdown]
# ## Effect of RU Spacing on Amplifier Saturation
#
# In the plots above you may have noticed that the booster amplifiers show
# **gain compression** (the AM/AM curve bending away from the ideal line) and
# **spectral regrowth** (the PSD widening after the amplifier). Both are
# signatures of power-amplifier **saturation**: the amplifier is driven close
# to its maximum output power.
#
# Why does this happen? Between every pair of Radio Units the signal travels
# through an **optical fibre** segment. Fibre introduces **attenuation** that
# grows with length. With the original RU spacing of 3 m, the fibre loss is
# substantial, so each booster amplifier must apply a **high gain** to
# compensate — pushing it into its nonlinear (saturated) region.
#
# **What if we reduce the RU spacing?**
#
# Shorter fibre → less attenuation → the booster amplifier needs **less gain**
# to maintain the signal level → it operates further from saturation → the
# AM/AM curve stays closer to the ideal straight line and spectral regrowth
# is reduced.
#
# To show this properly we need to redo the **entire signal path**: regenerate
# the deployment with denser RU placement, re-select the closest active RU
# (which will now be at a different position), rebuild the wireless channel to
# that new RU, and receive through the updated stripe. Simply swapping the
# fiber length alone is not enough — the RU positions and active-RU selection
# must also reflect the new spacing.
#
# Below we rebuild the entire deployment with **RU_SPACING = 1 m** (instead
# of 3 m), keeping the same room and the same number of RUs per stripe. We
# then redo the active-RU selection (which may pick a different, closer RU),
# rebuild the wireless channel to that RU, transmit the same UE signal, and
# receive through the new stripe. The AM/AM and PSD plots let you directly
# compare the two cases.

# %%
# ===========================================================================
# STEP 10 — Reduced RU spacing: impact on amplifier saturation
# ===========================================================================

RU_SPACING_REDUCED = 1.0  # metres (was 3.0 m in the original deployment)

# --- 1) Regenerate environment config with denser RU spacing ----------------
env_config_reduced = build_radio_stripe_config(
    room_x=ROOM_X,
    room_y=ROOM_Y,
    room_z=ROOM_Z,
    n_stripes=N_STRIPES,
    n_rus_per_stripe=N_RUS_PER_STRIPE,
    ru_spacing_m=RU_SPACING_REDUCED,
    stripe_spacing_m=STRIPE_SPACING_M,
)
env_config_reduced["ue_positions"] = env_config["ue_positions"]  # reuse same UEs

# --- 2) Update component config with the shorter fiber length ---------------
component_config_reduced = dict(component_config)
component_config_reduced["fiber"] = dict(component_config.get("fiber", {}))
component_config_reduced["fiber"]["length"] = float(RU_SPACING_REDUCED)

print(f"\n{'=' * 60}")
print(f"Rebuilding deployment with RU spacing = {RU_SPACING_REDUCED} m "
      f"(fiber length = {RU_SPACING_REDUCED} m)")
print(f"{'=' * 60}")

# --- 3) Re-select the active RU for the same UE ----------------------------
ACTIVE_STRIPE_IDX_R, ACTIVE_RU_IDX_R, active_ru_pos_r, active_ru_distance_r = select_active_radio_unit(
    ue_position=ue_pos,
    radio_stripes=env_config_reduced["radio_stripes"],
    mode=ACTIVE_RU_SELECTION_MODE,
)
print(
    f"Re-selected active RU (reduced spacing): "
    f"stripe={ACTIVE_STRIPE_IDX_R}, ru={ACTIVE_RU_IDX_R}, "
    f"distance={active_ru_distance_r:.3f} m, pos={active_ru_pos_r}"
)

# --- 4) Build stripes with reduced fiber length ----------------------------
stripes_reduced: list[RadioStripe] = []
for stripe_cfg in env_config_reduced["radio_stripes"]:
    stripes_reduced.append(
        RadioStripe.from_config_locations(
            stripe_cfg,
            component_config_reduced,
            env_config_reduced["antenna"]["N_antennas"],
            wf,
        )
    )

# --- 5) Build channel to the new active RU and transmit --------------------
channel_reduced = build_channel(
    channel_model="los",
    ue_coordinates=ue_pos,
    sim_env="office_space_perpendicular",
    component_config=component_config_reduced,
    stripe_positions=env_config_reduced["radio_stripes"],
    waveform=wf,
    Nr_ue_antennas=env_config_reduced["antenna"]["N_antennas"],
    Nr_ru_antennas=env_config_reduced["antenna"]["N_antennas"],
    debug=False,
    los_normalize_gain=False,
)

iq_data_rx_r = channel_reduced.transmit_ul_id(iq_data_tx, ACTIVE_STRIPE_IDX_R, ACTIVE_RU_IDX_R, wf)

# --- 6) Receive through the reduced-spacing stripe -------------------------
stripe_reduced = stripes_reduced[ACTIVE_STRIPE_IDX_R]
stripe_reduced.active_unit = ACTIVE_RU_IDX_R
ru_shift_r = stripe_reduced.radio_units[0].phase_shifter.get_phases(RU_BEAM)
y_reduced, imdata_reduced = stripe_reduced.receive(iq_data_rx_r, ru_shift_r)

# --- Build transition list for the reduced-spacing stripe -------------------
active_ru_label_r = f"Active RU {ACTIVE_RU_IDX_R} ({RU_SPACING_REDUCED} m spacing)"

transitions_reduced = []

# Active RU receive chain
transitions_reduced.append((imdata_reduced[0], imdata_reduced[1], f"{active_ru_label_r} — Antenna Amp"))
transitions_reduced.append((imdata_reduced[1], imdata_reduced[2], f"{active_ru_label_r} — Phase Shifter"))
transitions_reduced.append((imdata_reduced[2], imdata_reduced[3], f"{active_ru_label_r} — Combiner"))
transitions_reduced.append((imdata_reduced[3], imdata_reduced[4], f"{active_ru_label_r} — Coupler (in)"))
transitions_reduced.append((imdata_reduced[4], imdata_reduced[5], f"{active_ru_label_r} — Fiber"))

# Booster RUs
offset_r = 6
booster_rus_r = list(range(ACTIVE_RU_IDX_R - 1, -1, -1))
for booster_count, ru_idx in enumerate(booster_rus_r):
    base = offset_r + booster_count * 4
    if base + 3 >= len(imdata_reduced):
        break
    ru_label_r = f"Booster RU {ru_idx} ({RU_SPACING_REDUCED} m spacing)"
    transitions_reduced.append((imdata_reduced[base - 1], imdata_reduced[base],     f"{ru_label_r} — Coupler (in)"))
    transitions_reduced.append((imdata_reduced[base],     imdata_reduced[base + 1], f"{ru_label_r} — Boost Amp"))
    transitions_reduced.append((imdata_reduced[base + 1], imdata_reduced[base + 2], f"{ru_label_r} — Coupler (out)"))
    transitions_reduced.append((imdata_reduced[base + 2], imdata_reduced[base + 3], f"{ru_label_r} — Fiber"))

# --- AM/AM plots for reduced spacing ---------------------------------------
print(f"\nPlotting {len(transitions_reduced)} AM/AM transitions (RU spacing = {RU_SPACING_REDUCED} m) …")

n_plots_r = len(transitions_reduced)
n_cols_r = min(4, n_plots_r)
n_rows_r = int(np.ceil(n_plots_r / n_cols_r))

fig_amam_r, axes_amam_r = plt.subplots(n_rows_r, n_cols_r, figsize=(4.5 * n_cols_r, 3.5 * n_rows_r))
if n_plots_r == 1:
    axes_amam_r = np.array([axes_amam_r])
axes_amam_r = np.atleast_2d(axes_amam_r).reshape(-1)

for idx, (sig_in, sig_out, title) in enumerate(transitions_reduced):
    plotter.plot_amam_transition(sig_in, sig_out, title, ax=axes_amam_r[idx])

for idx in range(n_plots_r, len(axes_amam_r)):
    axes_amam_r[idx].set_visible(False)

fig_amam_r.suptitle(
    f"AM/AM per Component — RU spacing = {RU_SPACING_REDUCED} m (reduced fiber loss)",
    fontsize=12, y=1.02,
)
plt.tight_layout()
plt.show()

# --- PSD plots for reduced spacing ------------------------------------------
print(f"Plotting {len(transitions_reduced)} PSD transitions (RU spacing = {RU_SPACING_REDUCED} m) …")

fig_psd_r, axes_psd_r = plt.subplots(n_rows_r, n_cols_r, figsize=(4.5 * n_cols_r, 3.5 * n_rows_r))
if n_plots_r == 1:
    axes_psd_r = np.array([axes_psd_r])
axes_psd_r = np.atleast_2d(axes_psd_r).reshape(-1)

for idx, (sig_in, sig_out, title) in enumerate(transitions_reduced):
    plotter.plot_psd_transition(sig_in, sig_out, title, fs=wf.fs, ax=axes_psd_r[idx])

for idx in range(n_plots_r, len(axes_psd_r)):
    axes_psd_r[idx].set_visible(False)

fig_psd_r.suptitle(
    f"PSD per Component — RU spacing = {RU_SPACING_REDUCED} m (reduced fiber loss)",
    fontsize=12, y=1.02,
)
plt.tight_layout()
plt.show()

# --- Compare link quality ---------------------------------------------------
y_combined_freq_r = wf.ofdm_time_to_freq(y_reduced)
rx_subc_r = wf.extract_subcarriers(y_combined_freq_r)
H_est_r = wf.channel_estimate_ls(rx_subc_r, interp_mode="phase")
eq_subc_r = wf.equalize_one_tap(rx_subc_r, H_est_r)
y_qam_r = wf.demap_data_from_grid(eq_subc_r).flatten()
y_bits_r = wf.qam_to_bits(y_qam_r)
ber_r = wf.compute_ber(bits, y_bits_r)
evm_r = wf.compute_evm(xsubc, eq_subc_r)

# --- Room plot with the new active RU highlighted --------------------------
if plot_room:
    plotter.plot_room_with_active_selection(
        env_config_reduced,
        active_ue_pos=ue_pos,
        active_ru_pos=active_ru_pos_r,
    )
    plt.title(f"Reduced RU spacing ({RU_SPACING_REDUCED} m) — active RU selection")
    plt.show()

# --- Received constellation (before and after equalisation) -----------------
wf.plot_constellation(
    rx_subc_r[0, wf.pilot_indices],
    title=f"Received pilots (before EQ) — RU spacing {RU_SPACING_REDUCED} m",
)
wf.plot_constellation(
    eq_subc_r,
    symbols_tx=xsubc,
    title=f"Equalized received symbols — RU spacing {RU_SPACING_REDUCED} m",
)

print(f"\n{'=' * 50}")
print(f"  Comparison: RU spacing {RU_SPACING_M} m  vs  {RU_SPACING_REDUCED} m")
print(f"  Original active RU: stripe={ACTIVE_STRIPE_IDX}, ru={ACTIVE_RU_IDX}")
print(f"  Reduced  active RU: stripe={ACTIVE_STRIPE_IDX_R}, ru={ACTIVE_RU_IDX_R}")
print(f"  {'Metric':<10} {'Original':>12} {'Reduced':>12}")
print(f"  {'BER':<10} {ber:>12.3e} {ber_r:>12.3e}")
print(f"  {'EVM (%)':<10} {evm_percent:>12.2f} {evm_r:>12.2f}")
print(f"{'=' * 50}")

# %% [markdown]
# ## What's Next?
#
# You have now completed a full end-to-end walkthrough:
#
# 1. **Created** a RadioStripe deployment from scratch (room + stripes + UEs).
# 2. **Saved** the configuration to a reusable YAML file.
# 3. **Generated** an OFDM uplink signal (bits → QAM → OFDM).
# 4. **Transmitted** it through a LOS wireless channel to the nearest RU.
# 5. **Received** and equalised the signal, visualising the constellation.
#
# ### Suggested exercises
#
# - **Change the room geometry** — make the room larger or smaller and observe
#   how the automatic RU selection changes (Step 7).
# - **Add more stripes or RUs** — increase `N_STRIPES` or `N_RUS_PER_STRIPE`
#   and see how the deployment density affects coverage.
# - **Move the UE** — change `UE_INDEX_FOR_UPLINK` to pick a different user
#   and watch how the channel quality (constellation scatter) changes with
#   distance.
# - **Try beam steering** — set `UE_BEAM` or `RU_BEAM` to a non-zero value
#   to steer the antenna beam.
# - **Switch to the Sionna channel** — in `tutorial.py`, change
#   `channel_model = "sionna"` to use ray-traced multipath channels. This
#   requires the Sionna dataset files in `wireless_channel/sionna_dataset/`.
# - **Compare ideal vs. impaired hardware** — see `tutorial.py` (Sections 4
#   and 5) for a side-by-side comparison of ideal amplifiers vs. realistic
#   hardware models with noise, nonlinearity, and fiber/coupler effects.
#
# ### Further reading
#
# - `tutorial.py` — full uplink tutorial with ideal vs. impaired comparison.
# - `basic_ul_tutorial.ipynb` — interactive notebook version with Sionna
#   channel support.
# - `environments/` — pre-built deployment configs for different room types.

print("\nDone — see tutorial.py for the next steps.")