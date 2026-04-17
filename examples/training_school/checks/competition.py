# %% [markdown]
# # 🏆 RadioStripe Deployment Competition
#
# **Goal:** Optimize a sub-THz RadioStripe deployment to **maximize uplink
# coverage** across a room full of users.
#
# You just completed the tutorial and learned how the simulator works end-to-end.
# Now it's time to put that knowledge to work! Your task is to design the best
# possible RadioStripe deployment for a fixed room and a fixed set of UE
# positions.
#
# ## Rules
#
# 1. **Fixed room** — the room dimensions are given and **cannot be changed**.
# 2. **Fixed waveform** — modulation (QPSK), bandwidth (3 GHz), and number of
#    OFDM symbols are fixed and **cannot be changed**.
# 3. **Fixed UE positions** — you optimize on a *training set* of UEs (provided
#    below). The final ranking will be evaluated on a **hidden test set** with a
#    different random seed.
# 4. **RU budget** — you may use **at most `MAX_TOTAL_RUS` Radio Units** in
#    total across all stripes.
# 5. **Fiber length must match RU spacing** — the `fiber.length` in the
#    component config must equal the physical spacing between adjacent RUs on
#    each stripe. No cheating with zero-length fibers!
#
# ## What You Can Optimize
#
# | Knob | How | Trade-off |
# |---|---|---|
# | **Stripe / RU geometry** | Change `n_stripes`, `n_rus_per_stripe`, `ru_spacing_m`, `stripe_spacing_m` | Coverage area vs. component count |
# | **RU spacing** | Adjust `ru_spacing_m` | Shorter fiber → less PA saturation, but fewer RUs cover less area |
# | **TX power** | Change `TX_POWER_DBM` in the waveform config | Higher power → better SNR but more PA distortion |
# | **Beam steering** | Set `UE_BEAM` and `RU_BEAM` per UE | Directivity gain vs. alignment risk |
# | **Active RU selection** | Implement `custom_select_active_ru()` | Distance isn't always optimal — fewer booster hops = less distortion |
# | **Creative ideas** | Anything else within the rules! | Surprise us! |
#
# ## Scoring
#
# Your deployment is evaluated across **all UE positions**. For each UE the
# full uplink chain is simulated and the EVM (Error Vector Magnitude) is
# recorded. Two leaderboards rank participants:
#
# ### Leaderboard 1 — Signal Quality Score
# $$\text{Score}_1 = \text{mean}(\text{EVM}) + \lambda \cdot \max(\text{EVM})$$
# where $\lambda = 0.5$. This penalizes both poor average quality **and**
# leaving any single UE with terrible reception. **Lower is better.**
#
# ### Leaderboard 2 — Cost-Aware Score
# $$\text{Score}_2 = \text{mean}(\text{EVM}) + \lambda \cdot \max(\text{EVM}) + \mu \cdot N_{\text{RUs}}$$
# where $\lambda = 0.5$ and $\mu = 0.5$. This additionally penalizes using
# more hardware — rewarding efficient deployments. **Lower is better.**
#
# ## How to Submit
#
# 1. Edit **only** the `# === YOUR CODE HERE ===` sections below.
# 2. Run the full script to see your scores.
# 3. Report your two scores and your team name.
#
# Good luck!

# %%
# ===========================================================================
# Boilerplate — DO NOT MODIFY
# ===========================================================================
import os
import sys
from pathlib import Path
import copy

import matplotlib.pyplot as plt
import numpy as np
import yaml

# Resolve repository root.
root = Path(__file__).resolve().parent
while not (root / "wireless_channel").exists() and root != root.parent:
    root = root.parent
if not (root / "wireless_channel").exists():
    raise RuntimeError("Could not locate the repository root.")
os.chdir(root)
script_dir = str(Path(__file__).resolve().parent)
sys.path = [p for p in sys.path if p not in ("", script_dir)]
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

from examples.training_school.utils import (
    build_radio_stripe_config,
    generate_grid_ue_positions,
    generate_random_ue_positions,
    select_active_radio_unit,
)
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
# FIXED PARAMETERS — DO NOT MODIFY
# ===========================================================================

# --- Room (fixed) -----------------------------------------------------------
ROOM_X = 10.0   # metres
ROOM_Y = 20.0   # metres
ROOM_Z = 3.5    # metres

# --- RU budget (fixed) ------------------------------------------------------
MAX_TOTAL_RUS = 15  # maximum number of Radio Units across all stripes

# --- Scoring weights (fixed) ------------------------------------------------
LAMBDA = 0.5   # weight for max(EVM)
MU = 0.5       # weight for N_total_RUs (leaderboard 2 only)

# --- Training UE set (fixed) ------------------------------------------------
N_UES_TRAINING = 20
UE_Z = 1.5
UE_WALL_MARGIN = 0.1
UE_SEED_TRAINING = 31   # known seed — you can optimize on this
# UE_SEED_TEST = ???     # hidden seed — used for final ranking

# --- Grid UE set (fixed) ---------------------------------------------------
N_GRID_X = 5             # grid points along X
N_GRID_Y = 10            # grid points along Y

# --- Waveform (fixed) -------------------------------------------------------
configs_dir = Path(__file__).parent / "configs"
with open(configs_dir / "tutorial_waveform_config.yaml", "r", encoding="utf8") as _f:
    WAVEFORM_CONFIG_BASE = yaml.safe_load(_f)

# %% [markdown]
# ## ✏️ Your Solution
#
# Edit the sections below. You may change anything inside the
# `# === YOUR CODE HERE ===` blocks. Do **not** change the fixed parameters
# or the evaluation loop.

# %%
# ===========================================================================
# SECTION A — Deployment geometry
# ===========================================================================
# Choose how many stripes and RUs to use (within the budget), their spacing,
# and any other geometry parameters.
#
# Constraint: N_STRIPES * N_RUS_PER_STRIPE <= MAX_TOTAL_RUS
#
# === YOUR CODE HERE === (modify the values below) ===========================

TEAM_NAME = "Team Default"

N_STRIPES = 2
N_RUS_PER_STRIPE = 7          # 3 * 5 = 15 RUs total (max budget)
RU_SPACING_M = 2.0            # metres between adjacent RUs
STRIPE_SPACING_M = 4        # metres between stripes
TX_POWER_DBM = 20#30          # transmit power in dBm

# === END YOUR CODE ===========================================================

# --- Validate RU budget -----------------------------------------------------
N_TOTAL_RUS = N_STRIPES * N_RUS_PER_STRIPE
assert N_TOTAL_RUS <= MAX_TOTAL_RUS, (
    f"Budget exceeded! {N_TOTAL_RUS} RUs used, but maximum is {MAX_TOTAL_RUS}."
)

# %%
# ===========================================================================
# SECTION B — Active RU selection strategy
# ===========================================================================
# The default uses distance-based selection (closest RU to the UE).
# You may implement a smarter strategy here.
#
# Your function receives:
#   - ue_pos : dict with keys "x", "y", "z"
#   - radio_stripes : the stripe config (list of lists)
#   - component_config : the component config dict
#
# It must return:
#   - (stripe_idx, ru_idx)
#
# === YOUR CODE HERE === (edit the function body) =============================

def custom_select_active_ru(ue_pos, radio_stripes, component_config):
    """Select the active RU for a given UE position.

    Default: pick the closest RU by Euclidean distance.
    You can replace this with any strategy you like.
    """
    stripe_idx, ru_idx, _, _ = select_active_radio_unit(
        ue_position=ue_pos,
        radio_stripes=radio_stripes,
        mode="distance",
    )
    return stripe_idx, ru_idx

# === END YOUR CODE ===========================================================

# %%
# ===========================================================================
# SECTION C — Beam steering strategy
# ===========================================================================
# Return (ue_beam, ru_beam) for each UE. Default: broadside (0, 0).
#
# === YOUR CODE HERE === (edit the function body) =============================

BEAM_MODE = "zero"  # "zero", "random", or "custom"

def get_beam_indices(ue_pos, active_ru_pos):
    """Return (ue_beam_angle, ru_beam_angle) for a given UE-RU pair.

    Modes:
    - "zero":   broadside beam (0° steering) for both UE and RU.
    - "random": random beam angles in [-90, 90] degrees for both.
    - "custom": user-defined logic (edit the custom branch below).
    """
    if BEAM_MODE == "zero":
        return 0, 0
    elif BEAM_MODE == "random":
        ue_beam = np.random.uniform(-90, 90)
        ru_beam = np.random.uniform(-90, 90)
        return ue_beam, ru_beam
    elif BEAM_MODE == "custom":
        # === YOUR CUSTOM BEAM LOGIC HERE ===
        raise NotImplementedError("Implement your custom beam steering logic!")
    else:
        raise ValueError(f"Unknown BEAM_MODE: '{BEAM_MODE}'. Use 'zero', 'random', or 'custom'.")

# === END YOUR CODE ===========================================================

# %% [markdown]
# ## Evaluation Loop — DO NOT MODIFY
#
# The code below builds your deployment, loops over every UE position,
# simulates the full uplink chain, and computes the two scores.

# %%
# ===========================================================================
# BUILD DEPLOYMENT
# ===========================================================================

# --- Waveform config with your TX power ------------------------------------
waveform_config = copy.deepcopy(WAVEFORM_CONFIG_BASE)
waveform_config["tx_power"] = TX_POWER_DBM

# --- Environment config from your geometry ----------------------------------
env_config = build_radio_stripe_config(
    room_x=ROOM_X,
    room_y=ROOM_Y,
    room_z=ROOM_Z,
    n_stripes=N_STRIPES,
    n_rus_per_stripe=N_RUS_PER_STRIPE,
    ru_spacing_m=RU_SPACING_M,
    stripe_spacing_m=STRIPE_SPACING_M,
)

# --- UE positions (training set — random) ----------------------------------
ue_positions_random = generate_random_ue_positions(
    n_users=N_UES_TRAINING,
    room_x=ROOM_X,
    room_y=ROOM_Y,
    z_height=UE_Z,
    wall_margin=UE_WALL_MARGIN,
    seed=UE_SEED_TRAINING,
)

# --- UE positions (grid — uniformly sampled) --------------------------------
ue_positions_grid = generate_grid_ue_positions(
    n_x=N_GRID_X,
    n_y=N_GRID_Y,
    room_x=ROOM_X,
    room_y=ROOM_Y,
    z_height=UE_Z,
    wall_margin=UE_WALL_MARGIN,
)

N_RANDOM_UES = len(ue_positions_random)
ue_positions = ue_positions_random + ue_positions_grid
env_config["ue_positions"] = ue_positions

# --- Component config (fiber length must match RU spacing) ------------------
with open(configs_dir / "tutorial_component_config.yaml", "r", encoding="utf8") as _f:
    component_config = yaml.safe_load(_f)
component_config["fiber"]["length"] = float(RU_SPACING_M)

# --- Build waveform and generate TX signal once ----------------------------
freq_band_config = env_config["sub_thz"]
wf = Waveform.from_config(waveform_config, freq_band_config)

bits = wf.generate_bits()
qam = wf.qam_modulate()
ofdm_time = wf.ofdm_modulate()
x_combined_freq = wf.ofdm_time_to_freq(ofdm_time)
xsubc = wf.extract_subcarriers(x_combined_freq)

# --- Build UE transmit chain (fixed) ---------------------------------------
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
    x=0, y=0, z=0,
    boost_amp=amp, antenna_amp=amp,
    coup_in=coup, coup_out=coup,
    splitter=split, combiner=comb, pshift=ps,
)

# --- Build all stripes from your config ------------------------------------
stripes: list[RadioStripe] = []
for stripe_cfg in env_config["radio_stripes"]:
    stripes.append(
        RadioStripe.from_config_locations(
            stripe_cfg, component_config,
            env_config["antenna"]["N_antennas"], wf,
        )
    )

# --- Plot the room layout ---------------------------------------------------
plotter.plot_room(env_config)
plt.title(f"{TEAM_NAME} — Deployment ({N_TOTAL_RUS} RUs)")
plt.show()

print(f"\n{'=' * 60}")
print(f"  Team:     {TEAM_NAME}")
print(f"  Stripes:  {N_STRIPES}")
print(f"  RUs/stripe: {N_RUS_PER_STRIPE}  (total: {N_TOTAL_RUS})")
print(f"  RU spacing: {RU_SPACING_M} m")
print(f"  TX power:   {TX_POWER_DBM} dBm")
print(f"  UEs:        {len(ue_positions)} ({N_RANDOM_UES} random + {len(ue_positions_grid)} grid)")
print(f"{'=' * 60}")



# %%
# ===========================================================================
# EVALUATE ACROSS ALL UEs
# ===========================================================================

evm_per_ue = []
ber_per_ue = []

for ue_idx, ue_pos in enumerate(ue_positions):
    # --- Active RU selection (your function) --------------------------------
    stripe_idx, ru_idx = custom_select_active_ru(
        ue_pos, env_config["radio_stripes"], component_config,
    )

    # --- Beam steering (your function) --------------------------------------
    ru_entry = env_config["radio_stripes"][stripe_idx][ru_idx + 1]["radio_unit"]
    ue_beam, ru_beam = get_beam_indices(ue_pos, ru_entry)

    # --- UE transmit --------------------------------------------------------
    ue_shift = ue_ru.phase_shifter.get_phases(ue_beam)
    iq_data_tx, _ = ue_ru.transmit(ofdm_time_after_ue, ue_shift)

    # --- Wireless channel ---------------------------------------------------
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
    iq_data_rx = channel.transmit_ul_id(iq_data_tx, stripe_idx, ru_idx, wf)

    # --- Stripe receive -----------------------------------------------------
    stripe = stripes[stripe_idx]
    stripe.active_unit = ru_idx
    ru_shift = stripe.radio_units[0].phase_shifter.get_phases(ru_beam)
    y, _ = stripe.receive(iq_data_rx, ru_shift)

    # --- OFDM demodulation + equalization -----------------------------------
    y_combined_freq = wf.ofdm_time_to_freq(y)
    rx_subc = wf.extract_subcarriers(y_combined_freq)
    H_est = wf.channel_estimate_ls(rx_subc, interp_mode="phase")
    eq_subc = wf.equalize_one_tap(rx_subc, H_est)

    # --- Metrics ------------------------------------------------------------
    evm = wf.compute_evm(xsubc, eq_subc)
    y_qam = wf.demap_data_from_grid(eq_subc).flatten()
    y_bits = wf.qam_to_bits(y_qam)
    ber_val = wf.compute_ber(bits, y_bits)

    evm_per_ue.append(evm)
    ber_per_ue.append(ber_val)

    ue_type = "rand" if ue_idx < N_RANDOM_UES else "grid"
    print(f"  UE {ue_idx:2d} ({ue_type})  stripe={stripe_idx} ru={ru_idx:2d}  "
          f"EVM={evm:7.2f}%  BER={ber_val:.3e}")

evm_arr = np.array(evm_per_ue)
ber_arr = np.array(ber_per_ue)

# %%
# ===========================================================================
# COMPUTE SCORES
# ===========================================================================

mean_evm = np.mean(evm_arr)
max_evm = np.max(evm_arr)
score_1 = mean_evm + LAMBDA * max_evm
score_2 = mean_evm + LAMBDA * max_evm + MU * N_TOTAL_RUS

print(f"\n{'=' * 60}")
print(f"  RESULTS — {TEAM_NAME}")
print(f"{'=' * 60}")
print(f"  UEs evaluated     : {len(evm_arr)}")
print(f"  Mean EVM           : {mean_evm:.2f} %")
print(f"  Max  EVM           : {max_evm:.2f} %")
print(f"  Min  EVM           : {np.min(evm_arr):.2f} %")
print(f"  Std  EVM           : {np.std(evm_arr):.2f} %")
print(f"  Mean BER           : {np.mean(ber_arr):.3e}")
print(f"  Total RUs          : {N_TOTAL_RUS}")
print(f"{'─' * 60}")
print(f"  Leaderboard 1 (quality):    {score_1:.4f}")
print(f"    = mean({mean_evm:.2f}) + {LAMBDA}·max({max_evm:.2f})")
print(f"  Leaderboard 2 (cost-aware): {score_2:.4f}")
print(f"    = mean({mean_evm:.2f}) + {LAMBDA}·max({max_evm:.2f}) + {MU}·{N_TOTAL_RUS}")
print(f"{'=' * 60}")

# %%
# ===========================================================================
# VISUALIZATION
# ===========================================================================

# --- EVM per UE bar chart ---------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 4))
colors = plt.cm.RdYlGn_r(evm_arr / max(evm_arr.max(), 1))  # red=bad, green=good
ax.bar(range(len(evm_arr)), evm_arr, color=colors, edgecolor="gray", linewidth=0.5)
ax.axhline(mean_evm, color="blue", linestyle="--", linewidth=1.2, label=f"Mean EVM = {mean_evm:.2f}%")
ax.axhline(max_evm, color="red", linestyle="--", linewidth=1.2, label=f"Max EVM = {max_evm:.2f}%")
ax.set_xlabel("UE index")
ax.set_ylabel("EVM (%)")
ax.set_title(f"{TEAM_NAME} — EVM per UE  |  Score₁ = {score_1:.2f}  |  Score₂ = {score_2:.2f}")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# --- Room plot coloured by EVM per UE ---------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

# Plot stripes and RUs
for s_idx, stripe_cfg in enumerate(env_config["radio_stripes"]):
    cu = stripe_cfg[0]["central_unit"]
    ax.plot(cu["x"], cu["y"], "ks", markersize=8, markerfacecolor="none")
    for ru_entry in stripe_cfg[1:]:
        ru = ru_entry["radio_unit"]
        ax.plot(ru["x"], ru["y"], "b^", markersize=6)

# Plot UEs coloured by EVM — random (circles) and grid (squares)
sc_rand = ax.scatter(
    [u["x"] for u in ue_positions[:N_RANDOM_UES]],
    [u["y"] for u in ue_positions[:N_RANDOM_UES]],
    c=evm_arr[:N_RANDOM_UES], cmap="RdYlGn_r", s=80, marker="o",
    edgecolors="black", linewidths=0.5, zorder=5,
    vmin=evm_arr.min(), vmax=evm_arr.max(),
    label="Random UE",
)
sc_grid = ax.scatter(
    [u["x"] for u in ue_positions[N_RANDOM_UES:]],
    [u["y"] for u in ue_positions[N_RANDOM_UES:]],
    c=evm_arr[N_RANDOM_UES:], cmap="RdYlGn_r", s=80, marker="s",
    edgecolors="black", linewidths=0.5, zorder=5,
    vmin=evm_arr.min(), vmax=evm_arr.max(),
    label="Grid UE",
)
for i, ue in enumerate(ue_positions):
    ax.annotate(str(i), (ue["x"], ue["y"]), fontsize=7, ha="center", va="bottom",
                xytext=(0, 5), textcoords="offset points")

cbar = plt.colorbar(sc_rand, ax=ax)
cbar.set_label("EVM (%)")
margin = 0.5
ax.set_xlim(0 - margin, ROOM_X + margin)
ax.set_ylim(0 - margin, ROOM_Y + margin)
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_title(f"{TEAM_NAME} — Coverage Map  |  ▲ = RU, □ = CU")
ax.set_aspect("equal")
ax.legend(loc="upper right", fontsize=8)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# --- CDF of EVM ------------------------------------------------------------
evm_sorted = np.sort(evm_arr)
cdf = np.arange(1, len(evm_sorted) + 1) / len(evm_sorted)

fig, ax = plt.subplots(figsize=(8, 4))
ax.step(evm_sorted, cdf, where="post", linewidth=2)
ax.axvline(mean_evm, color="blue", linestyle="--", label=f"Mean = {mean_evm:.2f}%")
ax.set_xlabel("EVM (%)")
ax.set_ylabel("CDF")
ax.set_title(f"{TEAM_NAME} — EVM CDF")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print(f"\n🏆 Final Scores for {TEAM_NAME}:")
print(f"   Leaderboard 1 (quality):     {score_1:.4f}")
print(f"   Leaderboard 2 (cost-aware):  {score_2:.4f}")


# %%
# ===========================================================================
# INSPECT ACTIVE RU SELECTION FOR A SINGLE UE
# ===========================================================================
# Change UE_INDEX to inspect a different UE.

UE_INDEX = 1# <-- change this to inspect a different UE

_ue_pos = ue_positions[UE_INDEX]
_stripe_idx, _ru_idx = custom_select_active_ru(
    _ue_pos, env_config["radio_stripes"], component_config,
)
_ru_pos = env_config["radio_stripes"][_stripe_idx][_ru_idx + 1]["radio_unit"]

print(f"UE {UE_INDEX}: ({_ue_pos['x']:.2f}, {_ue_pos['y']:.2f}, {_ue_pos['z']:.2f})")
print(f"  → stripe {_stripe_idx}, RU {_ru_idx}: "
      f"({_ru_pos['x']:.2f}, {_ru_pos['y']:.2f}, {_ru_pos['z']:.2f})")

plotter.plot_room_with_active_selection(env_config, active_ue_pos=_ue_pos, active_ru_pos=_ru_pos)
plt.title(f"UE {UE_INDEX} → stripe {_stripe_idx}, RU {_ru_idx}")
plt.show()

# %%
# ===========================================================================
# INSPECT WIRELESS CHANNEL FOR A SINGLE UE
# ===========================================================================
# Uses the same UE_INDEX from above. Shows channel frequency response,
# impulse response, and path loss for the selected UE-RU link.

_channel = build_channel(
    channel_model="los",
    ue_coordinates=_ue_pos,
    sim_env="office_space_perpendicular",
    component_config=component_config,
    stripe_positions=env_config["radio_stripes"],
    waveform=wf,
    Nr_ue_antennas=env_config["antenna"]["N_antennas"],
    Nr_ru_antennas=env_config["antenna"]["N_antennas"],
    debug=False,
    los_normalize_gain=False,
)

# Get the channel frequency response for the active RU
H = _channel.get_csi(_stripe_idx, _ru_idx)  # (Nr_ue_ant, Nr_ru_ant, Nr_subcarriers)

# Use antenna 0 → antenna 0 for visualization
H_00 = H[0, 0, :]

# Frequency axis
fc = env_config["sub_thz"]["fc"]
bw = env_config["sub_thz"]["bw"]
n_sc = len(H_00)
freqs_ghz = (fc + np.linspace(-bw / 2, bw / 2, n_sc)) / 1e9

fig, axes = plt.subplots(2, 2, figsize=(14, 8))
fig.suptitle(f"Channel Inspection — UE {UE_INDEX} → stripe {_stripe_idx}, RU {_ru_idx}", fontsize=13)

# 1) Channel magnitude (frequency response)
ax = axes[0, 0]
ax.plot(freqs_ghz, 20 * np.log10(np.abs(H_00) + 1e-30), linewidth=0.8)
ax.set_xlabel("Frequency (GHz)")
ax.set_ylabel("|H| (dB)")
ax.set_title("Channel Frequency Response (magnitude)")
ax.grid(True, alpha=0.3)

# 2) Channel phase (frequency response)
ax = axes[0, 1]
ax.plot(freqs_ghz, np.unwrap(np.angle(H_00)), linewidth=0.8)
ax.set_xlabel("Frequency (GHz)")
ax.set_ylabel("Phase (rad)")
ax.set_title("Channel Frequency Response (phase)")
ax.grid(True, alpha=0.3)

# 3) Impulse response (IFFT of frequency response)
ax = axes[1, 0]
h_time = np.fft.ifft(H_00)
t_ns = np.arange(len(h_time)) / bw * 1e9
ax.stem(t_ns[:50], np.abs(h_time[:50]), markerfmt=".", basefmt="k-")
ax.set_xlabel("Delay (ns)")
ax.set_ylabel("|h(t)|")
ax.set_title("Channel Impulse Response (first 50 taps)")
ax.grid(True, alpha=0.3)

# 4) Path loss across all RUs on all stripes
ax = axes[1, 1]
for s_i, stripe_cfg in enumerate(env_config["radio_stripes"]):
    n_rus = len(stripe_cfg) - 1  # exclude CU entry
    pl_per_ru = []
    for r_i in range(n_rus):
        H_sr = _channel.get_csi(s_i, r_i)
        avg_gain = np.mean(np.abs(H_sr[0, 0, :]) ** 2)
        pl_db = -10 * np.log10(avg_gain + 1e-30)
        pl_per_ru.append(pl_db)
    ax.bar(
        [f"S{s_i}R{r}" for r in range(n_rus)],
        pl_per_ru,
        alpha=0.7,
        label=f"Stripe {s_i}",
    )
# Highlight the active RU
active_label = f"S{_stripe_idx}R{_ru_idx}"
xlabels = [t.get_text() for t in ax.get_xticklabels()]
ax.set_xlabel("Stripe / RU")
ax.set_ylabel("Path Loss (dB)")
ax.set_title("Path Loss to Each RU")
ax.tick_params(axis="x", rotation=45, labelsize=7)
ax.grid(True, alpha=0.3, axis="y")

plt.tight_layout()
plt.show()

# Print summary
avg_gain_active = np.mean(np.abs(H_00) ** 2)
pl_active = -10 * np.log10(avg_gain_active + 1e-30)
print(f"\nChannel summary for UE {UE_INDEX} → stripe {_stripe_idx}, RU {_ru_idx}:")
print(f"  Path loss:        {pl_active:.2f} dB")
print(f"  Mean |H|:         {np.mean(np.abs(H_00)):.6f}")
print(f"  |H| std:          {np.std(np.abs(H_00)):.6f}")
print(f"  3D distance:      {np.sqrt((float(_ue_pos['x'])-float(_ru_pos['x']))**2 + (float(_ue_pos['y'])-float(_ru_pos['y']))**2 + (float(_ue_pos['z'])-float(_ru_pos['z']))**2):.2f} m")
# %%
