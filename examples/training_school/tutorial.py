# %% [markdown]
# # RadioStripe Uplink Tutorial
#
# This tutorial walks through a simple uplink simulation using the 6GTandem simulator.
#
# You will:
# - configure a small RadioStripe deployment,
# - run an idealized baseline,
# - then enable realistic hardware impairments,
# - and compare BER, EVM, and constellation diagrams.
#
# The tutorial is intentionally compact so you can quickly understand how the
# simulator objects connect together.
#
# > **How to run:** open this file in VS Code and click *Run Cell* above each
# > `# %%` marker, or run the whole file with `python tutorial.py`.

# %%
import os
import sys
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

# Resolve repository root so this script can run from any working directory.
root = Path(__file__).resolve().parent
while not (root / "wireless_channel").exists() and root != root.parent:
    root = root.parent

if not (root / "wireless_channel").exists():
    raise RuntimeError("Could not locate repository root containing 'wireless_channel'.")

os.chdir(root)
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

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

print(f"Repository root: {root}")

# Plotting mode: set to True to save plots to files, False to display interactively
SAVE_PLOTS = True
if SAVE_PLOTS:
    plot_dir = Path(__file__).parent / "tutorial_plots"
    plot_dir.mkdir(exist_ok=True)
    print(f"Plots will be saved to: {plot_dir.resolve()}")
    plot_count = 0

def save_or_show_plot(title=None):
    """Save plot to file if SAVE_PLOTS=True, otherwise display interactively."""
    global plot_count
    if SAVE_PLOTS:
        plot_count += 1
        filename = plot_dir / f"{plot_count:02d}_{title.replace(' ', '_').lower()}.png" if title else plot_dir / f"{plot_count:02d}.png"
        plt.savefig(filename, dpi=100, bbox_inches="tight")
        print(f"  → Saved to {filename}")
        plt.close()
    else:
        plt.show()

# %% [markdown]
# ## 1) Configure A Small Scenario
#
# Three key architecture parameters are configurable here:
#
# | Variable | Meaning |
# |---|---|
# | `num_stripes` | How many RadioStripes to deploy |
# | `num_rus_per_stripe` | How many Radio Units to keep per stripe |
# | `active_ru_index` | Which RU within each stripe is the **active entry point** for uplink reception |
#
# Setting `stripe.active_unit = active_ru_index` before calling `stripe.receive()`
# tells the stripe to route the received signal starting from that specific RU
# along the chain towards the Central Unit.
#
# Local config files in `examples/training_school/configs/` keep this tutorial
# self-contained and separate from the shared environment configs.

# %%
# ---------------------------------------------------------------------------
# User-configurable knobs — change these to explore the architecture
# ---------------------------------------------------------------------------
num_stripes       = 3   # number of RadioStripes to use
num_rus_per_stripe = 5  # number of RUs per stripe (must be <= RUs in config)
ue_position_index = 0   # which UE position from the config to use
active_ru_index   = 0   # which RU within each stripe is the active reception entry point
plot_room         = True
channel_model     = "sionna"   # "sionna" or "los"
los_normalize_gain = False

# ---------------------------------------------------------------------------
# Load tutorial-local config files
# ---------------------------------------------------------------------------
environment = "office_space_perpendicular"
config_dir = root / "examples" / "training_school" / "configs"

with open(config_dir / "tutorial_environment.yaml", "r", encoding="utf8") as f:
    env_config = yaml.safe_load(f)

with open(config_dir / "tutorial_waveform_config.yaml", "r", encoding="utf8") as f:
    waveform_config = yaml.safe_load(f)

with open(config_dir / "tutorial_component_config.yaml", "r", encoding="utf8") as f:
    component_config_raw = yaml.safe_load(f)


def truncate_radio_stripes(radio_stripes, n_stripes, n_rus):
    """Keep only the first n_stripes stripes, each with at most n_rus RUs."""
    trimmed = []
    for stripe in radio_stripes[:n_stripes]:
        central_units = [entry for entry in stripe if "central_unit" in entry]
        rus = [entry for entry in stripe if "radio_unit" in entry][:n_rus]
        trimmed.append(central_units[:1] + rus)
    return trimmed


def make_ideal_component_config(component_cfg):
    """Return a copy of component_cfg with all hardware impairments turned off."""
    cfg = deepcopy(component_cfg)
    for amp_key in ("boost_amplifier", "antenna_amplifier"):
        cfg.setdefault(amp_key, {})
        cfg[amp_key]["mode"] = "ideal"
    return cfg


selected_stripes_cfg = truncate_radio_stripes(
    env_config["radio_stripes"],
    n_stripes=num_stripes,
    n_rus=num_rus_per_stripe,
)

component_config_ideal    = make_ideal_component_config(component_config_raw)
component_config_impaired = deepcopy(component_config_raw)

print(f"Environment      : {environment}")
print(f"Stripes          : {len(selected_stripes_cfg)}")
print(f"RUs per stripe   : {len(selected_stripes_cfg[0]) - 1}")
print(f"Active RU index  : {active_ru_index}")

# %% [markdown]
# ## 2) Build Waveform and RadioStripe Objects
#
# We instantiate the waveform and build **two** stripe sets:
# - **Ideal** — unity-gain amplifiers, no fiber/coupler distortion.
# - **Impaired** — full hardware models from `tutorial_component_config.yaml`.
#
# Both sets use exactly the same geometry so the later comparison is fair.

# %%
freq_band_config = env_config["sub_thz"]
wf = Waveform.from_config(waveform_config, freq_band_config)
n_antennas = env_config["antenna"]["N_antennas"]
ue_pos = env_config["ue_positions"][ue_position_index]


def build_stripes(stripe_cfgs, component_cfg, waveform, n_ant):
    return [
        RadioStripe.from_config_locations(cfg, component_cfg, antennas=n_ant, wf=waveform)
        for cfg in stripe_cfgs
    ]


stripes_ideal    = build_stripes(selected_stripes_cfg, component_config_ideal,    wf, n_antennas)
stripes_impaired = build_stripes(selected_stripes_cfg, component_config_impaired, wf, n_antennas)

if plot_room:
    plotter.plot_room(env_config)

print(f"Waveform : {wf.waveform_type}, carriers={wf.n_carriers}, OFDM symbols={wf.n_ofdm_symbols}")
for idx, stripe in enumerate(stripes_ideal):
    print(f"  Stripe {idx}: {len(stripe.radio_units)} RUs")
print(f"UE position [{ue_position_index}]: {ue_pos}")

# %% [markdown]
# ## 3) Step-By-Step Uplink Walkthrough
#
# This section goes through the full uplink chain one stage at a time and
# visualises the intermediate signals so you can see exactly what is happening
# at each step.
#
# | Step | What happens |
# |---|---|
# | 1 | Random bits are generated |
# | 2 | Bits → QAM symbols |
# | 3 | QAM → OFDM time-domain waveform |
# | 4 | UE front-end → wireless channel → RadioStripe reception |
# | 5 | Equalization → BER / EVM |
#
# We use **ideal hardware** here so the signal processing chain is easy to follow.
# The `run_uplink_scenario` helper defined at the end of this section is reused
# for the scenario comparisons in Sections 4 and 5.

# %% [markdown]
# ### Helper functions
#
# `build_ue_radio_unit` assembles the UE-side transmit hardware from the component
# config (amplifiers, coupler, splitter, combiner, and phase shifter).
#
# `equalize_and_measure` takes the received time-domain signal, runs OFDM
# demodulation, LS channel estimation, one-tap equalisation, and returns BER,
# EVM, and the raw/equalised QAM symbol arrays for plotting.

# %%
def build_ue_radio_unit(waveform, component_cfg, n_ant):
    coupler_cfg = component_cfg.get("coupler", {})
    if coupler_cfg:
        coup_in  = Coupler.from_config(coupler_cfg, waveform)
        coup_out = Coupler.from_config(coupler_cfg, waveform)
    else:
        coup_in  = Coupler(wf=waveform)
        coup_out = Coupler(wf=waveform)

    bw = waveform.bw * waveform.oversampling_factor
    boost_amp = Amplifier(**{**component_cfg.get("boost_amplifier", {}), "bw": bw})
    ant_amp   = Amplifier(**{**component_cfg.get("antenna_amplifier", {}), "bw": bw})

    splitter = Splitter(num_splits=n_ant)
    combiner = Combiner()

    phase_cfg    = component_cfg.get("phase_shifter", {})
    phase_shifter = PhaseShifter(num_shifters=n_ant, resolution=phase_cfg.get("resolution", 4))

    return RadioUnit(
        x=0.0, y=0.0, z=0.0,
        boost_amp=boost_amp, antenna_amp=ant_amp,
        coup_in=coup_in, coup_out=coup_out,
        splitter=splitter, combiner=combiner,
        pshift=phase_shifter,
    )


def equalize_and_measure(waveform, y_time, qam_tx, bits_tx):
    y_freq  = waveform.ofdm_time_to_freq(y_time)
    subc    = waveform.extract_subcarriers(y_freq)
    h_est   = waveform.channel_estimate_ls(subc)
    eq_subc = waveform.equalize_one_tap(subc, h_est)

    qam_raw = waveform.demap_data_from_grid(subc).flatten()
    qam_eq  = waveform.demap_data_from_grid(eq_subc).flatten()

    n = min(len(qam_tx), len(qam_eq))
    qam_tx, qam_raw, qam_eq = qam_tx[:n], qam_raw[:n], qam_eq[:n]

    # Normalise pre-EQ symbols for plotting so the scatter shape is visible.
    scale       = (np.sqrt(np.mean(np.abs(qam_tx) ** 2)) + 1e-12) / (np.sqrt(np.mean(np.abs(qam_raw) ** 2)) + 1e-12)
    qam_raw_plot = qam_raw * scale

    bits_rx = waveform.qam_to_bits(qam_eq)
    bits_tx = bits_tx[: len(bits_rx)]
    ber     = waveform.compute_ber(bits_tx, bits_rx)
    evm_rms = np.sqrt(np.mean(np.abs(qam_eq - qam_tx) ** 2) / np.mean(np.abs(qam_tx) ** 2))

    return {
        "ber": ber,
        "evm_percent": 100 * evm_rms,
        "qam_tx": qam_tx,
        "qam_raw": qam_raw,
        "qam_raw_plot": qam_raw_plot,
        "qam_eq": qam_eq,
    }

# %% [markdown]
# ### Step 1 — Generate transmit bits
#
# The waveform object generates a fresh block of random bits every time
# `generate_bits()` is called.  We plot a short window so the bit pattern is
# easy to inspect visually.

# %%
print("Step 1/5: Generate random TX bits")
bits_demo   = wf.generate_bits()
n_bits_plot = 200

fig, ax = plt.subplots(figsize=(12, 2.5))
ax.step(np.arange(n_bits_plot), bits_demo[:n_bits_plot], where="mid")
ax.set_title(f"First {n_bits_plot} generated TX bits")
ax.set_xlabel("Bit index")
ax.set_ylabel("Bit value")
ax.set_ylim(-0.2, 1.2)
ax.grid(True, alpha=0.3)
plt.tight_layout()
save_or_show_plot("tx_bits")

# %% [markdown]
# ### Step 2 — Map bits to QAM symbols
#
# Groups of log₂(M) bits are mapped to complex QAM constellation points.
# With QPSK (M=4) there are 4 distinct points; with 16-QAM there are 16, etc.
# The scatter plot should show the ideal grid — no noise yet.

# %%
print("Step 2/5: QAM mapping")
qam_demo    = wf.qam_modulate()
n_qam_plot  = 200

fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(np.real(qam_demo[:n_qam_plot]), np.imag(qam_demo[:n_qam_plot]),
           s=15, alpha=0.7, color="tab:blue")
ax.set_title(f"TX QAM constellation (first {n_qam_plot} symbols, {wf.qam_order}-QAM)")
ax.set_xlabel("In-phase (I)")
ax.set_ylabel("Quadrature (Q)")
ax.axhline(0, color="k", lw=0.8, alpha=0.5)
ax.axvline(0, color="k", lw=0.8, alpha=0.5)
ax.grid(True, alpha=0.25)
ax.set_aspect("equal", adjustable="box")
plt.tight_layout()
save_or_show_plot("tx_qam_constellation")

# %% [markdown]
# ### Step 3 — OFDM modulation
#
# The QAM symbols are mapped onto OFDM subcarriers (IFFT + cyclic prefix).
# The resulting time-domain signal looks noise-like — its power spectral density
# (PSD) should be flat within the occupied bandwidth and roll off sharply outside.

# %%
print("Step 3/5: OFDM modulation")
ofdm_time_demo = wf.ofdm_modulate()

# PSD
wf.plot_psd(ofdm_time_demo, title="TX OFDM power spectral density")

# Time-domain waveform
n_time_plot = 500
fig, ax = plt.subplots(figsize=(12, 3))
ax.plot(np.real(ofdm_time_demo[:n_time_plot]), label="Real (I)", lw=0.9)
ax.plot(np.imag(ofdm_time_demo[:n_time_plot]), label="Imag (Q)", lw=0.9, alpha=0.8)
ax.set_title(f"OFDM time-domain waveform — first {n_time_plot} samples")
ax.set_xlabel("Sample index")
ax.set_ylabel("Amplitude")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
save_or_show_plot("tx_ofdm_waveform")

# %% [markdown]
# ### Step 4 — UE front-end, wireless channel, and RadioStripe reception
#
# The signal now passes through:
# 1. **CentralUnit** (UE side) — models the UE baseband processing.
# 2. **UE RadioUnit** — applies amplifier gain and phase shifts before transmission.
# 3. **Wireless channel** — propagates the signal to each RU in the room
#    (`Channel.transmit_ul` returns one array per stripe × RU).
# 4. **RadioStripe reception** — for each stripe we set the active entry-point RU
#    via `stripe.active_unit` and call `stripe.receive()`.  The signal is then
#    routed from that RU along the daisy-chained fibre/coupler links towards the
#    Central Unit.

# %%
print("Step 4/5: UE front-end → channel → stripe reception")

cu_demo            = CentralUnit()
ofdm_after_cu_demo = cu_demo.run(ofdm_time_demo)

ue_ru_demo       = build_ue_radio_unit(wf, component_config_ideal, n_antennas)
phase_shifts_demo = np.zeros(n_antennas)   # broadside — no beam steering
iq_tx_demo, _    = ue_ru_demo.transmit(ofdm_after_cu_demo, phase_shifts_demo)

# Propagate through the wireless channel.
# debug=True uses an all-ones channel response (no fading) for clarity in this demo.
channel_demo           = build_channel(
    channel_model=channel_model,
    ue_coordinates=ue_pos,
    sim_env=environment,
    component_config=component_config_ideal,
    stripe_positions=selected_stripes_cfg,
    waveform=wf,
    Nr_ue_antennas=n_antennas,
    Nr_ru_antennas=n_antennas,
    debug=True,
    los_normalize_gain=los_normalize_gain,
)
channel_demo.Nr_stripes = len(stripes_ideal)
channel_demo.Nr_rus     = len(stripes_ideal[0].radio_units)
iq_rus_demo            = channel_demo.transmit_ul(iq_tx_demo, wf)

# Receive on each stripe using the selected active RU.
stripe_outputs_demo = []
for stripe_idx, stripe in enumerate(stripes_ideal):
    # Tell the stripe which RU the signal arrived at.
    stripe.active_unit = active_ru_index
    # Extract the data for the active RU: iq_rus_demo[stripe][ru][antennas x samples]
    # Route the signal from that RU to the Central Unit.
    # arg 3 = False  →  no windowing
    y_stripe_demo, _ = stripe.receive(iq_rus_demo[stripe_idx][active_ru_index], phase_shifts_demo, False)
    stripe_outputs_demo.append(y_stripe_demo)

# Combine outputs from all stripes (maximal-ratio combining approximation).
y_combined_demo = np.sum(np.array(stripe_outputs_demo), axis=0)
print(f"  Active RU index        : {active_ru_index}")
print(f"  Combined output length : {len(y_combined_demo)} samples")
save_or_show_plot("rx_psd_demo")

# %% [markdown]
# ### Step 5 — Equalisation and quality metrics
#
# The received time-domain block is OFDM-demodulated, LS channel estimates are
# computed from the pilots, and one-tap frequency-domain equalisation removes
# the channel response.  We then compute BER and EVM.
#
# The **before-equalisation** constellation shows the raw subcarrier symbols
# (normalised to TX power for a fair visual comparison).  The
# **after-equalisation** constellation should collapse back to the original
# QAM grid if the channel estimation worked well.

# %%
print("Step 5/5: Equalisation and quality metrics")
demo_meas = equalize_and_measure(wf, y_combined_demo, qam_demo, bits_demo)
print(f"  Demo BER : {demo_meas['ber']:.3e}")
print(f"  Demo EVM : {demo_meas['evm_percent']:.3f} %")

wf.plot_constellation(
    demo_meas["qam_raw_plot"][:2000],
    symbols_tx=demo_meas["qam_tx"][:2000],
    title="Before equalisation (normalised for display)",
)
save_or_show_plot("demo_before_eq")
wf.plot_constellation(
    demo_meas["qam_eq"][:2000],
    symbols_tx=demo_meas["qam_tx"][:2000],
    title="After equalisation",
)
save_or_show_plot("demo_after_eq")

# %% [markdown]
# ### Reusable scenario helper
#
# The function below bundles the full uplink chain so Sections 4 and 5 can run
# the same processing with a single call and different hardware configs.

# %%
def run_uplink_scenario(stripes, waveform, ue_position, sim_environment,
                        component_cfg, stripe_positions_cfg,
                        channel_model="sionna", los_normalize_gain=False,
                        debug_channel=False):
    bits     = waveform.generate_bits()
    qam      = waveform.qam_modulate()
    ofdm_time = waveform.ofdm_modulate()

    cu               = CentralUnit()
    ofdm_time_after_cu = cu.run(ofdm_time)

    ue_ru        = build_ue_radio_unit(waveform, component_cfg, n_antennas)
    phase_shifts = np.zeros(n_antennas)
    iq_data_tx, _ = ue_ru.transmit(ofdm_time_after_cu, phase_shifts)

    channel            = build_channel(
        channel_model=channel_model,
        ue_coordinates=ue_position,
        sim_env=sim_environment,
        component_config=component_cfg,
        stripe_positions=stripe_positions_cfg,
        waveform=waveform,
        Nr_ue_antennas=n_antennas,
        Nr_ru_antennas=n_antennas,
        debug=debug_channel,
        los_normalize_gain=los_normalize_gain,
    )
    channel.Nr_stripes = len(stripes)
    channel.Nr_rus     = len(stripes[0].radio_units)
    iq_data_rus        = channel.transmit_ul(iq_data_tx, waveform)

    # Set the active entry-point RU on each stripe, then receive.
    stripe_outputs = []
    for stripe_idx, stripe in enumerate(stripes):
        stripe.active_unit = active_ru_index
        # Extract data for the active RU from this stripe.
        y_stripe, _ = stripe.receive(iq_data_rus[stripe_idx][active_ru_index], phase_shifts, False)
        stripe_outputs.append(y_stripe)

    y_combined   = np.sum(np.array(stripe_outputs), axis=0)
    measurements = equalize_and_measure(waveform, y_combined, qam, bits)

    return {
        "tx_ofdm_time": ofdm_time,
        "tx_iq_data":   iq_data_tx,
        "y_combined":   y_combined,
        "measurements": measurements,
    }

# %% [markdown]
# ## 4) Scenario A — Idealized Baseline
#
# Hardware effects are switched off:
# - amplifiers are unity-gain and noiseless,
# - fiber and coupler responses are bypassed,
# - channel is set to debug mode (flat, all-ones CSI).
#
# This gives a clean lower-bound on BER/EVM so you have a reference point.

# %%
ideal_results = run_uplink_scenario(
    stripes=stripes_ideal,
    waveform=wf,
    ue_position=ue_pos,
    sim_environment=environment,
    component_cfg=component_config_ideal,
    stripe_positions_cfg=selected_stripes_cfg,
    channel_model=channel_model,
    los_normalize_gain=los_normalize_gain,
    debug_channel=True,
)

ideal_meas = ideal_results["measurements"]
print(f"Ideal BER : {ideal_meas['ber']:.3e}")
print(f"Ideal EVM : {ideal_meas['evm_percent']:.3f} %")

wf.plot_constellation(
    ideal_meas["qam_raw"][:2000],
    symbols_tx=ideal_meas["qam_tx"][:2000],
    title="Ideal: before equalisation",
)
save_or_show_plot("ideal_before_eq")
wf.plot_constellation(
    ideal_meas["qam_eq"][:2000],
    symbols_tx=ideal_meas["qam_tx"][:2000],
    title="Ideal: after equalisation",
)
save_or_show_plot("ideal_after_eq")
wf.plot_psd(ideal_results["tx_ofdm_time"], title="Ideal: UE transmit PSD")
save_or_show_plot("ideal_tx_psd")

# %% [markdown]
# ## 5) Scenario B — Realistic Hardware Impairments
#
# Now we enable the full hardware models from `tutorial_component_config.yaml`:
# - **Nonlinear PA** (5th-order polynomial model) — compresses large signal peaks.
# - **Amplifier noise figure** — adds thermal noise at each stage.
# - **Fiber dispersion** — inter-RU fibre links introduce frequency-dependent delay.
# - **Coupler insertion loss** — each tap on the daisy chain attenuates the signal.
# - **Real Sionna ray-traced channel** — realistic multipath fading.
#
# Compare the constellations and metrics against the ideal case above to see the
# combined impact of these effects.

# %%
impaired_results = run_uplink_scenario(
    stripes=stripes_impaired,
    waveform=wf,
    ue_position=ue_pos,
    sim_environment=environment,
    component_cfg=component_config_impaired,
    stripe_positions_cfg=selected_stripes_cfg,
    channel_model=channel_model,
    los_normalize_gain=los_normalize_gain,
    debug_channel=False,
)

impaired_meas = impaired_results["measurements"]
print(f"Impaired BER : {impaired_meas['ber']:.3e}")
print(f"Impaired EVM : {impaired_meas['evm_percent']:.3f} %")

wf.plot_constellation(
    impaired_meas["qam_raw"][:2000],
    symbols_tx=impaired_meas["qam_tx"][:2000],
    title="Impaired: before equalisation",
)
save_or_show_plot("impaired_before_eq")
wf.plot_constellation(
    impaired_meas["qam_eq"][:2000],
    symbols_tx=impaired_meas["qam_tx"][:2000],
    title="Impaired: after equalisation",
)
save_or_show_plot("impaired_after_eq")

# Side-by-side bar chart comparison
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
labels     = ["Ideal", "Impaired"]
ber_values = [ideal_meas["ber"],         impaired_meas["ber"]]
evm_values = [ideal_meas["evm_percent"], impaired_meas["evm_percent"]]

axes[0].bar(labels, ber_values, color=["tab:blue", "tab:orange"])
axes[0].set_yscale("log")
axes[0].set_title("BER comparison")
axes[0].set_ylabel("BER (log scale)")

axes[1].bar(labels, evm_values, color=["tab:blue", "tab:orange"])
axes[1].set_title("EVM comparison")
axes[1].set_ylabel("EVM (%)")

fig.tight_layout()
save_or_show_plot("metrics_comparison")

# %% [markdown]
# ## 6) What To Try Next
#
# You now have a compact, architecture-focused simulation flow. Some exercises:
#
# - **Scale the deployment** — increase `num_stripes` or `num_rus_per_stripe`
#   and watch BER/EVM improve as spatial diversity grows.
# - **Move the UE** — sweep `ue_position_index` to see how different channel
#   realisations affect performance.
# - **Move the active RU** — change `active_ru_index` from 0 to
#   `num_rus_per_stripe - 1`.  The signal travels further along the daisy chain
#   before reaching the CU, so you can observe the effect of additional
#   fiber/coupler hops.
# - **Add beam steering** — replace the all-zero phase shifts with
#   `stripe.radio_units[active_ru_index].phase_shifter.get_phases(beam_angle_deg)`
#   and pass the result to `stripe.receive()`.
# - **Compare hardware configs** — edit `tutorial_component_config.yaml` to
#   change noise figure, PA back-off, or fibre length and re-run Section 5.
