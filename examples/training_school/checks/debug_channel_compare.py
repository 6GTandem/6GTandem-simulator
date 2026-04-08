import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml


def resolve_repo_root() -> Path:
    """Resolve repository root so this script can run from any working directory."""
    root = Path(__file__).resolve().parent
    while not (root / "wireless_channel").exists() and root != root.parent:
        root = root.parent

    if not (root / "wireless_channel").exists():
        raise RuntimeError("Could not locate repository root containing 'wireless_channel'.")

    return root


ROOT = resolve_repo_root()
os.chdir(ROOT)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from wireless_channel.subTHz_channel import Channel, build_channel
from wireless_channel.waveforms import Waveform


# User settings: edit these values directly.
ENVIRONMENT = "office_space_perpendicular"
ENV_CONFIG_PATH = None
WAVEFORM_CONFIG_PATH = "examples/training_school/configs/tutorial_waveform_config.yaml"
COMPONENT_CONFIG_PATH = "examples/training_school/configs/tutorial_component_config.yaml"

UE_INDEX = 940
USE_MANUAL_UE_POSITION = False
UE_POSITION = None

STRIPE_IDX = 0
RU_IDX = 10
RX_ANT = 0
TX_ANT = 0

NUM_STRIPES = None
NUM_RUS_PER_STRIPE = None

LOS_NORMALIZE_GAIN = False
PLOT = True
PLOT_OUTPUT_PATH = "examples/training_school/channel_compare_plot.png"


def truncate_radio_stripes(radio_stripes, n_stripes, n_rus):
    """Keep only first n_stripes and first n_rus RUs per stripe."""
    if n_stripes is None:
        n_stripes = len(radio_stripes)

    trimmed = []
    for stripe in radio_stripes[:n_stripes]:
        central_units = [entry for entry in stripe if "central_unit" in entry]
        rus = [entry for entry in stripe if "radio_unit" in entry]
        if n_rus is not None:
            rus = rus[:n_rus]
        trimmed.append(central_units[:1] + rus)
    return trimmed


def summarize(name: str, values: np.ndarray) -> dict:
    return {
        "name": name,
        "min": float(np.min(values)),
        "p5": float(np.percentile(values, 5)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "p95": float(np.percentile(values, 95)),
        "max": float(np.max(values)),
    }


def print_summary_line(stats: dict, unit: str):
    print(
        f"{stats['name']:<24} "
        f"min={stats['min']:.4e}{unit} "
        f"p5={stats['p5']:.4e}{unit} "
        f"median={stats['median']:.4e}{unit} "
        f"mean={stats['mean']:.4e}{unit} "
        f"p95={stats['p95']:.4e}{unit} "
        f"max={stats['max']:.4e}{unit}"
    )


def ensure_index(name: str, value: int, size: int):
    if value < 0 or value >= size:
        raise ValueError(f"{name}={value} out of range [0, {size - 1}]")


def format_xyz(coords) -> str:
    if isinstance(coords, dict):
        x = coords["x"]
        y = coords["y"]
        z = coords["z"]
    else:
        x, y, z = coords
    return f"(x={x}, y={y}, z={z})"


def main():
    env_config_path = Path(ENV_CONFIG_PATH) if ENV_CONFIG_PATH else ROOT / "environments" / ENVIRONMENT / "config.yaml"
    waveform_config_path = Path(WAVEFORM_CONFIG_PATH)
    component_config_path = Path(COMPONENT_CONFIG_PATH)

    with open(env_config_path, "r", encoding="utf8") as f:
        env_config = yaml.safe_load(f)

    with open(waveform_config_path, "r", encoding="utf8") as f:
        waveform_config = yaml.safe_load(f)

    with open(component_config_path, "r", encoding="utf8") as f:
        component_config = yaml.safe_load(f)

    configured_pattern = str(component_config.get("antenna", {}).get("pattern", "")).strip().lower()
    los_pattern = configured_pattern if configured_pattern else "isotropic"

    stripes_cfg = truncate_radio_stripes(
        env_config["radio_stripes"],
        n_stripes=NUM_STRIPES,
        n_rus=NUM_RUS_PER_STRIPE,
    )

    ensure_index("stripe_idx", STRIPE_IDX, len(stripes_cfg))
    rus_in_stripe = [entry for entry in stripes_cfg[STRIPE_IDX] if "radio_unit" in entry]
    ensure_index("ru_idx", RU_IDX, len(rus_in_stripe))
    ru_pos = rus_in_stripe[RU_IDX]["radio_unit"]

    if USE_MANUAL_UE_POSITION:
        if UE_POSITION is None:
            raise ValueError("UE_POSITION must be set when USE_MANUAL_UE_POSITION=True.")
        ue_pos = UE_POSITION
    else:
        ensure_index("ue_index", UE_INDEX, len(env_config["ue_positions"]))
        ue_pos = env_config["ue_positions"][UE_INDEX]

    wf = Waveform.from_config(waveform_config, env_config["sub_thz"])
    n_antennas = int(env_config["antenna"]["N_antennas"])

    print(f"Repository root      : {ROOT}")
    print(f"Environment          : {ENVIRONMENT}")
    print(f"Environment config   : {env_config_path}")
    print(f"Waveform config      : {waveform_config_path}")
    print(f"Component config     : {component_config_path}")
    print(f"LOS antenna pattern  : {los_pattern}")
    print(f"Stripe/RU            : stripe={STRIPE_IDX}, ru={RU_IDX}")
    print(f"Selected RU coords   : {format_xyz(ru_pos)}")
    print(f"Selected UE coords   : {format_xyz(ue_pos)}")
    print(f"Antennas (rx, tx)    : rx={RX_ANT}, tx={TX_ANT}")
    print(f"LOS normalize gain   : {LOS_NORMALIZE_GAIN}")

    ru_xyz = np.array([ru_pos["x"], ru_pos["y"], ru_pos["z"]], dtype=float)
    ue_xyz = np.array([ue_pos["x"], ue_pos["y"], ue_pos["z"]], dtype=float)
    pattern_diag = Channel._compute_link_pattern_diagnostics(
        pattern=los_pattern,
        ue_xyz=ue_xyz,
        ru_xyz=ru_xyz,
        component_config=component_config,
    )
    print(
        "LOS boresights       : "
        f"UE(az={pattern_diag['ue_boresight_az_deg']}, el={pattern_diag['ue_boresight_el_deg']}), "
        f"RU(az={pattern_diag['ru_boresight_az_deg']}, el={pattern_diag['ru_boresight_el_deg']})"
    )
    print(
        "LOS local angles     : "
        f"tx(theta={pattern_diag['tx_theta_deg']:.2f}, phi={pattern_diag['tx_phi_deg']:.2f}), "
        f"rx(theta={pattern_diag['rx_theta_deg']:.2f}, phi={pattern_diag['rx_phi_deg']:.2f})"
    )
    print(
        "LOS pattern gains    : "
        f"tx={pattern_diag['tx_power_gain_db']:.2f} dBi, "
        f"rx={pattern_diag['rx_power_gain_db']:.2f} dBi, "
        f"total={pattern_diag['total_power_gain_db']:.2f} dB, "
        f"field={pattern_diag['field_gain']:.3f}"
    )

    channel_sionna = build_channel(
        channel_model="sionna",
        ue_coordinates=ue_pos,
        sim_env=ENVIRONMENT,
        stripe_positions=stripes_cfg,
        waveform=wf,
        Nr_ue_antennas=n_antennas,
        Nr_ru_antennas=n_antennas,
        debug=False,
        los_normalize_gain=False,
    )

    channel_los = build_channel(
        channel_model="los",
        ue_coordinates=ue_pos,
        sim_env=ENVIRONMENT,
        component_config=component_config,
        stripe_positions=stripes_cfg,
        waveform=wf,
        Nr_ue_antennas=n_antennas,
        Nr_ru_antennas=n_antennas,
        debug=False,
        los_normalize_gain=LOS_NORMALIZE_GAIN,
    )

    h_sionna = channel_sionna.get_csi(STRIPE_IDX, RU_IDX)
    h_los = channel_los.get_csi(STRIPE_IDX, RU_IDX)

    ensure_index("rx_ant", RX_ANT, h_sionna.shape[0])
    ensure_index("tx_ant", TX_ANT, h_sionna.shape[1])

    h_sionna_sel = h_sionna[RX_ANT, TX_ANT, :]
    h_los_sel = h_los[RX_ANT, TX_ANT, :]

    eps = 1e-18
    amp_sionna = np.abs(h_sionna_sel)
    amp_los = np.abs(h_los_sel)
    power_sionna = amp_sionna ** 2
    power_los = amp_los ** 2

    delta_amp_db = 20.0 * np.log10((amp_los + eps) / (amp_sionna + eps))
    delta_power_db = 10.0 * np.log10((power_los + eps) / (power_sionna + eps))

    link_power_sionna = np.mean(np.abs(h_sionna) ** 2)
    link_power_los = np.mean(np.abs(h_los) ** 2)
    link_power_delta_db = 10.0 * np.log10((link_power_los + eps) / (link_power_sionna + eps))

    print("\nPer-subcarrier stats for selected antenna pair")
    print_summary_line(summarize("|H| Sionna", amp_sionna), "")
    print_summary_line(summarize("|H| LOS", amp_los), "")
    print_summary_line(summarize("|H| delta (LOS/Sionna) dB", delta_amp_db), " dB")
    print_summary_line(summarize("|H|^2 Sionna", power_sionna), "")
    print_summary_line(summarize("|H|^2 LOS", power_los), "")
    print_summary_line(summarize("|H|^2 delta dB", delta_power_db), " dB")

    print("\nLink-average power across all rx/tx antennas and subcarriers")
    print(f"Sionna mean |H|^2     : {link_power_sionna:.6e}")
    print(f"LOS mean |H|^2        : {link_power_los:.6e}")
    print(f"Delta LOS-Sionna (dB) : {link_power_delta_db:.3f} dB")

    top_idx = np.argsort(np.abs(delta_power_db))[::-1][:10]
    print("\nTop 10 subcarriers by absolute power discrepancy")
    print("idx | |H|^2_sionna | |H|^2_los | delta_dB")
    for idx in top_idx:
        print(
            f"{int(idx):4d} | {power_sionna[idx]:.6e} | {power_los[idx]:.6e} | {delta_power_db[idx]:+.3f}"
        )

    if PLOT:
        subcarrier_idx = np.arange(h_sionna_sel.shape[0])
        fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

        axes[0].plot(subcarrier_idx, 20.0 * np.log10(amp_sionna + eps), label="Sionna")
        axes[0].plot(subcarrier_idx, 20.0 * np.log10(amp_los + eps), label="LOS", alpha=0.8)
        axes[0].set_ylabel("|H| (dB)")
        axes[0].set_title(
            f"Channel comparison: stripe={STRIPE_IDX}, ru={RU_IDX}, "
            f"rx={RX_ANT}, tx={TX_ANT}"
        )
        axes[0].grid(True, alpha=0.3)
        axes[0].legend()

        axes[1].plot(subcarrier_idx, 10.0 * np.log10(power_sionna + eps), label="Sionna")
        axes[1].plot(subcarrier_idx, 10.0 * np.log10(power_los + eps), label="LOS", alpha=0.8)
        axes[1].set_ylabel("|H|^2 (dB)")
        axes[1].grid(True, alpha=0.3)
        axes[1].legend()

        axes[2].plot(subcarrier_idx, delta_power_db, color="tab:red")
        axes[2].set_xlabel("Subcarrier index")
        axes[2].set_ylabel("LOS-Sionna (dB)")
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        output_path = ROOT / PLOT_OUTPUT_PATH
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved plot to: {output_path}")


if __name__ == "__main__":
    main()