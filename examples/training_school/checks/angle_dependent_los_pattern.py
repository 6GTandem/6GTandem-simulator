import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml


def resolve_repo_root() -> Path:
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

from wireless_channel.waveforms import Waveform


# User settings: edit these values directly.
ENVIRONMENT = "office_space_perpendicular"
ENV_CONFIG_PATH = None
WAVEFORM_CONFIG_PATH = "examples/training_school/configs/tutorial_waveform_config.yaml"

UE_INDEX = 0
USE_MANUAL_UE_POSITION = False
UE_POSITION = {"x": 2.0, "y": 5.0, "z": 1.5}

NUM_STRIPES = 2
NUM_RUS_PER_STRIPE = 20

# Antenna boresight orientation (global frame): azimuth and elevation in degrees.
# Azimuth: 0=+x, 90=+y. Elevation: 0=horizon, 90=+z.
# Default LOS deployment assumes RU arrays point downward and UE arrays upward.
UE_BORESIGHT_AZ_DEG = 0.0
UE_BORESIGHT_EL_DEG = 90.0
RU_BORESIGHT_AZ_DEG = 0.0
RU_BORESIGHT_EL_DEG = -90.0

MIN_DISTANCE_M = 1e-3
PLOT = True
PLOT_OUTPUT_PATH = "examples/training_school/angle_dependent_los_pattern.png"
CSV_OUTPUT_PATH = "examples/training_school/angle_dependent_los_pattern.csv"


C = 299792458.0
PI = np.pi


def truncate_radio_stripes(radio_stripes, n_stripes, n_rus):
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


def extract_ru_positions(stripe_positions):
    ru_positions = []
    for stripe_idx, stripe in enumerate(stripe_positions):
        ru_idx = 0
        for entry in stripe:
            ru = entry.get("radio_unit")
            if ru is None:
                continue
            ru_positions.append((stripe_idx, ru_idx, np.array([ru["x"], ru["y"], ru["z"]], dtype=float)))
            ru_idx += 1
    return ru_positions


def wrap_to_pi(phi_rad):
    return (phi_rad + PI) % (2.0 * PI) - PI


def tr38901_power_gain(theta_rad, phi_rad):
    """NumPy rewrite of the provided Dr.Jit/Mitsuba pattern (power gain, linear)."""
    phi = wrap_to_pi(phi_rad)

    theta_3db = np.deg2rad(65.0)
    phi_3db = np.deg2rad(65.0)
    a_max = 30.0
    sla_v = 30.0
    g_e_max = 20.0

    a_v = -np.minimum(12.0 * ((theta_rad - PI / 2.0) / theta_3db) ** 2, sla_v)
    a_h = -np.minimum(12.0 * (phi / phi_3db) ** 2, a_max)
    a_db = -np.minimum(-(a_v + a_h), a_max) + g_e_max
    return 10.0 ** (a_db / 10.0)


def unit_vector_from_az_el(az_deg, el_deg):
    az = np.deg2rad(az_deg)
    el = np.deg2rad(el_deg)
    x = np.cos(el) * np.cos(az)
    y = np.cos(el) * np.sin(az)
    z = np.sin(el)
    v = np.array([x, y, z], dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("Boresight vector norm too small.")
    return v / n


def local_frame_from_boresight(boresight_unit):
    """Build right-handed local frame where local +x is boresight."""
    x_hat = boresight_unit
    ref = np.array([0.0, 0.0, 1.0], dtype=float)
    if abs(np.dot(x_hat, ref)) > 0.95:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)

    y_hat = np.cross(ref, x_hat)
    y_hat /= np.linalg.norm(y_hat)
    z_hat = np.cross(x_hat, y_hat)
    z_hat /= np.linalg.norm(z_hat)
    return x_hat, y_hat, z_hat


def vector_to_local_spherical(direction_global, frame):
    """Return (theta, phi) in local spherical coordinates.

    theta is zenith angle from +z in [0, pi], phi is azimuth in [-pi, pi].
    """
    x_hat, y_hat, z_hat = frame
    x_local = np.dot(direction_global, x_hat)
    y_local = np.dot(direction_global, y_hat)
    z_local = np.dot(direction_global, z_hat)

    z_clamped = np.clip(z_local, -1.0, 1.0)
    theta = np.arccos(z_clamped)
    phi = np.arctan2(y_local, x_local)
    return theta, phi


def build_patterned_los_channels(ue_xyz, ru_positions, carrier_frequency_hz, n_subcarriers, subcarrier_spacing_hz):
    ue_boresight = unit_vector_from_az_el(UE_BORESIGHT_AZ_DEG, UE_BORESIGHT_EL_DEG)
    ru_boresight = unit_vector_from_az_el(RU_BORESIGHT_AZ_DEG, RU_BORESIGHT_EL_DEG)

    ue_frame = local_frame_from_boresight(ue_boresight)
    ru_frame = local_frame_from_boresight(ru_boresight)

    sampling_rate_hz = n_subcarriers * subcarrier_spacing_hz
    f_sub = np.fft.fftfreq(n_subcarriers, d=1.0 / sampling_rate_hz)
    wavelength = C / carrier_frequency_hz

    n_links = len(ru_positions)
    h_iso = np.zeros((n_links, n_subcarriers), dtype=complex)
    h_patterned = np.zeros((n_links, n_subcarriers), dtype=complex)

    rows = []
    for link_idx, (stripe_idx, ru_idx, ru_xyz) in enumerate(ru_positions):
        v_ue_to_ru = ru_xyz - ue_xyz
        distance = max(np.linalg.norm(v_ue_to_ru), MIN_DISTANCE_M)
        u_ue_to_ru = v_ue_to_ru / distance
        u_ru_to_ue = -u_ue_to_ru

        tx_theta, tx_phi = vector_to_local_spherical(u_ue_to_ru, ue_frame)
        rx_theta, rx_phi = vector_to_local_spherical(u_ru_to_ue, ru_frame)

        g_tx_lin = float(tr38901_power_gain(tx_theta, tx_phi))
        g_rx_lin = float(tr38901_power_gain(rx_theta, rx_phi))
        g_tot_lin = g_tx_lin * g_rx_lin

        tau = distance / C
        fspl_field = wavelength / (4.0 * PI * distance)
        phase = np.exp(-1j * 2.0 * PI * (carrier_frequency_hz + f_sub) * tau)

        h0 = fspl_field * phase
        hp = fspl_field * np.sqrt(g_tot_lin) * phase
        h_iso[link_idx, :] = h0
        h_patterned[link_idx, :] = hp

        p0 = float(np.mean(np.abs(h0) ** 2))
        pp = float(np.mean(np.abs(hp) ** 2))
        rows.append(
            {
                "stripe_idx": stripe_idx,
                "ru_idx": ru_idx,
                "distance_m": distance,
                "tx_theta_deg": np.rad2deg(tx_theta),
                "tx_phi_deg": np.rad2deg(tx_phi),
                "rx_theta_deg": np.rad2deg(rx_theta),
                "rx_phi_deg": np.rad2deg(rx_phi),
                "g_tx_dbi": 10.0 * np.log10(g_tx_lin + 1e-30),
                "g_rx_dbi": 10.0 * np.log10(g_rx_lin + 1e-30),
                "g_total_db": 10.0 * np.log10(g_tot_lin + 1e-30),
                "mean_power_iso": p0,
                "mean_power_patterned": pp,
                "pattern_vs_iso_db": 10.0 * np.log10((pp + 1e-30) / (p0 + 1e-30)),
            }
        )

    return h_iso, h_patterned, rows


def save_csv(rows, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    headers = list(rows[0].keys())
    lines = [",".join(headers)]
    for row in rows:
        values = [str(row[key]) for key in headers]
        lines.append(",".join(values))
    output_path.write_text("\n".join(lines) + "\n", encoding="utf8")


def plot_link_gains(rows, output_path):
    ru_labels = [f"s{int(r['stripe_idx'])}-ru{int(r['ru_idx'])}" for r in rows]
    gains_db = np.array([r["g_total_db"] for r in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(np.arange(len(gains_db)), gains_db, marker="o", linewidth=1.2)
    ax.set_title("TR38901 angle-dependent total antenna gain per link")
    ax.set_xlabel("Link index")
    ax.set_ylabel("G_tx + G_rx (dB)")
    ax.grid(True, alpha=0.3)

    if len(ru_labels) <= 30:
        ax.set_xticks(np.arange(len(ru_labels)))
        ax.set_xticklabels(ru_labels, rotation=60, ha="right", fontsize=8)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    env_config_path = Path(ENV_CONFIG_PATH) if ENV_CONFIG_PATH else ROOT / "environments" / ENVIRONMENT / "config.yaml"
    waveform_config_path = Path(WAVEFORM_CONFIG_PATH)

    with open(env_config_path, "r", encoding="utf8") as f:
        env_config = yaml.safe_load(f)

    with open(waveform_config_path, "r", encoding="utf8") as f:
        waveform_config = yaml.safe_load(f)

    selected_stripes = truncate_radio_stripes(
        env_config["radio_stripes"],
        n_stripes=NUM_STRIPES,
        n_rus=NUM_RUS_PER_STRIPE,
    )
    ru_positions = extract_ru_positions(selected_stripes)
    if len(ru_positions) == 0:
        raise ValueError("No radio units found after stripe/RU truncation.")

    if USE_MANUAL_UE_POSITION:
        ue_pos = UE_POSITION
    else:
        ue_pos = env_config["ue_positions"][UE_INDEX]

    ue_xyz = np.array([ue_pos["x"], ue_pos["y"], ue_pos["z"]], dtype=float)

    wf = Waveform.from_config(waveform_config, env_config["sub_thz"])
    n_subcarriers = int(wf.n_carriers)
    subcarrier_spacing_hz = float(wf.bw / wf.n_carriers)
    carrier_frequency_hz = float(wf.fc)

    h_iso, h_patterned, rows = build_patterned_los_channels(
        ue_xyz=ue_xyz,
        ru_positions=ru_positions,
        carrier_frequency_hz=carrier_frequency_hz,
        n_subcarriers=n_subcarriers,
        subcarrier_spacing_hz=subcarrier_spacing_hz,
    )

    p_iso = float(np.mean(np.abs(h_iso) ** 2))
    p_pat = float(np.mean(np.abs(h_patterned) ** 2))
    print(f"Environment          : {ENVIRONMENT}")
    print(f"UE position          : {ue_pos}")
    print(f"Num links            : {len(rows)}")
    print(f"UE boresight (az,el) : ({UE_BORESIGHT_AZ_DEG}, {UE_BORESIGHT_EL_DEG}) deg")
    print(f"RU boresight (az,el) : ({RU_BORESIGHT_AZ_DEG}, {RU_BORESIGHT_EL_DEG}) deg")
    print(f"Mean |H|^2 isotropic : {p_iso:.6e}")
    print(f"Mean |H|^2 patterned : {p_pat:.6e}")
    print(f"Pattern gain effect  : {10.0 * np.log10((p_pat + 1e-30) / (p_iso + 1e-30)):.3f} dB")

    rows_sorted = sorted(rows, key=lambda r: r["pattern_vs_iso_db"], reverse=True)
    print("\nTop 10 links by pattern boost (dB)")
    print("stripe ru distance_m tx_theta tx_phi rx_theta rx_phi g_total_db delta_db")
    for r in rows_sorted[:10]:
        print(
            f"{int(r['stripe_idx']):>4d} {int(r['ru_idx']):>2d} "
            f"{r['distance_m']:.3f} "
            f"{r['tx_theta_deg']:.1f} {r['tx_phi_deg']:.1f} "
            f"{r['rx_theta_deg']:.1f} {r['rx_phi_deg']:.1f} "
            f"{r['g_total_db']:.2f} {r['pattern_vs_iso_db']:.2f}"
        )

    csv_output = ROOT / CSV_OUTPUT_PATH
    save_csv(rows, csv_output)
    print(f"\nSaved link table to: {csv_output}")

    if PLOT:
        plot_output = ROOT / PLOT_OUTPUT_PATH
        plot_link_gains(rows, plot_output)
        print(f"Saved gain plot to: {plot_output}")


if __name__ == "__main__":
    main()