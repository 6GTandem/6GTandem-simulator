from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Config
ELEMENT_CSV = Path(__file__).with_name("models/antenna/element1.csv")
PHI_ANGLES_DEG = [0, 90]
CENTER_FREQ_HZ = 140e9
C = 299_792_458.0
WAVELENGTH = C / CENTER_FREQ_HZ


def parse_topology(topology: str) -> tuple[int, int]:
    parts = topology.lower().split("x")
    if len(parts) != 2:
        raise ValueError(f"Invalid topology '{topology}'. Use '1x4' or '2x2'.")
    nx, ny = (int(part) for part in parts)
    if nx < 1 or ny < 1:
        raise ValueError("Topology values must be positive integers.")
    return nx, ny


def array_factor_ula(
    theta_rad: np.ndarray,
    phi_rad: float,
    n_elements: int,
) -> np.ndarray:
    dx = WAVELENGTH / 2.0
    k = 2.0 * np.pi / WAVELENGTH
    af = np.zeros_like(theta_rad, dtype=np.complex128)
    sin_theta = np.sin(theta_rad)
    cos_phi = np.cos(phi_rad)
    for m in range(n_elements):
        phase = k * (m * dx) * sin_theta * cos_phi
        af += np.exp(1j * phase)
    return af


def main() -> None:
    df = pd.read_csv(ELEMENT_CSV)
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

    nx, ny = parse_topology(TOPOLOGY)

    dx = dy = WAVELENGTH / 2.0
    k = 2.0 * np.pi / WAVELENGTH

    for phi_angle in PHI_ANGLES_DEG:
        phi_cut = df[df["Phi[deg]"] == phi_angle]

        # Extract columns
        theta_rad = np.deg2rad(phi_cut["Theta[deg]"])
        mag = phi_cut["mag(rERHCP)[mV]"]

        # Convert magnitude to relative gain (dB)
        gain = 20 * np.log10(mag/np.max(mag))

        ax.plot(theta_rad, gain, label=f"Phi {phi_angle}")

        af =  20 * np.log10(np.abs(array_factor_ula(theta_rad, np.deg2rad(phi_angle), 4) ))

        ax.plot(theta_rad, gain + af - np.max(gain + af), label=f"Array {phi_angle}")

    ax.set_theta_zero_location("N")
    ax.set_rmin(-35)
    ax.set_rmax(5)
    ax.set_rlabel_position(135)
    ax.legend(loc="lower left", bbox_to_anchor=(1.05, 0.0))

    out_path = Path(__file__).with_name(f"antenna-polar-theta-af-1x4.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
