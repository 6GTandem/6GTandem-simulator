import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load measurement file
df = pd.read_csv("models/antenna/element1.csv")

fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

# Cut at a specific Phi angle
phi_angles = [0, 90]
for phi_angle in phi_angles:
    phi_cut = df[df["Phi[deg]"] == phi_angle]

    # Extract columns
    theta_rad = np.deg2rad(phi_cut["Theta[deg]"])
    mag = phi_cut["mag(rERHCP)[mV]"]
    deg = phi_cut["ang_deg(rERHCP)[deg]"]

    complex_mag = mag * np.sin(np.deg2rad(deg)) + 1j * mag * np.cos(np.deg2rad(deg))

    # Convert magnitude to relative gain (dB)
    gain_rel = 20 * np.log10(mag / np.max(mag))

    ax.plot(theta_rad, gain_rel, label=f"Phi {phi_angle}")

fig.savefig(f"antenna-polar-theta")
