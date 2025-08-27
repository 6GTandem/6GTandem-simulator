"""
% Matlab script developed using version R2023a for synthesis of array
% pattern using the radiaton pattern of a square patch antenna element in
% Phi and Theta polarized components, obtained from ANSYS HFSS version
% 2022 and saved to realized_gain.csv.
%
% The input impedance of the patch elements are matched to 50ohm.
%
% The configuration of the array is rectangular (along x-y plane).
% Main input parameters are:
% Nx = # array elements in the x dimension
% Ny = # array elements in the y dimension.
% dx = Center-to-center element spacing along the x dimensions
% dy = Center-to-center element spacing along the y dimensions
% PhiScan_deg = Scan angle in the Phi dimension (in degree)
% ThetaScan_deg = Scan angle in the Theta dimension (in degree)
%
% Output parameters are:
% Array pattern saved in the same format as the called element pattern
% with filename "arraypattern.mat"
% Plots of element/array patterns in 3D as well as Phi and Theta cuts.
%
% Antenna Toolbox is needed for the pattern plotting function patternCustom.m
%
% Developed for the HEU 6GTandem Project by
% Yuyan Cao and Buon Kiong Lau, Lund University, Sweden
%
% Current version: 20 Feb 2024
"""
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd


def polar_to_cartesian(gain, theta, phi):
    """Convert gain values in a polar grid to a cartesian grid.

    :param gain: Array containing the gain values.
    :param theta: Theta angle array in radians.
    :param phi: Phi angle array in radians.
    """
    # We have to convert everything to cartesian.
    gain_x = gain * np.sin(theta) * np.cos(phi)
    gain_y = gain * np.sin(theta) * np.sin(phi)
    gain_z = gain * np.cos(theta)

    return gain_x, gain_y, gain_z


m = pd.read_csv('antenna-model/realized_gain.csv')
m_db = m.copy()

# by definition, stored E-field data from Ansys HFSS is abs(Ephi)^2 abs(Etheta)^2
m_db_phi = 10 * np.log10(m_db['RealizedGainPhi'].to_numpy())
m_db_phi = np.reshape(m_db_phi, (len(m_db_phi) // 181, 181))
m_db_theta = 10 * np.log10(m_db['RealizedGainTheta'].to_numpy())
m_db_theta = np.reshape(m_db_theta, (len(m_db_theta) // 181, 181))

phi_deg = m['Phi[deg]'].to_numpy()
phi_deg = np.reshape(phi_deg, (len(phi_deg) // 181, 181))
phi_rad = np.deg2rad(m['Phi[deg]'].to_numpy())
phi_rad = np.reshape(phi_rad, (len(phi_rad) // 181, 181))

theta_deg = m['Theta[deg]'].to_numpy()
theta_deg = np.reshape(theta_deg, (len(theta_deg) // 181, 181))
theta_rad = np.deg2rad(m['Theta[deg]'].to_numpy())
theta_rad = np.reshape(theta_rad, (len(theta_rad) // 181, 181))

# converting the data back to abs(Ephi) and abs(theta) from squared quantities
field_phi = np.sqrt(m['RealizedGainPhi'].to_numpy())
field_phi = np.reshape(field_phi, (len(field_phi) // 181, 181))
field_theta = np.sqrt(m['RealizedGainTheta'].to_numpy())
field_theta = np.reshape(field_theta, (len(field_theta) // 181, 181))

c = 3e8  # speed of light in m/s
f = 145e9  # frequency in Hz
lbda = c/f  # vacuum wavelength in m
beta0 = 2 * np.pi / lbda

# input parameters
nx = 4  # number of antennas in the x dimension(rectangular array)
ny = 4  # number of antennas in the y dimension
dx = 0.55 * lbda  # inter-element spacing in x direction in m
dy = dx  # inter-element spacing in y direction in m
phi_scan_deg = 0  # scan angle \phi in degree
theta_scan_deg = 30  # scan angle \theta in degree

# Array factor in the look directions (φ and θ in degrees)
phi_x_fed_deg = np.rad2deg(
    beta0 * dx * np.sin(np.deg2rad(theta_scan_deg)) * np.cos(np.deg2rad(phi_scan_deg)))
phi_y_fed_deg = np.rad2deg(
    beta0 * dy * np.sin(np.deg2rad(theta_scan_deg)) * np.sin(np.deg2rad(phi_scan_deg)))

# Array Pattern synthesis
field_theta_syn = 0
field_phi_syn = 0
phi_xy = np.zeros(shape=(nx, ny))  # initial phase

for i in range(nx):
    for j in range(ny):
        # increased phase in x axis
        phi_x_array_rad = beta0 * dx * \
            np.sin(np.deg2rad(theta_deg)) * np.cos(np.deg2rad(phi_deg)) * (i-1)
        # increased phase in y axis
        phi_y_array_rad = beta0 * dy * \
            np.sin(np.deg2rad(theta_deg)) * np.sin(np.deg2rad(phi_deg)) * (j-1)
        # total increased phase at each position
        phi_array_rad = phi_x_array_rad + phi_y_array_rad

        # feeding phase at each position compared to reference antenna
        phi_xy[i, j] = phi_x_fed_deg * (i-1) + phi_y_fed_deg * (j-1)

        # element pattern * array factor
        # last exponential term applies conjugate beamforming
        field_theta_syn = field_theta_syn + field_theta * \
            np.exp(1j * phi_array_rad) * np.exp(-1j *
                                                np.deg2rad(phi_xy[i, j])) / np.sqrt(nx * ny)
        # array weight needs to be normalized by 1/sqrt(  # array elements)
        field_phi_syn = field_phi_syn + field_phi * \
            np.exp(1j * phi_array_rad) * np.exp(-1j *
                                                np.deg2rad(phi_xy[i, j])) / np.sqrt(nx * ny)


field_phi_syn_db = 20 * np.log10(abs(field_phi_syn))
field_theta_syn_db = 20 * np.log10(abs(field_theta_syn))

# Element and array pattern plots
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.set_title('Figure 1: Element Realized Gain (Phi)')
surf = ax.plot_surface(*polar_to_cartesian(field_phi, theta_rad, phi_rad),
                       cmap=cm.jet, linewidth=0, antialiased=True)
fig.colorbar(surf)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.set_title('Figure 2: Element Realized Gain (Theta)')
surf = ax.plot_surface(*polar_to_cartesian(field_theta, theta_rad,
                       phi_rad), cmap=cm.jet, linewidth=0, antialiased=True)
fig.colorbar(surf)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.set_title('Figure 3: Array Realized Gain (Phi)')
surf = ax.plot_surface(*polar_to_cartesian(np.abs(field_phi_syn), theta_rad,
                       phi_rad), cmap=cm.jet, linewidth=0, antialiased=True)
fig.colorbar(surf)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.set_title('Figure 3: Array Realized Gain (Theta)')
surf = ax.plot_surface(*polar_to_cartesian(np.abs(field_theta_syn),
                       theta_rad, phi_rad), cmap=cm.jet, linewidth=0, antialiased=True)
fig.colorbar(surf)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

total_gain = (np.abs(field_phi_syn) ** 2) + (np.abs(field_theta_syn) ** 2)

ax.set_title('Figure 7: Array Realized Gain')
surf = ax.plot_surface(*polar_to_cartesian(total_gain, theta_rad,
                       phi_rad), cmap=cm.jet, linewidth=0, antialiased=True)
fig.colorbar(surf)

fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': 'polar'})
fig.suptitle("Realized GainPhi comparison")
fig.tight_layout()
# Take a slice in the phi plane where phi=90deg
gain0 = field_phi[phi_rad == 0]
angle0 = theta_rad[phi_rad == 0]
gain90 = field_phi[phi_rad == (np.pi / 2)]
angle90 = theta_rad[phi_rad == (np.pi / 2)]
ax1.plot(angle0, gain0, label='ϕ=0°')
ax1.plot(angle90, gain90, label='ϕ=90°')
ax1.grid(True)
ax1.legend(loc='lower right')

gain0 = field_phi_syn[phi_rad == 0]
gain90 = field_phi_syn[phi_rad == (np.pi / 2)]
ax2.plot(angle0, gain0, label='ϕ=0°')
ax2.plot(angle90, gain90, label='ϕ=90°')
ax2.grid(True)
ax2.legend(loc='lower right')

fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': 'polar'})
fig.suptitle("Realized GainTheta comparison")
fig.tight_layout()
# Take a slice in the phi plane where phi=90deg
gain0 = field_theta[phi_rad == 0]
angle0 = theta_rad[phi_rad == 0]
gain90 = field_theta[phi_rad == (np.pi / 2)]
angle90 = theta_rad[phi_rad == (np.pi / 2)]
ax1.plot(angle0, gain0, label='ϕ=0°')
ax1.plot(angle90, gain90, label='ϕ=90°')
ax1.grid(True)
ax1.legend(loc='lower right')

gain0 = field_theta_syn[phi_rad == 0]
gain90 = field_theta_syn[phi_rad == (np.pi / 2)]
ax2.plot(angle0, gain0, label='ϕ=0°')
ax2.plot(angle90, gain90, label='ϕ=90°')
ax2.grid(True)
ax2.legend(loc='lower right')

plt.show()
