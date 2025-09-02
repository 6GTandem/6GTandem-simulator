import skrf as rf
import numpy as np
import pandas as pd

from sub_THz_stripe.coupler.coupler import Coupler
from sub_THz_stripe.fiber.fiber import Fiber
from scipy.signal import unit_impulse
import matplotlib.pyplot as pyplt


def plot_am_am(k: int, func, points):
    # Make AM/AM plot
    fig, ax = pyplt.subplots()

    # Generate linear input data.
    x = np.linspace(0, 1, num=points)
    ax.plot(np.abs(x), np.abs(x), 'o', label="ideal")

    y = np.array([x])
    # Run the k stages.
    for i in range(k):
        y = func(y)
        ax.plot(np.abs(x), np.abs(y[0]), 'o', label=f"K={i+1}")

    ax.set_xlabel("Input amplitude |x(0)|")
    ax.set_ylabel("Output amplitude |y(k)|")
    ax.legend()

    return fig


def plot_am_pm(k: int, func, points):
    # Make AM/PM plot
    fig, ax = pyplt.subplots()

    # Generate linear input data.
    x = np.linspace(0, 1, num=points)

    y = np.array([x])
    phase_in = np.angle(x)
    # Run the k stages.
    for i in range(k):
        y = func(y)

        phase_out = np.angle(y[0])
        phase_change = phase_out - phase_in

        ax.plot(np.abs(x), phase_change, 'o', label=f"K={i+1}")

    ax.set_xlabel("Input amplitude |x(0)|")
    ax.set_ylabel("Phase change [rad]")
    ax.legend()

    return fig


def plot_spars(freqs, spars):
    fig, ax = pyplt.subplots()

    ax.plot(freqs / 1e9, 20*np.log10(np.abs(spars)))
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("S-parameter [dB]")

    return fig


if __name__ == "__main__":
    # How many stages to plot?
    k = 8

    # Load the couplers S-parameter file
    coupler_spars = rf.Network('models/coupler/with_balun.s2p')

    # Extract frequency and S21 (transmission)
    coup_freqs = coupler_spars.f
    coup_s21 = coupler_spars.s[:, 1, 0]  # S21

    fig = plot_spars(coup_freqs, coup_s21)
    fig.savefig("coupler_balun_s21.pdf")

    # Compute impulse response
    impulse_response = np.fft.ifft(coup_s21)

    damping = 0  # in dB
    cp = Coupler(damping, impulse_response)

    # Confirm that the impulse response is correct by applying it to a unit impulse
    x = unit_impulse(len(coup_freqs))
    y = cp.run(np.array([x]))

    fft = np.fft.fft(y[0])
    fig = plot_spars(coup_freqs, fft)
    fig.savefig("coupler_ir_verification.pdf")

    fig = plot_am_am(8, cp.run, len(coup_freqs))
    fig.savefig("coupler_am_am.pdf")

    fig = plot_am_pm(8, cp.run, len(coup_freqs))
    fig.savefig("coupler_am_pm.pdf")

    fiber_spars = pd.read_csv('models/PMF/1m_taped.csv')

    fib_freqs = fiber_spars["freq[Hz]"]
    fib_ang_rad = np.deg2rad(fiber_spars["ang:Trc2_S21"])
    fib_mag = 10 ** (fiber_spars["db:Trc2_S21"] / 20.0)
    fib_s21 = fib_mag * np.cos(fib_ang_rad) + 1j * \
        fib_mag * np.sin(fib_ang_rad)

    fig = plot_spars(fib_freqs, fib_s21)
    fig.savefig("fiber_1m_taped_s21.pdf")

    # Compute impulse response
    impulse_response = np.fft.ifft(fib_s21)

    fib = Fiber(1, 0, filter=impulse_response)

    # Confirm that the impulse response is correct by applying it to a unit impulse
    x = unit_impulse(len(fib_freqs))
    y = fib.run(np.array([x]))

    fft = np.fft.fft(y[0])
    fig = plot_spars(fib_freqs, fft)
    fig.savefig("fiber_ir_verification.pdf")

    fig = plot_am_am(8, fib.run, len(fib_freqs))
    fig.savefig("fiber_am_am.pdf")

    fig = plot_am_pm(8, fib.run, len(fib_freqs))
    fig.savefig("fiber_am_pm.pdf")
