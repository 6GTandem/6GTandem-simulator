import numpy as np
import pandas as pd
import skrf as rf
import matplotlib.pyplot as plt

from sub_THz_stripe.amplifier.amplifier import Amplifier


def polynomial_fitting(order, gain, max_input_amplitude, filetype: str = "pdf"):
    """Plot how well the polynomial of the amplifier fits the model."""
    # Read out the measurement data.
    data = pd.read_csv("models/amplifier/Mag_150GHz_BalCascode2stage_V3.vcsv", skiprows=6, names=["x", "yre", "yim"])
    # Go from the input values in dBm to values in Watt.
    x_power = (10 ** (data["x"] / 10)) / 1000
    # Now convert the values in Watt to Voltages.
    x_voltage = np.sqrt(x_power * 50)
    # Calculate the magnitude of the output voltage based on the real and complex values.
    y_voltage = np.sqrt(data["yre"] ** 2 + data["yim"] ** 2)
    # What is the current gain of the amplifier?
    current_gain = (y_voltage[4] - y_voltage[0]) / (x_voltage[4] - x_voltage[0])
    scale = current_gain / gain
    # Re-scale the input and output to get the new gain.
    x_scaled = x_voltage * scale

    # Only perform fitting for the odd powers.
    odd_powers = list(range(1, order + 1, 2))
    x_powers = np.column_stack([x_scaled**p for p in odd_powers])

    coeffs, *_ = np.linalg.lstsq(x_powers, y_voltage, rcond=None)

    xlim = x_scaled.copy()
    xlim[np.abs(x_scaled) > max_input_amplitude] = max_input_amplitude
    xout = sum(c * xlim * np.abs(xlim) ** (p - 1) for c, p in zip(coeffs, range(1, order + 1, 2)))

    fig, ax = plt.subplots()

    ax.plot(x_scaled, y_voltage, "-o", label="model")
    ax.plot(x_scaled, xout, "-o", label="fitted")

    ax.set_ylabel("Output Amplitude |y|")
    ax.set_xlabel("Input Amplitude |x|")
    ax.grid(which="both")
    ax.legend()

    fig.savefig(f"models/amplifier/polynomial-fit.{filetype}")


def plot_amplifier(nf: int = 0, filetype: str = "pdf"):
    """Make plots of the amplifier model imported from Chalmers.

    Parameters
    ----------
        nf float
            The noise figure of the amplifier in dB.
    """
    mode = "poly5"
    smoothness = 1

    fig, ax = plt.subplots()

    for gain in range(2, 11, 2):
        x = np.random.rand(3000) * 0.5
        amp = Amplifier(bw=3e9, gain=gain, max_gain=30, noise_fig=nf, smoothness=smoothness, mode=mode)

        xout = amp.run(x)

        ax.plot(x, np.abs(xout), "o", label=gain, alpha=0.6)

    ax.set_ylabel("Output Amplitude |y|")
    ax.set_xlabel("Input Amplitude |x|")
    ax.grid(which="both")
    ax.legend()

    fig.savefig(f"models/amplifier/amam-plot-gains-nf{nf}db.{filetype}")


def plot_pmf(deembedding: float = 3.0, filetype: str = "pdf"):
    fig, ax = plt.subplots()

    for filename in ["pmf1m", "pmf2m"]:
        fiber_spars = rf.Network(f"models/PMF/{filename}.s2p")

        fib_freqs = fiber_spars.f
        fib_s21 = fiber_spars.s[:, 1, 0]

        fib_mag = np.abs(fib_s21)
        fib_phase = np.angle(fib_s21)

        ax.plot(fib_freqs / 1e9, fib_mag + deembedding, label=filename.lstrip("pmf"))

    ax.set_ylabel("Fiber Loss [dB]")
    ax.set_xlabel("Frequency [GHz]")
    ax.grid(which="both")
    ax.legend()

    fig.savefig(f"models/PMF/pmf.{filetype}")


def plot_coupler(filetype: str = "pdf"):
    fig, ax = plt.subplots()

    for filename in ["with_balun", "without_balun"]:
        coupler_spars = rf.Network(f"models/coupler/{filename}.s2p")

        coupler_freqs = coupler_spars.f
        coupler_s21 = coupler_spars.s[:, 1, 0]

        coupler_mag = np.abs(coupler_s21)
        coupler_phase = np.angle(coupler_s21)

        ax.plot(coupler_freqs / 1e9, coupler_mag, label=filename.replace("_", " "))

    ax.set_ylabel("Coupler Loss [dB]")
    ax.set_xlabel("Frequency [GHz]")
    ax.grid(which="both")
    ax.legend()

    fig.savefig(f"models/coupler/coupler.{filetype}")


def convert_files():
    files = ["PMF/pmf1m", "PMF/pmf2m", "coupler/with_balun", "coupler/without_balun"]
    for filename in files:
        spars = rf.Network(f"models/{filename}.s2p")

        freqs = spars.f
        s21 = spars.s[:, 1, 0]

        mag = np.abs(s21)
        phase = np.angle(s21)

        df = pd.DataFrame()
        df = df.assign(frequency=freqs, magnitude=mag, phase=phase)
        df.to_csv(f"models/{filename}.csv", index=False)


def process_amplifier(nf: int = 0):
    """Make csvs for the amplifier model imported from Chalmers."""
    mode = "poly5"
    smoothness = 1

    for gain in range(2, 11, 2):
        x = np.random.rand(3000) * 0.5
        amp = Amplifier(bw=3e9, gain=gain, max_gain=30, noise_fig=nf, smoothness=smoothness, mode=mode)

        data = amp.run(x)

        df = pd.DataFrame()
        df = df.assign(x=x, y=np.abs(data))
        df.to_csv(f"models/amplifier/amplifier-{gain}x-nf{nf}db.csv", index=False)

    # Read out the measurement data.
    data = pd.read_csv("models/amplifier/Mag_150GHz_BalCascode2stage_V3.vcsv", skiprows=6, names=["x", "yre", "yim"])
    # Go from the input values in dBm to values in Watt.
    x_power = (10 ** (data["x"] / 10)) / 1000
    # Now convert the values in Watt to Voltages.
    x_voltage = np.sqrt(x_power * 50)
    # Calculate the magnitude of the output voltage based on the real and complex values.
    y_voltage = np.sqrt(data["yre"] ** 2 + data["yim"] ** 2)

    # Only perform fitting for the odd powers.
    odd_powers = list(range(1, 5 + 1, 2))
    x_powers = np.column_stack([x_voltage**p for p in odd_powers])

    coeffs, *_ = np.linalg.lstsq(x_powers, y_voltage, rcond=None)

    xout = sum(c * x_voltage * np.abs(x_voltage) ** (p - 1) for c, p in zip(coeffs, range(1, 5 + 1, 2)))
    xout = np.abs(xout)

    df = pd.DataFrame()
    df = df.assign(x=x_voltage, y_model=y_voltage, y_fitted=xout)
    df.to_csv(f"models/amplifier/polynomial-fit.csv", index=False)


if __name__ == "__main__":
    plot_amplifier(nf=0)
    plot_amplifier(nf=10)
    plot_amplifier(nf=70)
    polynomial_fitting(5, 1, 5)
    plot_pmf()
    plot_coupler()
