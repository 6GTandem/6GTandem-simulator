import matplotlib.pyplot as plt
import numpy as np
import logging
from itertools import combinations
from scipy.signal import unit_impulse

from wireless_channel.waveforms import Waveform
from utils import remove_oversampling, cp_ofdm_to_freq, calculate_psd_per_symbol

logger = logging.getLogger(__name__)

booster_stages = ["fiber", "coupler", "amplifier", "coupler"]
tx_stages = ["fiber", "coupler", "splitter", "shifter", "amplifier"]
rx_stages = ["amplifier", "shifter", "combiner", "coupler", "fiber"]

def plot_im_data(imdata: list[np.ndarray], wf: Waveform):
    """Plot the intermediate data coming out of the stripe.

    Parameters
    ----------
        imdata (list[np.ndarray])
            aaa

    Returns
    -------
    matplotlib.Figure
        A matplotlib Figure object containing the created plot.
    """
    # Calculate how many booster stages are included in the intermediate data.
    stages = ((len(imdata) - 6) // 4) + 1
    # Make plot panes for all the stages
    plot_panes = stages
    if plot_panes % 2 != 0:
        plot_panes += 1

    # Make a figure for AM/AM plots and spec plots.
    fig, ax = plt.subplots(plot_panes // 2, 2)
    fig2, ax2 = plt.subplots(plot_panes // 2, 2)

    # Keep the array 1xm for simplicity.
    if len(ax.shape) == 1:
        ax = np.array([ax])
    if len(ax2.shape) == 1:
        ax2 = np.array([ax2])

    # Set the proper axis labels.
    ax[-1, 0].set_xlabel("Input amplitude |x|")
    ax[-1, -1].set_xlabel("Input amplitude |x|")
    ax[0, 0].set_ylabel("Output amplitude |y|")

    ax2[-1, 0].set_xlabel("Normalized Frequency")
    ax2[-1, -1].set_xlabel("Normalized Frequency")
    ax2[0, 0].set_ylabel("Power Spectral Density (dB/Hz)")

    # Loop over the data from all the stages. (Stages are RUs and BUs.)
    for i, (x, y) in enumerate(zip(imdata[:-1], imdata[1:])):
        # Calculate at which stage we are.
        # imdata is ((4 x boosters) + 5) + 1
        # Where 5 is the last transmit stage and the + 1 is the initial data being transmitted
        # at location 0.
        stage = (i // 4) + 1
        # The last stage is a RU containing 5 components instead of 4. We need to
        # compensate for this or we advance to a non-existing stage.
        if stage > stages:
            stage = stages

        # For the last three stages --- splitter, shifter and amplifier --- we select the data from
        # one single antenna to plot.
        component_in_last_stage = (i - ((stages - 1) * 4)) % 5
        if stage == stages:
            if component_in_last_stage >= 3:
                x = x[0]
            if component_in_last_stage >= 2:
                y = y[0]

        # Calculate the PSD.
        freq, psd = calculate_psd_per_symbol(y, wf.fs)

        if stage >= stages:
            label = tx_stages[component_in_last_stage]
            label += str(stages)
        else:
            label = booster_stages[i % 4]
            label += str(stage)
        column = 0
        if stage > (plot_panes // 2):
            column = 1
        row = (stage - 1) % (plot_panes // 2)

        # Select the data for one symbol only.
        ax[row, column].plot(np.abs(x[0]), np.abs(y[0]), "o", label=label)
        ax[row, column].legend()

        ax2[row, column].plot(freq, psd, label=label)
        ax2[row, column].legend()

    return fig, fig2

def plot_stripes(config: dict, stripes: list):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_aspect("auto")

    # Plot the room (corners and edges)
    x_size = config["room"]["x"]
    y_size = config["room"]["y"]
    z_size = config["room"]["z"]

    points = []
    for x in [0, x_size]:
        for y in [0, y_size]:
            for z in [0, z_size]:
                points.append([x, y, z])
    points = np.array(points)
    ax.scatter3D(points[:, 0], points[:, 1], points[:, 2], c="black")
    for s, e in combinations(points, 2):
        diff = list(s - e)
        if diff.count(0) == 2:
            ax.plot3D(*zip(s, e), color="k")

    ru_handle = None
    cu_handle = None

    # Plot the stripes
    for stripe in stripes:
        units = []
        for ru in stripe.radio_units:
            units.append([ru.x, ru.y, ru.z])
        units = np.array(units)
        ru_handle = ax.scatter3D(units[:, 0], units[:, 1], units[:, 2], c="red", label="Radio Unit")
        ax.plot3D(units[:, 0], units[:, 1], units[:, 2], color="red")

        cu_handle = ax.scatter3D(
            stripe.central_unit.x, stripe.central_unit.y, stripe.central_unit.z, c="yellow", label="Central Unit"
        )

    # Add XYZ coordinate system at the origin
    origin = np.array([[0, 0, 0]])
    axes = np.eye(3)
    ax.quiver(
        origin[:, 0],
        origin[:, 1],
        origin[:, 2],
        axes[:, 0],
        axes[:, 1],
        axes[:, 2],
        color=["k", "k", "k"],
        length=0.5,
        normalize=True,
    )
    ax.text(1, 0, 0, "X", color="k")
    ax.text(0, 1, 0, "Y", color="k")
    ax.text(0, 0, 1, "Z", color="k")

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())


def plot_room(config: dict):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_aspect("auto")

    points = []

    x_size = config["room"]["x"]
    y_size = config["room"]["y"]
    z_size = config["room"]["z"]

    ax.set_xlim([0, config["room"]["x"]])
    ax.set_ylim([0, config["room"]["y"]])
    ax.set_zlim([0, config["room"]["z"]])

    for x in [0, x_size]:
        for y in [0, y_size]:
            for z in [0, z_size]:
                points.append([x, y, z])

    points = np.array(points)
    ax.scatter3D(points[:, 0], points[:, 1], points[:, 2], c="black")
    for s, e in combinations(points, 2):
        diff = list(s - e)
        if diff.count(0) == 2:
            ax.plot3D(*zip(s, e), color="k")

    # Plot the stripes
    logger.debug("%d Stripes found", len(config["radio_stripes"]))
    for stripe in config["radio_stripes"]:
        units = []
        for unit in stripe:
            if "radio_unit" in unit:
                unit = unit["radio_unit"]
                if "x" in unit:
                    units.append([unit["x"], unit["y"], unit["z"]])
        units = np.array(units)
        logger.debug("Stripe units: %s", units)
        ax.scatter3D(units[:, 0], units[:, 1], units[:, 2], c="red")
        ax.plot3D(units[:, 0], units[:, 1], units[:, 2], color="red")

    # Plot the UEs
    ues = []
    logger.debug("%d UEs found", len(config["ue_positions"]))
    for ue_pos in config["ue_positions"]:
        if "x" in ue_pos:
            ues.append([ue_pos["x"], ue_pos["y"], ue_pos["z"]])
    ues = np.array(ues)
    logger.debug("UE positions: %s", ues)
    ax.scatter3D(ues[:, 0], ues[:, 1], ues[:, 2], c="blue", alpha=0.2)


def plot_constellation(iq_symbols: np.ndarray, labels: list[str] | None = None, title: str | None = None):
    """Plot a constellation of IQ symbols.

    Parameters
    ----------
    iq_symbols: np.ndarray
        Array containing the IQ symbols.
        Format nxm where n are different constellations plotted in different colors.
    labels: list[str] | None
        Labels for the different IQ constellations.
        No legend is plotted if labels is None.

    Returns
    -------
    matplotlib.Figure
        A matplotlib Figure object containing the created plot.
    """
    fig, ax = plt.subplots()

    const_labels = [""] * iq_symbols.shape[0]
    if labels is not None:
        const_labels = labels

    # Loop over all the provided constellations.
    for const, label in zip(iq_symbols, const_labels):
        ax.scatter(np.real(const), np.imag(const), label=label)

    # Provide axis labels and a legend.
    ax.set_xlabel("In-phase (I)")
    ax.set_ylabel("Quadrature (Q)")
    ax.axis("equal")
    if labels is not None:
        ax.legend()
    ax.set_title(title if title is not None else "IQ constellation")

    return fig


def plot_iq_per_symbol(
    cp_ofdm_time_data: np.ndarray, prefix_length: int, n_carriers: int, labels: list[str] | None = None
):
    """Plot the IQ constellation for a set of CP-OFDM symbols.

    Parameters
    ----------
    cpofdm_time_data: np.ndarray
        The CP-OFDM signal in the time domain.
    prefix_length: int
        See `utils.cp_ofdm_to_freq`.
    n_carriers: int
        See `remove_oversampling`.
    labels: list[str] | None
        See `plot_constellation`.

    Returns
    -------
        See `plot_constellation`.
    """
    # First transform the CP-OFDM signal from the time to the frequency domain using FFT.
    ofdm_freqs = cp_ofdm_to_freq(cp_ofdm_time_data, prefix_length)
    # Convert the frequency signal into QAM symbols which can be plotted in an IQ-plane.
    ofdm_us = remove_oversampling(ofdm_freqs, n_carriers)

    return plot_constellation(ofdm_us)


def plot_iq_per_carrier(
    cp_ofdm_time_data: np.ndarray, prefix_length: int, n_carriers: int, labels: list[str] | None = None
):
    """Plot the IQ constellation of all CP-OFDM signal for one symbol.

    Parameters
    ----------
    cpofdm_time_data: np.ndarray
        The CP-OFDM signal in the time domain containing a single symbol.
        Shape: 1xm
    prefix_length: int
        See `utils.cp_ofdm_to_freq`.
    n_carriers: int
        See `remove_oversampling`.
    labels: list[str] | None
        See `plot_constellation`.

    Returns
    -------
        See `plot_constellation`.
    """
    # First transform the CP-OFDM signal from the time to the frequency domain using FFT.
    ofdm_freqs = cp_ofdm_to_freq(cp_ofdm_time_data, prefix_length)
    # Convert the frequency signal into QAM symbols which can be plotted in an IQ-plane.
    ofdm_us = remove_oversampling(ofdm_freqs, n_carriers)

    # Put every subcarrier in a single row so that it is plotted in a different color.
    ofdm_carriers = ofdm_us.reshape(ofdm_us.shape[1], 1)

    return plot_constellation(ofdm_carriers)


def plot_psd_per_symbol(time_signal: np.ndarray, fs: float = 1, N: int = 1024):
    """Create a Power Spectral Density Plot (PSD) for every symbol.

    The PSD is calculated using Welch`s method with a 'hann' window.

    Parameters
    ----------
    time_signal: np.ndarray
        Time domain signal for which to plot the PSD.
    fs: float
        See `utils.calculate_psd_per_symbol`.
    N: int
        See `utils.calculate_psd_per_symbol`.

    Returns
    -------
    matplotlib.Figure
        Matplotlib Figure object containing the plot.
    """
    fig, ax = plt.subplots()

    # Loop over all the symbols and plot them.
    for n, symbol in enumerate(time_signal):
        f, s = calculate_psd_per_symbol(symbol, fs, N)

        ax.plot(f, s, label=f"Symbol{n+1}")

    ax.set_xlabel("Normalized Frequency")
    ax.set_ylabel("Power Spectral Density (dB/Hz)")
    ax.legend()
    ax.grid()

    return fig


def verify_impulse_response(func, freqs: np.ndarray):
    """Confirm that the impulse response for a certain component is correct.

    Verification is performed by applying the impulse response to a unit impulse.
    After performing an FFT on the output the result must be the original filter
    response.

    Parameters
    ----------
    func: Any
        Method for which to confirm that the filter works.
    freqs: np.ndarray
        Array containing the frequency points of the impulse response.

    Returns
    -------
    Figure
        Matplotlib Figure object containing the FFT plot of the filter output.
    """
    # Confirm that the impulse response is correct by applying it to a unit impulse.
    x = unit_impulse(len(freqs))
    x = np.concatenate([np.zeros(128), x])
    y = func(np.array([x]))

    fft = np.fft.fft(y[0][128:])

    fig, ax = plt.subplots()

    ax.plot(freqs, 20 * np.log10(np.abs(fft)))
    ax.set_xlabel("Normalized frequency [GHz]")
    ax.set_ylabel("S-parameter [dB]")

    return fig
