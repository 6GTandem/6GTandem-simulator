import logging
import numpy as np

from scipy.signal import welch, get_window

# Project-wide logger instance
logger = logging.getLogger("6GTandem")
logger.setLevel(logging.DEBUG)
if not logger.hasHandlers():
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

def remove_oversampling(ofdm_freqs: np.ndarray, n_carriers: int):
    """Remove the oversampling from an OFDM signal in the frequency domain.

    Oversampling adds zero symbols to the OFDM signal to increase the resolution.
    These symbols however contain no information and can be removed for plotting.
    
    Parameters
    ----------
    cpofdm_freqs: np.ndarray
        The OFDM signal in the frequency domain.
    n_carriers: int
        Number of carriers used for the original OFDM signal.

    Returns
    -------
    np.ndarray
        OFDM signal with the zeros removed.
    """
    qam_received = []

    for symbol in ofdm_freqs:
        half = (n_carriers) // 2
        row_symbols = np.concatenate([symbol[:half], symbol[-half:]])
        qam_received.append(row_symbols)
    return np.array(qam_received)

def cp_ofdm_to_freq(cp_ofdm_signal: np.ndarray, prefix_length: int):
    """Convert a CP-OFDM signal in the time domain to the frequency domain (FFT).

    The signal in the frequency domain no longer contains the cyclic prefix.

    Parameters
    ----------
    cp_ofdm_signal: np.ndarray
        The CP-OFDM signal in the time domain.
    prefix_length: int
        Length of the cyclic prefix in samples.

    Returns
    -------
    np.ndarray
        The CP-OFDM signal in the frequency domain.
        Shape: n_ofdm_symbols x (n_carriers * oversampling)
    """
    # Exclude the cyclic prefix from the FFT.
    n_fft = cp_ofdm_signal.shape[1] - prefix_length

    ofdm_freq = []
    # Apply an FFT to every symbol after removing the cyclic prefix.
    for symbol in cp_ofdm_signal:
        time_no_cp = symbol[prefix_length:]
        freq_domain = np.fft.fft(time_no_cp, n_fft)
        ofdm_freq.append(freq_domain)

    return np.array(ofdm_freq)

def ofdm_to_time(ofdm_freqs: np.ndarray, prefix_length: int):
    """Convert a OFDM signal in the frequency domain to the time domain.
    
    A cyclic prefix is added to the time domain signal.

    Parameters
    ----------
    ofdm_freqs: np.ndarray
        OFDM signal in the frequency domain.
    prefix_length: int
        See `cp_ofdm_to_freq`.

    Returns
    -------
    np.ndarray
        The time signal.
        Shape: n_ofdm_symbols x (n_carriers * oversampling + cp_length)
    """
    n_fft = ofdm_freqs.shape[1]

    cp_ofdm_time = []
    # Loop over every symbol, convert it to the time domain and add a cyclic prefix.
    for symbol in ofdm_freqs:
        time_domain = np.fft.ifft(symbol, n_fft)
        cp = time_domain[-prefix_length:]
        ofdm_symbol = np.concatenate([cp, time_domain])
        cp_ofdm_time.append(ofdm_symbol)

    return np.array(cp_ofdm_time)

def calculate_psd_per_symbol(time_signal: np.ndarray, fs: float = 1, N: int = 1024):
    """Calculate the Power Spectral Density Plot (PSD) for every symbol.

    The PSD is calculated using Welch`s method with a 'hann' window.
    
    Parameters
    ----------
    time_signal: np.ndarray
        Time domain signal for which to plot the PSD.
    fs: float
        Sampling frequency of time_signal.
    N: int
        Length of the segments used for the Welch method.
    
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Array containing the calculated PSD.
    """
    # Loop over all the symbols and calculate their PSD.
    psds = []
    freqs = []
    for symbol in time_signal:
        N = min(N, len(symbol) - 1)
        window = get_window("hann", N)
        f, Pxx = welch(
            symbol,
            fs=fs,
            window=window,
            nperseg=N,
            return_onesided=False,
            scaling="density",
        )
        psd = 10 * np.log10(Pxx)

        freqs.append(f)
        psds.append(psd)
    
    return np.array(freqs), np.array(psds)