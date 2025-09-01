import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch

class Waveform():
    """Definition of a waveform

    Components:
    Usage:
        >>> # Read the YAML file
        >>> with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
        >>>    config = yaml.safe_load(file)
        >>> # extract required params
        >>> waveform_config = config["waveform_config"]
        >>> freq_band_config = config['sub_thz']
        >>> # construct class instance
        >>> wf = Waveform.from_config(waveform_config, freq_band_config)
    """

    def __init__(self, waveform_type: str = 'cp-ofdm', **kwargs):
        self.waveform_type = waveform_type

        if self.waveform_type == "cp-ofdm":
            # OFDM-specific parameters with defaults
            self.n_ofdm_symbols = kwargs.get("n_ofdm_symbols", 2)
            self.n_carriers = kwargs.get("n_carriers", 1024)
            self.qam_order = kwargs.get("qam_order", 4)
            self.oversampling_factor = kwargs.get("oversampling_factor", 4)
            self.cp_length = kwargs.get("cp_length", 128)
            self.BW = kwargs.get("BW", 12.5e9)
            self.fc = kwargs.get("fc", 157.75e9)
            self.fs = self.BW * self.oversampling_factor
            self.fft_size = self.n_carriers * self.oversampling_factor

        else:
            raise ValueError(f"Unsupported waveform type: {self.waveform_type}")

    @classmethod
    def from_config(cls, waveform_config, freq_band_config):
        """
        Construct a Waveform from a configuration dictionary.

        :param waveform_config: dictionary containing waveform parameters
        :param freq_band_config: dictionary containing frequency band parameters (fc, bw, num_carriers)
        :return: Waveform instance
        """
        wf = cls(waveform_type=waveform_config['waveform_type'],
                 n_ofdm_symbols=waveform_config['n_ofdm_symbols'],
                 n_carriers=freq_band_config['num_subcarriers'],
                 qam_order=waveform_config['qam_order'],
                 oversampling_factor=waveform_config['oversampling_factor'],
                 cp_length=waveform_config['cp_length'],
                 BW=freq_band_config['bw'],
                 fc=freq_band_config['fc']
                 )
        return wf

    def __str__(self):
        """Return a human-readable summary of the waveform."""
        info = f"Waveform Type: {self.waveform_type}\n"
        if self.waveform_type == "cp-ofdm":
            info += (
                f"  Number of OFDM symbols: {self.n_ofdm_symbols}\n"
                f"  Number of carriers: {self.n_carriers}\n"
                f"  QAM order: {self.qam_order}\n"
                f"  Oversampling factor: {self.oversampling_factor}\n"
                f"  Cyclic prefix length: {self.cp_length}\n"
                f"  Bandwidth: {self.BW/1e9:.2f} GHz\n"
                f"  Carrier frequency: {self.fc/1e9:.2f} GHz\n"
                f"  Sampling frequency: {self.fs/1e9:.2f} GHz\n"
            )
        return info

    # ------------------------
    # OFDM methods
    # ------------------------
    def generate_bits(self):
        num_bits = self.n_ofdm_symbols * self.n_carriers * int(np.log2(self.qam_order))
        self.bits = np.random.randint(0, 2, num_bits)
        return self.bits

    def qam_modulate(self):
        bits = self.bits
        k = int(np.log2(self.qam_order))
        if len(bits) % k != 0:
            raise ValueError("Number of bits must be a multiple of log2(qam_order).")

        bit_groups = bits.reshape((-1, k))
        symbols_idx = np.array([int("".join(str(b) for b in grp), 2) for grp in bit_groups])
        m_side = int(np.sqrt(self.qam_order))

        def gray_map(x): return x ^ (x >> 1)

        I = gray_map(symbols_idx % m_side)
        Q = gray_map(symbols_idx // m_side)
        I = 2 * I - (m_side - 1)
        Q = 2 * Q - (m_side - 1)
        self.qam_symbols = (I + 1j * Q) / np.sqrt((2 / 3) * (self.qam_order - 1))
        return self.qam_symbols

    def plot_constellation(self):
        if self.qam_symbols is None:
            raise ValueError("QAM symbols not generated yet.")
        plt.scatter(np.real(self.qam_symbols), np.imag(self.qam_symbols))
        plt.title(f"{self.qam_order}-QAM Constellation")
        plt.xlabel("I")
        plt.ylabel("Q")
        plt.axis("equal")
        plt.show()

    def ofdm_modulate(self):
        qam_symbols = self.qam_symbols
        n_symbols = len(qam_symbols) // self.n_carriers
        qam_symbols = qam_symbols[:n_symbols * self.n_carriers]
        qam_matrix = qam_symbols.reshape((n_symbols, self.n_carriers))
        ofdm_time = []
        for row in qam_matrix:
            freq_domain = np.zeros(self.fft_size, dtype=complex)

            # Map QAM symbols to center of spectrum (baseband symmetric)
            half = self.n_carriers // 2
            freq_domain[:half] = row[:half]
            freq_domain[-half:] = row[half:]

            # IFFT to get time-domain signal
            time_domain = np.fft.ifft(freq_domain, self.fft_size)

            # Add cyclic prefix
            cp = time_domain[-self.cp_length:]
            ofdm_symbol = np.concatenate([cp, time_domain])
            ofdm_time.append(ofdm_symbol)

        self.ofdm_time = np.array(ofdm_time)
        return self.ofdm_time  # shape: n_ofdm_symbols x ( n_carriers * oversampling + cp_length)

    def ofdm_time_to_freq(self, received_time=None):
        if received_time is None:
            received_time = self.ofdm_time
        n_fft = received_time.shape[1] - self.cp_length
        ofdm_freq = []
        for symbol in received_time:
            time_no_cp = symbol[self.cp_length:]
            freq_domain = np.fft.fft(time_no_cp, n_fft)
            ofdm_freq.append(freq_domain)
        self.ofdm_freq = np.array(ofdm_freq)
        return self.ofdm_freq # shape: n_ofdm_symbols x (n_carriers * oversampling)

    def ofdm_freq_to_time(self, ofdm_freq=None):
        if ofdm_freq is None:
            ofdm_freq = self.ofdm_freq
        ofdm_time = []
        n_fft = ofdm_freq.shape[1]
        for symbol in ofdm_freq:
            time_domain = np.fft.ifft(symbol, n_fft)
            cp = time_domain[-self.cp_length:]
            ofdm_symbol = np.concatenate([cp, time_domain])
            ofdm_time.append(ofdm_symbol)
        return np.array(ofdm_time) # shape: n_ofdm_symbols x ( n_carriers * oversampling + cp_length)

    def pad_subcarriers(self, subcarriers):
        """
        Pad zeros in the middle of the subcarriers to restore the oversampled FFT size.

        Parameters
        ----------
        subcarriers : np.ndarray
            Frequency-domain signal with n_carriers points.
            Shape: (n_ofdm_symbols, n_carriers)

        Returns
        -------
        np.ndarray
            Frequency-domain signal with fft_size points (oversampled).
            Shape: (n_ofdm_symbols, fft_size)
        """
        n_symbols, n_carriers = subcarriers.shape
        if n_carriers != self.n_carriers:
            raise ValueError(f"Input must have {self.n_carriers} carriers, got {n_carriers}.")

        fft_size = self.fft_size
        half = n_carriers // 2
        freq_oversampled = np.zeros((n_symbols, fft_size), dtype=complex)

        # Copy the first half to the beginning
        freq_oversampled[:, :half] = subcarriers[:, :half]
        # Copy the second half to the end
        freq_oversampled[:, -half:] = subcarriers[:, half:]

        return freq_oversampled


    def extract_subcarriers(self, ofdm_freq_oversampled=None):
        """
        Extract the original subcarriers from the oversampled OFDM frequency-domain signal.

        Parameters
        ----------
        ofdm_freq_oversampled : np.ndarray, optional
            Frequency-domain OFDM signal with fft_size points.
            Shape: (n_ofdm_symbols, fft_size)
            If None, uses self.ofdm_freq.

        Returns
        -------
        np.ndarray
            Frequency-domain signal with n_carriers points only.
            Shape: (n_ofdm_symbols, n_carriers)
        """
        if ofdm_freq_oversampled is None:
            ofdm_freq_oversampled = self.ofdm_freq

        n_symbols, fft_size = ofdm_freq_oversampled.shape
        half = self.n_carriers // 2

        # Take the first half and last half of the original carriers
        subcarriers = np.zeros((n_symbols, self.n_carriers), dtype=complex)
        subcarriers[:, :half] = ofdm_freq_oversampled[:, :half]
        subcarriers[:, half:] = ofdm_freq_oversampled[:, -half:]

        return subcarriers

    def awgn(self, signal, snr_dB):
        snr_linear = 10 ** (snr_dB / 10)
        power_signal = np.mean(np.abs(signal) ** 2)
        noise_power = power_signal / snr_linear
        noise = np.sqrt(noise_power / 2) * (np.random.randn(*signal.shape) + 1j * np.random.randn(*signal.shape))
        return signal + noise

    def ofdm_to_qam(self, ofdm_freq_received=None):
        if ofdm_freq_received is None:
            ofdm_freq_received = self.ofdm_freq
        qam_received = []
        for sym in ofdm_freq_received:
            half = self.n_carriers // 2
            row_symbols = np.concatenate([sym[:half], sym[-half:]])
            qam_received.append(row_symbols)
        return np.array(qam_received).flatten()

    def qam_to_bits(self, qam_symbols):
        m_side = int(np.sqrt(self.qam_order))
        k = int(np.log2(self.qam_order))
        qam_symbols = qam_symbols * np.sqrt((2 / 3) * (self.qam_order - 1))
        I = np.clip(np.round((np.real(qam_symbols) + (m_side - 1) / 2)), 0, m_side - 1).astype(int)
        Q = np.clip(np.round((np.imag(qam_symbols) + (m_side - 1) / 2)), 0, m_side - 1).astype(int)

        def inv_gray(x):
            y = np.zeros_like(x)
            for shift in range(int(np.log2(m_side))):
                y ^= x >> shift
            return y

        I = inv_gray(I)
        Q = inv_gray(Q)
        symbols_idx = Q * m_side + I
        bits = (((symbols_idx[:, None] & (1 << np.arange(k)[::-1])) > 0)).astype(int)
        return bits.flatten()

    def compute_ber(self, bits_tx, bits_rx):
        if len(bits_tx) != len(bits_rx):
            raise ValueError("Transmitted and received bits must have the same length.")
        n_errors = np.sum(bits_tx != bits_rx)
        return n_errors / len(bits_tx)

    def plot_constellation(self, symbols_tx, symbols_rx):
        plt.scatter(np.real(symbols_rx), np.imag(symbols_rx), label='Rx symbols')
        plt.scatter(np.real(symbols_tx), np.imag(symbols_tx), label='Tx symbols')
        plt.title(f"{self.qam_order}-QAM Constellation")
        plt.xlabel("In-phase (I)")
        plt.ylabel("Quadrature (Q)")
        plt.axis("equal")
        plt.legend()
        plt.show()

    def plot_psd(self, ofdm_time):
        if ofdm_time is None:
            ofdm_time = self.ofdm_time

        ofdm_signal = ofdm_time.flatten()
        f, Pxx = welch(ofdm_signal, fs=self.fs, nperseg=1024, return_onesided=False)
        # Center around zero and shift
        Pxx_db = 10 * np.log10(np.fft.fftshift(Pxx))
        f_shifted = np.fft.fftshift(f)  # + fc  # shift to RF
        f_shifted_GHz = f_shifted / 1e9
        plt.figure(figsize=(10, 4))
        plt.plot(f_shifted_GHz, Pxx_db)
        plt.xlabel("Frequency (GHz)")
        plt.xlim([-10, 10])
        plt.ylabel("PSD (dB)")
        plt.title("OFDM PSD at RF")
        plt.grid(True)
        plt.show()