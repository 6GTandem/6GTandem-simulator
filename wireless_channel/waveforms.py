import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from scipy.interpolate import interp1d

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
            self.bw = kwargs.get("BW", 12.5e9)
            self.fc = kwargs.get("fc", 157.75e9)
            self.fs = self.bw * self.oversampling_factor
            self.fft_size = self.n_carriers * self.oversampling_factor

            # Pilot configuration
            self.pilot_spacing = kwargs.get("pilot_spacing", 16)  # every 16th carrier is a pilot by default
            self.pilot_mode = kwargs.get("pilot_mode", "interleaved") # "interleaved" = pilots interleaved every pilot_spacing
            # 'block' = first OFDM symbol is full pilot (all carriers), subsequent symbols are pure data



            self.tx_power = kwargs.get("tx_power", 10)

            # note that interleaved pilots might not work if the channel is very uncorrelated accross the carries
            # this is due to the ineffectiveness of the interpolation between the data and pilot carriers
            # as a solution block pilots can be used, where the first OFDM symbol is all pilots, and the remaining symbols are pure data

            if self.pilot_mode not in ("interleaved", "block"):
                raise ValueError("pilot_mode must be 'interleaved or 'block'")

            # derived
            if self.pilot_mode == "interleaved":
                self.pilot_indices = np.arange(0, self.n_carriers, self.pilot_spacing)
                self.n_pilots = len(self.pilot_indices)
                self.data_carriers = np.setdiff1d(np.arange(self.n_carriers), self.pilot_indices)
            elif self.pilot_mode == "block":
                self.n_pilots = self.n_carriers  # all carriers in first symbol are pilots
                self.pilot_indices = np.arange(0, self.n_carriers)

            # todo remove pilot_symbol from config
            self.pilot_symbols = self.set_pilots()
            # #self.pilot_symbol = kwargs.get("pilot_symbol", (1 + 1j))
            # self.pilot_symbol = kwargs.get("pilot_symbol", (1 + 1j) / np.sqrt((2 / 3) * (self.qam_order - 1)))

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
                 fc=freq_band_config['fc'],
                 pilot_spacing= waveform_config['pilot_spacing'],
                 pilot_symbol= waveform_config['pilot_symbol_re'] + 1j * waveform_config['pilot_symbol_imag'],
                 pilot_mode = waveform_config['pilot_mode'],
                 tx_power = waveform_config['tx_power']
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
                f"  Bandwidth: {self.bw/1e9:.2f} GHz\n"
                f"  Carrier frequency: {self.fc/1e9:.2f} GHz\n"
                f"  Sampling frequency: {self.fs/1e9:.2f} GHz\n"
                f" Pilot mode: {self.pilot_mode}"
                f" Pilot spacing: {self.pilot_spacing}\n"
                f" Number of pilots: {self.n_pilots}\n"
            )
        return info

    # ------------------------
    # OFDM methods
    # ------------------------

    def set_pilots(self):
        """
        generate pilot symbols as random QAM symbols
        """
        k = int(np.log2(self.qam_order))
        num_pilot_bits = self.n_pilots * k
        pilot_bits = np.random.randint(0, 2, num_pilot_bits)


        if len(pilot_bits) % k != 0:
            raise ValueError("Number of pilot bits must be a multiple of log2(qam_order).")

        bit_groups = pilot_bits.reshape((-1, k))
        symbols_idx = np.array([int("".join(str(b) for b in grp), 2) for grp in bit_groups])
        m_side = int(np.sqrt(self.qam_order))

        def gray_map(x): return x ^ (x >> 1)

        I = gray_map(symbols_idx % m_side)
        Q = gray_map(symbols_idx // m_side)
        I = 2 * I - (m_side - 1)
        Q = 2 * Q - (m_side - 1)
        pilot_symbols = (I + 1j * Q) / np.sqrt((2 / 3) * (self.qam_order - 1))
        return pilot_symbols

    def generate_bits(self):
        # number of data carriers per OFDM symbol depends on pilot_mode
        if self.pilot_mode == "interleaved":
            data_per_symbol = self.n_carriers - self.n_pilots
            n_data_symbols = self.n_ofdm_symbols
        else:  # block pilots: first symbol all pilots, remaining symbols are data
            if self.n_ofdm_symbols < 2:
                raise ValueError("For block pilot mode you need at least 2 OFDM symbols (one pilot + data).")
            data_per_symbol = self.n_carriers
            n_data_symbols = self.n_ofdm_symbols - 1

        num_bits = n_data_symbols * data_per_symbol * int(np.log2(self.qam_order))
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

    @staticmethod
    def generate_zadoff_chu_sequence(root, length):
        """
        Generate a Zadoff-Chu sequence.

        :param root: The root index of the sequence (integer, coprime with length).
        :param length: The length of the sequence (integer).
        :return: A numpy array containing the Zadoff-Chu sequence.
        """
        if root is None:
            # Deduce the root as the smallest integer coprime with the length
            root = next(r for r in range(1, length) if np.gcd(r, length) == 1)

        if np.gcd(root, length) != 1:
            raise ValueError("Root and length must be coprime.")

        n = np.arange(length)
        zc_sequence = np.exp(-1j * np.pi * root * n * (n + 1) / length)
        return zc_sequence

    def _map_data_and_pilots(self, data_symbols):
        """Given a 1D array of data QAM symbols of length data_per_symbol, return a full-length
        array of length n_carriers with pilots inserted at pilot_indices and data filled in the
        remaining carriers in order.
        """

        grid = np.zeros(self.n_carriers, dtype=complex)
        grid[self.pilot_indices] = self.pilot_symbols
        grid[self.data_carriers] = data_symbols
        return grid

    def ofdm_modulate(self):
        if not hasattr(self, 'qam_symbols') or self.qam_symbols is None:
            raise ValueError("QAM symbols not generated yet.")

        # compute how many data symbols we expect depending on pilot mode
        if self.pilot_mode == "interleaved":
            data_per_symbol = self.n_carriers - self.n_pilots
            n_symbols = len(self.qam_symbols) // data_per_symbol
            if n_symbols == 0:
                raise ValueError("Not enough QAM symbols for a single OFDM symbol with the current pilot spacing.")
            self.n_ofdm_symbols = n_symbols
            qam_symbols = self.qam_symbols[:n_symbols * data_per_symbol]
            qam_matrix = qam_symbols.reshape((n_symbols, data_per_symbol))
        else:  # block mode: first OFDM symbol is pilots, remaining are full-data symbols
            data_per_symbol = self.n_carriers
            n_data_symbols = len(self.qam_symbols) // data_per_symbol
            if n_data_symbols == 0:
                raise ValueError("Not enough QAM symbols for data symbols in block-pilot mode.")
            # total OFDM symbols = 1 pilot + n_data_symbols
            self.n_ofdm_symbols = 1 + n_data_symbols
            qam_symbols = self.qam_symbols[:n_data_symbols * data_per_symbol]
            qam_matrix = qam_symbols.reshape((n_data_symbols, data_per_symbol))

        ofdm_time = []

        if self.pilot_mode == "interleaved":
            # comb mode: each OFDM symbol contains pilots interleaved in frequency
            for row in qam_matrix:
                grid = self._map_data_and_pilots(row)  # length = n_carriers
                freq_oversampled = self.pad_subcarriers(grid[np.newaxis, :])[0]
                time_domain = np.fft.ifft(freq_oversampled, self.fft_size)
                cp = time_domain[-self.cp_length:]
                ofdm_symbol = np.concatenate([cp, time_domain])
                ofdm_time.append(ofdm_symbol)
        else:
            # block mode:
            # 1) first symbol: pilots on every carrier
            #pilot_grid = np.ones(self.n_carriers, dtype=complex) * self.pilot_symbol
            pilot_grid = self.pilot_symbols

            #todo change here
            #pilot_grid = self.generate_zadoff_chu_sequence(None, self.n_carriers)
            freq_oversampled = self.pad_subcarriers(pilot_grid[np.newaxis, :])[0]
            time_domain = np.fft.ifft(freq_oversampled, self.fft_size)
            cp = time_domain[-self.cp_length:]
            ofdm_time.append(np.concatenate([cp, time_domain]))

            # 2) subsequent symbols: pure data on all carriers (no pilots)
            for row in qam_matrix:
                # here row length == n_carriers (full-data symbol)
                # build grid directly (data on every carrier)
                grid = np.array(row, dtype=complex)
                freq_oversampled = self.pad_subcarriers(grid[np.newaxis, :])[0]
                time_domain = np.fft.ifft(freq_oversampled, self.fft_size)
                cp = time_domain[-self.cp_length:]
                ofdm_time.append(np.concatenate([cp, time_domain]))

        self.ofdm_time = np.array(ofdm_time)
        # Rescale the waveform to conform with the transmit power.
        pavg = np.mean(np.abs(self.ofdm_time) ** 2)
        txp = (10 ** (self.tx_power / 10)) / 1000
        alpha = np.sqrt(txp / pavg)
        self.ofdm_time *= alpha
        self.pilot_symbols *= alpha

        return self.ofdm_time


    # def ofdm_modulate(self):
    #     if self.qam_symbols is None:
    #         raise ValueError("QAM symbols not generated yet.")
    #
    #     data_per_symbol = self.n_carriers - self.n_pilots
    #     n_symbols = len(self.qam_symbols) // data_per_symbol
    #     if n_symbols == 0:
    #         raise ValueError("Not enough QAM symbols for a single OFDM symbol with the current pilot spacing.")
    #     self.n_ofdm_symbols = n_symbols  # update in case it was different
    #     qam_symbols = self.qam_symbols[:n_symbols * data_per_symbol]
    #     qam_matrix = qam_symbols.reshape((n_symbols, data_per_symbol))
    #     ofdm_time = []
    #     for row in qam_matrix:
    #         freq_domain = np.zeros(self.fft_size, dtype=complex)
    #
    #         # Map QAM symbols to center of spectrum (baseband symmetric)
    #         half = self.n_carriers // 2
    #         freq_domain[:half] = row[:half]
    #         freq_domain[-half:] = row[half:]
    #
    #         # IFFT to get time-domain signal
    #         time_domain = np.fft.ifft(freq_domain, self.fft_size)
    #
    #         # Add cyclic prefix
    #         cp = time_domain[-self.cp_length:]
    #         ofdm_symbol = np.concatenate([cp, time_domain])
    #         ofdm_time.append(ofdm_symbol)
    #
    #     self.ofdm_time = np.array(ofdm_time)
    #     return self.ofdm_time  # shape: n_ofdm_symbols x ( n_carriers * oversampling + cp_length)

    # ------------------------
    # Channel / RX side helpers
    # ------------------------
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

    def channel_estimate_ls(self, rx_subcarriers):
        n_sym = rx_subcarriers.shape[0]
        H_est = np.zeros_like(rx_subcarriers, dtype=complex)

        if self.pilot_mode == "interleaved":
            pilot_idx = self.pilot_indices
            for si in range(n_sym):
                y = rx_subcarriers[si]
                y_p = y[pilot_idx]
                h_p = y_p / self.pilot_symbols
                if len(pilot_idx) == 1:
                    H_est[si, :] = h_p[0]
                    continue
                x = pilot_idx
                xp = np.arange(self.n_carriers)
                f_re = interp1d(x, np.real(h_p), kind='linear', bounds_error=False, fill_value='extrapolate')
                f_im = interp1d(x, np.imag(h_p), kind='linear', bounds_error=False, fill_value='extrapolate')
                H_est[si, :] = f_re(xp) + 1j * f_im(xp)
        else:  # block pilot mode
            # Expect first symbol (index 0) to be full pilots
            # estimate H from first symbol and reuse for all OFDM symbols
            y0 = rx_subcarriers[0, :]
            H0 = y0 / self.pilot_symbols
            for si in range(n_sym):
                H_est[si, :] = H0

        return H_est


    def equalize_one_tap(self, rx_subcarriers, H_est, eps=1e-12):
        """One-tap equalizer (frequency domain). Avoids dividing by near-zero by thresholding.


        rx_subcarriers: (n_sym, n_carriers) received frequency samples
        H_est: (n_sym, n_carriers) estimated channel
        Returns: equalized symbols on every carrier (pilots too)."""

        # avoid division by zero
        H_safe = np.where(np.abs(H_est) < eps, eps, H_est)
        return rx_subcarriers / H_safe

    def demap_data_from_grid(self, symbol_grid):
        """Return data carriers only. In block mode (preamble), first symbol is pilots,
        so data are taken from symbols 1..end (all carriers)."""
        arr = np.atleast_2d(symbol_grid)
        if self.pilot_mode == "interleaved":
            data = arr[:, self.data_carriers]
        else:  # block mode: first row is pilot -> take subsequent rows and all carriers
            if arr.shape[0] < 2:
                return np.zeros((0, self.n_carriers), dtype=complex)
            data = arr[1:, :]  # shape (n_data_symbols, n_carriers)
        return data

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

        # Divide by oversampling factor for power scaling.
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

    def plot_constellation(self, symbols_rx, title=None, symbols_tx=None):
        fig, ax = plt.subplots()
        ax.scatter(np.real(symbols_rx), np.imag(symbols_rx), label='Rx symbols')
        if symbols_tx is not None:
            ax.scatter(np.real(symbols_tx), np.imag(symbols_tx), label='Tx symbols')
        ax.set_title(title if title is not None else f"{self.qam_order}-QAM Constellation")
        ax.set_xlabel("In-phase (I)")
        ax.set_ylabel("Quadrature (Q)")
        ax.axis("equal")
        ax.legend()

        return fig

    def plot_iq_time(self, iq, title=None):
        ywf = self.ofdm_time_to_freq(iq)
        ywf = ywf.reshape(self.n_ofdm_symbols, -1)
        qam = self.ofdm_to_qam(ywf)
        self.plot_constellation(qam, title=title)

    def plot_psd(self, ofdm_time=None, nperseg=1024, title=None):
        if ofdm_time is None:
            ofdm_time = self.ofdm_time

        ofdm_signal = ofdm_time.flatten()
        f, Pxx = welch(ofdm_signal, fs=self.fs, nperseg=nperseg, return_onesided=False)
        # Center around zero and shift
        Pxx_db = 10 * np.log10(np.fft.fftshift(Pxx))
        f_shifted = np.fft.fftshift(f)  # + fc  # shift to RF
        f_shifted_GHz = f_shifted / 1e9

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(f_shifted_GHz, Pxx_db)
        ax.set_xlabel("Frequency (GHz)")
        ax.set_xlim((-10, 10))
        ax.set_ylabel("PSD (dB)")
        ax.set_title(title if title is not None else "OFDM PSD Plot")
        ax.grid(True)

        return fig