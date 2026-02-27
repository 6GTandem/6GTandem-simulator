import numpy as np
import logging
import os
import xarray as xr
import matplotlib.pyplot as plt

from wireless_channel.waveforms import Waveform

logger = logging.getLogger(__name__)


class Channel:
    """
    Defines the wireless channel
    parameters:
    - Nr_ue_antennas
    - Nr_ru_antnennas
    -
    """

    # todo: check if default nr subcarriers is same as sionna data
    # todo: only works when num symples = num subcarriers
    # todo: account for oversampling, and multiple ofdm symbols => slice X_time into OFDM symbols of length Nr_subcarriers before FFT. Then apply channel per OFDM symbol.

    def __init__(
        self,
        channelmodel: str = "sionna",
        Nr_subcarriers: int = 1024,
        Nr_ue_antennas: int = 1,
        Nr_ru_antennas: int = 1,
        Nr_rus: int = 5,
        Nr_stripes: int = 2,
        ue_idx: int | None = None,
    ):
        # todo load all this based on csi
        self.channelmodel = channelmodel
        self.Nr_stripes = Nr_stripes
        self.Nr_rus = Nr_rus
        self.Nr_ue_antennas = Nr_ue_antennas
        self.Nr_ru_antennas = Nr_ru_antennas
        self.Nr_subcarriers = Nr_subcarriers
        self.ue_idx = ue_idx

        if channelmodel == "subTHz-Rayleigh":
            self.csi = self.subTHz_Rayleigh()

    @property
    def csi(self):
        """Return the channel state information (CSI)."""
        return self._csi

    @csi.setter
    def csi(self, value):
        """Set the CSI (Channel State Information) which is an xarray"""
        self._csi = value

    def get_csi(self, stripe_idx: int, ru_idx: int):
        """
        Return the CSI for a specific stripe and RU as a numpy array.

        Parameters
        ----------
        stripe_idx : int
            Index of the stripe
        ru_idx : int
            Index of the RU

        Returns
        -------
        xr.DataArray
            CSI with dimensions (rx_ant, tx_ant, subcarrier)
        """
        if not hasattr(self, "_csi"):
            raise AttributeError("CSI has not been set yet.")

        match = (self.csi["stripe_idx"] == stripe_idx) & (self.csi["RU_idx"] == ru_idx)
        if match.any():
            tx_index = match.argmax().item()  # first match
            channel = self.csi["channel"].isel(tx_pair=tx_index).values
            return channel

        raise ValueError(f"Stripe/RU combination does not exist: Stripe: {stripe_idx}, RU: {ru_idx}")

    @classmethod
    def from_sionna(cls, ue_coordinates: dict, sim_env: str, debug: bool = False):
        """Load the channel state information (CSI) for a specific UE.

        Parameters
        ----------
        ue_cooridnates : dict
            Coordinates of the UE from which to load the CSI. Dictionary containing the coordinates under the
            x, y and z keys.
        sim_env : str
            Simulation environment to use.
        debug : bool
            Use a channel only containing 1s for debugging when True.

        Returns
        -------
        Channel
            `Channel` object for the specified UE.
        """
        # todo load locations metadata => ue_idx
        dir_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = os.path.join(dir_path, "sionna_dataset", sim_env)
        ue_ds = xr.load_dataset(os.path.join(dir_path, "ue_locations", "ue_locations.nc"))

        # based on coordinates get UE idx
        x, y, z = ue_coordinates["x"], ue_coordinates["y"], ue_coordinates["z"]
        matched_user = ue_ds.where((ue_ds["x"] == x) & (ue_ds["y"] == y) & (ue_ds["z"] == z), drop=True)
        ue_idx = int(matched_user["user_id"].values.item())

        # based on UE idx load CSI
        csi_file = f"channels_thz_ue_{ue_idx}.nc"
        ds_sub_thz = xr.load_dataset(os.path.join(dir_path, "sub_thz_channels", csi_file))

        # extract needed params for Channel class
        channelmodel = "sionna"
        Nr_subcarriers = ds_sub_thz.sizes["subcarrier"]
        Nr_ue_antennas = ds_sub_thz.sizes["rx_ant"]
        Nr_ru_antennas = ds_sub_thz.sizes["tx_ant"]
        Nr_rus = ds_sub_thz["RU_idx"].max().item() + 1
        Nr_stripes = ds_sub_thz["stripe_idx"].max().item() + 1

        channel = cls(channelmodel, Nr_subcarriers, Nr_ue_antennas, Nr_ru_antennas, Nr_rus, Nr_stripes, ue_idx)

        csi_channel = ds_sub_thz["channel"]
        if csi_channel.dtype.fields is not None and "r" in csi_channel.dtype.fields and "i" in csi_channel.dtype.fields:
            # Convert structured array to complex
            csi_complex = csi_channel.values["r"] + 1j * csi_channel.values["i"]
            # Put back as a DataArray, preserving dims and coords
            ds_sub_thz["channel"] = xr.DataArray(
                csi_complex,
                dims=csi_channel.dims,
                coords=csi_channel.coords,
                name=csi_channel.name,
                attrs=csi_channel.attrs,
            )

        channel.csi = ds_sub_thz

        if debug:
            # all ones channel
            channel.csi["channel"] = xr.ones_like(channel.csi["channel"])

            # rayleigh channel (frequency uncorrelated)
            # channel.csi["channel"] = xr.ones_like(channel.csi["channel"]) * channel.subTHz_Rayleigh()

            # rayleigh channel (frequency correlated)
            # channel.csi["channel"] = xr.ones_like(channel.csi["channel"]) * channel.correlated_freq_channel()

        return channel

    def correlated_freq_channel(
        self, n_taps=8, max_delay_samples=None, power_profile="exponential", normalize=True, seed=None
    ):
        """
        Create correlated frequency-selective Rayleigh channels.

        Returns H with shape (n_links, n_rx, n_tx, n_subcarriers), complex.

        Parameters
        ----------
        n_links : int
            Number of (stripe x RU) links (your previous self.Nr_stripes*self.Nr_rus).
        n_rx : int
            Number of UE antennas.
        n_tx : int
            Number of RU antennas.
        n_subcarriers : int
            Number of OFDM subcarriers (your self.n_carriers).
        n_taps : int
            Number of time-domain multipath taps (L). Typical small values: 4..16.
        max_delay_samples : int or None
            Maximum delay (in samples) for taps. If None, use n_taps (dense).
            Larger max_delay_samples -> larger delay spread -> more frequency selectivity.
        power_profile : {'exponential', 'uniform'}
            Power delay profile shape.
        normalize : bool
            If True, normalize the per-link average power to 1.
        seed : int or None
            RNG seed for reproducibility.
        """
        rng = np.random.RandomState(seed)

        if max_delay_samples is None:
            max_delay_samples = n_taps

        # build PDP (power for each tap index 0..n_taps-1)
        if power_profile == "exponential":
            # exponential decay across taps; shape parameter controls spread
            delays = np.arange(n_taps)
            # choose a decay constant so later taps still contribute; tweak as needed
            decay = 1.0  # larger -> faster decay (less delay spread)
            pdp = np.exp(-decay * delays)
        elif power_profile == "uniform":
            pdp = np.ones(n_taps)
        else:
            raise ValueError("unknown power_profile")

        pdp = pdp / np.sum(pdp)  # normalize so total power = 1

        # allocate output
        n_links = self.Nr_stripes * self.Nr_rus
        n_rx = self.Nr_ue_antennas
        n_tx = self.Nr_ru_antennas
        n_subcarriers = self.Nr_subcarriers
        H_freq = np.zeros((n_links, n_rx, n_tx, n_subcarriers), dtype=complex)

        # for each link / rx / tx generate time-domain taps and FFT to get freq response
        for link in range(n_links):
            for rx in range(n_rx):
                for tx in range(n_tx):
                    # place the taps at random small delays within max_delay_samples
                    # simple model: taps at integer delays 0..(n_taps-1) (you can randomize if desired)
                    # Create an impulse response of length `n_subcarriers` (zero-padded)
                    h_time = np.zeros(n_subcarriers, dtype=complex)

                    # generate complex Gaussian taps scaled by PDP sqrt
                    # (independent Rayleigh fading per tap)
                    tap_amps = (rng.normal(size=n_taps) + 1j * rng.normal(size=n_taps)) / np.sqrt(2.0)
                    tap_amps *= np.sqrt(pdp)  # scale by sqrt(power per tap)

                    # optionally randomize tap positions within the first max_delay_samples
                    # Here we place taps at 0,1,2,... by default, but you may choose random delays:
                    # delays_idx = rng.randint(0, max_delay_samples, size=n_taps)
                    delays_idx = np.arange(n_taps)  # simple, contiguous delays

                    # build time-domain CIR
                    for amp, d in zip(tap_amps, delays_idx):
                        if d < n_subcarriers:
                            h_time[d] += amp

                    # compute frequency response (n_subcarriers)
                    Hf = np.fft.fft(h_time, n_subcarriers)

                    if normalize:
                        # normalize so average power across subcarriers = 1
                        Hf = Hf / np.sqrt(np.mean(np.abs(Hf) ** 2) + 1e-16)

                    H_freq[link, rx, tx, :] = Hf

        return H_freq

    def subTHz_Rayleigh(self, p: int = 1):
        """
        CSI dimension: Nr_stripes x Nr_rus x Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers
        """
        variance = p / 2
        stdev = np.sqrt(variance)
        H = np.random.normal(
            0, stdev, (self.Nr_stripes * self.Nr_rus, self.Nr_ue_antennas, self.Nr_ru_antennas, self.Nr_subcarriers)
        ) + 1j * np.random.normal(
            0, stdev, (self.Nr_stripes * self.Nr_rus, self.Nr_ue_antennas, self.Nr_ru_antennas, self.Nr_subcarriers)
        )
        return H

    def transmit_ul(self, X, wf: Waveform):
        """
        transmit from UE to all RUs
        Y_ul = H^T X_f
        shapes:
        Y_ul: Nr_stripes x Nr_rus x Nr_ru_antennas x Nr_subcarriers
        H: Nr_stripes x Nr_rus x Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers
        X: Nr_stripes x Nr_ue_antennas x Nr_samples
        X_f: fft of the time domain samples X
        """
        Y_f = np.zeros((wf.n_ofdm_symbols, self.Nr_ue_antennas, self.Nr_subcarriers), dtype=complex)
        Y_f_padded = np.zeros((wf.n_ofdm_symbols, self.Nr_ru_antennas, wf.fft_size), dtype=complex)
        Y_time = np.zeros(
            (self.Nr_stripes, self.Nr_rus, self.Nr_ru_antennas, wf.n_ofdm_symbols, wf.fft_size + wf.cp_length),
            dtype=complex,
        )

        assert self.Nr_subcarriers == wf.n_carriers, "Nr_subcarriers in Channel must match n_carriers in Waveform."

        for s in range(self.Nr_stripes):
            for r in range(self.Nr_rus):
                X_f = np.zeros((self.Nr_ue_antennas, wf.n_ofdm_symbols, wf.n_carriers), dtype=complex)
                # Convert time-domain signal to frequency domain.
                for m in range(self.Nr_ue_antennas):
                    X_f[m, :, :] = wf.extract_subcarriers(wf.ofdm_time_to_freq(X[m, :, :]))

                # select channel
                H = self.get_csi(s, r)  # [Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers]
                H = np.transpose(H, (1, 0, 2))
                logger.debug(f"channel shape: {H.shape}")

                # Apply the channel.
                for sym in range(wf.n_ofdm_symbols):
                    for k in range(self.Nr_subcarriers):
                        Y_f[sym, :, k] = H[:, :, k] @ X_f[:, sym, k]

                # pad zeros
                for m in range(self.Nr_ru_antennas):
                    Y_f_padded[:, m, :] = wf.pad_subcarriers(
                        Y_f[:, m, :]
                    )  # (n_ofdm_symbols x n_carriers) => (n_ofdm_symbosl x fftsize)

                    # IFFT back to time domain (should error => fft size + cp length
                    Y_time[s, r, m, :, :] = wf.ofdm_freq_to_time(Y_f_padded[:, m, :])

        return Y_time

    def transmit_ul_id(self, X, stripe_idx, ru_idx, wf: Waveform):
        """
        transmit from UE to all RUs
        Y_ul = H^T X_f
        shapes:
        Y_ul: Nr_stripes x Nr_rus x Nr_ru_antennas x Nr_subcarriers
        H: Nr_stripes x Nr_rus x Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers
        X: Nr_stripes x Nr_ue_antennas x Nr_samples
        X_f: fft of the time domain samples X
        """
        Y_f = np.zeros((wf.n_ofdm_symbols, self.Nr_ue_antennas, self.Nr_subcarriers), dtype=complex)
        Y_f_padded = np.zeros((wf.n_ofdm_symbols, self.Nr_ru_antennas, wf.fft_size), dtype=complex)
        Y_time = np.zeros(
            (self.Nr_ru_antennas, wf.n_ofdm_symbols, wf.fft_size + wf.cp_length),
            dtype=complex,
        )

        assert self.Nr_subcarriers == wf.n_carriers, "Nr_subcarriers in Channel must match n_carriers in Waveform."

        X_f = np.zeros((self.Nr_ue_antennas, wf.n_ofdm_symbols, wf.n_carriers), dtype=complex)
        # Convert time-domain signal to frequency domain.
        for m in range(self.Nr_ue_antennas):
            X_f[m, :, :] = wf.extract_subcarriers(wf.ofdm_time_to_freq(X[m, :, :]))

        # select channel
        H = self.get_csi(stripe_idx, ru_idx)  # [Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers]
        H = np.transpose(H, (1, 0, 2))
        logger.debug(f"channel shape: {H.shape}")

        # Apply the channel.
        for sym in range(wf.n_ofdm_symbols):
            for k in range(self.Nr_subcarriers):
                Y_f[sym, :, k] = H[:, :, k] @ X_f[:, sym, k]

        # pad zeros
        for m in range(self.Nr_ru_antennas):
            Y_f_padded[:, m, :] = wf.pad_subcarriers(
                Y_f[:, m, :]
            )  # (n_ofdm_symbols x n_carriers) => (n_ofdm_symbosl x fftsize)

            # IFFT back to time domain (should error => fft size + cp length
            Y_time[m, :, :] = wf.ofdm_freq_to_time(Y_f_padded[:, m, :])

        return Y_time

    def transmit_dl(self, X_list: list[np.ndarray], active_ru_idxes: list[int], waveform: Waveform):

        # todo debug and see if makes sense!!!
        """ "
        transmit per stripe from a RU to UE

         X_list (list of np.ndarray):
            List of transmit signals, one per stripe.
            Each element is either:
                - None  (if no RU in that stripe transmits)
                - np.ndarray of shape [Nr_ru_antennas x Nr_ofdm_symbols x fft_size+cp_length] for the chosen RU
         Returns:
                np.ndarray: Received signal at UE of shape
                            [Nr_ue_antennas x Nr_subcarriers]
        """

        Y_f = np.zeros((waveform.n_ofdm_symbols, self.Nr_ue_antennas, self.Nr_subcarriers), dtype=complex)

        for stripe_idx, active_ru_idx in enumerate(active_ru_idxes):
            if active_ru_idx is None or X_list[stripe_idx] is None:
                continue  # no transmission from this stripe

            # Time-domain transmit signal for this RU
            X_time = X_list[
                stripe_idx
            ]  # shape: nr_RU_antennas x n_ofdm_symbols x ( n_carriers * oversampling + cp_length)

            assert (
                self.Nr_subcarriers == waveform.n_carriers
            ), "Nr_subcarriers in Channel must match n_carriers in Waveform"

            # FFT to frequency domain
            X_f = np.zeros((self.Nr_ru_antennas, waveform.n_ofdm_symbols, waveform.n_carriers), dtype=complex)
            for m in range(self.Nr_ru_antennas):
                # n_ofdm_symbols x (n_carriers)
                X_f[m, :, :] = waveform.extract_subcarriers(waveform.ofdm_time_to_freq(X_time[m, :, :]))

            # select channel
            H = self.get_csi(stripe_idx, active_ru_idx)  # [Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers]
            # H = np.ones((4, 4, 1024))  # debug with all ones channel
            logger.debug(f"channel shape: {H.shape}")

            # plot channel
            plt.stem(np.abs(H[0, 0, :]) ** 2)
            plt.xlabel("SUBCARRIER INDEX")
            plt.ylabel("CHANNEL GAIN |H|^2")
            plt.show()

            # apply the channel
            for sym in range(waveform.n_ofdm_symbols):
                for k in range(self.Nr_subcarriers):
                    # X_f[:, sym, k] shape: [Nr_ru_antennas]
                    # H[:, :, k] shape: [Nr_ue_antennas x Nr_ru_antennas]
                    Y_f[sym, :, k] += H[:, :, k] @ X_f[:, sym, k]  # sum because signals come from multiple stripes

        # pad zeros
        Y_f_padded = np.zeros((waveform.n_ofdm_symbols, self.Nr_ue_antennas, waveform.fft_size), dtype=complex)
        for m in range(self.Nr_ue_antennas):
            Y_f_padded[:, m, :] = waveform.pad_subcarriers(
                Y_f[:, m, :]
            )  # (n_ofdm_symbols x n_carriers) => (n_ofdm_symbosl x fftsize)

        # IFFT back to time domain (should error => fft size + cp length
        Y_time = np.zeros(
            (self.Nr_ue_antennas, waveform.n_ofdm_symbols, waveform.fft_size + waveform.cp_length), dtype=complex
        )
        for m in range(self.Nr_ue_antennas):
            Y_time[m, :, :] = waveform.ofdm_freq_to_time(Y_f_padded[:, m, :])

        # # move back to time domain
        # Y_time = waveform.ofdm_freq_to_time(Y_f) # todo fix

        return Y_time

    def __str__(self):
        """
        Return a human-readable string summary of the Channel object.
        This is called when you do `print(channel)`.
        """
        info = (
            f"Channel model          : {self.channelmodel}\n"
            f"Number of stripes      : {self.Nr_stripes}\n"
            f"Number of RUs          : {self.Nr_rus}\n"
            f"Number of UE antennas  : {self.Nr_ue_antennas}\n"
            f"Number of RU antennas  : {self.Nr_ru_antennas}\n"
            f"Number of subcarriers  : {self.Nr_subcarriers}\n"
        )
        if hasattr(self, "csi"):
            info += f"CSI shape              : {self.csi.coords}\n"
        return info
