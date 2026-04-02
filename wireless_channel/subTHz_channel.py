import numpy as np
import logging
import os
import xarray as xr
import matplotlib.pyplot as plt
from collections.abc import Sequence

from wireless_channel.waveforms import Waveform

logger = logging.getLogger(__name__)
C = 299792458.0


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
        print(f"ue_ds: {ue_ds.where((ue_ds['user_id'] == 2), drop=True)}")
        print(f' matched user: {matched_user}')
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
            #channel.csi["channel"] = xr.ones_like(channel.csi["channel"])

            # rayleigh channel (frequency uncorrelated)
            # channel.csi["channel"] = xr.ones_like(channel.csi["channel"]) * channel.subTHz_Rayleigh()

            # rayleigh channel (frequency correlated)
            channel.csi["channel"] = xr.ones_like(channel.csi["channel"]) * channel.correlated_freq_channel()

        return channel

    @staticmethod
    def _extract_ru_positions(stripe_positions: list[list[dict]]) -> list[tuple[int, int, dict]]:
        """Extract RU positions and preserve stripe/RU indices.

        Expected stripe_positions format is a list of stripes where each stripe is a list
        containing entries like {"radio_unit": {"x": ..., "y": ..., "z": ...}}.
        """
        ru_positions = []
        for stripe_idx, stripe in enumerate(stripe_positions):
            ru_idx = 0
            for entry in stripe:
                ru = entry.get("radio_unit")
                if ru is None:
                    continue
                ru_positions.append((stripe_idx, ru_idx, ru))
                ru_idx += 1
        return ru_positions

    @staticmethod
    def _build_csi_dataset(
        channel_values: np.ndarray,
        stripe_ids: np.ndarray,
        ru_ids: np.ndarray,
    ) -> xr.Dataset:
        """Build an xarray dataset compatible with get_csi/transmit methods."""
        return xr.Dataset(
            data_vars={
                "channel": (("tx_pair", "rx_ant", "tx_ant", "subcarrier"), channel_values),
            },
            coords={
                "tx_pair": np.arange(channel_values.shape[0]),
                "rx_ant": np.arange(channel_values.shape[1]),
                "tx_ant": np.arange(channel_values.shape[2]),
                "subcarrier": np.arange(channel_values.shape[3]),
                "stripe_idx": ("tx_pair", stripe_ids.astype(int)),
                "RU_idx": ("tx_pair", ru_ids.astype(int)),
            },
        )

    @staticmethod
    def _to_linear_gains(gain_dbi, n_antennas: int, label: str) -> np.ndarray:
        """Convert antenna gain(s) in dBi to linear power gain per antenna."""
        if np.isscalar(gain_dbi):
            gains_dbi = np.full(n_antennas, float(gain_dbi), dtype=float)
        elif isinstance(gain_dbi, Sequence):
            gains_dbi = np.asarray(gain_dbi, dtype=float)
            if gains_dbi.ndim != 1 or gains_dbi.size != n_antennas:
                raise ValueError(
                    f"{label} must be a scalar or a sequence with length {n_antennas}."
                )
        else:
            raise ValueError(
                f"{label} must be a scalar or a sequence with length {n_antennas}."
            )

        return 10.0 ** (gains_dbi / 10.0)

    @staticmethod
    def _wrap_to_pi(phi_rad: float) -> float:
        return (phi_rad + np.pi) % (2.0 * np.pi) - np.pi

    @staticmethod
    def _tr38901_power_gain(theta_rad: float, phi_rad: float) -> float:
        """TR 38.901 antenna power gain pattern (linear)."""
        phi = Channel._wrap_to_pi(phi_rad)

        theta_3db = np.deg2rad(65.0)
        phi_3db = np.deg2rad(65.0)
        a_max = 30.0
        sla_v = 30.0
        g_e_max = 20.0

        a_v = -min(12.0 * ((theta_rad - np.pi / 2.0) / theta_3db) ** 2, sla_v)
        a_h = -min(12.0 * (phi / phi_3db) ** 2, a_max)
        a_db = -min(-(a_v + a_h), a_max) + g_e_max
        return float(10.0 ** (a_db / 10.0))

    @staticmethod
    def _unit_vector_from_az_el(az_deg: float, el_deg: float) -> np.ndarray:
        az = np.deg2rad(az_deg)
        el = np.deg2rad(el_deg)
        vec = np.array(
            [
                np.cos(el) * np.cos(az),
                np.cos(el) * np.sin(az),
                np.sin(el),
            ],
            dtype=float,
        )
        return vec / (np.linalg.norm(vec) + 1e-16)

    @staticmethod
    def _local_frame_from_boresight(boresight_unit: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build a right-handed local frame with local +x aligned to boresight."""
        x_hat = boresight_unit
        ref = np.array([0.0, 0.0, 1.0], dtype=float)
        if abs(float(np.dot(x_hat, ref))) > 0.95:
            ref = np.array([0.0, 1.0, 0.0], dtype=float)

        y_hat = np.cross(ref, x_hat)
        y_hat = y_hat / (np.linalg.norm(y_hat) + 1e-16)
        z_hat = np.cross(x_hat, y_hat)
        z_hat = z_hat / (np.linalg.norm(z_hat) + 1e-16)
        return x_hat, y_hat, z_hat

    @staticmethod
    def _vector_to_local_spherical(
        direction_global: np.ndarray,
        frame: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> tuple[float, float]:
        """Return local spherical (theta, phi). theta from +z and phi from +x in local frame."""
        x_hat, y_hat, z_hat = frame
        x_local = float(np.dot(direction_global, x_hat))
        y_local = float(np.dot(direction_global, y_hat))
        z_local = float(np.dot(direction_global, z_hat))

        theta = float(np.arccos(np.clip(z_local, -1.0, 1.0)))
        phi = float(np.arctan2(y_local, x_local))
        return theta, phi

    @staticmethod
    def _resolve_los_pattern(component_config: dict | None, antenna_pattern: str | None) -> str:
        if antenna_pattern is not None and str(antenna_pattern).strip() != "":
            pattern = str(antenna_pattern).strip().lower()
        else:
            cfg_pattern = ""
            if component_config is not None:
                cfg_pattern = str(component_config.get("antenna", {}).get("pattern", "")).strip().lower()
            pattern = cfg_pattern if cfg_pattern else "isotropic"

        normalized = pattern.replace("-", "").replace("_", "").replace(".", "")
        if normalized in ("isotropic", "iso"):
            return "isotropic"
        if normalized in ("tr38901", "tr38901pattern"):
            return "tr38901"
        raise ValueError(f"Unsupported LOS antenna pattern '{pattern}'. Supported: isotropic, tr38901.")

    @classmethod
    def _compute_link_pattern_field_gain(
        cls,
        pattern: str,
        ue_xyz: np.ndarray,
        ru_xyz: np.ndarray,
        component_config: dict | None,
    ) -> float:
        if pattern == "isotropic":
            return 1.0

        if pattern != "tr38901":
            raise ValueError(f"Unsupported LOS antenna pattern '{pattern}'.")

        antenna_cfg = (component_config or {}).get("antenna", {})
        ue_az = float(antenna_cfg.get("ue_boresight_az_deg", 0.0))
        ue_el = float(antenna_cfg.get("ue_boresight_el_deg", 0.0))
        ru_az = float(antenna_cfg.get("ru_boresight_az_deg", 180.0))
        ru_el = float(antenna_cfg.get("ru_boresight_el_deg", 0.0))

        ue_frame = cls._local_frame_from_boresight(cls._unit_vector_from_az_el(ue_az, ue_el))
        ru_frame = cls._local_frame_from_boresight(cls._unit_vector_from_az_el(ru_az, ru_el))

        direction_ue_to_ru = ru_xyz - ue_xyz
        direction_ue_to_ru = direction_ue_to_ru / (np.linalg.norm(direction_ue_to_ru) + 1e-16)
        direction_ru_to_ue = -direction_ue_to_ru

        tx_theta, tx_phi = cls._vector_to_local_spherical(direction_ue_to_ru, ue_frame)
        rx_theta, rx_phi = cls._vector_to_local_spherical(direction_ru_to_ue, ru_frame)

        g_tx_lin = cls._tr38901_power_gain(tx_theta, tx_phi)
        g_rx_lin = cls._tr38901_power_gain(rx_theta, rx_phi)
        return float(np.sqrt(g_tx_lin * g_rx_lin))

    @classmethod
    def from_los(
        cls,
        ue_coordinates: dict,
        stripe_positions: list[list[dict]],
        n_subcarriers: int,
        carrier_frequency_hz: float,
        subcarrier_spacing_hz: float,
        Nr_ue_antennas: int = 1,
        Nr_ru_antennas: int = 1,
        ue_antenna_gain_dbi: float | Sequence[float] = 0.0,
        ru_antenna_gain_dbi: float | Sequence[float] = 0.0,
        component_config: dict | None = None,
        antenna_pattern: str | None = None,
        normalize_gain: bool = False,
        min_distance_m: float = 1e-3,
    ):
        """Build a deterministic LOS channel with OFDM phase slope across subcarriers."""
        ru_positions = cls._extract_ru_positions(stripe_positions)
        if len(ru_positions) == 0:
            raise ValueError("No radio_unit entries found in stripe_positions.")

        stripe_ids = np.array([s for s, _, _ in ru_positions], dtype=int)
        ru_ids = np.array([r for _, r, _ in ru_positions], dtype=int)
        n_links = len(ru_positions)

        # Match modem subcarrier ordering used by extract_subcarriers/pad_subcarriers:
        # [0, 1, ..., N/2-1, -N/2, ..., -1].
        sampling_rate_hz = n_subcarriers * subcarrier_spacing_hz
        f_sub = np.fft.fftfreq(n_subcarriers, d=1.0 / sampling_rate_hz)

        ue_xyz = np.array([ue_coordinates["x"], ue_coordinates["y"], ue_coordinates["z"]], dtype=float)
        channel_values = np.zeros((n_links, Nr_ue_antennas, Nr_ru_antennas, n_subcarriers), dtype=complex)
        distances = []

        # Friis in complex baseband: |h| = lambda/(4*pi*d)*sqrt(G_rx*G_tx)
        # where G_rx and G_tx are linear power gains.
        ue_gains_lin = cls._to_linear_gains(ue_antenna_gain_dbi, Nr_ue_antennas, "ue_antenna_gain_dbi")
        ru_gains_lin = cls._to_linear_gains(ru_antenna_gain_dbi, Nr_ru_antennas, "ru_antenna_gain_dbi")
        antenna_field_gain = np.sqrt(np.outer(ue_gains_lin, ru_gains_lin))
        pattern_name = cls._resolve_los_pattern(component_config, antenna_pattern)

        wavelength = C / carrier_frequency_hz
        for link_idx, (_, _, ru_pos) in enumerate(ru_positions):
            ru_xyz = np.array([ru_pos["x"], ru_pos["y"], ru_pos["z"]], dtype=float)
            distance = max(np.linalg.norm(ue_xyz - ru_xyz), min_distance_m)
            distances.append(distance)

            pattern_field_gain = cls._compute_link_pattern_field_gain(
                pattern_name,
                ue_xyz,
                ru_xyz,
                component_config,
            )

            tau = distance / C
            gain = wavelength / (4.0 * np.pi * distance)
            phase = np.exp(-1j * 2.0 * np.pi * (carrier_frequency_hz + f_sub) * tau)
            h = gain * pattern_field_gain * phase

            #print(f"Link {link_idx}: distance={distance:.3f}m, gain={gain:.3e}, pattern_gain={pattern_field_gain:.3f}")

            for rx_ant in range(Nr_ue_antennas):
                for tx_ant in range(Nr_ru_antennas):
                    channel_values[link_idx, rx_ant, tx_ant, :] = h * antenna_field_gain[rx_ant, tx_ant]

        if normalize_gain:
            power = np.mean(np.abs(channel_values) ** 2)
            channel_values = channel_values / np.sqrt(power + 1e-16)

        Nr_stripes = stripe_ids.max().item() + 1
        Nr_rus = ru_ids.max().item() + 1
        channel = cls(
            channelmodel="los",
            Nr_subcarriers=n_subcarriers,
            Nr_ue_antennas=Nr_ue_antennas,
            Nr_ru_antennas=Nr_ru_antennas,
            Nr_rus=Nr_rus,
            Nr_stripes=Nr_stripes,
            ue_idx=None,
        )
        channel.csi = cls._build_csi_dataset(channel_values, stripe_ids, ru_ids)

        logger.info(
            "LOS channel created: pattern=%s, links=%d, d_min=%.3fm, d_max=%.3fm, fc=%.3fGHz, df=%.3fMHz, G_ue=%.2fdBi, G_ru=%.2fdBi",
            pattern_name,
            n_links,
            min(distances),
            max(distances),
            carrier_frequency_hz / 1e9,
            subcarrier_spacing_hz / 1e6,
            float(np.mean(10.0 * np.log10(ue_gains_lin))),
            float(np.mean(10.0 * np.log10(ru_gains_lin))),
        )
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


def build_channel(
    channel_model: str,
    ue_coordinates: dict,
    *,
    sim_env: str | None = None,
    component_config: dict | None = None,
    stripe_positions: list[list[dict]] | None = None,
    waveform: Waveform | None = None,
    Nr_ue_antennas: int = 1,
    Nr_ru_antennas: int = 1,
    debug: bool = False,
    los_normalize_gain: bool = False,
    los_ue_antenna_gain_dbi: float | Sequence[float] = 0.0,
    los_ru_antenna_gain_dbi: float | Sequence[float] = 0.0,
    los_antenna_pattern: str | None = None,
) -> Channel:
    """Factory helper for selecting channel model at runtime."""
    model = channel_model.lower()

    if model == "sionna":
        if sim_env is None:
            raise ValueError("sim_env is required when channel_model='sionna'.")
        return Channel.from_sionna(ue_coordinates, sim_env, debug=debug)

    if model == "los":
        if stripe_positions is None:
            raise ValueError("stripe_positions is required when channel_model='los'.")
        if waveform is None:
            raise ValueError("waveform is required when channel_model='los'.")
        subcarrier_spacing_hz = waveform.bw / waveform.n_carriers
        return Channel.from_los(
            ue_coordinates=ue_coordinates,
            stripe_positions=stripe_positions,
            n_subcarriers=waveform.n_carriers,
            carrier_frequency_hz=waveform.fc,
            subcarrier_spacing_hz=subcarrier_spacing_hz,
            Nr_ue_antennas=Nr_ue_antennas,
            Nr_ru_antennas=Nr_ru_antennas,
            ue_antenna_gain_dbi=los_ue_antenna_gain_dbi,
            ru_antenna_gain_dbi=los_ru_antenna_gain_dbi,
            component_config=component_config,
            antenna_pattern=los_antenna_pattern,
            normalize_gain=los_normalize_gain,
        )

    raise ValueError(f"Unknown channel model '{channel_model}'. Supported: sionna, los.")
