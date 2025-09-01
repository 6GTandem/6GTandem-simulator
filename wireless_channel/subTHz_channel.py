import numpy as np
import os
import xarray as xr


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

    def __init__(self, channelmodel: str = 'sionna',
                 Nr_subcarriers: int = 1024, Nr_ue_antennas: int = 1, Nr_ru_antennas: int =1,
                 Nr_rus: int=5, Nr_stripes: int=2):
        # todo load all this based on csi
        self.channelmodel = channelmodel
        self.Nr_stripes = Nr_stripes
        self.Nr_rus = Nr_rus
        self.Nr_ue_antennas = Nr_ue_antennas
        self.Nr_ru_antennas = Nr_ru_antennas
        self.Nr_subcarriers = Nr_subcarriers

        if channelmodel == 'subTHz-Rayleigh':
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

        match = ((self.csi["stripe_idx"] == stripe_idx) & (self.csi["RU_idx"] == ru_idx))
        if match.any():
            tx_index = match.argmax().item()  # first match
            channel = self.csi["channel"].isel(tx_pair=tx_index).values
            return channel
        else:
            print("No matching stripe/RU combination found.")
            return None

    @classmethod
    def from_sionna(cls, ue_coordinates):
        """
        Construct a Waveform from a configuration dictionary.

        :param waveform_config: dictionary containing waveform parameters
        :param freq_band_config: dictionary containing frequency band parameters (fc, bw, num_carriers)
        :return: Waveform instance
        """
        # todo load locations metadata => ue_idx
        config_file = "ue_locations_5681.nc"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        config_path = os.path.join(dir_path, "..", "configurations")
        ue_ds = xr.load_dataset(os.path.join(config_path, config_file))

        # based on coordinates get UE idx
        x, y, z = ue_coordinates['x'], ue_coordinates['y'], ue_coordinates['z']
        matched_user = ue_ds.where(
            (ue_ds['x'] == x) &
            (ue_ds['y'] == y) &
            (ue_ds['z'] == z),
            drop=True
        )
        ue_idx = int(matched_user['user_id'].values.item())

        # based on UE idx load CSI
        csi_file = f'channels_thz_ue_{ue_idx}.nc'
        ds_sub_thz = xr.load_dataset(os.path.join(dir_path, "sionna_dataset", "sub_thz_channels", csi_file))

        # extract needed params for Channel class
        channelmodel = 'sionna'
        Nr_subcarriers = ds_sub_thz.sizes['subcarrier']
        Nr_ue_antennas = ds_sub_thz.sizes['rx_ant']
        Nr_ru_antennas = ds_sub_thz.sizes['tx_ant']
        Nr_rus = ds_sub_thz["RU_idx"].max().item() + 1
        Nr_stripes = ds_sub_thz["stripe_idx"].max().item() + 1

        channel = Channel(channelmodel, Nr_subcarriers, Nr_ue_antennas, Nr_ru_antennas, Nr_rus, Nr_stripes)
        channel.csi = ds_sub_thz

        return channel

    def subTHz_Rayleigh(self, p: int = 1):
        """
        CSI dimension: Nr_stripes x Nr_rus x Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers
        """
        variance = p / 2
        stdev = np.sqrt(variance)
        H = (np.random.normal(0, stdev, (self.Nr_stripes, self.Nr_rus, self.Nr_ue_antennas,
                                         self.Nr_ru_antennas, self.Nr_subcarriers)) +
             1j * np.random.normal(0, stdev, (self.Nr_stripes, self.Nr_rus, self.Nr_ue_antennas,
                                              self.Nr_ru_antennas, self.Nr_subcarriers)))
        return H

    def transmit_ul(self, X):
        """
        transmit from UE to all RUs
        Y_ul = H^T X_f
        shapes:
        Y_ul: Nr_stripes x Nr_rus x Nr_ru_antennas x Nr_subcarriers
        H: Nr_stripes x Nr_rus x Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers
        X: Nr_stripes x Nr_ue_antennas x Nr_samples
        X_f: fft of the time domain samples X
        """
        Nr_samples = X.shape[-1]
        if Nr_samples < self.Nr_subcarriers:
            raise ValueError("Nr_samples must be >= Nr_subcarriers for FFT processing")

        # Convert time-domain signal to frequency domain
        X_f = np.fft.fft(X, n=self.Nr_subcarriers, axis=-1)  # shape: [Nr_stripes x Nr_ue_antennas x Nr_subcarriers]

        Y_ul = np.zeros((self.Nr_stripes, self.Nr_rus, self.Nr_ru_antennas, self.Nr_subcarriers), dtype=complex)

        # Loop over Stripes, RUs and subcarriers
        for s in range(self.Nr_stripes):
            for r in range(self.Nr_rus): # loop over RUs
                H_r = self.csi[s, r, :, :, :]  # channel at RU r: shape: [Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers]
                for k in range(self.Nr_subcarriers): # loop over carriers
                    Y_ul[s, r, :, k] = H_r[s, :, :, k].T @ X_f[s, :, k]


        # move back to time domain
        Y_time = np.fft.ifft(Y_ul, n=self.Nr_subcarriers, axis=-1)

        return Y_time


    def transmit_dl(self, X_list, active_ru_idxes, waveform):

        #todo debug and see if makes sens!!!
        """"
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
            X_time = X_list[stripe_idx]  # shape: nr_RU_antennas x n_ofdm_symbols x ( n_carriers * oversampling + cp_length)

            # FFT to frequency domain
            X_f = np.zeros((self.Nr_ru_antennas, waveform.n_ofdm_symbols, waveform.n_carriers) ,dtype=complex)
            for m in range(self.Nr_ru_antennas):
                #n_ofdm_symbols x (n_carriers)
                X_f[m, :, :] = waveform.extract_subcarriers(waveform.ofdm_time_to_freq(X_time[m, :, :]))

            # select channel
            H = self.get_csi(stripe_idx, active_ru_idx)  # [Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers]
            #H = np.ones((4, 4, 1024)) # debug with all ones channel
            print(f'channel shape: {H.shape}')

            # apply the channel
            for sym in range(waveform.n_ofdm_symbols):
                for k in range(self.Nr_subcarriers):
                    # X_f[:, sym, k] shape: [Nr_ru_antennas]
                    # H[:, :, k] shape: [Nr_ue_antennas x Nr_ru_antennas]
                    Y_f[sym, :, k] += H[:, :, k] @ X_f[:, sym, k] # sum because signals come from multiple stripes

            # pad zeros
            Y_f_padded = np.zeros((waveform.n_ofdm_symbols, self.Nr_ue_antennas, waveform.fft_size), dtype=complex)
            for m in range(self.Nr_ue_antennas):
                Y_f_padded[:, m, :] = waveform.pad_subcarriers(Y_f[:, m, :]) # (n_ofdm_symbols x n_carriers) => (n_ofdm_symbosl x fftsize)

            # IFFT back to time domain (should error => fft size + cp length
            Y_time = np.zeros((self.Nr_ue_antennas, waveform.n_ofdm_symbols, waveform.fft_size+waveform.cp_length), dtype=complex)
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








