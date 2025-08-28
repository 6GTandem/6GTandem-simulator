import numpy as np
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

    def __init__(self, channelmodel: str = 'subTHz-Rayleigh',
                 Nr_subcarriers: int = 1024, Nr_ue_antennas: int = 1, Nr_ru_antennas: int =1,
                 Nr_rus: int=5, Nr_stripes: int=2):
        self.Nr_stripes = Nr_stripes
        self.Nr_rus = Nr_rus
        self.Nr_ue_antennas = Nr_ue_antennas
        self.Nr_ru_antennas = Nr_ru_antennas
        self.Nr_subcarriers = Nr_subcarriers

        if channelmodel == 'subTHz-Rayleigh':
            self.csi = self.subTHz_Rayleigh()
        elif channelmodel == 'sionna':
            raise NotImplementedError("sionna channel model is not implemented yet.")


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


    def transmit_dl(self, X_list, active_ru_per_stripe):
        """"
        transmit per stripe from a RU to UE

         X_list (list of np.ndarray):
            List of transmit signals, one per stripe.
            Each element is either:
                - None  (if no RU in that stripe transmits)
                - np.ndarray of shape [Nr_ru_antennas x Nr_samples] for the chosen RU
         active_ru_per_stripe (list of int or None):
                    For each stripe, the index of the transmitting RU.
                    If None, that stripe is silent.
         Returns:
                np.ndarray: Received signal at UE of shape
                            [Nr_ue_antennas x Nr_subcarriers]
        """

        Y_dl = np.zeros((self.Nr_ue_antennas, self.Nr_subcarriers), dtype=complex)

        for stripe_idx, ru_idx in enumerate(active_ru_per_stripe):
            if ru_idx is None or X_list[stripe_idx] is None:
                continue  # no transmission from this stripe

            # Time-domain transmit signal for this RU
            X_time = X_list[stripe_idx]  # shape: [Nr_ru_antennas x Nr_samples]

            # FFT to frequency domain
            X_f = np.fft.fft(X_time, n=self.Nr_subcarriers,
                             axis=-1)  # [Nr_ru_antennas x Nr_subcarriers]

            # select channel
            H = self.csi[stripe_idx, ru_idx, :, :, :]  # [Nr_ue_antennas x Nr_ru_antennas x Nr_subcarriers]

            for k in range(self.Nr_subcarriers):
                Y_dl[:, k] += H[:, :, k] @ X_f[:, k]

        # move back to time domain
        Y_time = np.fft.ifft(Y_dl, n=self.Nr_subcarriers, axis=-1)
        return Y_time








