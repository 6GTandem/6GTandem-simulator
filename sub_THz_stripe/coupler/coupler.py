import os
import numpy as np
import logging
from scipy.signal import lfilter
import skrf as rf

from wireless_channel.waveforms import Waveform
from ..component.component import Component
from ..utils import db_to_magnitude

logger = logging.getLogger(__name__)

class Coupler(Component):
    def __init__(self, wf: Waveform, damping: float = 0, filter: np.ndarray = np.array([1]), filter_mode='time_domain', *args, **kwargs):
        """Initialize a coupler instance.

        :param damping: The couplers damping in dB.

        When the filter is used, set the damping to 0 as it is included in the filter itself.
        """
        self.damping = damping
        self.filter = filter
        self.filter_mode = filter_mode
        self.wf = wf

        super().__init__(*args, **kwargs)

    def run(self, x):
        if self.filter_mode == 'time_domain':
            logger.debug("Coupler model is being applied in the time domain.")
            taps = np.fft.ifft(np.fft.ifftshift(self.filter))
            logger.debug(f"Filter taps used: {taps}")
            #x_filt = delay(lfilter(taps, [1.0], x.flatten(), [self.delay]) #todo delay needed or not???
            x_filt = lfilter(taps, [1.0], x.flatten())
            xout = np.reshape(x_filt, (self.wf.n_ofdm_symbols, -1))

        elif self.filter_mode == 'freq_domain':
            logger.debug("Coupler model being applied in the frequency domain.")
            logger.debug(f'Filter shape: {self.filter.shape}')
            logger.debug(f'Filter: {self.filter}')

            x_freq = self.wf.ofdm_time_to_freq(x)
            logger.debug(f'xfreq shape: {x_freq.shape}')

            f_shift = np.fft.fftshift(self.filter)
            x_freq_filtered = x_freq * f_shift

            logger.debug(f'xfreq filtered shape: {x_freq_filtered.shape}')
            xout = self.wf.ofdm_freq_to_time(x_freq_filtered) #back to time

        return xout * db_to_magnitude(self.damping)
    
    @classmethod
    def from_config(cls, config: dict, wf: Waveform):
        damping = config.get("damping", 0)  # in dB

        if "model" in config:
            # Load the couplers S-parameter file
            base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            coupler_spars = rf.Network(os.path.join(base_path, f"../models/coupler/{config['model']}"))

            # Subcarrier frequencies
            ofdm_freqs = np.linspace(
                wf.fc - wf.bw / 2, wf.fc + wf.bw / 2, wf.n_carriers * wf.oversampling_factor
            )

            # Extract frequency and S21 (transmission)
            coup_freqs = coupler_spars.f
            coup_s21 = coupler_spars.s[:, 1, 0]  # S21

            # Interpolate magnitude and phase separately for better accuracy
            coup_s21_mag = np.abs(coup_s21)
            coup_s21_phase = np.angle(coup_s21)

            interp_mag = np.interp(ofdm_freqs, coup_freqs, coup_s21_mag)
            interp_phase = np.interp(ofdm_freqs, coup_freqs, coup_s21_phase)
            coup_s21_ofdm = interp_mag * np.exp(1j * interp_phase)

            filter_mode = "freq_domain"
            coup = Coupler(damping=damping, filter=coup_s21_ofdm, filter_mode=filter_mode, wf=wf)
        else:
            coup = Coupler(damping=damping, wf=wf)

        return coup