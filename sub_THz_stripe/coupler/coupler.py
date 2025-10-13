import numpy as np
import logging
from scipy.signal import lfilter

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
