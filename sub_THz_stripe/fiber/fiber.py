import numpy as np
import logging

from ..component.component import Component
from wireless_channel.waveforms import Waveform
from ..utils import db_to_magnitude, delay
from scipy.signal import lfilter

logger = logging.getLogger(__name__)


class Fiber(Component):
    def __init__(self, wf: Waveform, length: float = 1, damping_per_meter: float = 0, fs: float = 15e9, filter_mode='time_domain',
                 filter: np.ndarray = np.array([1]), *args, **kwargs):
        """Initialize a fiber component.

        :param length: Length of the fiber in meter.
        :param damping_per_meter: Damping in dB per meter.
        :param filter: The impulse response of the fiber.


        If filter is used, the damping and delay are included in the filter. In this case set length and
        damping_per_meter to zero.
        """
        self.length = length
        self.damping_per_meter = damping_per_meter
        self.fs = fs
        self.filter = filter # passed in frequency domain
        self.filter_mode = filter_mode
        self.wf = wf

        super().__init__(*args, **kwargs)

    def run(self, x):
        if self.filter_mode == 'time_domain':
            logger.debug("Fiber model is being applied in the time domain.")
            taps = np.fft.ifft(np.fft.ifftshift(self.filter))
            logger.debug(f"Filter taps used: {taps}")
            #x_filt = delay(lfilter(taps, [1.0], x.flatten(), [self.delay]) #todo delay needed or not???
            x_filt = lfilter(taps, [1.0], x.flatten())
            xout = np.reshape(x_filt, (self.wf.n_ofdm_symbols, -1))

        elif self.filter_mode == 'freq_domain':
            logger.debug("Fiber model being applied in the frequency domain.")
            logger.debug(f'Filter shape: {self.filter.shape}')
            logger.debug(f'Filter: {self.filter}')

            x_freq = self.wf.ofdm_time_to_freq(x)
            logger.debug(f'xfreq shape: {x_freq.shape}')

            f_shift = np.fft.fftshift(self.filter)
            x_freq_filtered = x_freq * f_shift

            logger.debug(f'xfreq filtered shape: {x_freq_filtered.shape}')
            xout = self.wf.ofdm_freq_to_time(x_freq_filtered) #back to time

        return xout * db_to_magnitude(self.damping)

    @property
    def delay(self):
        return self.length * 1.5 / 3e8 * self.fs  # 1.5 / 3e8 speed of EM waves in fiber

    @property
    def damping(self):
        return self.length * self.damping_per_meter
