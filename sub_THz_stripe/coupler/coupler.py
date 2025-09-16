import numpy as np
from scipy.signal import lfilter, filtfilt

from wireless_channel.waveforms import Waveform
from ..component.component import Component
from ..utils import db_to_magnitude


class Coupler(Component):
    def __init__(self, wf: Waveform, damping: float = 0, filter: np.ndarray = np.array([1]), filter_mode='time_domain', *args, **kwargs):
        """Initialize a coupler instance.

        :param damping: The couplers damping in dB.
        """
        self.damping = damping
        self.filter = filter
        self.filter_mode = filter_mode
        self.wf = wf

        if self.filter == 'freq_domain' and wf is None:
            raise ValueError("wf should not be None when using frequency domain filtering.")

        self.wf = wf

        super().__init__(*args, **kwargs)

    def run(self, x):
        if self.filter_mode == 'time_domain':
            taps = np.fft.ifft(np.fft.ifftshift(self.filter))
            #x_filt = delay(lfilter(taps, [1.0], x.flatten(), [self.delay]) #todo delay needed or not???
            x_filt = lfilter(taps, [1.0], x.flatten())
            xout = x_filt.reshape(self.wf.n_ofdm_symbols, -1)

        elif self.filter_mode == 'freq_domain':
            print(f' i a m in freq domaiiinnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn')
            print(f'filter : {self.filter.shape}')
            print(f'filter: {self.filter}')
            x_freq = self.wf.ofdm_time_to_freq(x)
            print(f'xfreq : {x_freq.shape}')

            x_freq_filtered = x_freq * self.filter

            print(f'xfreq filtered : {x_freq_filtered.shape}')
            xout = self.wf.ofdm_freq_to_time(x_freq_filtered) #back to time

        return xout * db_to_magnitude(self.damping)
