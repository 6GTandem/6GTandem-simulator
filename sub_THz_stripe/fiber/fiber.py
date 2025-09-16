from ..component.component import Component
from ..utils import db_to_magnitude, delay
from scipy.signal import lfilter
import numpy as np


class Fiber(Component):
    def __init__(self, length: float = 1, damping_per_meter: float = 0, fs: float = 15e9,
                 filter: np.ndarray = np.array([1]), *args, **kwargs):
        """Initialize a fiber component.

        :param length: Length of the fiber in meter.
        :param damping_per_meter: Damping in dB per meter.
        :param filter: The impulse response of the fiber.


        If filter is used, the damping is included in the filter. In this case
        set the damping_per_meter to zero.
        """
        self.length = length
        self.damping_per_meter = damping_per_meter
        self.fs = fs
        self.filter = filter

        super().__init__(*args, **kwargs)

    def run(self, x):
        #xout = delay(lfilter(self.filter, [1.0], x), [self.delay])

        #xout = lfilter(self.filter, [1.0], x)

        cp_length = 128
        n_fft = x.shape[1] - cp_length
        ofdm_freq = []
        for symbol in x:
            time_no_cp = symbol[cp_length:]
            freq_domain = np.fft.fft(time_no_cp, n_fft)
            ofdm_freq.append(freq_domain)
        ofdm_freq = np.array(ofdm_freq)

        xout = ofdm_freq * self.filter

        ofdm_time = []
        n_fft = xout.shape[1]
        for symbol in xout:
            time_domain = np.fft.ifft(symbol, n_fft)
            cp = time_domain[-cp_length:]
            ofdm_symbol = np.concatenate([cp, time_domain])
            ofdm_time.append(ofdm_symbol)
        xout = np.array(ofdm_time) # shape: n_ofdm_symbols x ( n_carriers * oversampling + cp_length)


        return xout * db_to_magnitude(self.damping)

    @property
    def delay(self):
        return self.length * 1.5 / 3e8 * self.fs  # 1.5 / 3e8 speed of EM waves in fiber

    @property
    def damping(self):
        return self.length * self.damping_per_meter
