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
        """
        self.length = length
        self.damping_per_meter = damping_per_meter
        self.fs = fs
        self.filter = filter

        super().__init__(*args, **kwargs)

    def run(self, x):
        xout = delay(lfilter(self.filter, [1.0], x), [self.delay])

        xout = lfilter(self.filter, [1.0], x)

        return xout * db_to_magnitude(self.damping)

    @property
    def delay(self):
        return self.length * 1.5 / 3e8 * self.fs  # 1.5 / 3e8 speed of EM waves in fiber

    @property
    def damping(self):
        return self.length * self.damping_per_meter
