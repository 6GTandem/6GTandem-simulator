import numpy as np
from scipy.signal import lfilter

from ..component.component import Component
from ..utils import db_to_magnitude


class Coupler(Component):
    def __init__(self, damping: float = 0, filter: np.ndarray = np.array([1]), *args, **kwargs):
        """Initialize a coupler instance.

        :param damping: The couplers damping in dB.
        """
        self.damping = damping
        self.filter = filter

        super().__init__(*args, **kwargs)

    def run(self, x):
        xout = lfilter(self.filter, [1.0], x)

        return xout * db_to_magnitude(self.damping)
