import numpy as np

from .component import Component
from .utils import db_to_magnitude


class Coupler(Component):
    def __init__(self, damping: float = 0, *args, **kwargs):
        """Initialize a coupler instance.

        :param damping: The couplers damping in dB.
        """
        self.damping = damping

        super().__init__(*args, **kwargs)

    def run(self, x):
        return np.array(x) * db_to_magnitude(self.damping)
