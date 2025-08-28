import numpy as np
from ..component.component import Component


class Combiner(Component):
    def run(self, x):
        """Combines incoming data into one stream by adding everything up.
        """
        return np.array([np.sum(x, axis=0)])
