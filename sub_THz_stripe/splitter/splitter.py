import numpy as np
from ..component.component import Component


class Splitter(Component):
    def __init__(self, num_splits: int, *args, **kwargs):
        """
        :param num_splits: How many branches the splitter splits into.
        """
        self.splits = num_splits
        super().__init__(*args, **kwargs)

    def run(self, x):
        """Returns the data split into `num_splits` branches.

        A fixed attenuation is applied based on the number of branches: 1 / num_splits
        """
        attenuation = 1 / np.sqrt(self.splits)

        x = np.broadcast_to(x, (self.splits,) + x.shape).copy()
        return x * attenuation
