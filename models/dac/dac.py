import numpy as np

from .component import Component
from .utils import limiter


class Dac(Component):
    def __init__(self, trunc_level: float = np.inf, nobits: float = np.inf, *args, **kwargs):
        self.trunc_level = trunc_level
        self.nobits = nobits

        super().__init__(*args, **kwargs)

    def run(self, yin):
        yin = np.array(yin)

        if self.trunc_level > 1e98:
            return yin
        if self.nobits > 20:
            yout, p = limiter(yin, self.trunc_level)
        else:
            step = 2 * self.trunc_level / (2 ** self.nobits)
            yout, p = limiter(
                np.round(yin / step) * step, self.trunc_level)

        return yout
