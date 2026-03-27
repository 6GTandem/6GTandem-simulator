import numpy as np

from ..component.component import Component

# TODO REIMPLEMENT


class PhaseShifter(Component):
    def __init__(self, num_shifters: int, resolution: int, in_degrees=False, *args, **kwargs):
        """

        :param num_shifters: Number of phase shifters.
        :param resolution: Resolution in number of bits.
        """
        self.num_shifters = num_shifters
        self.resolution = resolution
        self.in_degrees = in_degrees

        super().__init__(*args, **kwargs)

    def run(self, x, shifts: list[int] | np.ndarray):
        """Apply the phase shift to the input data.

        :param x: Array with the input data.
        :param shifts: List of phase shift indexes to apply to the input data.
                       The phase shift index is an integer between 0 - (2 ^ resolution) - 1
        """

        def shift(x, k):
            return x * np.exp(1j * k)

        y = []
        assert self.num_shifters == x.shape[0], f"Input data shape ({x.shape[0]}) does not match number of phase shifters ({self.num_shifters})."
        for r, k in zip(x, shifts):
            if self.in_degrees:
                k = np.deg2rad(k)
            y.append(shift(r, k))

        return np.array(y)

    def get_phases(self, beam_angle: float):
        assert (
            beam_angle <= 90 and beam_angle >= -90
        ), f"The requested beam angle ({beam_angle}) doesn't lie between 90 and -90deg."

        spacing = np.arange(1, self.num_shifters + 1) / 2
        phi = 2 * np.pi * spacing * np.sin(np.deg2rad(beam_angle))

        return phi
