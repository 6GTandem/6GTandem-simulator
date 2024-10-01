import math
import models.utils
import numpy as np
import numpy.random

from abc import ABC, abstractmethod


def db_to_power(db_value: float | int):
    """Convert a value in decibels (dB) to a power value.

    The power/decibel relationship is as follows:
    db_value = 10 * log10(power_value)
    """
    power = 10.0 ** (db_value / 10.0)

    return power


def db_to_magnitude(db_value: float | int):
    """Convert a value in decibels (dB) to a magnitude value.

    The power/decibel relationship is as follows:
    db_value = 20 * log10(power_value)
    """
    power = 10.0 ** (db_value / 20.0)

    return power


class Component(ABC):
    def __init__(self, mode: str = 'ideal'):
        self._mode = mode

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, mode):
        """Set a new mode.

        Valid options are:
            'ideal','linear','atan', 'tanh', 'poly3','poly3_pm',
            'poly5','limiter','softlimiter','weblab','6gtandem'
        """
        if mode not in ('ideal', 'linear', 'atan', 'tanh', 'poly3', 'poly3_pm', 'poly5', 'limiter',
                        'softlimiter', 'weblab', '6gtandem'):
            raise ValueError(f"Invalid mode: {mode}")

        self._mode = mode

    @abstractmethod
    def run(self, x):
        ...


class Amplifier(Component):
    """Amplifier class.

    This Amplifier class implements a nonlinear function with a low-signal gain of 1 and a max
    output amplitude of 1 (13 dBm).

    Y = G * f(X + W)
    where X is the input signal and W is noise. f() is a nonlinear function with max amplitude of
    1 V. G is the output gain.

    The gain, nonlinearity function and noise power level can be controlled.
    Example usage:
        >>> pa = Amplifier()
        >>> pa.mode = 'tanh'
        >>> pa.gain = 3.2
        >>> y = pa.run(x)
    """

    def __init__(self, gain=1, max_out_amp=1, noise_var=0, smoothness=1):
        """Initialize the amplifier object with the given parameters.

        :param gain: The low signal-gain of the amplifier.
        :param max_out_amp:
        :param noise_var: The variance of the AWGN [V^2] per channel.
        :param smoothness: Used in Mode 'softlimiter'.
        """
        self.gain = gain
        self._max_output_amplitude = max_out_amp
        self._noise_var = noise_var
        self.smoothness = smoothness

        super().__init__()

    @property
    def noise_var(self):
        return self._noise_var

    @noise_var.setter
    def noise_var(self, t, b, nfdb):
        k = 1.3806E-23
        noise_density = k * t * db_to_power(nfdb)
        noise_power = noise_density * b
        # P = U^2/R
        voltage_power = noise_power * 50
        # Noise variance per channel.
        self._noise_var = voltage_power / 2

    @property
    def max_output_amplitude(self):
        return self._max_output_amplitude

    @max_output_amplitude.setter
    def max_output_amplitude(self, max_power_dbm):
        """ Set the maximum output amplitude.

        :param max_power_dbm: The maximum output power in dBm.
        """
        self._max_output_amplitude = math.sqrt(
            50 * db_to_power(max_power_dbm) * 1e-3)

    def set_average_power(self, pin, pout):
        """The correct gain for a desired output power is calculated.

        Only holds true if there is no non-linearity. There is some discrepancy with NL.
        """
        self.gain = db_to_magnitude(pout - pin)

    def run(self, x):
        x = self.gain / self._max_output_amplitude * (x + math.sqrt(self._noise_var) * (
            numpy.random.normal(size=np.shape(x)) + 1j * numpy.random.normal(size=np.shape(x))))

        match self.mode:
            case 'ideal' | 'linear':
                xout = x
            case 'atan':
                # alpha=0.6340 # This factor gives 1dB compression at x=1
                alpha = 2 / np.pi  # This factor gives Amax=1;
                xout = 1 / alpha * \
                    np.atan(alpha * abs(x)) * np.exp(1j * np.angle(x))
            case 'tanh':
                # alpha=0.6125# This factor gives 1dB compression at x=1
                alpha = 1  # This factor gives Amax=1;
                xout = 1 / alpha * \
                    np.tanh(alpha * abs(x)) * np.exp(1j * np.angle(x))
            case 'poly3':
                # alpha=1-10^(-0.05*1) # This factor gives 1dB compression at x=1
                alpha = 4/27  # This factor gives Amax=1;
                xout = x * (1 - alpha * abs(x) ** 2)
                # dont allow for negative gain
                peakx = math.sqrt(1 / 3 / alpha)
                xout[abs(x) > peakx] = peakx * (1 - alpha * abs(peakx)
                                                ** 2) * np.exp(1j * np.angle(xout[abs(x) > peakx]))
            case 'poly3_pm':
                # alpha=(1-10^(-0.05*1))*exp(1i*0.2) # This factor gives 1dB compression at x=1
                alpha = (4 / 27) * np.exp(1j * 0.2)  # This factor Amax=1
                xout = x * (1 - alpha * abs(x) ** 2)
                # dont allow for negative gain
                peakx = math.sqrt(1 / 3 / abs(alpha))
                xout[abs(x) > peakx] = peakx * (1 - alpha * abs(peakx)
                                                ** 2) * np.exp(1j * np.angle(xout[abs(x) > peakx]))
            case 'poly5':
                # This factor gives 1dB compression at x=1
                alpha = (10 / 9) - np.sqrt((20 / 9 * 10) ** (-0.05 - 80 / 81))
                # gives minimum derivative=0, at x=sqrt(2/3/alpha)
                beta = 9 * alpha ** (2 / 20)
                xout = self._max_output_amplitude * x * \
                    (1 - alpha * abs(x) ** 2 + beta * abs(x) ** 4)
            case 'limiter':
                xout = x
                # abs(xout) > 1 gives saturation at 1 V.
                xout[abs(xout) > 1] = xout[abs(xout) > 1] / \
                    abs(xout[abs(xout) > 1])
            case 'softlimiter' | '6gtandem':
                xout = models.utils.softlimiter(x, self.smoothness)

        return self._max_output_amplitude * xout
