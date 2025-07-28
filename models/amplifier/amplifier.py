import math
import numpy as np
import numpy.random

from .component import Component
from .utils import db_to_power, db_to_magnitude, softlimiter


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

    def __init__(self, gain=1, max_out_amp=1, noise_var=0, smoothness=1, *args, **kwargs):
        """Initialize the amplifier object with the given parameters.

        :param gain: The low signal-gain of the amplifier.
        :param max_out_amp: The maximum output amplitude in Volts.
        :param noise_var: The variance of the AWGN [V^2] per channel.
        :param smoothness: Used in Mode 'softlimiter'.
        """
        self.gain = gain
        self.max_output_amplitude = max_out_amp
        self.noise_var = noise_var
        self.smoothness = smoothness

        super().__init__(*args, **kwargs)

    def set_maximum_output_power(self, max_power_dbm):
        """ Set the maximum output amplitude.

        :param max_power_dbm: The maximum output power in dBm.
        """
        self.max_output_amplitude = math.sqrt(
            50 * db_to_power(max_power_dbm) * 1e-3)

    def set_gain(self, pin, pout):
        """The correct gain for a desired output power is calculated.

        Only holds true if there is no non-linearity. There is some discrepancy with NL.

        <Matlab function> set_average_power
        """
        self.gain = db_to_magnitude(pout - pin)

    def set_noise_var(self, t, b, nfdb):
        k = 1.3806e-23
        noise_density = k * t * db_to_power(nfdb)
        noise_power = noise_density * b
        # P = U^2/R
        voltage_power = noise_power * 50
        # Noise variance per channel.
        self.noise_var = voltage_power / 2

    def run(self, x):
        x = np.array(x)
        x = self.gain / self.max_output_amplitude * (x + math.sqrt(self.noise_var) * (
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
                alpha = ((10 / 9) -
                         np.sqrt(((20 / 9) * (10 ** -0.05)) - (80 / 81)))
                # gives minimum derivative=0, at x=sqrt(2/3/alpha)
                beta = (9 * (alpha ** 2)) / 20
                xout = (x * (1 - (alpha * (abs(x) ** 2)) + (beta * (abs(x) ** 4))))
            case 'limiter':
                xout = x
                # abs(xout) > 1 gives saturation at 1 V.
                xout[abs(xout) > 1] = xout[abs(xout) > 1] / \
                    abs(xout[abs(xout) > 1])
            case 'softlimiter' | '6gtandem':
                xout = softlimiter(x, self.smoothness)

        return self.max_output_amplitude * xout
