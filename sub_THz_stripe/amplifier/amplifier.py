import math
import numpy as np
import pandas as pd
import numpy.random

from ..component.component import Component
from ..utils import db_to_power, db_to_magnitude, softlimiter


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

    def __init__(self, bw, gain=1, max_gain=1, noise_fig=0, smoothness=1, *args, **kwargs):
        """Initialize the amplifier object with the given parameters.

        :param gain: The low signal-gain of the amplifier.
        :param max_gain: The maximum gain that can be physically achieved by the amplifier.
        :param max_out_amp: The maximum output amplitude in Volts.
        :param noise_var: The variance of the AWGN [V^2] per channel.
        :param smoothness: Used in Mode 'softlimiter'.
        """
        self.coeffs = None
        self.max_input_amplitude = None
        self.polynomial = None
        self._gain = 1
        self.max_gain = max_gain
        self.gain = gain
        self.noise_fig = noise_fig
        self.set_noise_var(25 + 273.15, bw, self.noise_fig)
        self.smoothness = smoothness

        super().__init__(*args, **kwargs)

        if self.mode.startswith('poly'):
            self.polynomial = int(self.mode[-1])
            self.coeffs, self.max_input_amplitude = self.polynomial_fitting(self.polynomial, self.gain)
    
    @property
    def gain(self):
        return self._gain
    
    @gain.setter
    def gain(self, new_gain):
        assert new_gain <= self.max_gain, f"Gain is higher than the maximum achievable gain: {new_gain} > {self.max_gain}"
        self._gain = new_gain

        # The polynomial coefficients have to be recalculated in case the gain changes.
        if self.polynomial is not None:
            self.coeffs, self.max_input_amplitude = self.polynomial_fitting(self.polynomial, self.gain)

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
        noise_power_watts =  10 ** (nfdb / 10) * k * t * b

        self.noise_var = noise_power_watts

    
    def polynomial_fitting(self, order, gain):
        # Read out the measurement data.
        data = pd.read_csv('models/amplifier/Mag_150GHz_BalCascode2stage_V3.vcsv', skiprows=6, names=["x", "yre", "yim"])
        # Go from the input values in dBm to values in Watt.
        x_power = (10 ** (data["x"] / 10)) / 1000
        # Now convert the values in Watt to Voltages.
        x_voltage = np.sqrt(x_power * 50)
        # Calculate the magnitude of the output voltage based on the real and complex values.
        y_voltage = np.sqrt(data["yre"] ** 2 + data["yim"] ** 2)
        # What is the current gain of the amplifier?
        current_gain = ((y_voltage[4] - y_voltage[0]) / (x_voltage[4] - x_voltage[0]))
        scale = current_gain / gain
        # Re-scale the input and output to get the new gain.
        x_scaled = x_voltage * scale

        # Only perform fitting for the odd powers.
        odd_powers = list(range(1, order + 1, 2))
        x_powers = np.column_stack([x_scaled ** p for p in odd_powers])
    
        coeffs, *_ = np.linalg.lstsq(x_powers, y_voltage, rcond=None)

        return coeffs, np.max(x_scaled)

    def run(self, x: np.ndarray):
        noise = np.random.normal(0, np.sqrt(self.noise_var / 2), size=np.shape(x)) + 1j * numpy.random.normal(0, np.sqrt(self.noise_var / 2), size=np.shape(x))
        x_noise = x + noise

        match self.mode:
            case 'ideal' | 'linear':
                xout = x_noise * self.gain
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
                if self.polynomial is None or self.max_input_amplitude is None or self.coeffs is None:
                    raise ValueError("Selected polynomial model but no coefficients configured")
                xlim = x_noise.copy()
                xlim[np.abs(x_noise) > self.max_input_amplitude] = self.max_input_amplitude
                
                xout = sum(c * xlim * np.abs(xlim) ** (p-1) for c, p in zip(self.coeffs, range(1, self.polynomial + 1, 2)))
            case 'limiter':
                xout = x
                # abs(xout) > 1 gives saturation at 1 V.
                xout[abs(xout) > 1] = xout[abs(xout) > 1] / \
                    abs(xout[abs(xout) > 1])
            case 'softlimiter' | '6gtandem':
                xout = softlimiter(x, self.smoothness)

        return xout
