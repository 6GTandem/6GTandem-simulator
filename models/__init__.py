import math
import models.utils
import numpy as np
import numpy.fft
import numpy.random

from abc import ABC, abstractmethod
from scipy.signal import lfilter


def db_to_power(db_value: float | int):
    """Convert a value in decibels (dB) to a power value.

    The power/decibel relationship is as follows:
    db_value = 10 * log10(power_value)
    """
    return 10.0 ** (db_value / 10.0)


def power_to_db(power_value: float | int):
    """Convert a power value to decibels (dB).

    The power/decibel relationship is as follows:
    db_value = 10 * log10(power_value)
    """
    return 10 * np.log10(power_value)


def db_to_magnitude(db_value: float | int):
    """Convert a value in decibels (dB) to a magnitude value.

    The power/decibel relationship is as follows:
    db_value = 20 * log10(power_value)
    """
    return 10.0 ** (db_value / 20.0)


class Component(ABC):
    modes = ['ideal', 'linear', 'atan', 'tanh', 'poly3', 'poly3_pm',
             'poly5', 'limiter', 'softlimiter', 'weblab', '6gtandem']

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
        if mode not in self.modes:
            raise ValueError(f"Invalid mode: {mode}")

        self._mode = mode

    @abstractmethod
    def run(self, x):
        ...


class Coupler(Component):
    def __init__(self, damping: float = 0, *args, **kwargs):
        """Initialize a coupler instance.

        :param damping: The couplers damping in dB.
        """
        self.damping = damping

        super().__init__(*args, **kwargs)

    def run(self, x):
        return np.array(x) * db_to_magnitude(self.damping)


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

    def set_average_power(self, pin, pout):
        """The correct gain for a desired output power is calculated.

        Only holds true if there is no non-linearity. There is some discrepancy with NL.
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
                xout = (self.max_output_amplitude * x *
                        (1 - (alpha * (abs(x) ** 2)) + (beta * (abs(x) ** 4))))
            case 'limiter':
                xout = x
                # abs(xout) > 1 gives saturation at 1 V.
                xout[abs(xout) > 1] = xout[abs(xout) > 1] / \
                    abs(xout[abs(xout) > 1])
            case 'softlimiter' | '6gtandem':
                xout = models.utils.softlimiter(x, self.smoothness)

        return self.max_output_amplitude * xout


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
            yout, p = models.utils.limiter(yin, self.trunc_level)
        else:
            step = 2 * self.trunc_level / (2 ** self.nobits)
            yout, p = models.utils.limiter(
                np.round(yin / step) * step, self.trunc_level)

        return yout


class Fiber(Component):
    def __init__(self, length: float = 5, damping_per_meter: float = 5, fs: float = 15e9,
                 filter: float = 1, *args, **kwargs):
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
        xout = models.utils.delay(lfilter(self.filter, 1, x), [self.delay])

        return models.utils.setdbm(xout, models.utils.getdbm(x) - self.damping)

    @property
    def delay(self):
        return self.length * 1.5 / 3e8 * self.fs

    @property
    def damping(self):
        return self.length * self.damping_per_meter


class IQModem(Component):
    def __init__(self, iqi_coef: int = 0, iqi_filter: int = 1, iqi_delay_imbalance: int = 0,
                 dc_offset: int = 0, *args, **kwargs):
        self.iqi_coef = iqi_coef
        self.iqi_filter = iqi_filter
        self.iqi_delay_imbalance = iqi_delay_imbalance
        self.dc_offset = dc_offset
        self.modes.extend(['filter', 'static'])

        super().__init__(*args, **kwargs)

    def run(self, yin, phasor):
        yin = np.array(yin)
        phasor = np.array(phasor)
        assert (len(yin) == len(phasor)
                ), "Yin and Phasor data lengths do not match."

        match self.mode:
            case 'ideal':
                yout = yin
            case 'filter':
                d = (self.iqi_filter - 1) // 2
                data = list(np.conj(yin)) + list(np.zeros((d, 1)))
                xc = lfilter(self.iqi_filter, 1, data)
                xc = xc[d:]
                yout = yin + xc + self.dc_offset
            case 'static':
                yout = (models.utils.delay(np.real(yin), [-self.iqi_delay_imbalance / 2]) +
                        1j * models.utils.delay(np.imag(yin), [self.iqi_delay_imbalance / 2]))
                yout = yout + self.iqi_coef * np.conj(yout) + self.dc_offset

        return yout * phasor


class Link(Component):
    def __init__(self, *args, **kwargs):
        self.amp = Amplifier()
        self.fiber = Fiber()
        self.coupler_in = Coupler()
        self.coupler_out = Coupler()

        super().__init__(*args, **kwargs)

    def run(self, y):
        z = self.fiber.run(y)
        x1 = self.coupler_in.run(z)
        x = self.amp.run(x1)

        return self.coupler_out.run(x)


class Oscillator(Component):
    """Oscillator class.

    Example usage:
        >>> x = Oscillator()
        >>> pn = x.run_phase(1000) # Generates 1000 phase noise samples.

    # modes of pn generation 'ideal', 'cfo', 'model', 'spectrum'
    """

    def __init__(self, fs: float = 1.5e9, cfo: float = 0, cfo_std: float = 0, l100_db: float = -200,
                 linf_db: float = -300, f3db: float = 0, freq: list = [], spec: list = [], *args,
                 **kwargs):
        """

        :param fs: Sampling frequency in Hertz.
        :param cfo: Carrier frequency offset (receiver only) in Hertz.
        :param cfo_std: Carrier frequency random offset (receiver only) in Hertz
        :param l100_db: Pn level at 100 kHz in dB.
        :param linf_db: With pn level in dB.
        :param f3db: Cut-off frequency in Hertz.
        :param freq: Definition of spectrum response in Hertz.
        :param spec: Definition of spectrum response in dB.
        """
        self.fs = fs
        self.cfo = cfo
        self.cfo_std = cfo_std
        self.l100_db = l100_db
        self.linf_db = linf_db
        self.f3db = f3db
        self.freq = freq
        self.spec = spec

        self.last_phasor = 0
        self.current_phase = numpy.random.normal() * 2 * np.pi

        super().__init__(*args, **kwargs)

    @property
    def l0_db(self):
        return power_to_db(1e10 * (db_to_power(self.l100_db) - db_to_power(self.linf_db)) / self.f3db ** 2)

    @property
    def a(self):
        return np.exp(-2 * np.pi * self.f3db / self.fs)

    @property
    def pnvar(self):
        a2 = self.a ** 2
        if (1 - a2) < 1e-10:
            pnvar = (4 * (np.pi ** 2) * 1e10 *
                     db_to_power(self.l100_db) / self.fs)
        else:
            pnvar = ((1 - a2) * np.pi * 1e10 *
                     db_to_power(self.l100_db) / self.f3db)

        return pnvar

    @property
    def pnmean(self):
        # Extra variance to the CFO offset.
        return 2 * np.pi * (self.cfo + self.cfo_std * np.random.normal()) / self.fs

    @property
    def whitepnvar(self):
        return db_to_power(self.linf_db) * self.fs

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, mode):
        """Set a new mode.

        Valid options are:
            'ideal','cfo','model', 'spectrum', 'iid'
        """
        if mode not in ('ideal', 'cfo', 'model', 'spectrum', 'iid'):
            raise ValueError(f"Invalid mode: {mode}")

        self._mode = mode

    def run_phase(self, nosamples):
        match self.mode:
            case 'ideal':
                fi = np.zeros((nosamples, 1))
            case 'cfo':
                fi = self.run_phase_cfo(nosamples)
            case 'model':
                fi = self.run_phase_model(nosamples)
            case 'spectrum':
                fi = self.run_phase_spectrum(nosamples)

        return fi

    def run(self, nosamples):
        x = np.exp(1j * self.run_phase(nosamples))
        # Add AWGN
        if self.mode == 'model':
            x = x + np.sqrt(self.whitepnvar) * (numpy.random.normal(
                size=np.shape(x)) + 1j * numpy.random.normal(size=np.shape(x)))

        self.last_phasor = x

        return x

    def variance_phase_spectrum(self):
        n = 2 ** 23
        f = np.array([np.finfo(float).eps] + list(self.freq) + [self.fs])
        s = np.array([self.spec[0]] + list(self.spec) + [self.spec[-1]])
        flin = np.transpose(np.linspace(0, self.fs / 2, int(np.floor(n + 1))))
        x = db_to_power(np.interp(np.log(flin), np.log(f), s) / 2)
        x = x * np.sqrt(flin[1] - flin[0]) * n * 2

        return sum(abs(x[1:]) ** 2) / 2 / len(x) ** 2

    def run_phase_spectrum(self, nosamples):
        n = nosamples / 2
        f = np.array([np.finfo(float).eps] + list(self.freq) + [self.fs])
        s = np.array([self.spec[0]] + list(self.spec) + [self.spec[-1]])
        flin = np.linspace(0, self.fs / 2, int(np.floor(n + 1)))
        flin = flin[1:]
        x = db_to_power(np.interp(np.log(flin), np.log(f), s) / 2)
        x = x * np.sqrt(flin[1] - flin[0]) * n * 2
        x = x * np.sqrt(0.5) * (models.utils.randn_c(cols=np.shape(x)[0]))

        b = np.array([0] + list(x[1:]) + [0] + list(np.flipud(np.conj(x[1:]))))
        fi = np.real(numpy.fft.ifft(b))
        fi = fi + self.run_phase_cfo(np.shape(fi))
        self.current_phase = fi[-1]

        return fi

    def run_phase_model(self, nosamples):
        u = np.sqrt(self.pnvar) * \
            numpy.random.normal(size=nosamples) + self.pnmean
        fi = lfilter(1, [1, -self.a], np.insert(u, 0, self.current_phase))
        fi = fi[1:]
        self.current_phase = fi[-1]

        return fi

    def run_phase_cfo(self, nosamples):
        fi = np.cumsum(self.pnmean * np.ones(nosamples)) + self.current_phase
        self.current_phase = fi[-1]

        return fi


class Transmitter(Component):
    """Definition of a transmitter

    Components: c_oscillator, c_iqmod, c_pa, c_dac, c_antenna
    Usage:
        >>> tx = Transmitter()
        >>> y = tx.run(x)
    """

    def __init__(self, delay: float = 0, *args, **kwargs):
        self.oscillator = Oscillator()
        self.iqmodem = IQModem()
        self.amplifier = Amplifier()
        self.dac = Dac()
        self.delay = delay

        super().__init__(*args, **kwargs)

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
        if mode not in self.modes:
            raise ValueError(f"Invalid mode: {mode}")

        self.amplifier.mode = mode

        self._mode = mode

    def run(self, x, phasor=None):
        if phasor is None:
            phasor = self.oscillator.run(len(x))

        x1 = self.dac.run(x)
        x2 = self.iqmodem.run(x1, phasor)
        xout = self.amplifier.run(x2)

        if self.delay != 0:
            xout = models.utils.delay(xout, self.delay)

        return xout


class RadioStripe(Component):
    """Class representing a radio stripe.

    This class implements a radiostripe in the 6GTANDEM project.
    The stripe consists of a transmitter component, and a set of links.
    Usage:
        >>> rs = RadioStripe()
        >>> y = rs.run(x) # Y is a matrix
    """

    def __init__(self, nolinks: int = 3, bandwith: float = 5e9, os=5, *args, **kwargs):
        """
        :param nolinks: Number of links in the vector of links.
        """
        self.transmitter = Transmitter()
        self.bandwidth = bandwith
        self.os = os

        self.max_power = 10
        self.average_power = 5
        self.transmitter.mode = '6gtandem'
        self.transmitter.amplifier.set_maximum_output_power(self.max_power)
        self.transmitter.amplifier.set_average_power(15, self.average_power)

        self.links = []
        for l in range(nolinks):
            link = Link()
            link.amp.mode = '6gtandem'
            link.amp.set_maximum_output_power(self.max_power)
            link.amp.set_average_power(self.average_power - link.fiber.damping -
                                       link.coupler_in.damping - link.coupler_out.damping, self.average_power)
            link.amp.set_noise_var(300, self.bandwidth * self.os, 10)

            self.links.append(link)

        super().__init__(*args, **kwargs)

    def run(self, x):
        y = np.zeros((len(x), 1 + len(self.links)), dtype=np.complex128)
        y[:, 0] = self.transmitter.run(x)[:, 0]

        for i, link in enumerate(self.links):
            y[:, i+1] = link.run(np.transpose([y[:, i]]))[:, 0]

        return y

    def calibrate(self, x, desired_amplifier_dbm):
        """This function sets the small-signal gain of the link amplifiers.

        To give an approximate constant power(DesiredAmplifierDBM) at the output of each link.
        The calibration is valid for a given input signal, and recalibration must be performed if
        the signal statistics changes.
        """
        for j in range(3):
            z = self.transmitter.run(x)
            scale = db_to_magnitude(
                desired_amplifier_dbm - models.utils.getdbm(z))
            self.transmitter.amplifier.gain = self.transmitter.amplifier.gain * scale

        z = self.transmitter.run(x)
        for link in self.links:
            for j in range(3):
                z2 = link.run(z)
                scale = db_to_magnitude(
                    desired_amplifier_dbm - models.utils.getdbm(z2))
                link.amp.gain = link.amp.gain * scale

            z = link.run(z)
