import numpy as np
import numpy.random

from ..component.component import Component
from scipy.signal import lfilter
from ..utils import power_to_db, db_to_power, randn_c


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
                fi = np.zeros((1, nosamples))
            case 'cfo':
                fi = self.run_phase_cfo(nosamples)
            case 'model':
                fi = self.run_phase_model(nosamples)
            case 'spectrum':
                fi = self.run_phase_spectrum(nosamples)
            case _:
                fi = np.zeros((1, nosamples))

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
        x = x * np.sqrt(0.5) * (randn_c(cols=np.shape(x)[0]))

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
