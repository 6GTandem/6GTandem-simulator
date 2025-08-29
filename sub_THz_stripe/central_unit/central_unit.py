import numpy as np

from ..oscillator.oscillator import Oscillator
from ..iqmodem.iqmodem import IQModem
from ..amplifier.amplifier import Amplifier
from ..dac.dac import Dac
from ..component.component import Component
from ..utils import delay, randconst, pulseshape


class CentralUnit(Component):
    """Definition of a transmitter

    Components: c_oscillator, c_iqmod, c_pa, c_dac, c_antenna
    Usage:
        >>> tx = Transmitter()
        >>> y = tx.run(x)
    """

    def __init__(self, oscillator: Oscillator, iqmodem: IQModem, amplifier: Amplifier, dac: Dac, delay: float = 0, nosamples: int = 1000, phasor: np.ndarray | None = None, waveform: str = "Gaussian", *args, **kwargs):
        """Instantiate a transmitter based on an oscillator, iqmodem, amplifier and dac.
        n : number of IQ samples
        waveform: wave form can be selected from the list: ["Gaussian-ideal", "Gaussian-imparied", "OFDM", "CP-OFDM"]
        """
        self.waveform = waveform
        self.delay = delay

        self.x = x
        self.y = y
        self.z = z

        self.oscillator = oscillator
        self.iqmodem = iqmodem
        self.amplifier = amplifier
        self.dac = dac
        self.nosamples = nosamples
        self.phasor = phasor

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

    def getGaussianSymbols(self, K=1, Ndata=1024, p=1):
        """Doc
        :param Ndata: number of symbols to generate
        :param p: signal variance
        :return: K x Ndata: symbols sampled from a complex gaussian with variance p

        Note that this generates a variable which is drawn from a complex gaussian distribution with variance p
        which is equivalent to a + bj with a and b sampled from a gaussian distribution with variance p/2.
        Here we first sample a and b from a gaussian with mean 0 and variance 1, by multiplying with sqrt(p)/sqrt(2)
        we obtain variance p/2 for both a and b, given that var(constant * X) = constant^2 var(X)
        """
        s = np.sqrt(p) / np.sqrt(2) * (np.random.randn(K, Ndata) +
                                       1j * np.random.randn(K, Ndata))
        return s.astype(np.complex64)

    def run(self, x=None, phasor=None):
        """
        Generates and processes waveform symbols based on the selected waveform type.
        Parameters:
            x (optional): Input symbols. For most waveform types, these should be RRC pulse-shaped QAM symbols.
                          For "Gaussian-ideal" and "Gaussian-impaired", symbols are generated internally.
            phasor (optional): Phasor values for modulation. If not provided, generated using the oscillator.
        Returns:
            np.ndarray: The processed output waveform after DAC, IQ modulation, amplification, and optional delay.
        Notes:
            - For "Gaussian-ideal", Gaussian symbols are generated and returned.
            - For "Gaussian-impaired", Gaussian symbols are generated, processed through DAC, IQ modem, and amplifier.
            - For other waveform types, input symbols are processed through DAC, IQ modem, amplifier, and optional delay.
            - The function is designed to eventually generate RRC pulse-shaped QAM symbols internally, removing the need for 'x' as an input.
        """
        # todo: code of Thomas E. required RRC pulse shaped qam symbols to be input here (x)
        # todo: however this should be generated here on the fly! and removed as input parameters!
        match self.waveform:
            case "Gaussian":
                xout = generate_gaussian_symbols(ndata=self.nosamples)
            case "QAM":
                # generate a QAM signal, and use RRC pulse shaping
                os = 5
                # number of symbols, and oversampling factor
                xout, c = randconst(1, self.nosamples)
                xout = pulseshape(xout, os, 0.1)
            case _:
                raise ValueError(f"Unknown waveform: {self.waveform}")

        if self.mode != "ideal":
            xout = generate_gaussian_symbols()

            phasor = self.phasor
            if phasor is None:
                phasor = self.oscillator.run(xout.shape[1])

            xout = self.dac.run(xout)
            xout = self.iqmodem.run(xout, phasor)
            xout = self.amplifier.run(xout)

            if self.delay != 0:
                xout = delay(xout, [self.delay])

        return xout
