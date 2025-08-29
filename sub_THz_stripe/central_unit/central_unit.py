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

    def __init__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0, oscillator: Oscillator = None, iqmodem: IQModem = None, amplifier: Amplifier = None,
                 dac: Dac = None, delay: float = 0, nosamples: int = 1000, phasor: np.ndarray | None = None
                 , *args, **kwargs):
        """Instantiate a transmitter based on an oscillator, iqmodem, amplifier and dac.
        n : number of IQ samples
        waveform: wave form can be selected from the list: ["Gaussian-ideal", "Gaussian-imparied", "OFDM", "CP-OFDM"]
        """
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

    def run(self, x, phasor=None):
        """
        Generates and processes waveform symbols based on the selected waveform type.
        Parameters:
            x : Input symbols.
            phasor (optional): Phasor values for modulation. If not provided, generated using the oscillator.
        Returns:
            np.ndarray: The processed output waveform after DAC, IQ modulation, amplification, and optional delay.
        Notes:
            - For "Gaussian-ideal", Gaussian symbols are generated and returned.
            - For "Gaussian-impaired", Gaussian symbols are generated, processed through DAC, IQ modem, and amplifier.
            - For other waveform types, input symbols are processed through DAC, IQ modem, amplifier, and optional delay.
            - The function is designed to eventually generate RRC pulse-shaped QAM symbols internally, removing the need for 'x' as an input.
        """
        xout = x

        if self.mode != "ideal":
            phasor = self.phasor
            if phasor is None:
                phasor = self.oscillator.run(xout.shape[1])

            xout = self.dac.run(xout)
            xout = self.iqmodem.run(xout, phasor)
            xout = self.amplifier.run(xout)

            if self.delay != 0:
                xout = delay(xout, [self.delay])


        return xout
