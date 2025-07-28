from .oscillator import Oscillator
from .iqmodem import IQModem
from .amplifier import Amplifier
from .dac import Dac
from .component import Component
from .utils import delay


class Transmitter(Component):
    """Definition of a transmitter

    Components: c_oscillator, c_iqmod, c_pa, c_dac, c_antenna
    Usage:
        >>> tx = Transmitter()
        >>> y = tx.run(x)
    """

    def __init__(self, osc: Oscillator | None = None, iqmodem: IQModem | None = None, amp: Amplifier | None = None,
                 dac: Dac | None = None, delay: float = 0, *args, **kwargs):
        """Instantiate an transmitter based on an oscillator, iqmodem, amplifier and dac.
        """
        if osc is not None:
            self.oscillator = Oscillator()
        if iqmodem is not None:
            self.iqmodem = IQModem()
        if amp is not None:
            self.amplifier = Amplifier()
        if dac is not None:
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
            xout = delay(xout, self.delay)

        return xout
