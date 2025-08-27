import numpy as np

from ..radio_unit.radio_unit import RadioUnit
from ..fiber.fiber import Fiber
from ..component.component import Component
from ..utils import db_to_magnitude, getdbm


class RadioStripe(Component):
    """Class representing a radio stripe.

    This class implements a radiostripe in the 6GTANDEM project.
    The stripe consists of a transmitter component, and a set of links.
    Usage:
        >>> rs = RadioStripe()
        >>> y = rs.run(x) # Y is a matrix

    Current assumptions:
        - Only one radio unit is active at a time.
        - All radiostripes have one common central unit. This is instantiated outside of this class.
        - All components can have different settings but are initialized with the same parameters if none are given.

    RadioStripe:
        Fiber -> RadioUnit -> Fiber -> RadioUnit -> ...
    """

    def __init__(self, radio_units: int | list[RadioUnit] = 3, fibers: int | list[Fiber] = 3, active_unit: int = 0, *args, **kwargs):
        """
        :param : .
        """
        self._active_unit = 0

        self.radio_units = []
        if type(radio_units) is int:
            for n in range(radio_units):
                self.radio_units.append(RadioUnit())
        elif type(radio_units) is list:
            self.radio_units = radio_units

        self.fibers = []
        if type(fibers) is int:
            for n in range(fibers):
                self.fibers.append(Fiber())
        elif type(fibers) is list:
            self.fibers = fibers

        assert len(self.radio_units) == len(
            self.fibers), "Amount of radio units/fibers don't match."

        self.active_unit = active_unit

        super().__init__(*args, **kwargs)

    @property
    def active_unit(self):
        return self._active_unit

    @active_unit.setter
    def active_unit(self, nactive_unit: int):
        """Set a new radio unit to be the active one being used.

        :param nactive_unit: 0-based index for the new active unit. [0, len(radio_units)[
        """
        self._active_unit = nactive_unit

    def transmit(self, x: np.ndarray, shifts: list[int]):
        """IQ-data coming from the central unit is passed through the radio stripe and transmitted by the active RU.

        :param x: IQ-data in the form of a 1 dimensional array.
        :param shifts: See `PhaseShifter`.

        :returns: A (n x m) matrix with the amount of rows n equal to the amount of splits (See `Splitter`).
        """
        # Data comes from the central unit and first passes through the chain of RUs.
        y = x
        for ru, fib in zip(self.radio_units[:self.active_unit], self.fibers[:self.active_unit]):
            y = ru.boost(y)
            y = fib.run(y)
            yield y

        # Data transmitted by the active radio unit.
        y = self.radio_units[self.active_unit].transmit(y, shifts)

        yield y

    def receive(self, x: np.ndarray, shifts: list[int]):
        """Takes incoming IQ-data on the antennas and runs it along the stripe towards the central unit.

        :param x: IQ-data in the form of an (n x m) array with n the amount of rows being equal to the amount of splits.
        :param shifts: See `PhaseShifter`.

        :returns: A single 1-dimensional IQ-data array.
        """
        # Data is received by the active radio unit.
        y = self.radio_units[self.active_unit].receive(x, shifts)

        # Data passes through the radio units between the active unit and the central unit.
        for ru in self.radio_units[:self.active_unit][::-1]:
            y = ru.boost(y)

        return y

    def calibrate(self, x, desired_amplifier_dbm):
        """This function sets the small-signal gain of the link amplifiers.

        To give an approximate constant power(DesiredAmplifierDBM) at the output of each link.
        The calibration is valid for a given input signal, and recalibration must be performed if
        the signal statistics changes.
        """
        z = self.radio_units[0].boost(x)
        for ru in self.radio_units:
            for j in range(3):
                z2 = ru.boost(z)
                scale = db_to_magnitude(
                    desired_amplifier_dbm - getdbm(z2))
                ru.amp.gain = ru.amp.gain * scale

            z = ru.boost(z)

    @classmethod
    def from_configuration(cls, config):
        """Construct a RadioStripe from a configuration dictionary.

        :param config_file: Dictionary containing all the necessary configuration information.
        """
        bw = config["sub_THz"]["bw"]
        rs = cls()

        # Configure the amplifier
        amp_cfg = config["sub_THz"]["amplifier"]
        rs.transmitter.amplifier.set_maximum_output_power(amp_cfg["max_power"])
        rs.transmitter.amplifier.mode = amp_cfg["mode"]
        rs.transmitter.amplifier.gain = amp_cfg["gain"]
        rs.transmitter.amplifier.set_noise_var(300, bw, amp_cfg["noise_var"])
        rs.transmitter.amplifier.smoothness = amp_cfg["smoothness"]

        # Configure the limk amplifiers

        # Configure the fiber
        fiber_cfg = config["sub_THz"]["fiber"]
        # Configure the coupler
        coupler_cfg = config["sub_THz"]["coupler"]
        for link in rs.links:
            link.coupler_in.damping = coupler_cfg["damping"]
            link.coupler_in.mode = coupler_cfg["mode"]
            link.coupler_out.damping = coupler_cfg["damping"]
            link.coupler_out.mode = coupler_cfg["mode"]

            link.fiber.damping_per_meter = fiber_cfg["damping_per_meter"]
            link.fiber.filter = fiber_cfg["filter"]

        # Configure the dac
        dac_cfg = config["sub_THz"]["dac"]
        rs.transmitter.dac.trunc_level = dac_cfg["trunc_level"]
        rs.transmitter.dac.nobits = dac_cfg["num_bits"]

        # Configure the iqmodem
        iqmodem_cfg = config["sub_THz"]["iqmodem"]
        rs.transmitter.iqmodem.mode = iqmodem_cfg["mode"]
        rs.transmitter.iqmodem.iqi_coef = iqmodem_cfg["iqi_coef"]
        rs.transmitter.iqmodem.iqi_filter = iqmodem_cfg["iqi_num_filter_taps"]
        rs.transmitter.iqmodem.iqi_delay_imbalance = iqmodem_cfg["iqi_delay_imbalance"]
        rs.transmitter.iqmodem.dc_offset = iqmodem_cfg["dc_offset"]

        # Configure the oscillator
        osc_cfg = config["sub_THz"]["oscillator"]
        rs.transmitter.oscillator.cfo = osc_cfg["cfo"]
