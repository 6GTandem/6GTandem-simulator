import numpy as np

from .transmitter import Transmitter
from .link import Link
from .component import Component
from .utils import db_to_magnitude, getdbm


class RadioStripe(Component):
    """Class representing a radio stripe.

    This class implements a radiostripe in the 6GTANDEM project.
    The stripe consists of a transmitter component, and a set of links.
    Usage:
        >>> rs = RadioStripe()
        >>> y = rs.run(x) # Y is a matrix
    """

    def __init__(self, links: int | list[Link] = 3, bandwith: float = 5e9, os=5, *args, **kwargs):
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
        self.transmitter.amplifier.set_gain(15, self.average_power)

        if type(links) is int:
            self.links = [Link() for l in range(links)]
            for link in self.links:
                link.amp.mode = '6gtandem'
                link.amp.set_maximum_output_power(self.max_power)
                link.amp.set_gain(self.average_power - link.fiber.damping -
                                link.coupler_in.damping - link.coupler_out.damping, self.average_power)
                link.amp.set_noise_var(300, self.bandwidth * self.os, 10)
        else:
            self.links = links

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
                desired_amplifier_dbm - getdbm(z))
            self.transmitter.amplifier.gain = self.transmitter.amplifier.gain * scale

        z = self.transmitter.run(x)
        for link in self.links:
            for j in range(3):
                z2 = link.run(z)
                scale = db_to_magnitude(
                    desired_amplifier_dbm - getdbm(z2))
                link.amp.gain = link.amp.gain * scale

            z = link.run(z)

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
