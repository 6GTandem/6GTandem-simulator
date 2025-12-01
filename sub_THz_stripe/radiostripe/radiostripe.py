import os
import numpy as np
import pandas as pd
import logging

from sub_THz_stripe.central_unit.central_unit import CentralUnit
from wireless_channel.waveforms import Waveform
from typing import Any
import skrf as rf

from ..radio_unit.radio_unit import RadioUnit
from ..fiber.fiber import Fiber
from ..amplifier.amplifier import Amplifier
from ..coupler.coupler import Coupler
from ..component.component import Component
from ..utils import db_to_magnitude, getdbm

logger = logging.getLogger(__name__)


class RadioStripe(Component):
    """Class representing a radio stripe.

    This class implements a radiostripe in the 6GTANDEM project.
    The stripe consists of a transmitter component, and a set of links.
    Usage:
        >>> rs = RadioStripe()
        >>> y = rs.run(x) # Y is a matrix

    Current assumptions:
        - Only one radio unit is active at a time.
        - All radiostripes have one common central unit, each stripe has a dedicated port on the CU. This is instantiated outside of this class.
        - All components can have different settings but are initialized with the same parameters if none are given.

    RadioStripe:
        Fiber -> RadioUnit -> Fiber -> RadioUnit -> ...
    """

    def __init__(
        self,
        radio_units: int | list[RadioUnit] = 3,
        fibers: int | list[Fiber] = 3,
        active_unit: int = 0,
        central_unit: CentralUnit | None = None,
        waveform: Waveform | None = None,
        *args,
        **kwargs,
    ):
        """
        :param : .
        """
        self._active_unit = 0
        self.central_unit = central_unit
        self.wf = waveform

        self.radio_units = []
        if type(radio_units) is int:
            for n in range(radio_units):
                self.radio_units.append(RadioUnit(n, 0, 0, self.wf))
        elif type(radio_units) is list:
            self.radio_units = radio_units

        self.fibers = []
        if type(fibers) is int:
            for n in range(fibers):
                self.fibers.append(Fiber(self.wf))
        elif type(fibers) is list:
            self.fibers = fibers

        assert len(self.radio_units) == len(self.fibers), "Amount of radio units/fibers don't match."

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
        assert nactive_unit >= 0 and nactive_unit < len(
            self.radio_units
        ), f"Radio unit {nactive_unit} does not exist. [0, {len(self.radio_units)}["
        self._active_unit = nactive_unit

    def transmit(self, x: np.ndarray, shifts: list[int]):
        """IQ-data coming from the central unit is passed through the radio stripe and transmitted by the active RU.

        :param x: IQ-data in the form (1 x n_iq_samples).
        :param shifts: See `PhaseShifter`.

        :returns: A (n x m) matrix with the amount of rows n equal to the amount of splits (See `Splitter`).
        """
        # Data comes from the central unit and first passes through the chain of RUs.
        y = x
        imdata = [x]
        for i, (ru, fib) in enumerate(
            zip(self.radio_units[: self.active_unit + 1], self.fibers[: self.active_unit + 1])
        ):
            y = fib.run(y)
            imdata.append(y)
            if i == self.active_unit:
                y, imd = ru.transmit(y, shifts)
                imdata.extend(imd)
            else:
                y, imd = ru.boost(y)
                imdata.extend(imd)

        return y, imdata

    def receive(self, x: np.ndarray, shifts: list[int], delay: bool = False):
        """Receives incoming IQ-data on the antennas of the active unit and runs it along the stripe.

        :param x: IQ-data in the form of an (n x m) array with n the amount of rows being equal to the amount of splits.
        :param shifts: See `PhaseShifter`.

        :returns: A 1xm IQ-data array arriving at the central unit.
        """
        imdata = [x]
        # Data is received by the active radio unit.
        y, im = self.radio_units[self.active_unit].receive(x[self.active_unit], shifts)
        imdata.extend(im)
        y = self.fibers[self.active_unit].run(y, delay=delay)
        imdata.append(y)

        # Data passes through the radio units between the active unit and the central unit.
        if self.active_unit != 0:
            for ru, fib in zip(
                self.radio_units[: self.active_unit - 1][::-1],
                self.fibers[: self.active_unit - 1][::-1],
            ):
                y, im = ru.boost(y)
                imdata.extend(im)
                y = fib.run(y, delay=delay)
                imdata.append(y)

        return y, imdata

    def receive_all(self, x: np.ndarray, shifts: list[int], delay: bool = False, window: str = "fixed"):
        """Takes incoming IQ-data on the antennas and runs it along the stripe towards the central unit.

        :param x: Array containing all the IQ data being received on all the radio units. The data at x[0] is the unit
                  closest to the central unit.
        :param shifts: See `PhaseShifter`.

        :returns: IQ-data array arriving at the end of the stripe.
        """
        y_in = np.zeros(x[0][0].shape, dtype=np.complex128)
        imdata = []

        for ru, fib, data in zip(self.radio_units[::-1], self.fibers[::-1], x[::-1]):
            imdata.append(data)
            rdata, im = ru.receive(data, shifts)
            imdata.extend(im)
            imdata.append(y_in)
            y_in, im = ru.boost(y_in)
            imdata.extend(im)
            # If the window is dynamic we need to append zeros to the received data to make
            # its length match the data coming from the fiber.
            if window != "fixed":
                rdata = np.pad(
                    rdata,
                    (
                        (0, 0),
                        (0, y_in.shape[1] - rdata.shape[1]),
                    ),
                )
            y_combined = np.sum((rdata, y_in), axis=0)
            y_out = fib.run(y_combined, delay=delay, window=window)
            y_in = y_out  # becomes new in

        return y_in, imdata

    def calibrate(self, x, desired_signal_dbm):
        """This function sets the small-signal gain of the link amplifiers.

        To give an approximate constant power(DesiredAmplifierDBM) at the output of each link.
        The calibration is valid for a given input signal, and recalibration must be performed if
        the signal statistics changes.
        """
        for ru, fib in zip(self.radio_units, self.fibers):
            z2 = x
            for j in range(10):
                z, _ = ru.boost(x)
                z2 = fib.run(z)
                current_amp_dbm = getdbm(z2)
                scale = db_to_magnitude(desired_signal_dbm - current_amp_dbm)
                ru.amp.gain = ru.amp.gain * scale

            z = z2

    def calibrate_losses(self, x):
        """This function sets the small-signal gain of the link amplifiers to compensate the stripe losses.

        The small signal gain of the amplifiers is set so that after all the losses from the Fiber and Couplers
        the same amplitude is achieved at the output of the RadioUnit as received on its input.
        """
        for ru, fib in zip(self.radio_units, self.fibers):
            ru.amp.gain = 1
            z, _ = ru.boost(x)
            z2 = fib.run(z)
            pavg_out = np.mean(np.abs(z2) ** 2)
            pavg_in = np.mean(np.abs(x) ** 2)
            scale = np.sqrt(pavg_in / pavg_out)
            ru.amp.gain = scale

    @classmethod
    def from_config_locations(
        cls, stripe_config: list[dict[str, dict[str, float]]], component_config: dict[str, dict[str, Any]], wf: Waveform
    ):
        """Initialize a RadioStripe from a stripe configuration containing radio unit locations.

        Only location is considered; all other parameters are default.

        :param stripe_config: List of dicts, each with a 'radio_unit' key containing 'x', 'y', 'z'.
        :param component_config: Configuration of the different components used in the RadioStripe. Components
                                 which are not specified here are given default values.

        :return: RadioStripe instance with radio units at specified locations.
        """
        radio_units = []
        units = []

        central_unit = None
        amplifier_config = component_config.get("amplifier", None)
        fiber_config = component_config.get("fiber", None)
        coupler_config = component_config.get("coupler", None)

        for unit_cfg in stripe_config:
            if "radio_unit" in unit_cfg:
                amp = None
                coup_in = None
                coup_out = None

                if amplifier_config is not None:
                    amp = Amplifier(**amplifier_config)
                if coupler_config is not None:
                    coup_in = Coupler.from_config(coupler_config, wf)
                    coup_out = Coupler.from_config(coupler_config, wf)

                loc = unit_cfg["radio_unit"]
                ru = RadioUnit(
                    x=loc.get("x", 0),
                    y=loc.get("y", 0),
                    z=loc.get("z", 0),
                    wf=wf,
                    amp=amp,
                    coup_in=coup_in,
                    coup_out=coup_out,
                )
                radio_units.append(ru)
                unit = ru
            else:
                loc = unit_cfg["central_unit"]
                central_unit = CentralUnit(x=loc.get("x", 0), y=loc.get("y", 0), z=loc.get("z", 0))
                unit = central_unit

            units.append(unit)

        fibers = []
        for i in range(len(radio_units)):
            if fiber_config is not None:
                fiber = Fiber.from_config(fiber_config, wf)
            else:
                p1 = np.array([units[i].x, units[i].y, units[i].z])
                p2 = np.array([units[i + 1].x, units[i + 1].y, units[i + 1].z])
                length = np.linalg.norm(p2 - p1)
                fiber = Fiber(length=length, wf=wf)

            fibers.append(fiber)

        logger.info(f"Constructed {len(radio_units)} radio units and {len(fibers)} fibers.")
        return cls(radio_units=radio_units, fibers=fibers, central_unit=central_unit, waveform=wf)

    def __str__(self):
        ru_strs = []
        for idx, ru in enumerate(self.radio_units):
            ru_strs.append(f"RU{idx}: ({ru.x:.2f}, {ru.y:.2f}, {ru.z:.2f})")
        fiber_strs = []
        for idx, fiber in enumerate(self.fibers):
            fiber_strs.append(f"Fiber{idx}: length={getattr(fiber, 'length', 'N/A'):.2f}m")
        return (
            f"RadioStripe(\n"
            f"  Central Unit: {self.central_unit if self.central_unit else 'None'},\n"
            f"  Active Unit Index: {self.active_unit},\n"
            f"  Radio Units:\n    " + "\n    ".join(ru_strs) + "\n"
            f"  Fibers:\n    " + "\n    ".join(fiber_strs) + "\n"
            f")"
        )
