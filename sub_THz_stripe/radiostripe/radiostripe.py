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
            **kwargs):
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

        assert len(self.radio_units) == len(
            self.fibers
        ), "Amount of radio units/fibers don't match."

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
        assert nactive_unit >= 0 and nactive_unit < len(self.radio_units), f"Radio unit {nactive_unit} does not exist \
                                                                             . [0, {len(self.radio_units)}["
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
        for i, (ru, fib) in enumerate(zip(self.radio_units[: self.active_unit+1], self.fibers[: self.active_unit+1])):
            y = fib.run(y)
            imdata.append(y)
            if i == self.active_unit:
                y, imd = ru.transmit(y, shifts)
                imdata.extend(imd)
            else:
                y, imd = ru.boost(y)
                imdata.extend(imd)

        return y, imdata

    def receive(self, x: np.ndarray, shifts: list[int]):
        """Receives incoming IQ-data on the antennas of the acitve unit and runs it along the stripe.

        :param x: IQ-data in the form of an (n x m) array with n the amount of rows being equal to the amount of splits.
        :param shifts: See `PhaseShifter`.

        :returns: A 1xm IQ-data array arriving at the central unit.
        """
        # Data is received by the active radio unit.
        y = self.radio_units[self.active_unit].receive(x, shifts)
        y = self.fibers[self.active_unit].run(y)
        yield y

        # Data passes through the radio units between the active unit and the central unit.
        for ru, fib in zip(
            self.radio_units[: self.active_unit - 1][::-1],
            self.fibers[: self.active_unit - 1][::-1],
        ):
            y = ru.boost(y)
            y = fib.run(y)
            yield y

    def receive_all(self, x: list[np.ndarray], shifts: list[int]):
        """Takes incoming IQ-data on the antennas and runs it along the stripe towards the central unit.

        :param x: List containing all the IQ data being received on all the radio units. The data at x[0] is the unit
                  closest to the central unit.
        :param shifts: See `PhaseShifter`.

        :returns: IQ-data array arriving at the end of the stripe.
        """
        y_in = np.zeros(x[0].shape, dtype=np.complex128)
        for ru, fib, data in zip(self.radio_units[::-1], self.fibers[::-1], x):
            rdata = ru.receive(data, shifts)
            y_in = ru.boost(y_in)
            y_combined = np.sum((rdata, y_in))
            y_out = fib.run(y_combined)
            y_in = y_out # becomes new in
            yield y_out

    def calibrate(self, x, desired_signal_db):
        """This function sets the small-signal gain of the link amplifiers.

        To give an approximate constant power(DesiredAmplifierDBM) at the output of each link.
        The calibration is valid for a given input signal, and recalibration must be performed if
        the signal statistics changes.
        """
        for ru, fib in zip(self.radio_units, self.fibers):
            z2 = x
            for j in range(3):
                z, _ = ru.boost(x)
                z2 = fib.run(z)
                current_amp_db = getdbm(z2)
                scale = db_to_magnitude(desired_signal_db - current_amp_db)
                ru.amp.gain = ru.amp.gain * scale
            
            z = z2

    @classmethod
    def from_config_locations(
            cls,
            stripe_config: list[dict[str, dict[str, float]]],
            component_config: dict[str, dict[str, Any]],
            wf: Waveform):
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
                    if 'model' in coupler_config:
                        # Load the couplers S-parameter file
                        base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                        coupler_spars = rf.Network(os.path.join(base_path, f'../models/{coupler_config['model']}'))

                        # Subcarrier frequencies
                        ofdm_freqs = np.linspace(wf.fc - wf.bw/2, wf.fc + wf.bw/2, wf.n_carriers * wf.oversampling_factor)

                        # Extract frequency and S21 (transmission)
                        coup_freqs = coupler_spars.f
                        coup_s21 = coupler_spars.s[:, 1, 0]  # S21

                        # Interpolate magnitude and phase separately for better accuracy
                        coup_s21_mag = np.abs(coup_s21)
                        coup_s21_phase = np.angle(coup_s21)

                        interp_mag = np.interp(ofdm_freqs, coup_freqs, coup_s21_mag)
                        interp_phase = np.interp(ofdm_freqs, coup_freqs, coup_s21_phase)
                        coup_s21_ofdm = interp_mag * np.exp(1j * interp_phase)

                        damping = 0  # in dB
                        filter_mode = 'freq_domain'
                        coup_in = Coupler(damping=damping, filter=coup_s21_ofdm, filter_mode=filter_mode, wf=wf)
                        coup_out = Coupler(damping=damping, filter=coup_s21_ofdm, filter_mode=filter_mode, wf=wf)
                    else:
                        coup_in = Coupler(**coupler_config)
                        coup_out = Coupler(**coupler_config)

                loc = unit_cfg["radio_unit"]
                ru = RadioUnit(
                    x=loc.get("x", 0),
                    y=loc.get("y", 0),
                    z=loc.get("z", 0),
                    wf=wf,
                    amp=amp, coup_in=coup_in, coup_out=coup_out)
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
                if 'model' in fiber_config:
                    # Load the fiber model from the given file.
                    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                    fiber_spars = pd.read_csv(os.path.join(base_path, f'../models/{fiber_config['model']}'))

                    # Subcarrier frequencies
                    ofdm_freqs = np.linspace(wf.fc - wf.bw/2, wf.fc + wf.bw/2, wf.n_carriers * wf.oversampling_factor)

                    fib_freqs = fiber_spars["freq[Hz]"]
                    fib_phase = np.unwrap(np.deg2rad(fiber_spars["ang:Trc2_S21"]))
                    fib_mag = 10 ** (fiber_spars["db:Trc2_S21"] / 20.0)
                    
                    interp_mag = np.interp(ofdm_freqs, fib_freqs, fib_mag)
                    interp_phase = np.interp(ofdm_freqs, fib_freqs, fib_phase)
                    fib_s21_ofdm = interp_mag * np.exp(1j * interp_phase)

                    filter_mode = 'freq_domain'
                    fiber = Fiber(damping_per_meter=0, length=0, filter=fib_s21_ofdm, filter_mode=filter_mode, wf=wf)
                else:
                    fiber = Fiber(**fiber_config)
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
