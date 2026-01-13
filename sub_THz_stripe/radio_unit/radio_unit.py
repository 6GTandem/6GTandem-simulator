import numpy as np

from ..coupler.coupler import Coupler
from ..splitter.splitter import Splitter
from ..combiner.combiner import Combiner
from ..component.component import Component
from ..amplifier.amplifier import Amplifier
from ..phase_shifter.phase_shifter import PhaseShifter
from wireless_channel.waveforms import Waveform


class RadioUnit(Component):
    def __init__(
        self,
        x,
        y,
        z,
        boost_amp: Amplifier,
        antenna_amp: Amplifier,
        coup_in: Coupler,
        coup_out: Coupler,
        splitter: Splitter,
        combiner: Combiner,
        pshift: PhaseShifter,
        *args,
        **kwargs,
    ):
        """Instantiate a Radio Unit.

        Radio unit:
        Description here.
        By default the radio unit is in receiving mode.
        In a stripe only one radio unit is active at once. All the other radio units operate in booster mode and just propagate the signal over the fiber with an amplification. Only when the switch within the unit is active the signal propagates to the antennas.

        Acting as a booster:
        --------------------
        Coupler In -> Amplifier -> Coupler out

        Operating as radio unit:
        ------------------------
        Coupler In -> Splitter -> Phase shifter -> Amplifier -> Antenna

        """
        # Instantiate the components if they are not given.
        self.boost_amp = boost_amp
        self.antenna_amp = antenna_amp
        self.coupler_in = coup_in
        self.coupler_out = coup_out
        self.splitter = splitter
        self.combiner = combiner
        self.phase_shifter = pshift

        self.x = x
        self.y = y
        self.z = z

        super().__init__(*args, **kwargs)

    def transmit(self, idata, shifts: list[int] | np.ndarray):
        # From the input coupler to the antennas.
        c1data = self.coupler_in.run(idata)
        sdata = self.splitter.run(c1data)
        psdata = self.phase_shifter.run(sdata, shifts)
        adata = self.antenna_amp.run(psdata)
        imdata = [c1data, sdata, psdata, adata]
        odata = adata

        return odata, imdata

    def receive(self, idata, shifts: list[int] | np.ndarray):
        # From the antennas to the input coupler.
        adata = self.antenna_amp.run(idata)
        psdata = self.phase_shifter.run(adata, shifts)
        cdata = self.combiner.run(psdata)
        odata = self.coupler_in.run(cdata)
        imdata = [adata, psdata, cdata, odata]

        return odata, imdata

    def boost(self, idata):
        # From the input coupler to the output coupler.
        c1data = self.coupler_in.run(idata)
        adata = self.boost_amp.run(c1data)
        odata = self.coupler_out.run(adata)
        imdata = [c1data, adata, odata]

        return odata, imdata

    def __str__(self):
        """Human-readable summary of the Radio Unit."""
        info = (
            f"RadioUnit @ ({self.x}, {self.y}, {self.z})\n"
            f"  Amplifier:      {self.boost_amp.__class__.__name__}\n"
            f"  Coupler In:     {self.coupler_in.__class__.__name__}\n"
            f"  Coupler Out:    {self.coupler_out.__class__.__name__}\n"
            f"  Splitter:       {self.splitter.__class__.__name__} "
            f"(outputs={getattr(self.splitter, 'n_outputs', 'N/A')})\n"
            f"  Combiner:       {self.combiner.__class__.__name__}\n"
            f"  Phase Shifter:  {self.phase_shifter.__class__.__name__} "
            f"(n={getattr(self.phase_shifter, 'n_elements', 'N/A')}, "
            f"bits={getattr(self.phase_shifter, 'n_bits', 'N/A')})\n"
        )
        return info
