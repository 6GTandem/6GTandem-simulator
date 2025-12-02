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
        wf: Waveform,
        boost_amp: Amplifier | None = None,
        antenna_amp: Amplifier | None = None,
        coup_in: Coupler | None = None,
        coup_out: Coupler | None = None,
        splitter: Splitter | None = None,
        combiner: Combiner | None = None,
        pshift: PhaseShifter | None = None,
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
        if boost_amp is None:
            self.boost_amp = Amplifier()
        else:
            self.boost_amp = boost_amp

        if antenna_amp is None:
            self.antenna_amp = Amplifier()
        else:
            self.antenna_amp = antenna_amp

        if coup_in is None:
            self.coupler_in = Coupler(wf)
        else:
            self.coupler_in = coup_in

        if coup_out is None:
            self.coupler_out = Coupler(wf)
        else:
            self.coupler_out = coup_out

        if splitter is None:
            self.splitter = Splitter(num_splits=4)
        else:
            self.splitter = splitter

        if combiner is None:
            self.combiner = Combiner()
        else:
            self.combiner = combiner

        if pshift is None:
            self.phase_shifter = PhaseShifter(num_shifters=4, resolution=2)
        else:
            self.phase_shifter = pshift

        self.x = x
        self.y = y
        self.z = z

        super().__init__(*args, **kwargs)

    def transmit(self, idata, shifts: list[int]):
        # From the input coupler to the antennas.
        c1data = self.coupler_in.run(idata)
        sdata = self.splitter.run(c1data)
        psdata = self.phase_shifter.run(sdata, shifts)
        adata = self.antenna_amp.run(psdata)
        imdata = [c1data, sdata, psdata, adata]
        odata = adata

        return odata, imdata

    def receive(self, idata, shifts: list[int]):
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
