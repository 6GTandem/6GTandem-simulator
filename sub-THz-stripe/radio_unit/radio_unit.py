from ..coupler.coupler import Coupler
from ..splitter.splitter import Splitter
from ..combiner.combiner import Combiner
from ..component.component import Component
from ..amplifier.amplifier import Amplifier
from ..phase_shifter.phase_shifter import PhaseShifter


class RadioUnit(Component):
    def __init__(self, amp: Amplifier | None = None, coup_in: Coupler | None = None,
                 coup_out: Coupler | None = None, splitter: Splitter | None = None, combiner: Combiner | None = None, pshift: PhaseShifter | None = None, state: str = "tx", *args, **kwargs):
        """Instantiate a Radio Unit.

        TODO: Add Phase shifters and antennas to the RadioUnit.

        Radio unit:
        Description here.
        By default the radio unit is in receiving mode.
        In a stripe only one radio unit is active at once. All the other radio units operate in booster mode and just propagate the signal over the fiber with an amplification. Only when the switch within the unit is active the signal propagates to the antennas.

        Acting as a booster:
        --------------------
        Coupler In -> Amplifier -> Coupler out

        Operating as radio unit:
        ------------------------
        Coupler In -> Amplifier -> Coupler out
                      Switch    -> Splitter -> Phase shifter -> Antenna

        """
        # Keep track of what mode the radio unit is in. Is it receiving or transmitting?
        self.state = state

        # Instantiate the components if they are not given.
        if amp is None:
            self.amp = Amplifier()
        if coup_in is None:
            self.coupler_in = Coupler()
        if coup_out is None:
            self.coupler_out = Coupler()
        if splitter is None:
            self.splitter = Splitter(4)
        if combiner is None:
            self.combiner = Combiner()
        if pshift is None:
            self.phase_shifter = PhaseShifter(4, 2)

        super().__init__(*args, **kwargs)

    def run(self, idata, shifts):
        match self.state:
            case "tx":
                # From the input coupler to the antennas.
                c1data = self.coupler_in.run(idata)
                sdata = self.splitter.run(c1data)
                psdata = self.phase_shifter.run(sdata, shifts)
                adata = self.amp.run(psdata)
                # TODO: Add Antenna
                odata = adata
            case "rx":
                # From the antennas to the input coupler.
                # TODO: Add Antenna
                adata = self.amp.run(idata)
                psdata = self.phase_shifter.run(adata, shifts)
                cdata = self.combiner.run(psdata)
                odata = self.coupler_in.run(cdata)
            case "boost":
                # From the input coupler to the output coupler.
                c1data = self.coupler_in.run(idata)
                adata = self.amp.run(c1data)
                odata = self.coupler_out.run(adata)

        return odata
