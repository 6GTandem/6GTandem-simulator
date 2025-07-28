from ..coupler.coupler import Coupler
from ..component.component import Component
from ..amplifier.amplifier import Amplifier


class RadioUnit(Component):
    def __init__(self, amp: Amplifier | None = None, coup_in: Coupler | None = None,
                 coup_out: Coupler | None = None, mode: str = "rx", *args, **kwargs):
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
        self.mode = mode

        # Instantiate the components if they are not given.
        if amp is not None:
            self.amp = Amplifier()
        if coup_in is not None:
            self.coupler_in = Coupler()
        if coup_out is not None:
            self.coupler_out = Coupler()

        super().__init__(*args, **kwargs)

        raise NotImplementedError(
            "Phase shifters and antenna not implemented.")

    def run(self, idata):
        match self.mode:
            case "tx":
                # From the input coupler to the antennas.
                c1data = self.coupler_in.run(idata)
                # TODO: Add splitter -> phase shifter -> PA -> Antenna
            case "rx":
                # From the antennas to the input coupler.
                # TODO: Add Antenna -> PA -> phase shifter -> splitter
                odata = self.coupler_out.run(idata)
            case "boost":
                # From the input coupler to the output coupler.
                c1data = self.coupler_in.run(idata)
                adata = self.amp.run(c1data)
                odata = self.coupler_out.run(adata)

        return odata
