from .amplifier import Amplifier
from .fiber import Fiber
from .coupler import Coupler
from .component import Component


class Link(Component):
    def __init__(self, amp: Amplifier | None = None, fiber: Fiber | None = None, coup_in: Coupler | None = None,
                 coup_out: Coupler | None = None, *args, **kwargs):
        """Instantiate a link component.

        Coupler In -> Amplifier -> Fiber -> Coupler out
        """
        if amp is not None:
            self.amp = Amplifier()
        if fiber is not None:
            self.fiber = Fiber()
        if coup_in is not None:
            self.coupler_in = Coupler()
        if coup_out is not None:
            self.coupler_out = Coupler()

        super().__init__(*args, **kwargs)

    def run(self, y):
        z = self.fiber.run(y)
        x1 = self.coupler_in.run(z)
        x = self.amp.run(x1)

        return self.coupler_out.run(x)
