import numpy as np
from .radio_unit import RadioUnit

if __name__ == "__main__":
    ru = RadioUnit()
    data = np.ones(20)
    shifts = [0, 1, 2, 3]

    odata = ru.transmit(data, shifts)

    print(data)
    print(np.abs(odata))
    print(np.angle(odata))

    rxdata = ru.receive(odata, [0, 3, 2, 1])

    print(np.abs(rxdata))
    print(np.angle(rxdata))
