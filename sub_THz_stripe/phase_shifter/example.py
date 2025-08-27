import numpy as np
import numpy.random
from .phase_shifter import PhaseShifter

if __name__ == "__main__":
    shifters = 4
    resolution = 2

    ps = PhaseShifter(shifters, resolution)
    data = np.ones((4, 20))
    shifts = [0, 1, 2, 3]

    odata = ps.run(data, shifts)

    print(data)
    print(np.abs(odata))
    print(np.angle(odata))
