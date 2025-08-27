import numpy as np
import numpy.random
from .splitter import Splitter

if __name__ == "__main__":
    splits = 4

    sp = Splitter(splits)
    data = np.arange(1, 11)

    odata = sp.run(data)

    print(data)
    print(odata)
