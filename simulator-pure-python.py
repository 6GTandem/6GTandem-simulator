# if any troubles with the engine, please consult https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html and our readme
import os

import numpy as np
import yaml
import models
import models.utils

from utils import spec
from plotter import plotter


config_file = "example.yml"

dir_path = os.path.dirname(os.path.realpath(__file__))

config_path = os.path.join(dir_path, "configurations")

# Read the YAML file
with open(os.path.join(dir_path, "configurations", config_file), "r", encoding="utf8") as file:
    data = yaml.safe_load(file)


plotter.plot_room(data)


if __name__ == "__main__":
    # Define a radiostripe with a transmitter and nolinks links/repeaters.
    nolinks = 5
    rs = models.RadioStripe(nolinks)

    # generate a QAM signal, and use RRC pulse shaping
    n = 1e4
    os = 5
    # number of symbols, and oversampling factor
    x = models.utils.randconst(n, 1)
    x2 = models.utils.pulseshape(x, os, 0.1)

    # t = np.linspace(0,25, 1/OS)
    # x = eng.sin(np.pi*t*(1-0.1))

    # Run signal over stripe. The input is a vector N*1. The output is a matrix
    # N*(nolinks+1), where the first column is the output signal of the
    # transmitter, and the rest of the columns are the outputs of each link.
    Y = np.asarray(rs.run(x2))

    spec(Y, plot=True)

    rs.calibrate(x2, 0)

    Y2 = np.asarray(rs.run(x2))

    spec(Y2, plot=True)
