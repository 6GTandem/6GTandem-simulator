# if any troubles with the engine, please consult https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html and our readme
import os

import matlab.engine
import numpy as np
import yaml

from utils import *
from plotter import plotter


config_file = "example.yml"

eng = matlab.engine.start_matlab() # ensure MATLAB is on


dir_path = os.path.dirname(os.path.realpath(__file__))

config_path = os.path.join(dir_path,"configurations")

# Read the YAML file
with open(os.path.join(dir_path, "configurations", config_file), "r", encoding="utf8") as file:
    data = yaml.safe_load(file)


plotter.plot_room(data)


if __name__ == "__main__":
    eng.cd(dir_path, nargout=0)
    # To add a path to all subfolders
    s = eng.genpath(os.path.join(dir_path, "sub-THz-link-model"))
    eng.addpath(s, nargout=0)

    # Define a radiostripe with a transmitter and nolinks links/repeaters.
    nolinks = 5
    RS = eng.c_radiostripe(nolinks)

    # generate a QAM signal, and use RRC pulse shaping
    N = 1e4
    OS = 5
    # number of symbols, and oversampling factor
    x = eng.Usefulfunctions.randconst(N, 1)
    x2 = eng.Usefulfunctions.pulseshape(x, matlab.double(OS), 0.1)

    # t = np.linspace(0,25, 1/OS)
    # x = eng.sin(np.pi*t*(1-0.1))

    # Run signal over stripe. The input is a vector N*1. The output is a matrix
    # N*(nolinks+1), where the first column is the output signal of the
    # transmitter, and the rest of the columns are the outputs of each link.
    Y = np.asarray(eng.run(RS, x2))

    spec(Y, plot=True)

    eng.calibrate(RS, x2, matlab.double(0, is_complex=True), nargout=0)

    Y2 = np.asarray(eng.run(RS, x2))

    spec(Y2, plot=True)
