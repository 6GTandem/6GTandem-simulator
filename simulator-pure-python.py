# if any troubles with the engine, please consult https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html and our readme
import os

import numpy as np
import yaml
from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.central_unit.central_unit import CentralUnit
from sub_THz_stripe import utils

from utils import spec
from plotter import plotter


config_file = "example.yml"

dir_path = os.path.dirname(os.path.realpath(__file__))

config_path = os.path.join(dir_path, "configurations")

# Read the YAML file
with open(os.path.join(dir_path, "configurations", config_file), "r", encoding="utf8") as file:
    data = yaml.safe_load(file)


# plotter.plot_room(data)


def getGaussianSymbols(K=1, Q=1, Ndata=5000, p=1):
    """Doc

    :param K: number of users
    :param Q: number of subcarriers
    :param Ndata: number of symbols to generate
    :param p: signal variance
    :return: K x Q x Ndata: symbols sampled from a complex gaussian with variance p

    Note that this generates a variable which is drawn from a complex gaussian distribution with variance p
    which is equivalent to a + bj with a and b sampled from a gaussian distribution with variance p/2.
    Here we first sample a and b from a gaussian with mean 0 and variance 1, by multiplying with sqrt(p)/sqrt(2)
    we obtain variance p/2 for both a and b, given that var(constant * X) = constant^2 var(X)
    """
    s = np.sqrt(p) / np.sqrt(2) * (np.random.randn(K, Q, Ndata) +
                                   1j * np.random.randn(K, Q, Ndata))
    return s.astype(np.complex64)


if __name__ == "__main__":
    # Define a radiostripe with a transmitter and nolinks links/repeaters.
    nolinks = 5
    shifts = [0, 0, 0, 0]
    rs = RadioStripe(nolinks, nolinks, active_unit=4)

    # generate a QAM signal, and use RRC pulse shaping
    n = 10000
    os = 5
    # number of symbols, and oversampling factor
    x, c = utils.randconst(1, n)
    x2 = utils.pulseshape(x, os, 0.1)

    Y = []
    Y.append(x2[0])

    # Generate random data
    # x2 = getGaussianSymbols()
    # x2 = x2[0][0]

    # t = np.linspace(0,25, 1/OS)
    # x = eng.sin(np.pi*t*(1-0.1))

    # Run signal over stripe. The input is a vector N*1.
    for cdata in rs.transmit(x2, shifts):
        Y.append(cdata[0])

    spec(np.array(Y[:-1]), plot=True)

    rs.calibrate(x2, 0)

    Y2 = rs.transmit(x2, shifts)

    spec(Y2, plot=True)
