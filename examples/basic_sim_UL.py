# Put sugar in Gilles coffee.
import os

import numpy as np
import yaml
from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.central_unit.central_unit import CentralUnit

from utils import spec
from plotter import plotter


config_file = "example.yml"

dir_path = os.path.dirname(os.path.realpath(__file__))

parent_dir = os.path.dirname(dir_path)

config_path = os.path.join(parent_dir, "configurations")

# Read the YAML file
with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
    data = yaml.safe_load(file)


def getGaussianSymbols(K=1, Ndata=5000, p=1):
    """Doc
    :param Ndata: number of symbols to generate
    :param p: signal variance
    :return: K x Ndata: symbols sampled from a complex gaussian with variance p

    Note that this generates a variable which is drawn from a complex gaussian distribution with variance p
    which is equivalent to a + bj with a and b sampled from a gaussian distribution with variance p/2.
    Here we first sample a and b from a gaussian with mean 0 and variance 1, by multiplying with sqrt(p)/sqrt(2)
    we obtain variance p/2 for both a and b, given that var(constant * X) = constant^2 var(X)
    """
    s = np.sqrt(p) / np.sqrt(2) * (np.random.randn(K, Ndata) +
                                   1j * np.random.randn(K, Ndata))
    return s.astype(np.complex64)


if __name__ == "__main__":
    # Define a radiostripe with a transmitter and nolinks links/repeaters. todo this should be loaded from config file
    nolinks = 5
    shifts = [0, 0, 0, 0]
    rs = RadioStripe(nolinks, nolinks, active_unit=4)

    nr_samples = 1000
    waveform = "Gaussian-ideal"
    ue = CentralUnit(nr_samples, waveform)
    x = ue.run()

    # Split x in 4 equal streams.
    x = x / 4
    x = np.tile(x, (4, 1))

    print(f'x shape: {x.shape}')
    print(f'{np.mean(np.abs(x**2))}')


    Y = []
    Y.append(x[0])

    # Run signal over stripe. The input is a vector N*4.
    for cdata in rs.receive(x, shifts):
        Y.append(cdata[0])

    spec(np.array(Y), plot=True)