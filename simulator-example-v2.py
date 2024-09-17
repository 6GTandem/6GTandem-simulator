import os

import matlab.engine
import numpy as np
import yaml

from utils import *
from simulator.simulator import Simulator as sim


CONFIG_FILE = "example.yml"

sim = sim()
sim.load_env(CONFIG_FILE)
sim.load_config(CONFIG_FILE)
sim.plot_env()

