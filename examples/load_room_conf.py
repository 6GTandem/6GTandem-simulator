import os
import sys
import numpy as np
import yaml

# Add project root to sys.path for utils import
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from utils import spec
from plotter import plotter

config_file = "office_config.yml"

dir_path = os.path.dirname(os.path.realpath(__file__))

config_path = os.path.join(dir_path, "..", "configurations")

# Read the YAML file
with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
    data = yaml.safe_load(file)

plotter.plot_room(data)
