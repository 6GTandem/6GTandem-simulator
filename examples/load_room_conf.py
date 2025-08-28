import os
import sys
import numpy as np
import yaml

# Add project root to sys.path for utils import
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from utils import spec
from plotter import plotter

config_file = "office_config.yml"

dir_path = os.path.dirname(os.path.realpath(__file__))

config_path = os.path.join(dir_path, "..", "configurations")

# Read the YAML file
with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
    config = yaml.safe_load(file)

# plotter.plot_room(config)
print(f"{len(config['radio_stripes'])} stripes in the room")

stripes = []
for stripe_cfg in config["radio_stripes"]:
    stripes.append(RadioStripe.from_config_locations(stripe_cfg))

# stripes = stripes[0:3]
plotter.plot_stripes(config, stripes)
