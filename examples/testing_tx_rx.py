import os

import numpy as np
import yaml
from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.central_unit.central_unit import CentralUnit
from sub_THz_stripe import utils
import sys
from utils import spec
from plotter import plotter

# Add project root to sys.path for utils import
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from utils import spec
from plotter import plotter
from wireless_channel.subTHz_channel import Channel


from sub_THz_stripe.oscillator.oscillator import Oscillator
from sub_THz_stripe.iqmodem.iqmodem import IQModem
from sub_THz_stripe.amplifier.amplifier import Amplifier
from sub_THz_stripe.dac.dac import Dac

if __name__ == "__main__":

    # step 1 create stripes
    config_file = "office_config.yml"

    dir_path = os.path.dirname(os.path.realpath(__file__))

    config_path = os.path.join(dir_path, "..", "configurations")

    # Read the YAML file
    with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
        config = yaml.safe_load(file)

    print(f"{len(config['radio_stripes'])} stripes in the room")

    stripes = []
    for stripe_cfg in config["radio_stripes"]:
        stripes.append(RadioStripe.from_config_locations(stripe_cfg))

    stripes = stripes[0:1]
    nr_stripes = len(stripes)
    # set active unit
    RadioStripe.active_unit = 3
    #plotter.plot_stripes(config, stripes)

    # step 2 create CU
    nr_samples = 1000 # todo add
    waveform = "Gaussian-ideal"
    os = Oscillator()
    iqmod = IQModem()
    amp = Amplifier()
    dac = Dac()
    cu = CentralUnit(os, iqmod, amp, dac, waveform=waveform)
    x = cu.run() # 1 x Nrsamples

    # step 3 create or load channel
    num_carriers = 1024
    channel = Channel("subTHz-Rayleigh", num_carriers, Nr_ue_antennas=4,
                      Nr_ru_antennas=4, Nr_rus=42, Nr_stripes=nr_stripes)

    Y = [] # list with element per ru
    Y.append(x)

    # step 4 Run signal over stripe. The input is a vector N*1.
    # todo how to deal with multiple carriers?
    shifts = [0, 0, 0, 0] # shifts for phase shifter
    for stripe in stripes:
        for cdata in stripe.transmit(x, shifts):
            Y.append(cdata)

    # todo make compatible with muliple stripes

    # setp 5 transmit over the air
    X_list = Y
    active_ru_per_stripe = [3]
    Rx_sig = channel.transmit_dl(X_list, active_ru_per_stripe)
    print(f'done')