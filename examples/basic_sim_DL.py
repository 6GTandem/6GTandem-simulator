import os
import sys
import numpy as np
import yaml
from wireless_channel.waveforms import Waveform

# Add project root to sys.path for utils import
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.central_unit.central_unit import CentralUnit
from wireless_channel.subTHz_channel import Channel
from utils import spec
from plotter import plotter


if __name__ == "__main__":
    # read config file
    config_file = "office_config.yml"
    dir_path = os.path.dirname(os.path.realpath(__file__))
    config_path = os.path.join(dir_path, "..", "configurations")
    with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
        config = yaml.safe_load(file)

    # plot full room with all stripes and all possible ue locations
    plotter.plot_room(config)
    print(f"{len(config['radio_stripes'])} stripes in the room")

    # build all radio stripes
    stripes = []
    for stripe_cfg in config["radio_stripes"]:
        stripes.append(RadioStripe.from_config_locations(stripe_cfg))

    # continue with 3 stripes, separated by 1m
    stripes = stripes[5:11:2]
    plotter.plot_stripes(config, stripes)

    # construct waveform class
    waveform_config = config["waveform_config"]
    freq_band_config = config['sub_thz']
    wf = Waveform.from_config(waveform_config, freq_band_config)
    print(wf)

    # generate ofdm waveform in time domain
    bits = wf.generate_bits()
    qam = wf.qam_modulate()
    ofdm_time = wf.ofdm_modulate()
    wf.plot_psd(ofdm_time)

    # load channels
    # select first ue position:
    # todo this might be a bit backward, you might first need to interact with the
    # channels dataset meta data to figure out which UE you want... but ok for now
    ue_pos = config["ue_positions"][0]
    channel = Channel.from_sionna(ue_pos)
    print(channel)
    # todo this loads the full channel to all RUs and Stripes, we might need a function
    # to only select the ones we need

    # todo build CU
    cu = CentralUnit() # todo are these configs loadable?
    ofdm_time_after_cu = cu.run(ofdm_time) # shape: nr_ofdm_symbols x fft_size
    print(f' shape of ofdm timee: {ofdm_time_after_cu.shape}')
    print(f'np alike: {np.allclose(ofdm_time, ofdm_time_after_cu)}')
    # todo check with impairments if something changes

    # todo send over stripe

    # todo send over channel
    # todo when sent over stripe can be a list with multiple signals form multiple stripe
    # also now only 1 antenna
    active_ru_per_stripe = [4] # todo get from config (now just dummy get forth RU)
    channel.transmit_dl([ofdm_time_after_cu], active_ru_per_stripe, wf)
    # todo debug this

    # todo combine at UE (phase shifters, combiners?)

