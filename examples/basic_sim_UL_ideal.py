import os
import sys
# Add project root to sys.path for local imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
###########################################
# DO NOT MOVE ANY IMPORTS ABOVE THIS LINE #
###########################################

# todo add a RU to act as UE transmitter (splitter, phase shifters, PAs, disable coupler in!!)

# todo transmit using RU

# todo send over channel

# todo test receive

# todo test receive all

from matplotlib import pyplot as plt
import numpy as np
import logging
import yaml

from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.central_unit.central_unit import CentralUnit
from sub_THz_stripe.radio_unit.radio_unit import RadioUnit
from wireless_channel.subTHz_channel import Channel
from wireless_channel.waveforms import Waveform
from utils import setup_logging
from plotter import plotter

booster_stages = ["fiber", "coupler", "amplifier", "coupler"]
tx_stages = ["fiber", "coupler", "splitter", "shifter", "amplifier"]


if __name__ == "__main__":
    # Logfile name the same as the current script name.
    script_name = __file__.split(".")[0]
    logfile = f"{script_name}.log"

    setup_logging(logfile)
    logger = logging.getLogger("6GTandemBasicDL")

    # read config file
    config_file = "office_config.yml"
    dir_path = os.path.dirname(os.path.realpath(__file__))
    config_path = os.path.join(dir_path, "..", "configurations")
    with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
        config = yaml.safe_load(file)

    # todo load from config
    nr_antennas = 4

    # plot full room with all stripes and all possible ue locations
    plotter.plot_room(config)
    logger.debug("%d stripes in the room", len(config['radio_stripes']))

    # construct waveform class
    waveform_config = config["waveform_config"]
    freq_band_config = config["sub_thz"]
    wf = Waveform.from_config(waveform_config, freq_band_config)
    logger.debug("%s", wf)

    # generate ofdm waveform in time domain
    bits = wf.generate_bits()
    qam = wf.qam_modulate()
    ofdm_time = wf.ofdm_modulate()  # shape: nr_ofdm_symb x (fftsize + cp length)
    wf.plot_psd(ofdm_time)

    # build all radio stripes
    stripes = []
    component_config = config["component_config"]
    for stripe_cfg in config["radio_stripes"]:
        stripes.append(RadioStripe.from_config_locations(stripe_cfg, component_config, wf))

    # continue with 3 stripes, separated by 1m
    stripes = [stripes[5]]  # stripes[5:11:2]
    active_ru_idxes = [4]  # , 4, 6
    plotter.plot_stripes(config, stripes)
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug("stripe: %d: %s", stripe_idx, stripe)

    # load channels
    ue_pos = config["ue_positions"][0]
    channel = Channel.from_sionna(ue_pos, debug=True)
    logger.debug("%s", channel)

    # We use a single radio unit to represent the UE.
    ue = RadioUnit(ue_pos['x'], ue_pos['y'], ue_pos['z'], wf)
    logger.debug('ue RU: %s', ue)
    shifts = [0, 0, 0, 0]
    iq_data_tx, imdata = ue.transmit(ofdm_time, shifts)

    wf.plot_iq_time(iq_data_tx[0], title="Transmitted at UE")

    iq_data_rus = channel.transmit_ul(iq_data_tx, wf)
    wf.plot_iq_time(iq_data_rus[0][0][0], title="After wireless channel")

    active_ru_idxes = [2]  # , 4, 6
    logger.debug('Receiving data on all stripes.')
    iq_at_last_rus = []
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug('stripe: %d - active ru %d', stripe_idx, active_ru_idxes[stripe_idx])
        stripe.active_unit = active_ru_idxes[stripe_idx]

        logger.debug('shape of iq data: %s', iq_data_rus.shape)
        phase_shifts = [0, 0, 0, 0]
        iq_out, imdata = stripe.receive_all(iq_data_rus, phase_shifts)

        iq_out_reshaped = iq_out.reshape(nr_antennas, wf.n_ofdm_symbols, -1)
        logger.debug('reshaped after stripe: %s', iq_out_reshaped.shape)
        iq_at_last_rus.append(iq_out_reshaped)

    cu = CentralUnit()  # todo are these configs loadable?
    logger.debug(f"CU: {cu}")
    ofdm_time_after_cu = cu.run(iq_at_last_rus)  # shape: nr_ofdm_symbols x (fft_size + cp length)
    logger.debug("shape of ofdm timee: %s", ofdm_time_after_cu.shape)
    logger.debug("np alike: %s", np.allclose(ofdm_time, ofdm_time_after_cu))

    y_combined_freq = wf.ofdm_time_to_freq(ofdm_time_after_cu)
    # # todo do we need equalization?

    y_qam = wf.ofdm_to_qam(y_combined_freq)
    # wf.plot_constellation(y_qam, symbols_tx=qam, title="RX'ed symbols")

    y_bits = wf.qam_to_bits(y_qam)

    ber = wf.compute_ber(bits, y_bits)
    # logger.debug('BER: %f', ber)
    plt.show()
