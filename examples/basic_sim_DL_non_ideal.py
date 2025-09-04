import os
import sys
# Add project root to sys.path for local imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
###########################################
# DO NOT MOVE ANY IMPORTS ABOVE THIS LINE #
###########################################

from matplotlib import pyplot as plt
import numpy as np
import yaml

from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.central_unit.central_unit import CentralUnit
from sub_THz_stripe.radio_unit.radio_unit import RadioUnit
from wireless_channel.subTHz_channel import Channel
from wireless_channel.waveforms import Waveform
from utils import spec, logger
from plotter import plotter

booster_stages = ["fiber", "coupler", "amplifier", "coupler"]
tx_stages = ["fiber", "coupler", "splitter", "shifter", "amplifier"]


if __name__ == "__main__":
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
    for stripe_cfg in config["radio_stripes"]:
        stripes.append(RadioStripe.from_config_locations(stripe_cfg, wf))

    # continue with 3 stripes, separated by 1m
    stripes = [stripes[5], stripes[6]]  # stripes[5:11:2]
    plotter.plot_stripes(config, stripes)
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug("stripe: %d: %s", stripe_idx, stripe)

    # load channels
    ue_pos = config["ue_positions"][0]
    channel = Channel.from_sionna(ue_pos, debug=True)
    logger.debug("%s", channel)

    cu = CentralUnit() # todo are these configs loadable?
    print(f'CU: {cu}')
    ofdm_time_after_cu = cu.run(ofdm_time) # shape: nr_ofdm_symbols x (fft_size + cp length)
    logger.debug('shape of ofdm timee: %s', ofdm_time_after_cu.shape)
    logger.debug('np alike: %s', np.allclose(ofdm_time, ofdm_time_after_cu))

    active_ru_idxes = [2, 2]  # , 4, 6
    logger.debug('transmitting over the stripe...')
    iq_at_last_rus = []
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug('stripe: %d - active ru %d', stripe_idx, active_ru_idxes[stripe_idx])
        stripe.active_unit = active_ru_idxes[stripe_idx]

        iq_data = ofdm_time_after_cu.reshape(1, -1) # flatten to (1 x nr_iq_symbols)

        iq_data = iq_data.reshape(wf.n_ofdm_symbols, -1)  # flatten to (1 x nr_iq_symbols)

        logger.debug('shape of iq data: %s', iq_data.shape) # 1 d array
        phase_shifts = [0, 0, 0, 0]
        iq_out, imdata = stripe.transmit(iq_data, phase_shifts)

        fig, ax = plt.subplots()
        ax.set_xlabel("Input amplitude |x|")
        ax.set_ylabel("Output amplitude |y|")
        for i in range(len(imdata)-1):
            # Make AM/AM plots
            x = imdata[i][0]
            y = imdata[i+1][0]
            if len(imdata[i].shape) >= 3:
                x = imdata[i][0][0]
            if len(imdata[i+1].shape) >= 3:
                y = imdata[i+1][0][0]
            ax.plot(np.abs(x), np.abs(y), 'o', label=f"stage{i}")
        ax.legend()
        fig.savefig(f"am_am_plot_stripe{stripe_idx}.pdf")

        iq_out_reshaped = iq_out.reshape(nr_antennas, wf.n_ofdm_symbols, -1)
        logger.debug('reshaped after stripe: %s', iq_out_reshaped.shape)
        iq_at_last_rus.append(iq_out_reshaped)

    wf.plot_psd(iq_at_last_rus.flatten())

    y_ue = channel.transmit_dl(iq_at_last_rus, active_ru_idxes, wf)
    logger.debug('received signal at ue: %s', y_ue.shape)

    wf.plot_iq_time(y_ue[0], title="After wireless channel")

    ue = RadioUnit(ue_pos['x'], ue_pos['y'], ue_pos['z'])
    logger.debug('ue RU: %s', ue)
    shifts = [0, 0, 0, 0]
    y_combined_time = ue.receive(y_ue, shifts)
    logger.debug('y combined shape: %s', y_combined_time.shape)
    y_combined_time = np.squeeze(y_combined_time, axis=0)

    wf.plot_iq_time(y_combined_time, title="At UE")

    y_combined_freq = wf.ofdm_time_to_freq(y_combined_time)
    # # todo do we need equalization?

    y_qam = wf.ofdm_to_qam(y_combined_freq)
    # wf.plot_constellation(y_qam, symbols_tx=qam, title="RX'ed symbols")

    y_bits = wf.qam_to_bits(y_qam)

    ber = wf.compute_ber(bits, y_bits)
    # logger.debug('BER: %f', ber)
    plt.show()
