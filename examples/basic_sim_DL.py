import os
import sys
# Add project root to sys.path for local imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
###########################################
# DO NOT MOVE ANY IMPORTS ABOVE THIS LINE #
###########################################


import numpy as np
import yaml


from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.central_unit.central_unit import CentralUnit
from sub_THz_stripe.radio_unit.radio_unit import RadioUnit
from wireless_channel.subTHz_channel import Channel
from wireless_channel.waveforms import Waveform
from utils import spec
from plotter import plotter


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
    print(f"{len(config['radio_stripes'])} stripes in the room")

    # build all radio stripes
    stripes = []
    for stripe_cfg in config["radio_stripes"]:
        stripes.append(RadioStripe.from_config_locations(stripe_cfg))

    # continue with 3 stripes, separated by 1m
    stripes = stripes[5:11:2]
    plotter.plot_stripes(config, stripes)
    for stripe_idx, stripe in enumerate(stripes):
        print(f'stripe: {stripe_idx}: {stripe}')

    # construct waveform class
    waveform_config = config["waveform_config"]
    freq_band_config = config['sub_thz']
    wf = Waveform.from_config(waveform_config, freq_band_config)
    print(wf)

    # generate ofdm waveform in time domain
    bits = wf.generate_bits()
    qam = wf.qam_modulate()
    ofdm_time = wf.ofdm_modulate() # shape: nr_ofdm_symb x (fftsize + cp length)
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
    ofdm_time_after_cu = cu.run(ofdm_time) # shape: nr_ofdm_symbols x (fft_size + cp length)
    print(f' shape of ofdm timee: {ofdm_time_after_cu.shape}')
    print(f'np alike: {np.allclose(ofdm_time, ofdm_time_after_cu)}')
    # todo check with impairments if something changes

    # send over stripes
    active_ru_idxes = [2, 4, 6]
    print(f'transmitting over the stripe...')
    iq_at_last_rus = []
    for stripe_idx, stripe in enumerate(stripes):
        print(f'stripe: {stripe_idx} - active ru {active_ru_idxes[stripe_idx]}')
        # set active units
        stripe.active_unit = active_ru_idxes[stripe_idx]

        # transmit over stripe
        iq_data = ofdm_time_after_cu.reshape(1, -1) # flatten to (1 x nr_iq_symbols)
        #todo for now just same IQ data over all stripes, change if we want multiplexing
        print(f'shape of iq data: {iq_data.shape}') # 1 d array
        phase_shifts = [0, 0, 0, 0]
        # loop over RUs
        for ru_idx, iq_out in enumerate(stripe.transmit(iq_data, phase_shifts)):
            print(f'ru {ru_idx}: iq out shape: {iq_out.shape}')
            # store iq at active RU
            if ru_idx == active_ru_idxes[stripe_idx]:
                # reshape back to shape: nr_ofdm_symbols x (fft_size + cp_length)
                iq_out_reshaped = iq_out.reshape(nr_antennas, wf.n_ofdm_symbols, -1)
                print(f'reshaped after stripe: {iq_out_reshaped.shape}')
                iq_at_last_rus.append(iq_out_reshaped)


    # send over channel
    y_ue = channel.transmit_dl(iq_at_last_rus, active_ru_idxes, wf)
    # expected shape nr_ue_antennas x nr_ofdm_symb x fft_length+cp_length
    print(f'received signal at ue: {y_ue.shape}')

    # RU that acts as UE => no couplers!
    ue = RadioUnit(ue_pos['x'], ue_pos['y'], ue_pos['z'])
    print(f'ue RU: {ue}')
    shifts = [0, 0, 0, ]
    y_combined_time = ue.receive(y_ue, shifts)
    print(f'y combined shape: {y_combined_time.shape}')
    y_combined_time = np.squeeze(y_combined_time, axis=0)

    # convert back to f domain
    y_combined_freq = wf.ofdm_time_to_freq(y_combined_time)

    # ofdm to qam
    y_qam = wf.ofdm_to_qam(y_combined_freq)
    wf.plot_constellation(qam, y_qam)

    # qam to bits
    y_bits = wf.qam_to_bits(y_qam)

    # compute ber
    ber = wf.compute_ber(bits, y_bits)
    print(f'BER: {ber}')
