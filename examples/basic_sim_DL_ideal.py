import os
"""
This script simulates a basic downlink transmission in a sub-THz radio stripe system using OFDM waveforms,
assuming only IDEAL hardware (no hardware impairments are modeled).

It loads configuration parameters from a YAML file, generates OFDM signals, models radio stripes and their transmission,
applies a wireless channel, and evaluates the received signal at a user equipment (UE) location.

Main steps:
1. Loads simulation and waveform configuration from a YAML file.
2. Plots the room layout and radio stripes.
3. Generates OFDM waveform and visualizes its power spectral density.
4. Constructs radio stripes and selects active radio units for transmission.
5. Processes the OFDM signal through a central unit and radio stripes.
6. Simulates transmission over a wireless channel to the UE.
7. Receives and processes the signal at the UE, including reshaping and combining.
8. Converts the received signal back to frequency domain, demodulates QAM symbols, and computes bit error rate (BER).
9. Visualizes various stages of the signal (time domain, after channel, at UE).

Dependencies:
- numpy
- matplotlib
- yaml
- Custom modules: sub_THz_stripe, wireless_channel, utils, plotter

Note:
- The script is intended to be run as a standalone module.
- Only IDEAL hardware is simulated (no hardware impairments).
- Some configuration loading and equalization steps are marked as TODO.
"""
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
from utils import logger, calculate_psd_per_symbol
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
    wf.plot_psd(ofdm_time, nperseg=wf.fft_size)


    # sanity check
    tx_freq_oversampled = wf.ofdm_time_to_freq(ofdm_time)  # uses your method
    tx_subc = wf.extract_subcarriers(tx_freq_oversampled)  # shape (n_sym, n_carriers)

    # check pilot magnitudes for the first OFDM symbol
    print("pilot indices:", wf.pilot_indices)
    print("pilot values (TX) for first symbol:", tx_subc[0, wf.pilot_indices])
    print("pilot magnitudes mean:", np.mean(np.abs(tx_subc[:, wf.pilot_indices])))

    # build all radio stripes
    stripes = []
    for stripe_cfg in config["radio_stripes"]:
        stripes.append(RadioStripe.from_config_locations(stripe_cfg, wf))

    # continue with 3 stripes, separated by 1m
    stripes = [stripes[5]]#, stripes[6]]  # stripes[5:11:2]
    active_ru_idxes = [1]#, 2]  # , 4, 6

    plotter.plot_stripes(config, stripes)
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug("stripe: %d: %s", stripe_idx, stripe)

    # load channels
    ue_pos = config["ue_positions"][0]
    channel = Channel.from_sionna(ue_pos, debug=True) # debug=True enables a dummy channel of all ones
    logger.debug("%s", channel)

    cu = CentralUnit() # todo are these configs loadable?
    print(f'CU: {cu}')
    ofdm_time_after_cu = cu.run(ofdm_time) # shape: nr_ofdm_symbols x (fft_size + cp length)
    logger.debug('shape of ofdm timee: %s', ofdm_time_after_cu.shape)
    logger.debug('np alike: %s', np.allclose(ofdm_time, ofdm_time_after_cu))

    logger.debug('transmitting over the stripe...')
    iq_at_last_rus = []
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug('stripe: %d - active ru %d', stripe_idx, active_ru_idxes[stripe_idx])
        stripe.active_unit = active_ru_idxes[stripe_idx]

        wf.plot_iq_time(ofdm_time_after_cu, title="Before reshape")

        iq_data = ofdm_time_after_cu.reshape(1, -1) # flatten to (1 x nr_iq_symbols)

        iq_data = iq_data.reshape(wf.n_ofdm_symbols, -1)  # flatten to (1 x nr_iq_symbols)
        wf.plot_iq_time(iq_data, title="After reshape")

        logger.debug('shape of iq data: %s', iq_data.shape) # 1 d array
        phase_shifts = [0, 0, 0, 0]
        iq_out, imdata = stripe.transmit(iq_data, phase_shifts)
        logger.debug('iq out shape: %s', iq_out.shape)
        iq_out_reshaped = iq_out.reshape(nr_antennas, wf.n_ofdm_symbols, -1)
        logger.debug('reshaped after stripe: %s', iq_out_reshaped.shape)
        iq_at_last_rus.append(iq_out_reshaped)

    y_ue = channel.transmit_dl(iq_at_last_rus, active_ru_idxes, wf)
    logger.debug('received signal at ue: %s', y_ue.shape)

    wf.plot_iq_time(y_ue[0], title="After wireless channel")

    ue = RadioUnit(ue_pos['x'], ue_pos['y'], ue_pos['z'])
    logger.debug('ue RU: %s', ue)
    shifts = [0, 0, 0, 0]
    y_combined_time, imdata = ue.receive(y_ue, shifts)
    logger.debug('y combined shape: %s', y_combined_time.shape)

    wf.plot_iq_time(y_combined_time, title="At UE")

    y_combined_freq = wf.ofdm_time_to_freq(y_combined_time)

    # sanity check
    rx_freq_oversampled = wf.ofdm_time_to_freq(y_combined_time)  # your received time -> freq
    rx_subc = wf.extract_subcarriers(rx_freq_oversampled)
    # look at a few carriers around pilots and data
    print("RX pilot bins first symbol:", rx_subc[0, wf.pilot_indices])
    print("RX some data bins first symbol (first 10):", rx_subc[0, wf.data_carriers[:10]])
    wf.plot_constellation(rx_subc[0, wf.pilot_indices], title="received pilots")

    # channel estimation
    subc = wf.extract_subcarriers(y_combined_freq)
    H_est = wf.channel_estimate_ls(subc)

    # sanity check
    print("H_est shape:", H_est.shape)
    # show a summary for first symbol
    print("H_est at pilot bins:", H_est[0, wf.pilot_indices])
    print("H_est magnitude stats:", np.min(np.abs(H_est)), np.median(np.abs(H_est)), np.max(np.abs(H_est)))

    # equalization
    eq_subc = wf.equalize_one_tap(subc, H_est)

    # sanity check
    i = wf.data_carriers[0]
    print("raw rx on that carrier (first sym):", rx_subc[0, i])
    print("H_est there:", H_est[0, i])
    print("after equalize:", eq_subc[0, i])

    # Demap data carriers and rebuild stream
    data_symbols = wf.demap_data_from_grid(eq_subc).flatten()  # these are the received QAM symbols

    y_qam = data_symbols

    # quick sanity
    assert y_qam.shape[0] == qam.shape[0], f"Lengths differ: rx {y_qam.shape[0]} tx {qam.shape[0]}"

    wf.plot_constellation(y_qam, symbols_tx=qam, title="equalized symbols")

    y_bits = wf.qam_to_bits(y_qam)

    ber = wf.compute_ber(bits, y_bits)
    logger.debug('BER: %f', ber)
    plt.show()
