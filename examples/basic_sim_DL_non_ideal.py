import os
import sys
import logging

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
from utils import calculate_psd_per_symbol, setup_logging
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
    logger.debug("%d stripes in the room", len(config["radio_stripes"]))

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
    stripes = [stripes[5]]
    active_ru_idxes = [4]  # , 4, 6
    # , stripes[6]]  # stripes[5:11:2]
    plotter.plot_stripes(config, stripes)
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug("stripe: %d: %s", stripe_idx, stripe)

    # load channels
    ue_pos = config["ue_positions"][0]
    channel = Channel.from_sionna(ue_pos)  # , debug=True)
    logger.debug("%s", channel)

    cu = CentralUnit()  # todo are these configs loadable?
    logger.debug(f"CU: {cu}")
    ofdm_time_after_cu = cu.run(ofdm_time)  # shape: nr_ofdm_symbols x (fft_size + cp length)
    logger.debug("shape of ofdm timee: %s", ofdm_time_after_cu.shape)
    logger.debug("np alike: %s", np.allclose(ofdm_time, ofdm_time_after_cu))

    logger.debug("transmitting over the stripe...")
    iq_at_last_rus = []
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug("stripe: %d - active ru %d", stripe_idx, active_ru_idxes[stripe_idx])
        stripe.active_unit = active_ru_idxes[stripe_idx]

        iq_data = ofdm_time_after_cu.reshape(1, -1)  # flatten to (1 x nr_iq_symbols)

        iq_data = iq_data.reshape(wf.n_ofdm_symbols, -1)  # flatten to (1 x nr_iq_symbols)
        stripe.calibrate(iq_data, -30)

        logger.debug("shape of iq data: %s", iq_data.shape)  # 1 d array
        phase_shifts = [0, 0, 0, 0]
        iq_out, imdata = stripe.transmit(iq_data, phase_shifts)

        fig, fig2 = plotter.plot_im_data_tx(imdata, wf)

        fig.savefig(f"am_am_plot_stripe{stripe_idx}.pdf")
        fig2.savefig(f"spec_plot_stripe{stripe_idx}.pdf")

        iq_out_reshaped = iq_out.reshape(nr_antennas, wf.n_ofdm_symbols, -1)
        logger.debug("reshaped after stripe: %s", iq_out_reshaped.shape)
        iq_at_last_rus.append(iq_out_reshaped)

    iq_at_last_rus = np.array(iq_at_last_rus)
    wf.plot_psd(iq_at_last_rus)

    y_ue = channel.transmit_dl(iq_at_last_rus, active_ru_idxes, wf)
    logger.debug(f"received signal at ue: {y_ue.shape}")

    wf.plot_iq_time(y_ue[0], title="After wireless channel")
    wf.plot_psd(y_ue[0])

    ue = RadioUnit(ue_pos["x"], ue_pos["y"], ue_pos["z"], wf)
    logger.debug(f"ue RU: {ue}")
    shifts = [0, 0, 0, 0]
    y_combined_time, imdata = ue.receive(y_ue, shifts)
    logger.debug(f"y combined shape: {y_combined_time.shape}")

    y_combined_freq = wf.ofdm_time_to_freq(y_combined_time)

    # sanity check
    rx_freq_oversampled = wf.ofdm_time_to_freq(y_combined_time)  # your received time -> freq
    rx_subc = wf.extract_subcarriers(rx_freq_oversampled)
    # look at a few carriers around pilots and data
    logger.debug(f"RX pilot bins first symbol: {rx_subc[0, wf.pilot_indices]}")
    logger.debug(f"RX some data bins first symbol (first 10): {rx_subc[0, wf.data_carriers[:10]]}")
    wf.plot_constellation(rx_subc[0, wf.pilot_indices], title="received pilots")

    # channel estimation
    subc = wf.extract_subcarriers(y_combined_freq)
    H_est = wf.channel_estimate_ls(subc)

    # sanity check
    logger.debug(f"H_est shape: {H_est.shape}")
    # show a summary for first symbol
    logger.debug(f"H_est at pilot bins: {H_est[0, wf.pilot_indices]}")
    logger.debug(f"H_est magnitude stats: {np.min(np.abs(H_est))}, {np.median(np.abs(H_est))}, {np.max(np.abs(H_est))}")

    # equalization
    eq_subc = wf.equalize_one_tap(subc, H_est)

    # sanity check
    i = wf.data_carriers[0]
    logger.debug(f"raw rx on that carrier (first sym): {rx_subc[0, i]}")
    logger.debug(f"H_est there: {H_est[0, i]}")
    logger.debug(f"After equalize: {eq_subc[0, i]}")

    # Demap data carriers and rebuild stream
    data_symbols = wf.demap_data_from_grid(eq_subc).flatten()  # these are the received QAM symbols

    y_qam = data_symbols

    # quick sanity
    assert y_qam.shape[0] == qam.shape[0], f"Lengths differ: rx {y_qam.shape[0]} tx {qam.shape[0]}"

    wf.plot_constellation(y_qam, symbols_tx=qam, title="equalized symbols")

    y_bits = wf.qam_to_bits(y_qam)

    ber = wf.compute_ber(bits, y_bits)
    logger.debug(f"BER: {ber}")
    plt.show()

    labels = [f"Symbol{n}" for n in range(y_combined_time.shape[0])]
    fout = plotter.plot_iq_per_symbol(y_combined_time, wf.cp_length, wf.n_carriers, labels)
    fout.savefig("iq_per_symbol.pdf")

    fout = plotter.plot_iq_per_carrier(np.delete(y_combined_time, np.s_[1:], axis=0), wf.cp_length, wf.n_carriers)
    fout.savefig("iq_per_carrier.pdf")

    fout = plotter.plot_psd_per_symbol(y_combined_time, wf.fs)
    fout.savefig("psd_per_symbol.pdf")
