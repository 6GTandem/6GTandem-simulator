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

from PySide6 import QtWidgets
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
    app = QtWidgets.QApplication([])
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

    # Construct waveform class.
    waveform_config = config["waveform_config"]
    freq_band_config = config["sub_thz"]
    wf = Waveform.from_config(waveform_config, freq_band_config)
    logger.debug(f"{wf}")

    # Generate a CP-OFDM waveform in time domain.
    bits = wf.generate_bits()
    qam = wf.qam_modulate()
    ofdm_time = wf.ofdm_modulate()  # shape: nr_ofdm_symb x (fftsize + cp length)
    wf.plot_psd(ofdm_time)

    # Build all radio stripes.
    stripes = []
    component_config = config["component_config"]
    for stripe_cfg in config["radio_stripes"]:
        stripes.append(RadioStripe.from_config_locations(stripe_cfg, component_config, wf))

    # As a simplified example only simulate the first two stripes.
    # We cannot select stripes arbitrarily. Otherwise the stripe indexes no longer match the
    # ones used by the channel. For example if we select only stripe 6 it gets index 0 in the
    # stripes list. The channel will select the channel for stripe index 0 which does not
    # match with the actual selected stripe. Same goes for selecting radio units.
    stripes = stripes[0:2]
    plotter.plot_stripes(config, stripes)
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug("stripe: %d: %s", stripe_idx, stripe)

    # load channels. Select only a single UE position to start with. Same indexing remark
    # applies here as for the stripes.
    ue_pos = config["ue_positions"][3017]
    channel = Channel.from_sionna(ue_pos, debug=True)
    channel.Nr_stripes = 2
    logger.debug(f"{channel}")

    # Currently we use the central unit at the UE side to transmit the OFDM signal.
    cu = CentralUnit()  # TODO: are these configs loadable?
    logger.debug(f"CU: {cu}")
    ofdm_time_after_cu = cu.run(ofdm_time)  # shape: nr_ofdm_symbols x (fft_size + cp length)
    logger.debug(f"shape of ofdm timee: {ofdm_time_after_cu.shape}")
    logger.debug(f"np alike: {np.allclose(ofdm_time, ofdm_time_after_cu)}")

    # We use a single radio unit to represent the UE.
    ue = RadioUnit(ue_pos["x"], ue_pos["y"], ue_pos["z"], wf)
    logger.debug(f"ue RU: {ue}")
    shifts = [0, 0, 0, 0]
    iq_data_tx, imdata = ue.transmit(ofdm_time_after_cu, shifts)
    # Make an IQ-plot and PSD-plot of the data transmitted at the UE (one antenna).
    wf.plot_iq_time(iq_data_tx[0], title="Transmitted at UE")
    wf.plot_psd(iq_data_tx[0], title="Transmitted at UE")

    iq_data_rus = channel.transmit_ul(iq_data_tx, wf)
    wf.plot_psd(iq_data_rus[0][0][0], title="After wireless channel")

    # Plot the received IQ-data for several RUs of the first stripe.
    ru_ofdm_freqs = []
    for i, ru_data in enumerate(iq_data_rus[0][0:10]):
        ru_ofdm_freqs.append(wf.ofdm_time_to_freq(ru_data[0])[0])

    labels = [f"RU{i+1}" for i in range(20)]
    fig = plotter.plot_constellation(np.array(ru_ofdm_freqs), labels=labels, title="IQ received at the different RUs.")

    logger.debug("Receiving data on all stripes.")
    iq_stripes = []
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug(f"stripe: {stripe_idx}")
        logger.debug("shape of iq data: %s", iq_data_rus.shape)
        phase_shifts = [0, 0, 0, 0]
        iq_out, imdata = stripe.receive_all(iq_data_rus[stripe_idx], phase_shifts)

        fig, fig2 = plotter.plot_im_data_rx_all(imdata[:9*5], wf)

        fig.savefig(f"am_am_plot_stripe{stripe_idx}.pdf")
        fig2.savefig(f"spec_plot_stripe{stripe_idx}.pdf")

        iq_stripes.append(iq_out)

    y_combined_freq = wf.ofdm_time_to_freq(iq_stripes[0])

    # sanity check
    rx_freq_oversampled = wf.ofdm_time_to_freq(iq_stripes[0])  # your received time -> freq
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
    #plt.show(block=True)

    widget = plotter.DynamicPlotter(imdata, len(stripes[0].radio_units), "rx", wf)
    widget.show()

    sys.exit(app.exec())
