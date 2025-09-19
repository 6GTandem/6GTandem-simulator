import os
import sys
# Add project root to sys.path for local imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
###########################################
# DO NOT MOVE ANY IMPORTS ABOVE THIS LINE #
###########################################

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import skrf as rf
import yaml

from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.central_unit.central_unit import CentralUnit
from sub_THz_stripe.radio_unit.radio_unit import RadioUnit
from sub_THz_stripe.coupler.coupler import Coupler
from sub_THz_stripe.fiber.fiber import Fiber
from sub_THz_stripe.amplifier.amplifier import Amplifier
from wireless_channel.subTHz_channel import Channel
from wireless_channel.waveforms import Waveform
from utils import logger, calculate_psd_per_symbol
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

    # Load the couplers S-parameter file
    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # todo check if still works in vscode
    coupler_spars = rf.Network(os.path.join(base_path, 'models/coupler/with_balun.s2p')) #rf.Network('models/coupler/with_balun.s2p')


    # OFDM parameters
    fc = 157.75e9  # center frequency
    bw = 12.5e9   # bandwidth
    num_subcarriers = wf.fft_size

    # Subcarrier frequencies
    ofdm_freqs = np.linspace(fc - bw/2, fc + bw/2, num_subcarriers)

    # Extract frequency and S21 (transmission)
    coup_freqs = coupler_spars.f
    coup_s21 = coupler_spars.s[:, 1, 0]  # S21

    # Interpolate magnitude and phase separately for better accuracy
    coup_s21_mag = np.abs(coup_s21)
    coup_s21_phase = np.angle(coup_s21)

    interp_mag = np.interp(ofdm_freqs, coup_freqs, coup_s21_mag)
    interp_phase = np.interp(ofdm_freqs, coup_freqs, coup_s21_phase)
    coup_s21_ofdm = interp_mag * np.exp(1j * interp_phase)

    # Compute impulse response
    impulse_response = np.fft.ifft(coup_s21_ofdm)

    damping = 0  # in dB
    filter_mode = 'freq_domain'
    cp = Coupler(damping=damping, filter=coup_s21_ofdm, filter_mode=filter_mode, wf=wf)

    #fig = plotter.verify_impulse_response(cp.run, ofdm_freqs)
    #fig.savefig("coupler_interpolated.pdf")

    fiber_spars = pd.read_csv(os.path.join(base_path, 'models/PMF/with_tape/1m_taped.csv'))

    fib_freqs = fiber_spars["freq[Hz]"]
    fib_phase = np.unwrap(np.deg2rad(fiber_spars["ang:Trc2_S21"]))
    fib_mag = 10 ** (fiber_spars["db:Trc2_S21"] / 20.0)
    fib_s21 = fib_mag * np.exp(1j * fib_phase)
    
    interp_mag = np.interp(ofdm_freqs, fib_freqs, fib_mag)
    interp_phase = np.interp(ofdm_freqs, fib_freqs, fib_phase)
    fib_s21_ofdm = interp_mag * np.exp(1j * interp_phase)
    print(f' fiber filter taps: {fib_s21_ofdm.shape} - {fib_s21_ofdm}')

    fib_group_delay = -np.gradient(interp_phase, ofdm_freqs)

    # Compute impulse response
    impulse_response = np.fft.ifft(fib_s21_ofdm)

    fib = Fiber(damping_per_meter=0, length=0, filter=fib_s21_ofdm, filter_mode=filter_mode, wf=wf)

    #fig = plotter.verify_impulse_response(fib.run, ofdm_freqs)
    #fig.savefig("fiber_interpolated.pdf")

    fig, ax = plt.subplots()
    ax.plot(ofdm_freqs, fib_group_delay)
    ax.set_xlabel("Frequency")
    ax.set_ylabel("Group Delay")
    fig.savefig("fiber_group_delay.pdf")

    # build all radio stripes
    stripes = []
    for stripe_cfg in config["radio_stripes"]:
        stripes.append(RadioStripe.from_config_locations(stripe_cfg, wf))
    
    amp = Amplifier(8, max_out_amp=0.3, mode='poly3')
    amp.set_noise_var(273.5 + 30, 160e9, 0)
    
    for i in range(5):
        stripes[5].fibers[i] = fib
        stripes[5].radio_units[i].amp = amp
        stripes[5].radio_units[i].coupler_in = cp
        stripes[5].radio_units[i].coupler_out = cp
    
    # continue with 3 stripes, separated by 1m
    stripes = [stripes[5]]
    active_ru_idxes = [4]  # , 4, 6
    #, stripes[6]]  # stripes[5:11:2]
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

    logger.debug('transmitting over the stripe...')
    iq_at_last_rus = []
    for stripe_idx, stripe in enumerate(stripes):
        logger.debug('stripe: %d - active ru %d', stripe_idx, active_ru_idxes[stripe_idx])
        stripe.active_unit = active_ru_idxes[stripe_idx]

        iq_data = ofdm_time_after_cu.reshape(1, -1) # flatten to (1 x nr_iq_symbols)

        iq_data = iq_data.reshape(wf.n_ofdm_symbols, -1)  # flatten to (1 x nr_iq_symbols)
        stripe.calibrate(iq_data, -30)

        logger.debug('shape of iq data: %s', iq_data.shape) # 1 d array
        phase_shifts = [0, 0, 0, 0]
        iq_out, imdata = stripe.transmit(iq_data, phase_shifts)

        # Make AM/AM plots for every stage in the stripe.
        stages = ((len(imdata) - 6) // 4) + 1
        plot_panes = stages
        if stages % 2 != 0:
            plot_panes += 1

        # Make spec plots for every stage in the stripe.
        fig, ax = plt.subplots(plot_panes // 2, 2)
        fig2, ax2 = plt.subplots(plot_panes // 2, 2)

        # Keep the array 1xm for simplicity.
        if len(ax.shape) == 1:
            ax = np.array([ax])
        if len(ax2.shape) == 1:
            ax2 = np.array([ax2])

        # Set the proper axis labels.
        ax[-1, 0].set_xlabel("Input amplitude |x|")
        ax[-1, -1].set_xlabel("Input amplitude |x|")
        ax[0, 0].set_ylabel("Output amplitude |y|")

        ax2[-1, 0].set_xlabel("Normalized Frequency")
        ax2[-1, -1].set_xlabel("Normalized Frequency")
        ax2[0, 0].set_ylabel("Power Spectral Density (dB/Hz)")

        # Loop over the data from all the stages. (Stages are RUs and BUs.)
        for i in range(len(imdata)-1):
            # Calculate at which stage we are.
            stage = (i // 4) + 1
            # The last stage is a RU containing 5 components instead of 4. We need to
            # compensate for this or we advance to an non-existing stage.
            if stage > stages:
                stage = stages

            # Extract the data going into the stage and coming out of it.
            x = imdata[i][0]
            y = imdata[i+1]
            # If the matrix is three dimensional we need to extract one level deeper.
            if len(imdata[i].shape) >= 3:
                x = imdata[i][0][0]
            if len(imdata[i+1].shape) >= 3:
                y = imdata[i+1][0]

            # Calculate the PSD.
            freq, psd = calculate_psd_per_symbol(y, wf.fs)

            if stage >= stages:
                label = tx_stages[i%5]
                label += str(stages)
            else:
                label = booster_stages[i%4]
                label += str((i // 4) + 1)
            column = 0
            if stage > (plot_panes // 2):
                column = 1
            row = (stage-1) % (plot_panes // 2)

            ax[row, column].plot(np.abs(x), np.abs(y[0]), 'o', label=label)
            ax[row, column].legend()

            ax2[row, column].plot(freq, psd, label=label)
            ax2[row, column].legend()

        fig.savefig(f"am_am_plot_stripe{stripe_idx}.pdf")
        fig2.savefig(f"spec_plot_stripe{stripe_idx}.pdf")

        iq_out_reshaped = iq_out.reshape(nr_antennas, wf.n_ofdm_symbols, -1)
        logger.debug('reshaped after stripe: %s', iq_out_reshaped.shape)
        iq_at_last_rus.append(iq_out_reshaped)

    iq_at_last_rus = np.array(iq_at_last_rus)
    wf.plot_psd(iq_at_last_rus)

    y_ue = channel.transmit_dl(iq_at_last_rus, active_ru_idxes, wf)
    logger.debug('received signal at ue: %s', y_ue.shape)

    wf.plot_iq_time(y_ue[0], title="After wireless channel")
    wf.plot_psd(y_ue[0])

    ue = RadioUnit(ue_pos['x'], ue_pos['y'], ue_pos['z'], wf)
    logger.debug('ue RU: %s', ue)
    shifts = [0, 0, 0, 0]
    y_combined_time, imdata = ue.receive(y_ue, shifts)
    logger.debug('y combined shape: %s', y_combined_time.shape)

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

    labels = [f"Symbol{n}" for n in range(y_combined_time.shape[0])]
    fout = plotter.plot_iq_per_symbol(y_combined_time, wf.cp_length, wf.n_carriers, labels)
    fout.savefig("iq_per_symbol.pdf")

    fout = plotter.plot_iq_per_carrier(np.delete(y_combined_time, np.s_[1:], axis=0), wf.cp_length, wf.n_carriers)
    fout.savefig("iq_per_carrier.pdf")

    fout = plotter.plot_psd_per_symbol(y_combined_time, wf.fs)
    fout.savefig("psd_per_symbol.pdf")