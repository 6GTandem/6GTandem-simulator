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
from wireless_channel.waveforms import Waveform

if __name__ == "__main__":

    # step 1 create stripes
    config_file = "office_config.yml"

    dir_path = os.path.dirname(os.path.realpath(__file__))

    config_path = os.path.join(dir_path, "..", "configurations")

    # Read the YAML file
    with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
        config = yaml.safe_load(file)

    # construct object
    waveform_config = config["waveform_config"]
    freq_band_config = config['sub_thz']
    wf = Waveform.from_config(waveform_config, freq_band_config)
    print(wf)

    # generate bits, convert to QAM symbols, convert to ofdm symbols
    bits = wf.generate_bits()
    qam = wf.qam_modulate()
    print(f'qam pwr: {np.mean(np.abs(qam)**2)}')
    ofdm_time = wf.ofdm_modulate()
    wf.plot_psd(ofdm_time)

    # Send over noisy channel
    ofdm_noisy = wf.awgn(ofdm_time.flatten(), snr_dB=1).reshape(ofdm_time.shape)

    # Back to frequency
    ofdm_freq_rx = wf.ofdm_time_to_freq(ofdm_noisy)

    # Back to QAM symbols
    qam_rx = wf.ofdm_to_qam(ofdm_freq_rx)

    # plot constellation
    wf.plot_constellation(qam, qam_rx)

    # Back to bits
    bits_rx = wf.qam_to_bits(qam_rx)

    ber = wf.compute_ber(bits, bits_rx)
    print(f"BER: {ber:.6f}")

    print(f'qam pwr: {np.mean(np.abs(qam)**2)}')
    print(f'ofdm pwr: {np.mean(np.abs(ofdm_time)**2) }')
    print(f'theoretical pwr: {(wf.n_carriers/(wf.fft_size+ wf.cp_length)) * (1/wf.fft_size)}')
    # note on theoretical pwr:
    # n_carriers each have QAM symbol with unit power,
    # this get's divided over (fft_size + cp_length) time domain samples,
    # on top of this numpy scales the power in the ifft with 1/fftsize
    # todo rescaling etc should probably be handled by the calibrate function
    print(f'ofdm time shape: {ofdm_time.shape}')