import yaml
import math
import numpy as np
import pandas as pd
import concurrent.futures
import logging
from datetime import datetime

from pathlib import Path, PurePath

# Configure logging to file and console
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%H:%M:%S',
    handlers=[
        logging.FileHandler('flickering.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)
from wireless_channel.waveforms import Waveform
from sub_THz_stripe.radiostripe.radiostripe import RadioStripe
from sub_THz_stripe.radio_unit.radio_unit import RadioUnit
from wireless_channel.subTHz_channel import Channel
from sub_THz_stripe.central_unit.central_unit import CentralUnit
from sub_THz_stripe.utils import getdbm, calculate_sndr
from sub_THz_stripe.amplifier.amplifier import Amplifier
from sub_THz_stripe.combiner.combiner import Combiner
from sub_THz_stripe.coupler.coupler import Coupler
from sub_THz_stripe.phase_shifter.phase_shifter import PhaseShifter
from sub_THz_stripe.splitter.splitter import Splitter

start_time = datetime.now()
logger.info("Starting flickering.py")

environment = "office_space_reduced_inline"

logger.info("Loading configuration files...")
# Read out the config file.
config_file = f"environments/{environment}/config.yaml"
with open(config_file, "r", encoding="utf8") as file:
    config = yaml.safe_load(file)

stripe_config = config["stripe_config"]

# Read out the Waveform config file.
config_file = "environments/waveform_config.yaml"
with open(config_file, "r", encoding="utf8") as file:
    waveform_config = yaml.safe_load(file)

# Construct the waveform class based on the loaded configuration.
freq_band_config = config["sub_thz"]
wf = Waveform.from_config(waveform_config, freq_band_config)

# Generate an CP-OFDM waveform in the time domain.
bits = wf.generate_bits()
qam = wf.qam_modulate()
ofdm_time = wf.ofdm_modulate()  # shape: nr_ofdm_symb x (fftsize + cp length)
x_combined_freq = wf.ofdm_time_to_freq(ofdm_time)
xsubc = wf.extract_subcarriers(x_combined_freq)

# The UE has the same hardware as a central unit.
ue = CentralUnit()  # todo are these configs loadable?
ofdm_time_after_ue = ue.run(ofdm_time)  # shape: nr_ofdm_symbols x (fft_size + cp length)

# A single radio unit is used to represent the transmit stage of the ue. The UE position does not influence the
# transmission so we can set it arbitrarily.
amp = Amplifier(bw=config["sub_thz"]["bw"])
coup = Coupler(wf=wf)
split = Splitter(config["antenna"]["N_antennas"])
comb = Combiner()
ps = PhaseShifter(num_shifters=config["antenna"]["N_antennas"], resolution=4)
ue_ru = RadioUnit(
    x=0, y=0, z=0, boost_amp=amp, antenna_amp=amp, coup_in=coup, coup_out=coup, splitter=split, combiner=comb, pshift=ps
)

# Read out the Component config file.
config_file = "environments/component_config.yaml"
with open(config_file, "r", encoding="utf8") as file:
    component_config = yaml.safe_load(file)

# Build all the RadioStripe objects based on the configuration.
logger.info(f"Building {len(config['radio_stripes'])} radio stripes...")
stripes: list[RadioStripe] = []
for stripe_cfg in config["radio_stripes"]:
    stripes.append(RadioStripe.from_config_locations(stripe_cfg, component_config, config["antenna"]["N_antennas"], wf))
logger.info(f"Setup complete. Processing {len(config['ue_positions'])} UE positions...")


def receive_worker(stripe: RadioStripe, data, phase_shifts, ue_idx, stripe_idx, ru_idx):
    stripe.active_unit = ru_idx
    y, imdata = stripe.receive(data, phase_shifts, False)

    out = [ue_idx, stripe_idx, ru_idx, y]
    return out


dataset_path = Path("wireless_channel/sionna_dataset/", environment, "flickering/")
dataset_path.mkdir(parents=True, exist_ok=True)


def find_closest_ru(data, ue, n_stripes, n_rus):
    """
    data: nested list of dictionaries
    ue: dict with 'x', 'y', 'z' of the UE
    """

    ue_x, ue_y, ue_z = ue["x"], ue["y"], ue["z"]

    closest_distance = float("inf")
    closest_index = None
    closest_ru_coords = None

    # Loop over outer list (groups)
    ru_index = 0
    for group in data:

        # Loop over elements inside each group
        for entry in group:
            if "radio_unit" in entry:  # skip central_unit
                ru = entry["radio_unit"]
                dx = ru["x"] - ue_x
                dy = ru["y"] - ue_y
                dz = ru["z"] - ue_z
                dist = math.sqrt(dx**2 + dy**2 + dz**2)

                if dist < closest_distance:
                    closest_distance = dist
                    closest_index = ru_index
                    closest_ru_coords = (ru["x"], ru["y"], ru["z"])

                ru_index += 1
            
    closest_stripe = closest_index // n_rus
    closest_ru = closest_index % n_rus

    return closest_distance, closest_ru_coords, closest_stripe, closest_ru


def ue_worker(ue_pos: dict, wf: Waveform, stripes: list[RadioStripe]):
    ue_data = []
    beam_angles = [-30, -10, 0, 10, 30]
    for ue_beam in beam_angles:
        for ru_beam in beam_angles:
            # Transmit the data towards the stripe. We only perform this task once. Every UE in the room will transmit the same
            # signal towards the stripe. Only its position and hence the channel will change.
            ue_shift = ue_ru.phase_shifter.get_phases(ue_beam)
            iq_data_tx, imdata = ue_ru.transmit(ofdm_time_after_ue, ue_shift)

            # Construct the wireless channel for this UE.
            channel = Channel.from_sionna(ue_pos, environment)

            # First determine which RUs are closest to the user as we will only analyze the 9 closest units.
            n_stripes = stripe_config["N_stripes"]
            n_rus = stripe_config["N_RUs"]
            dist, _, stripe_idx, ru_idx = find_closest_ru(config["radio_stripes"], ue_pos, n_stripes, n_rus)
            # print(f"UE coords: {ue_pos}")
            # print(f"Distance: {dist}")
            # print(f"RU coords: {_}")
            # print(f"RU index: {stripe_idx}, {ru_idx}")

            if stripe_idx == n_stripes - 1:
                stripe_idx = n_stripes - 2
            if ru_idx == n_rus - 1:
                ru_idx = n_rus - 3
            if stripe_idx == 0:
                stripe_idx = 1
            if ru_idx == 0:
                ru_idx = 2

            xgrid, ygrid = np.meshgrid((stripe_idx - 1, stripe_idx, stripe_idx + 1), (ru_idx - 2, ru_idx - 1, ru_idx, ru_idx + 1, ru_idx + 2))

            # Now receive this data on all the stripes.
            for row, col in zip(xgrid, ygrid):
                for rux, ruy in zip(row, col):
                    # for stripe_idx, (stripe, data) in enumerate(zip(stripes, iq_data_rx)):
                    if rux >= 0 and ruy >= 0 and rux < n_stripes and ruy < n_rus:
                        stripe = stripes[rux]
                        # Transmit the data over the channel and get the data at all the radio units.
                        selected_ru = ruy
                        if stripe_config.get('invert_stripe_dir', False):
                            selected_ru = n_rus - 1 - ruy
                        data = channel.transmit_ul_id(iq_data_tx, rux, selected_ru, wf)
                        ru_shift = stripe.radio_units[0].phase_shifter.get_phases(ru_beam)
                        # for ru_idx in range(len(stripe.radio_units)):
                        # Receive the data over the stripe coming from ru_idx.
                        stripe.active_unit = ruy
                        y, imdata = stripe.receive(data, ru_shift, False)
                        # print(f"Stripe {rux}")
                        # print(f"Antenna gain: {stripe.radio_units[ruy].antenna_amp.gain}")
                        # print(f"Antenna coeffs: {stripe.radio_units[ruy].antenna_amp.coeffs}")
                        # for ruid, ru in enumerate(stripe.radio_units[:ruy]):
                        #     meta_data.append([rux, ruid, ru.boost_amp.gain, ru.boost_amp.coeffs])
                        #     print(f"RU id: {ruid}")
                        #     print(f"Boost gain: {ru.boost_amp.gain}")
                        #     print(f"Boost coeffs: {ru.boost_amp.coeffs}")

                        y_combined_freq = wf.ofdm_time_to_freq(y)

                        # channel estimation
                        subc = wf.extract_subcarriers(y_combined_freq)
                        H_est = wf.channel_estimate_ls(subc)

                        # equalization
                        eq_subc = wf.equalize_one_tap(subc, H_est)

                        # Demap data carriers and rebuild stream
                        data_symbols = wf.demap_data_from_grid(eq_subc).flatten()  # these are the received QAM symbols

                        y_qam = data_symbols

                        # quick sanity
                        assert y_qam.shape[0] == qam.shape[0], f"Lengths differ: rx {y_qam.shape[0]} tx {qam.shape[0]}"

                        y_bits = wf.qam_to_bits(y_qam)

                        ber = wf.compute_ber(bits, y_bits)

                        nmse = np.sum(np.abs(eq_subc - xsubc) ** 2) / np.sum(np.abs(xsubc) ** 2)
                        nmse = 10 * np.log10(nmse)

                        # Calculate the sndr and average power received on the antennas of the RU.
                        pavg_ant = np.mean(getdbm(data))
                        pavg_ru = np.mean(getdbm(imdata[3]))
                        sndr_ru, _ = calculate_sndr(ofdm_time, imdata[4], wf)
                        sndr_ru = np.mean(sndr_ru)
                        # Calculate the sndr and average power received at the CU.
                        pavg_cu = np.mean(getdbm(y))
                        sndr_cu, _ = calculate_sndr(ofdm_time, y, wf)
                        sndr_cu = np.mean(sndr_cu)
                        # Save the data.
                        ue_data.append(
                            [
                                channel.ue_idx,
                                rux,
                                ruy,
                                ue_beam,
                                ru_beam,
                                pavg_ru,
                                pavg_ant,
                                pavg_cu,
                                sndr_ru,
                                sndr_cu,
                                nmse,
                                ber,
                            ]
                        )
                        # time_data.append([ofdm_time, y, x_combined_freq, y_combined_freq])

            # Save the data to disk.
            df = pd.DataFrame(
                ue_data,
                columns=[
                    "ue_id",
                    "stripe_id",
                    "ru_id",
                    "ue_beam_id",
                    "ru_beam_id",
                    "pavg_ru",
                    "pavg_ant",
                    "pavg_cu",
                    "sndr_ru",
                    "sndr_cu",
                    "nmse",
                    "ber",
                ],
            )
            file_path = PurePath(dataset_path, f"flickering_data_{channel.ue_idx}.pkl")
            df.to_pickle(file_path)
            # df = pd.DataFrame(time_data, columns=["x", "y", "x_freq", "y_freq"])
            # meta_df = pd.DataFrame(meta_data, columns=["stripe_idx", "ru_idx", "gain", "coeffs"])
            # file_path = PurePath(dataset_path, f"time_data_{channel.ue_idx}.hdf5")
            # df.to_hdf(file_path, key="data")
            # meta_df.to_hdf(file_path, key="meta")


# Loop over all the UE positions and perform an exhaustive search.
workers = []
total_ues = len(config["ue_positions"])
processed_count = 0
# for ue in config["ue_positions"][:10]:
#     print(f"UE Pos: {ue}")
#     ue_worker(ue, wf, stripes)
#ue_worker(config["ue_positions"][-1], wf, stripes)
with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
    for ue_idx, ue_pos in enumerate(config["ue_positions"]):
        worker = executor.submit(ue_worker, ue_pos, wf, stripes)
        workers.append(worker)

    # Monitor progress with timing
    for future in concurrent.futures.as_completed(workers):
        processed_count += 1
        elapsed = datetime.now() - start_time
        rate = processed_count / elapsed.total_seconds() if elapsed.total_seconds() > 0 else 0
        remaining = (total_ues - processed_count) / rate if rate > 0 else 0
        logger.info(f"UE progress: {processed_count}/{total_ues} | Elapsed: {elapsed.seconds//60}m {elapsed.seconds%60}s | Est. remaining: {int(remaining)//60}m {int(remaining)%60}s")

total_elapsed = datetime.now() - start_time
logger.info("Flickering analysis complete!")
logger.info(f"Total time: {total_elapsed.seconds//60}m {total_elapsed.seconds%60}s ({total_elapsed.total_seconds():.1f}s)")
