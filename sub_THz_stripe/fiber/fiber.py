import os
import numpy as np
import logging
import skrf as rf

from ..component.component import Component
from wireless_channel.waveforms import Waveform
from ..utils import db_to_magnitude
from ..utils import delay as delay_signal
from scipy.signal import lfilter

logger = logging.getLogger(__name__)


class Fiber(Component):
    def __init__(
        self,
        wf: Waveform,
        length: float = 1,
        damping_per_meter: float = 0,
        filter_mode="time_domain",
        filter: np.ndarray = np.array([1]),
        deembedding: float = 0,
        *args,
        **kwargs,
    ):
        """Initialize a fiber component.

        :param length: Length of the fiber in meter.
        :param damping_per_meter: Damping in dB per meter.
        :param filter: The impulse response of the fiber.


        If filter is used, the damping and delay are included in the filter. In this case set length and
        damping_per_meter to zero.
        """
        self.length = length
        self.damping_per_meter = damping_per_meter
        self.fs = wf.fs
        self.filter = filter  # passed in frequency domain
        self.filter_mode = filter_mode
        self.wf = wf
        self.deembedding = deembedding

        super().__init__(*args, **kwargs)

    def run(self, x, delay: bool = False, window: str = "fixed"):
        match self.filter_mode:
            case "time_domain":
                logger.debug("Fiber model is being applied in the time domain.")
                taps = np.fft.ifft(np.fft.ifftshift(self.filter))
                logger.debug(f"Filter taps used: {taps}")
                # x_filt = delay(lfilter(taps, [1.0], x.flatten(), [self.delay]) #todo delay needed or not???
                x_filt = lfilter(taps, [1.0], x.flatten())
                if delay:
                    x_filt = delay_signal(x_filt, self.delay, window=window)
                xout = np.reshape(x_filt, (self.wf.n_ofdm_symbols, -1))
            case "freq_domain":
                logger.debug("Fiber model being applied in the frequency domain.")
                logger.debug(f"Filter shape: {self.filter.shape}")
                logger.debug(f"Filter: {self.filter}")

                x_freq = self.wf.ofdm_time_to_freq(x)
                logger.debug(f"xfreq shape: {x_freq.shape}")

                f_shift = np.fft.fftshift(self.filter)
                x_freq_filtered = x_freq * f_shift

                logger.debug(f"xfreq filtered shape: {x_freq_filtered.shape}")
                xout = self.wf.ofdm_freq_to_time(x_freq_filtered)  # back to time
                xout = xout.flatten()
                if delay:
                    xout = delay_signal(xout, self.delay, window=window)
                xout = np.reshape(xout, (self.wf.n_ofdm_symbols, -1))
            case _:
                xout = x.flatten()
                if delay:
                    xout = delay_signal(xout, self.delay, window=window)
                xout = np.reshape(xout, (self.wf.n_ofdm_symbols, -1))

        return xout * db_to_magnitude(self.damping) * db_to_magnitude(self.deembedding)

    @property
    def delay(self):
        return self.length * 1.5 / 3e8 * self.fs  # 1.5 / 3e8 speed of EM waves in fiber

    @property
    def damping(self):
        return self.length * self.damping_per_meter

    @classmethod
    def from_config(cls, config: dict, wf: Waveform):
        """Construct a Fiber object based on a configuration dictionary.

        Parameters
        ----------
        config : dict
            Dictionary containing the following keys:
            'model': A model file located in the models/PMF folder.
            'length': Physical length of the fiber in meter. (optional)
            'damping_per_meter': Additional damping per meter in dB. (optional)
        wf : Waveform
            Waveform object used to extract the right wavefrom characteristics.

        Returns
        -------
        fib : Fiber
            Fiber object with the requested parameters.
        """
        if "model" in config:
            # Load the fiber model from the given file.
            base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            fiber_spars = rf.Network(os.path.join(base_path, f"../models/PMF/{config['model']}"))

            # Subcarrier frequencies
            ofdm_freqs = np.linspace(wf.fc - wf.bw / 2, wf.fc + wf.bw / 2, wf.n_carriers * wf.oversampling_factor)

            fib_freqs = fiber_spars.f
            fib_s21 = fiber_spars.s[:, 1, 0]
            
            fib_mag = np.abs(fib_s21)
            fib_phase = np.angle(fib_s21)

            interp_mag = np.interp(ofdm_freqs, fib_freqs, fib_mag)
            interp_phase = np.interp(ofdm_freqs, fib_freqs, fib_phase)
            fib_s21_ofdm = interp_mag * np.exp(1j * interp_phase)

            filter_mode = "freq_domain"
            damping = config.get("damping_per_meter", 0)
            length = config.get("length", 0)
            deembedding = config.get("deembedding", 0)

            fiber = cls(damping_per_meter=damping, length=length, filter=fib_s21_ofdm, filter_mode=filter_mode, wf=wf, deembedding=deembedding)
        else:
            fiber = cls(**config, filter_mode="no_filter", wf=wf)

        return fiber
