import numpy as np

from ..component.component import Component
from ..utils import delay
from scipy.signal import lfilter


class IQModem(Component):
    def __init__(self, iqi_coef: int = 0, iqi_filter: int = 1, iqi_delay_imbalance: int = 0,
                 dc_offset: int = 0, *args, **kwargs):
        self.iqi_coef = iqi_coef
        self.iqi_filter = iqi_filter
        self.iqi_delay_imbalance = iqi_delay_imbalance
        self.dc_offset = dc_offset
        self.modes.extend(['filter', 'static'])

        super().__init__(*args, **kwargs)

    def run(self, yin, phasor):
        phasor = np.array(phasor)
        assert (len(yin) == len(phasor)
                ), "Yin and Phasor data lengths do not match."

        match self.mode:
            case 'ideal':
                yout = yin
            case 'filter':
                d = (self.iqi_filter - 1) // 2
                data = np.concatenate((np.conj(yin), np.zeros((1, d))), axis=1)
                xc = lfilter(self.iqi_filter, 1, data)
                xc = xc[d:]
                yout = yin + xc + self.dc_offset
            case 'static':
                yout = (delay(np.real(yin), [-self.iqi_delay_imbalance / 2]) +
                        1j * delay(np.imag(yin), [self.iqi_delay_imbalance / 2]))
                yout = yout + self.iqi_coef * np.conj(yout) + self.dc_offset

        return yout * phasor
