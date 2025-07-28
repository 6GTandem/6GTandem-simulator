from abc import ABC, abstractmethod

class Component(ABC):
    modes = ['ideal', 'linear', 'atan', 'tanh', 'poly3', 'poly3_pm',
             'poly5', 'limiter', 'softlimiter', 'weblab', '6gtandem']

    def __init__(self, mode: str = 'ideal'):
        self._mode = mode

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, mode):
        """Set a new mode.

        Valid options are:
            'ideal','linear','atan', 'tanh', 'poly3','poly3_pm',
            'poly5','limiter','softlimiter','weblab','6gtandem'
        """
        if mode not in self.modes:
            raise ValueError(f"Invalid mode: {mode}")

        self._mode = mode

    @abstractmethod
    def run(self, x):
        ...