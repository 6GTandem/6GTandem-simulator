import models
import numpy as np
import matplotlib.pyplot as pyplt


if __name__ == '__main__':
    gain = 1
    amax = 1
    nv = 0
    smoothness = 1
    modes = ['ideal', 'linear', 'atan', 'tanh', 'poly3',
             'poly3_pm', 'poly5', 'limiter', 'softlimiter', '6gtandem']

    fig, ax = pyplt.subplots()

    for mode in modes:
        x = np.random.rand(100)
        amp = models.Amplifier(gain, amax, nv, smoothness, mode)

        xout = amp.run(x)

        ax.plot(x, xout, 'o', label=mode, alpha=0.5)
        ax.set_ylabel("xout")
        ax.set_xlabel("x")
        ax.grid(which='both')
        ax.legend()

        fig.savefig(f"amplifiers.pdf")
