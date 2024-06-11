from scipy.signal import welch, get_window
import matplotlib.pyplot as plt
import numpy as np

def spec(inputsignal, fs=1, N=1024, plot=False):
    

    num_signals = inputsignal.shape[1]

    for k in range(num_signals):
        x = inputsignal[:, k]
        N = min(N, len(x) - 1)
        window = get_window("hann", N)
        f, Pxx = welch(
            x,
            fs=fs,
            window=window,
            nperseg=N,
            return_onesided=False,
            scaling="density",
        )

        s = 10 * np.log10(Pxx)

        if plot:
            plt.plot(
                np.fft.fftshift(f), np.fft.fftshift(s), linewidth=2, label=str(k + 1)
            )

    if plot:
        plt.xlabel("Normalized Frequency")
        plt.ylabel("Power Spectral Density (dB/Hz)")
        plt.legend()
        plt.grid()
        plt.show()

    return s
