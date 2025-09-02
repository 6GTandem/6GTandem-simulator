import pandas as pd
import numpy as np
from scipy.fft import ifft, ifftshift, fft, fftshift, irfft, rfft, fftfreq
from matplotlib import pyplot as plt

from scipy.constants import c

import os

n_fiber = 1 / 0.67
c_fiber = c / n_fiber

file_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

fig, axes = plt.subplots(nrows=3,ncols=1)
corrs = [[]]*2
for i, f in enumerate(["1m_not_taped", "1m_taped"]):
    print(F"Processing: {f}")
    data = pd.read_csv(
        f+".csv",
        sep=r",",
        dtype={
            "freq[Hz]": np.float64,
            "db:Trc2_S21": np.float64,
            "ang:Trc2_S21": np.float64,
        },
    )

    # Extract the relevant columns
    freq = data['freq[Hz]'].values
    db = data['db:Trc2_S21'].values
    ang = np.deg2rad(data["ang:Trc2_S21"].values)

    fs = np.abs(freq[-1] - freq[0])
    deltaT = 1 / fs
    N = len(freq)
    tau = np.arange(0, N / fs, 1 / fs)

    # Convert dB to linear scale
    linear_magnitude = 10 ** (db / 20.0)  # TODO assume S21 is in power scale

    # Convert polar coordinates (magnitude and angle) to complex numbers
    complex_response = linear_magnitude * np.exp(1j * ang)

    print(f"Average S21: {20*np.log10(np.mean(np.abs(complex_response))):.2f} dB")
    print(f"Min S21: {20*np.log10(np.min(np.abs(complex_response))):.2f} dB")

    # Compute the Power Delay Profile (PDP) using the inverse FFT
    # Given that there is only one measurement, we cannot take the average

    # first shift so DC is in bin 0
    pdp = np.abs(ifft(fftshift(complex_response))) ** 2

    axes[0].set_title("Spectrum Magnitude")
    axes[0].plot(freq / 1e9, db, label=f, alpha=0.5)
    axes[0].set_xlabel("Frequency [GHz]")

    axes[1].set_title("PDP")

    print(f"Estimated fiber length: {tau[np.argmax(pdp)] * c_fiber}")

    print(deltaT * c_fiber)

    axes[1].plot(
        tau * 1e9, 10 * np.log10(pdp), "-o", label=f, markersize=0.1, alpha=0.5
    )
    axes[1].set_xlabel("Delay [ns]")

    # pd.DataFrame({"t": tau * 1e9, "pdp": 10 * np.log10(pdp)}).to_csv(
    #     os.path.join(file_dir, f"TXT/pdp_{f}.txt"), index=False
    # )

    # pdp_filterd = np.abs(pdp)
    # # pdp_filterd[10*np.log10(pdp)<-50.0] = 1.0e-10
    # pdp_filterd[(tau*1e6<5) | (tau*1e6>10)] = 1.0e-10
    # axes[1].plot(tau*1e6, 10*np.log10(pdp_filterd), '-o', label=f+" filtered", markersize = 0.1, alpha=0.5)

    # pdp_complex = ifft(fftshift(complex_response))
    # pdp_complex[(tau*1e6<5) | (tau*1e6>10)] = 1.0e-10

    # axes[1].plot(tau*1e9, fft(fft(pdp_complex)), label=f, alpha=0.5)

    # # Calculate the delay spread
    Pm = np.sum(pdp)
    tau_mean = np.sum(pdp * tau) / Pm
    tau_sq = np.sum(pdp * tau**2) / Pm
    tau_rms = np.sqrt(tau_sq - tau_mean**2)
    print(f"tau mean for {f} is {tau_mean*1e9:.4f} ns")
    print(f"tau RMS for {f} is {tau_rms*1e9:.4f} ns")

    print(f"fiber_length {c*tau_mean:.2f}m")

    print(f"Coherecen B (0.9) {f} is {1/(50*tau_rms*1e9):.4f} GHz")
    print(f"Coherence B (0.5) {f} is {1/(5*tau_rms*1e9):.4f} GHz")

    axes[2].set_title("frequency-domain correlation function")
    xf = fftfreq(len(pdp), deltaT)
    xf = fftshift(xf)
    y = np.abs(fftshift(fft(pdp)))
    idx1 = np.argmin(np.abs((y / np.max(y)) - 0.5))
    idx2 = np.argmin(np.abs((y / np.max(y)) - 0.7))
    idx3 = np.argmin(np.abs((y / np.max(y)) - 0.9))
    print(f"Coherence B (0.5) {f} is {np.abs(xf[idx1]/1e9)*2}")
    print(f"Coherence B (0.7) {f} is {np.abs(xf[idx2]/1e9)*2}")
    print(f"Coherence B (0.9) {f} is {np.abs(xf[idx3]/1e9)*2}")

    axes[2].plot(xf / 1e9,y / np.max(y), "-o", label=f, markersize=0.1)
    axes[2].set_xlabel("$\Delta F$ [GHz]")

    # pd.DataFrame({"xf": xf / 1e9, "corr":y / np.max(y)}).to_csv(
    #     os.path.join(file_dir, f"TXT/corr_{f}.txt"), index=False
    # )

    corrs[i] = np.abs(y / np.max(y))

axes[2].axhline(y=0.5, color="r", linestyle="--", alpha=0.5)
axes[2].axhline(y=0.9, color="r", linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show(block=True)


# only look at the positive frequencies
corrs0 = corrs[0][len(corrs[0])//2:]
corrs1 = corrs[1][len(corrs[1]) // 2 :]
xf = xf[len(xf)//2 :]


coh_diff = np.zeros_like(corrs0)

for j, c in enumerate(corrs1):
    # find index which is closest

    closest_idx = np.argmin(np.abs(corrs0-c))

    coh_diff[j] = (xf[j] - xf[closest_idx])*2

plt.figure()
# plt.plot(xf / 1e9, corrs0, "-o", label=f, markersize=0.1)
# plt.plot(xf / 1e9, corrs1, "-o", label=f, markersize=0.1)
plt.plot(corrs1, coh_diff / 1e9)
plt.tight_layout()
plt.show(block=True)

# pd.DataFrame({"x": corrs1, "y": coh_diff / 1e9}).to_csv(
#     os.path.join(file_dir, "TXT/corr_diff.txt"), index=False
# )
