import pandas as pd
import numpy as np
from scipy.fft import ifft, ifftshift, fft, fftshift, irfft, rfft, fftfreq
from matplotlib import pyplot as plt


fig, axes = plt.subplots(nrows=3,ncols=1)
for f in ['No tape.csv', 'With tape.csv']:
    data = pd.read_csv(f, sep=r'\s+')

    # Extract the relevant columns
    freq = data['freq[Hz]'].values 
    db = data['db:Trc1_S21'].values
    ang = data['ang:Trc1_S21'].values

    bw = freq[-1]-freq[0]
    deltaT = 1/bw
    N = len(freq)

    tau = np.arange(N)*deltaT

    # Convert dB to linear scale
    linear_magnitude = 10 ** (db / 20.0)

    # Convert polar coordinates (magnitude and angle) to complex numbers
    complex_response = linear_magnitude * np.exp(1j * ang)



    # Compute the Power Delay Profile (PDP) using the inverse FFT
    # Given that there is only one measurement, we cannot take the average

    # first shift so DC is in bin 0
    pdp =np.abs(ifft(fftshift(complex_response)))
    
    axes[0].set_title("Spectrum Magnitude")
    axes[0].plot(freq/1e9, db, label=f, alpha=0.5)
    axes[0].set_xlabel("Frequency [GHz]")

    axes[1].set_title("PDP")

    axes[1].plot(tau*1e9, 20*np.log10(pdp), '-o', label=f, markersize = 0.1, alpha=0.5)
    axes[1].set_xlabel("Delay [ns]")

    # # Calculate the delay spread
    Pm = np.sum(pdp)
    tau_mean = np.sum(pdp*tau) / Pm
    tau_sq = np.sum(pdp*tau**2) / Pm
    tau_rms = np.sqrt(tau_sq-tau_mean**2)
    print(f"tau RMS for {f} is {tau_rms*1e9} ns")

    axes[2].set_title("frequency-domain correlation function")
    xf = fftfreq(len(pdp), deltaT)
    xf = fftshift(xf)
    axes[2].plot(xf/1e9, np.abs(fftshift(fft(pdp))), '-o', label=f, markersize = 0.1, alpha=0.5)
    axes[2].set_xlabel("$\Delta F$ [GHz]")


plt.legend()
plt.tight_layout()
plt.show()



# Calculate the coherence bandwidth