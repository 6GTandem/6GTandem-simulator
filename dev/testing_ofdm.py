import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch

def get_bits(n_carriers, qam_order, n_ofdm_symbols):
    num_bits = n_ofdm_symbols * n_carriers * int(np.log2(qam_order))
    bits = np.random.randint(0, 2, num_bits)
    return bits

def get_qam_symbols(bits, qam_order):
    """Map bits to QAM symbols using Gray coding."""
    k = int(np.log2(qam_order))  # bits per symbol
    if len(bits) % k != 0:
        raise ValueError("Number of bits must be a multiple of log2(qam_order).")

    # Reshape bits into groups
    bit_groups = bits.reshape((-1, k))

    # Convert bit groups into integers
    symbols_idx = np.array([int("".join(str(b) for b in grp), 2) for grp in bit_groups])

    # Square QAM constellation
    m_side = int(np.sqrt(qam_order))
    if m_side ** 2 != qam_order:
        raise ValueError("qam_order must be a square number (e.g., 4, 16, 64).")

    # Gray coding mapping for I and Q separately
    def gray_map(x):
        return x ^ (x >> 1)

    I = gray_map(symbols_idx % m_side)
    Q = gray_map(symbols_idx // m_side)

    # Center constellation around zero
    I = 2 * I - (m_side - 1)
    Q = 2 * Q - (m_side - 1)

    # Normalize average power to 1
    symbols = (I + 1j * Q) / np.sqrt((2 / 3) * (qam_order - 1))

    return symbols

def plot_constellation(qam_symbols, qam_order):
    """Plot QAM constellation scatter diagram."""
    plt.scatter(np.real(qam_symbols), np.imag(qam_symbols))
    plt.title(f"{qam_order}-QAM Constellation")
    plt.xlabel("In-phase (I)")
    plt.ylabel("Quadrature (Q)")
    plt.axis("equal")
    plt.show()

def get_ofdm_symbols(qam_symbols, n_carriers=1024, cp_length=16, oversampling=4):
    """Generate OFDM symbols with cyclic prefix and oversampling."""

    # Number of OFDM symbols we can form
    n_symbols = len(qam_symbols) // n_carriers
    print(f'n_symbols:{n_symbols}')
    qam_symbols = qam_symbols[: n_symbols * n_carriers]  # truncate extra

    # Reshape into (n_symbols, n_carriers)
    qam_matrix = qam_symbols.reshape((n_symbols, n_carriers))

    # Set fft size: oversampling by zero-padding in frequency domain
    n_fft = n_carriers * oversampling

    ofdm_time = []
    for row in qam_matrix:
        freq_domain = np.zeros(n_fft, dtype=complex)

        # Map QAM symbols to center of spectrum (baseband symmetric)
        half = n_carriers // 2
        freq_domain[:half] = row[:half]
        freq_domain[-half:] = row[half:]

        # IFFT to get time-domain signal
        time_domain = np.fft.ifft(freq_domain, n_fft)

        # Add cyclic prefix
        cp = time_domain[-cp_length:]
        ofdm_symbol = np.concatenate([cp, time_domain])

        ofdm_time.append(ofdm_symbol)

    return np.array(ofdm_time) # shape: n_ofdm_symbols x ( n_carriers * oversampling + cp_length)

def ofdm_time_to_freq(ofdm_time, cp_length=16):
    """Remove cyclic prefix and perform FFT to get OFDM symbols in frequency domain."""
    n_fft = ofdm_time.shape[1] - cp_length
    ofdm_freq = []
    for symbol in ofdm_time:
        time_no_cp = symbol[cp_length:]
        freq_domain = np.fft.fft(time_no_cp, n_fft)
        ofdm_freq.append(freq_domain)
    return np.array(ofdm_freq) # shape: n_ofdm_symbols x (n_carriers * oversampling)


def ofdm_freq_to_time(ofdm_freq, cp_length=16):
    """Convert OFDM frequency-domain symbols back to time-domain with cyclic prefix."""
    n_fft = ofdm_freq.shape[1]
    ofdm_time = []
    for symbol in ofdm_freq:
        time_domain = np.fft.ifft(symbol, n_fft)
        cp = time_domain[-cp_length:]
        ofdm_symbol = np.concatenate([cp, time_domain])
        ofdm_time.append(ofdm_symbol)
    return np.array(ofdm_time) # shape: n_ofdm_symbols x ( n_carriers * oversampling + cp_length)

def awgn(signal, snr_dB):
    """Add AWGN noise to a signal given SNR in dB."""
    snr_linear = 10 ** (snr_dB / 10)
    power_signal = np.mean(np.abs(signal) ** 2)
    noise_power = power_signal / snr_linear
    noise = np.sqrt(noise_power / 2) * (np.random.randn(*signal.shape) + 1j * np.random.randn(*signal.shape))
    return signal + noise

def ofdm_to_qam(ofdm_freq_received, n_carriers, oversampling_factor):
    n_fft = n_carriers * oversampling_factor
    qam_received = []
    for sym in ofdm_freq_received:
        half = n_carriers // 2
        row_symbols = np.concatenate([sym[:half], sym[-half:]])
        qam_received.append(row_symbols)
    qam_received = np.array(qam_received).flatten()
    return qam_received
def qam_to_bits(qam_symbols, qam_order):
    """Fast demapping of QAM symbols to bits using vectorized operations."""
    m_side = int(np.sqrt(qam_order))
    k = int(np.log2(qam_order))

    # Undo normalization
    qam_symbols = qam_symbols * np.sqrt((2 / 3) * (qam_order - 1))

    # Map to nearest I and Q points
    I = np.clip(np.round((np.real(qam_symbols) + (m_side - 1) / 2)), 0, m_side-1).astype(int)
    Q = np.clip(np.round((np.imag(qam_symbols) + (m_side - 1) / 2)), 0, m_side-1).astype(int)

    # Inverse Gray code
    def inv_gray(x):
        y = np.zeros_like(x)
        for shift in range(int(np.log2(m_side))):
            y ^= x >> shift
        return y

    I = inv_gray(I)
    Q = inv_gray(Q)

    # Combine to symbol indices
    symbols_idx = Q * m_side + I

    # Convert indices to bits
    bits = (((symbols_idx[:, None] & (1 << np.arange(k)[::-1])) > 0)).astype(int)
    return bits.flatten()

def compute_ber(bits_tx, bits_rx):
    """Compute Bit Error Rate between transmitted and received bits."""
    if len(bits_tx) != len(bits_rx):
        raise ValueError("Transmitted and received bits must have the same length.")
    n_errors = np.sum(bits_tx != bits_rx)
    ber = n_errors / len(bits_tx)
    return ber

if __name__ == '__main__':
    # config
    n_ofdm_symbols = 2
    n_carriers = 1024
    qam_order = 16 # qam order
    oversampling_factor = 4
    cp_length = 128
    BW = 12.5e9
    fc = 157.75e9
    fs = BW * oversampling_factor

    # generate bit sequence
    bits = get_bits(n_carriers, qam_order, n_ofdm_symbols)

    # generate QAM symbols
    qam_symbols = get_qam_symbols(bits, qam_order)

    print(f"Generated {len(bits)} bits and {len(qam_symbols)} QAM symbols")
    print(f"Avg constellation pwr: {np.mean(np.abs(qam_symbols)**2)}")

    # get ofdm symbol
    ofdm_time = get_ofdm_symbols(qam_symbols, n_carriers=n_carriers, cp_length=cp_length, oversampling=oversampling_factor)
    print(f"Generated {ofdm_time.shape[0]} OFDM symbols of length {ofdm_time.shape[1]}")

    # move to freq domain
    ofdm_freq = ofdm_time_to_freq(ofdm_time, cp_length=cp_length)
    print(f"OFDM symbols in frequency domain: {ofdm_freq.shape}")

    # move to time domain
    ofdm_time_check = ofdm_freq_to_time(ofdm_freq, cp_length=cp_length)
    print(f"Generated {ofdm_time_check.shape[0]} OFDM symbols of length {ofdm_time_check.shape[1]}")

    # Compute PSD
    ofdm_signal = ofdm_time.flatten()
    f, Pxx = welch(ofdm_signal, fs=fs, nperseg=1024, return_onesided=False)
    # Center around zero and shift
    Pxx_db = 10 * np.log10(np.fft.fftshift(Pxx))
    f_shifted = np.fft.fftshift(f) #+ fc  # shift to RF
    f_shifted_GHz = f_shifted / 1e9
    plt.figure(figsize=(10, 4))
    plt.plot(f_shifted_GHz, Pxx_db)
    plt.xlabel("Frequency (GHz)")
    plt.xlim([-10, 10])
    plt.ylabel("PSD (dB)")
    plt.title("OFDM PSD at RF")
    plt.grid(True)
    plt.show()
    # note dip at DC is because there is no QAM symbol mapped to the DC carrier

    # send over awgn channel as test
    snr_dB = 10
    ofdm_time_noisy = awgn(ofdm_time.flatten(), snr_dB)
    ofdm_time_noisy = ofdm_time_noisy.reshape(ofdm_time.shape)  # reshape back to OFDM symbols

    # back to f domain
    ofdm_freq_received = ofdm_time_to_freq(ofdm_time_noisy, cp_length=cp_length)

    # back to qam symbols
    qam_received = ofdm_to_qam(ofdm_freq_received, n_carriers, oversampling_factor)

    # back to bits
    received_bits = qam_to_bits(qam_received, qam_order)
    print(f'bits all close: {np.allclose(bits, received_bits)}')

    ber = compute_ber(bits, received_bits)
    print(f"BER at {snr_dB} dB SNR: {ber:.6f}")

    plt.scatter(np.real(qam_received), np.imag(qam_received), label='Rx symbols')
    plt.scatter(np.real(qam_symbols), np.imag(qam_symbols), label='Tx symbols')
    plt.title(f"{qam_order}-QAM Constellation - BER: {ber}")
    plt.xlabel("In-phase (I)")
    plt.ylabel("Quadrature (Q)")
    plt.axis("equal")
    plt.legend()
    plt.show()