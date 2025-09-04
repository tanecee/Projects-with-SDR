#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch

def read_dot_file(filename: str) -> np.ndarray:
    raw = np.fromfile(filename, dtype=np.float32)
    if raw.size % 2:
        raw = raw[:-1]
    I = raw[0::2]
    Q = raw[1::2]
    return I.astype(np.float32) + 1j * Q.astype(np.float32)

def build_indices(nfft: int):
    act_scs = nfft // 2
    data_scs = act_scs // 2
    pilot_scs = act_scs // 2
    base = nfft // 4  # 64 when nfft=256

    # 1-based MATLAB style
    data_1b = np.arange(1, 2*data_scs+1, 2) + base         # 65,67,... for nfft=256
    pilot_1b = np.arange(2, 2*pilot_scs+1, 2) + base        # 66,68,...

    # Convert to 0-based (correct for Python/C++)
    data_0b = data_1b - 1                                   # 64,66,...
    pilot_0b = pilot_1b - 1

    # Also provide +1-shifted variant to catch TX code that used 1-based indices directly
    data_0b_shiftp1 = data_0b + 1                           # 65,67,...
    pilot_0b_shiftp1 = pilot_0b + 1

    return {
        "data_0b": data_0b.astype(int),
        "pilot_0b": pilot_0b.astype(int),
        "data_0b_shiftp1": data_0b_shiftp1.astype(int),
        "pilot_0b_shiftp1": pilot_0b_shiftp1.astype(int),
        "data_count": data_scs,
        "pilot_count": pilot_scs,
    }

def extract_constellation(x, fs, nfft, cp_len, nsym, zc_len, use_shift=False):
    if x.size < zc_len + (nfft + cp_len):
        raise ValueError("File too short for given parameters.")

    x_ofdm = x[zc_len:]
    sym_len = nfft + cp_len
    total_needed = nsym * sym_len
    if x_ofdm.size < total_needed:
        total_needed = (x_ofdm.size // sym_len) * sym_len
    x_ofdm = x_ofdm[:total_needed]
    K = x_ofdm.size // sym_len

    X = x_ofdm.reshape(K, sym_len)
    X_no_cp = X[:, cp_len:]
    F = np.fft.fft(X_no_cp, nfft, axis=1)
    F = np.fft.fftshift(F, axes=1)

    idx = build_indices(nfft)
    if use_shift:
        data_idx = idx["data_0b_shiftp1"]
    else:
        data_idx = idx["data_0b"]

    data_eq = F[:, data_idx]
    data_flat = data_eq.reshape(-1)
    return data_flat

def analyze_tx_dot(path, fs=20e6, nfft=256, cp_len=64, nsym=14, zc_len=255, show_time=True):
    x = read_dot_file(path)

    # Try both index hypotheses and choose the one with more energy
    d0 = extract_constellation(x, fs, nfft, cp_len, nsym, zc_len, use_shift=False)
    d1 = extract_constellation(x, fs, nfft, cp_len, nsym, zc_len, use_shift=True)

    pwr0 = np.mean(np.abs(d0)**2)
    pwr1 = np.mean(np.abs(d1)**2)
    use_shift = pwr1 > pwr0
    data_flat = d1 if use_shift else d0

    # Plots
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle("OFDM TX Analysis", fontsize=14, fontweight="bold")

    ax1 = fig.add_subplot(2, 2, 1)
    if show_time:
        Nshow = min(4000, x.size)
        t_ms = np.arange(Nshow)/fs*1e3
        ax1.plot(t_ms, np.real(x[:Nshow]))
        ax1.plot(t_ms, np.imag(x[:Nshow]))
        ax1.set_title("Time Domain (first {} samples)".format(Nshow))
        ax1.set_xlabel("Time (ms)")
        ax1.set_ylabel("Amplitude")
        ax1.grid(True)
    else:
        ax1.axis("off")
        ax1.set_title("Time Domain (disabled)")

    ax2 = fig.add_subplot(2, 2, 2)
    nperseg = min(4096, x.size)
    noverlap = nperseg // 2
    f, pxx = welch(x, fs=fs, nperseg=nperseg, noverlap=noverlap, return_onesided=False)
    ax2.plot(np.fft.fftshift(f)/1e6, 10*np.log10(np.fft.fftshift(pxx) + 1e-20))
    ax2.set_title("Power Spectrum (Welch)")
    ax2.set_xlabel("Frequency (MHz)")
    ax2.set_ylabel("PSD (dB/Hz)")
    ax2.grid(True)

    ax3 = fig.add_subplot(2, 2, (3, 4))
    ax3.plot(np.real(data_flat), np.imag(data_flat), ".", markersize=4)
    ax3.set_aspect("equal", adjustable="box")
    tag = "+1 shift" if use_shift else "0-based"
    ax3.set_title(f"Constellation (Data SCs) [{tag}]")
    ax3.set_xlabel("I")
    ax3.set_ylabel("Q")
    ax3.grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True, help="Path to .dot file (float32 interleaved IQ)")
    ap.add_argument("--fs", type=float, default=20e6, help="Sample rate (Hz)")
    ap.add_argument("--nfft", type=int, default=256, help="OFDM FFT size")
    ap.add_argument("--cp", type=int, default=64, help="Cyclic prefix length (samples)")
    ap.add_argument("--nsym", type=int, default=14, help="Number of OFDM symbols in the frame")
    ap.add_argument("--zc_len", type=int, default=255, help="ZC preamble length at the start of the frame")
    ap.add_argument("--no_time", action="store_true", help="Disable time-domain plot")
    args = ap.parse_args()

    analyze_tx_dot(
        path=args.file,
        fs=args.fs,
        nfft=args.nfft,
        cp_len=args.cp,
        nsym=args.nsym,
        zc_len=args.zc_len,
        show_time=not args.no_time,
    )