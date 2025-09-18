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
    base = nfft // 4  # 64 when nfft=256

    # 1-based MATLAB style - C++ kodundaki gibi
    data_1b = np.arange(1, 2*data_scs+1, 2) + base         # 65,67,69,...,191
    
    # DC'yi çıkar (129)
    data_1b = data_1b[data_1b != 129]
    
    # Convert to 0-based (Python için)
    data_0b = data_1b - 1

    return {
        "data_0b": data_0b.astype(int),
        "data_1b": data_1b.astype(int),
        "data_count": len(data_1b),
    }

def extract_constellation(x, fs, nfft, cp_len, nsym, zc_len):
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
    
    # FFT ve ölçekleme (C++ kodundaki gibi)
    F = np.fft.fft(X_no_cp, nfft, axis=1)
    F = np.fft.fftshift(F, axes=1)
    F = F / np.sqrt(nfft)  # C++'daki sqrt(Nfft) ölçeklemesi

    idx = build_indices(nfft)
    data_idx = idx["data_0b"]

    data_eq = F[:, data_idx]
    data_flat = data_eq.reshape(-1)
    
    # Normalize et
    if len(data_flat) > 0:
        max_val = np.max(np.abs(data_flat))
        if max_val > 0:
            data_flat = data_flat / max_val
    
    return data_flat

def analyze_tx_dot(path, fs=20e6, nfft=256, cp_len=64, nsym=14, zc_len=255, show_time=True):
    x = read_dot_file(path)
    print(f"File loaded: {len(x)} samples")
    
    # Debug info
    idx_info = build_indices(nfft)
    print(f"Data indices (1-based): {idx_info['data_1b']}")
    print(f"Data indices (0-based): {idx_info['data_0b']}")
    print(f"Expected data symbols: {nsym * idx_info['data_count']}")
    
    try:
        data_flat = extract_constellation(x, fs, nfft, cp_len, nsym, zc_len)
        print(f"Extracted {len(data_flat)} data symbols")
        print(f"Data power: {10*np.log10(np.mean(np.abs(data_flat)**2) + 1e-20):.2f} dB")
    except Exception as e:
        print(f"Error extracting constellation: {e}")
        return

    # Plots
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle("OFDM TX Analysis", fontsize=14, fontweight="bold")

    # Time domain
    ax1 = fig.add_subplot(2, 2, 1)
    if show_time:
        Nshow = min(4000, x.size)
        t_ms = np.arange(Nshow)/fs*1e3
        ax1.plot(t_ms, np.real(x[:Nshow]), label='Real')
        ax1.plot(t_ms, np.imag(x[:Nshow]), label='Imag')
        ax1.set_title("Time Domain (first {} samples)".format(Nshow))
        ax1.set_xlabel("Time (ms)")
        ax1.set_ylabel("Amplitude")
        ax1.legend()
        ax1.grid(True)
    else:
        ax1.axis("off")

    # Spectrum
    ax2 = fig.add_subplot(2, 2, 2)
    nperseg = min(4096, x.size)
    f, pxx = welch(x, fs=fs, nperseg=nperseg, return_onesided=False)
    ax2.plot(np.fft.fftshift(f)/1e6, 10*np.log10(np.fft.fftshift(pxx) + 1e-20))
    ax2.set_title("Power Spectrum (Welch)")
    ax2.set_xlabel("Frequency (MHz)")
    ax2.set_ylabel("PSD (dB/Hz)")
    ax2.grid(True)

    # Constellation
    ax3 = fig.add_subplot(2, 2, (3, 4))
    if len(data_flat) > 0:
        ax3.scatter(np.real(data_flat), np.imag(data_flat), s=4, alpha=0.6)
        ax3.set_aspect("equal", adjustable="box")
        ax3.set_title("Constellation (Data Subcarriers)")
        ax3.set_xlabel("I")
        ax3.set_ylabel("Q")
        ax3.grid(True)
        

    else:
        ax3.text(0.5, 0.5, "No data symbols extracted", ha='center', va='center')
        ax3.set_title("Constellation - No Data")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True, help="Path to .dot file")
    ap.add_argument("--fs", type=float, default=20e6, help="Sample rate (Hz)")
    ap.add_argument("--nfft", type=int, default=256, help="OFDM FFT size")
    ap.add_argument("--cp", type=int, default=64, help="Cyclic prefix length")
    ap.add_argument("--nsym", type=int, default=14, help="Number of OFDM symbols")
    ap.add_argument("--zc_len", type=int, default=255, help="ZC preamble length")
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

#Bu kodu command window ile python3 ofdm_tx_analyzer.py --file tx_ofdm_suan.dot --fs 20e6 --nfft 256 --cp 64 --nsym 14 --zc_len 255 ile cagirip kaydettiğiniz sinyalin analizini yapabilirsiniz.
