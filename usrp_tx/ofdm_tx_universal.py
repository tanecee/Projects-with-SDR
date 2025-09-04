#!/usr/bin/env python3
import argparse
import os
import numpy as np

# Switch to non-interactive backend if DISPLAY is not set
_headless = os.environ.get("DISPLAY", "") == ""
if _headless:
    import matplotlib
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from scipy.signal import welch
from datetime import datetime

def read_dot_file(filename: str) -> np.ndarray:
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    raw = np.fromfile(filename, dtype=np.float32)
    if raw.size % 2:
        raw = raw[:-1]
    if raw.size == 0:
        raise ValueError("Empty or unreadable .dot file")
    I = raw[0::2]
    Q = raw[1::2]
    return I.astype(np.float32) + 1j * Q.astype(np.float32)

def build_indices(nfft: int):
    act_scs = nfft // 2
    data_scs = act_scs // 2
    pilot_scs = act_scs // 2
    base = nfft // 4

    data_1b = np.arange(1, 2*data_scs+1, 2) + base
    pilot_1b = np.arange(2, 2*pilot_scs+1, 2) + base

    data_0b = data_1b - 1
    pilot_0b = pilot_1b - 1

    return {
        "data_0b": data_0b.astype(int),
        "pilot_0b": pilot_0b.astype(int),
        "data_0b_shiftp1": (data_0b + 1).astype(int),
        "pilot_0b_shiftp1": (pilot_0b + 1).astype(int),
        "data_count": data_scs,
        "pilot_count": pilot_scs,
    }

def ideal_constellation(mod: str) -> np.ndarray:
    mod = mod.lower()
    if mod in ["bpsk","2","m2"]:
        pts = np.array([-1+0j, 1+0j], dtype=np.complex64)
        return pts / np.sqrt(np.mean(np.abs(pts)**2))
    if mod in ["qpsk","4","m4"]:
        re = np.array([-1, 1], dtype=np.float32)/np.sqrt(2.0)
        pts = np.array([r + 1j*i for r in re for i in re], dtype=np.complex64)
        return pts / np.sqrt(np.mean(np.abs(pts)**2))
    if mod in ["16qam","qam16","16"]:
        m_side = 4
    elif mod in ["64qam","qam64","64"]:
        m_side = 8
    else:
        raise ValueError("Unsupported modulation: " + mod)
    levels = np.arange(-(m_side-1), (m_side), 2, dtype=np.float32)
    xx, yy = np.meshgrid(levels, levels)
    pts = (xx + 1j*yy).ravel().astype(np.complex64)
    pts = pts / np.sqrt(np.mean(np.abs(pts)**2))
    return pts

def hard_decision_evmMER(samples: np.ndarray, mod_set: np.ndarray):
    if samples.size == 0:
        return np.nan, np.nan
    diffs = samples.reshape(-1,1) - mod_set.reshape(1,-1)
    d2 = np.abs(diffs)**2
    idx_min = np.argmin(d2, axis=1)
    ref = mod_set[idx_min]
    err = samples - ref
    evm_rms = np.sqrt(np.mean(np.abs(err)**2) / (np.mean(np.abs(ref)**2) + 1e-12))
    evm_pct = 100.0 * evm_rms
    mer_db = 20.0 * np.log10(1.0 / (evm_rms + 1e-12))
    return evm_pct, mer_db

def extract_freq_grid(x, fs, nfft, cp_len, nsym, zc_len, debug=False):
    if x.size < zc_len + (nfft + cp_len):
        raise ValueError("File too short for given parameters.")
    x_ofdm = x[zc_len:]
    sym_len = nfft + cp_len
    total_needed = nsym * sym_len
    if x_ofdm.size < total_needed:
        total_needed = (x_ofdm.size // sym_len) * sym_len
    x_ofdm = x_ofdm[:total_needed]
    K = x_ofdm.size // sym_len
    if K == 0:
        raise ValueError("No full OFDM symbol in file with given parameters.")
    X = x_ofdm.reshape(K, sym_len)
    X_no_cp = X[:, cp_len:]
    F = np.fft.fft(X_no_cp, nfft, axis=1)
    F = np.fft.fftshift(F, axes=1)
    if debug:
        print(f"[DEBUG] K(sym)={K}, nfft={nfft}, cp={cp_len}, total_samples={x.size}")
    return F

def pick_best_bins(F, nfft, debug=False):
    idx = build_indices(nfft)
    candidates = [
        ("data",  "0b",  idx["data_0b"]),
        ("data",  "+1",  idx["data_0b_shiftp1"]),
        ("pilot", "0b",  idx["pilot_0b"]),
        ("pilot", "+1",  idx["pilot_0b_shiftp1"]),
    ]
    best = None
    best_pwr = -1.0
    for role, tag, inds in candidates:
        inds = np.asarray(inds)
        inds = inds[(inds >= 0) & (inds < nfft)]
        if inds.size == 0:
            continue
        v = F[:, inds].reshape(-1)
        pwr = float(np.mean(np.abs(v)**2))
        if debug:
            print(f"[DEBUG] {role} [{tag}] bins={inds.size}, pwr={pwr:.4e}")
        if pwr > best_pwr:
            best_pwr = pwr
            best = (role, tag, v)
    if best is None:
        return "none", "none", np.array([], dtype=np.complex64)
    return best

def auto_detect_mod(samples: np.ndarray, mods_try=("bpsk","qpsk","16qam","64qam")):
    best_mod = None
    best_evm = np.inf
    best_mer = -np.inf
    norm = np.sqrt(np.mean(np.abs(samples)**2)) + 1e-12
    s = samples / norm
    for m in mods_try:
        ref = ideal_constellation(m)
        evm, mer = hard_decision_evmMER(s, ref)
        if evm < best_evm:
            best_mod = m
            best_evm = evm
            best_mer = mer
    return best_mod, best_evm, best_mer

def analyze(path, fs=20e6, nfft=256, cp_len=64, nsym=14, zc_len=255, mod="auto", show_time=True, save_path=None, debug=False):
    x = read_dot_file(path)
    F = extract_freq_grid(x, fs, nfft, cp_len, nsym, zc_len, debug=debug)
    role, tag, sc = pick_best_bins(F, nfft, debug=debug)

    if sc.size > 0:
        norm = np.sqrt(np.mean(np.abs(sc)**2)) + 1e-12
        sc_norm = sc / norm
    else:
        sc_norm = sc

    if mod.lower() == "auto":
        chosen_mod, evm_pct, mer_db = auto_detect_mod(sc_norm)
    else:
        chosen_mod = mod.lower()
        ref = ideal_constellation(chosen_mod)
        evm_pct, mer_db = hard_decision_evmMER(sc_norm, ref)

    # Plot
    import matplotlib.pyplot as plt
    from scipy.signal import welch
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle("OFDM TX Analysis", fontsize=14, fontweight="bold")

    ax1 = fig.add_subplot(2, 2, 1)
    if show_time:
        Nshow = min(4000, x.size)
        t_ms = np.arange(Nshow)/fs*1e3
        ax1.plot(t_ms, np.real(x[:Nshow]))
        ax1.plot(t_ms, np.imag(x[:Nshow]))
        ax1.set_title(f"Time Domain (first {Nshow} samples)")
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
    ax2.plot(np.fft.fftshift(f)/1e6, 10*np.log10(np.fft.fftshift(pxx)+1e-20))
    ax2.set_title("Power Spectrum (Welch)")
    ax2.set_xlabel("Frequency (MHz)")
    ax2.set_ylabel("PSD (dB/Hz)")
    ax2.grid(True)

    ax3 = fig.add_subplot(2, 2, (3, 4))
    if sc_norm.size > 0:
        ax3.plot(np.real(sc_norm), np.imag(sc_norm), ".", markersize=5)
        ax3.set_aspect("equal", adjustable="box")
        lim = 1.6 if chosen_mod in ["bpsk","qpsk"] else 2.5 if chosen_mod=="16qam" else 4.0
        ax3.set_xlim(-lim, lim); ax3.set_ylim(-lim, lim)
        ax3.grid(True)
        ax3.set_title(f"Constellation ({role} [{tag}]) • Mod={chosen_mod.upper()} • EVM={evm_pct:.2f}% • MER={mer_db:.2f} dB")
        ax3.set_xlabel("I"); ax3.set_ylabel("Q")
    else:
        ax3.set_title("Constellation (no bins found)")
        ax3.axis("off")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Save/Show logic
    auto_saved = None
    if save_path:
        plt.savefig(save_path, dpi=150)
    elif _headless:
        # Auto-save when headless and no --save provided
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        auto_saved = f"ofdm_tx_{stamp}.png"
        plt.savefig(auto_saved, dpi=150)
    else:
        plt.show()

    # Console summary
    print("=== Summary ===")
    print(f"Best bin set : {role} [{tag}]")
    print(f"Modulation   : {chosen_mod.upper()}")
    print(f"EVM (RMS %)  : {evm_pct:.3f}")
    print(f"MER (dB)     : {mer_db:.3f}")
    if save_path:
        print(f"Figure saved : {save_path}")
    if auto_saved:
        print(f"Figure saved (auto): {auto_saved}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True)
    parser.add_argument("--fs", type=float, default=20e6)
    parser.add_argument("--nfft", type=int, default=256)
    parser.add_argument("--cp", type=int, default=64)
    parser.add_argument("--nsym", type=int, default=14)
    parser.add_argument("--zc_len", type=int, default=255)
    parser.add_argument("--mod", type=str, default="auto", help="auto|bpsk|qpsk|16qam|64qam")
    parser.add_argument("--no_time", action="store_true")
    parser.add_argument("--save", type=str, default="")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    analyze(
        path=args.file,
        fs=args.fs,
        nfft=args.nfft,
        cp_len=args.cp,
        nsym=args.nsym,
        zc_len=args.zc_len,
        mod=args.mod,
        show_time=not args.no_time,
        save_path=(args.save if args.save else None),
        debug=args.debug,
    )
