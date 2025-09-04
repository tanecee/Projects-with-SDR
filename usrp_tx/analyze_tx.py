import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import argparse
import os

def read_dot_file(filename):
    """ .dot dosyasını oku ve complex array'e dönüştür """
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    real_part = data[0::2]
    imag_part = data[1::2]
    return real_part + 1j * imag_part



def welch_or_fft_psd(sig, fs, nperseg=1024):
    try:
        f, Pxx = signal.welch(sig, fs, nperseg=nperseg, return_onesided=False)
        return f, Pxx
    except Exception:
        # NumPy FFT fallback
        N = len(sig)
        win = np.hanning(N)
        S = np.fft.fftshift(np.fft.fft(sig * win))
        Pxx = (np.abs(S) ** 2) / (fs * np.sum(win**2))
        f = np.fft.fftshift(np.fft.fftfreq(N, d=1.0/fs))
        return f, Pxx

def analyze_signal(sig, fs, title="Sinyal Analizi"):
    """ Sinyali detaylı analiz et (yalnızca TX/tek dosya) """

    # Plot ayarları
    plt.style.use('default')
    fig = plt.figure(figsize=(15, 12))
    fig.suptitle(title, fontsize=16, fontweight='bold')

    # 1) Zaman domain - İlk 1000 örnek
    ax1 = plt.subplot(3, 3, 1)
    Ntime = min(1000, len(sig))
    t = np.arange(Ntime) / fs * 1000.0  # ms
    ax1.plot(t, np.real(sig[:Ntime]), label='I', alpha=0.8)
    ax1.plot(t, np.imag(sig[:Ntime]), label='Q', alpha=0.8)
    ax1.set_xlabel('Zaman (ms)')
    ax1.set_ylabel('Genlik')
    ax1.set_title('Zaman Domain - İlk {} Örnek'.format(Ntime))
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2) Spektrum (PSD)
    ax2 = plt.subplot(3, 3, 2)
    f, Pxx = welch_or_fft_psd(sig, fs, nperseg=1024)
    ax2.plot(np.fft.fftshift(f) / 1e3, 10 * np.log10(np.fft.fftshift(Pxx) + 1e-20))
    ax2.set_xlabel('Frekans (kHz)')
    ax2.set_ylabel('Güç (dB/Hz)')
    ax2.set_title('Güç Spektral Yoğunluğu')
    ax2.grid(True, alpha=0.3)

    # 3) Constellation (ham)
    ax3 = plt.subplot(3, 3, 3)
    ax3.scatter(np.real(sig), np.imag(sig), s=1, alpha=0.5)
    ax3.set_xlabel('I')
    ax3.set_ylabel('Q')
    ax3.set_title('Constellation Diagram')
    ax3.grid(True, alpha=0.3)
    ax3.axis('equal')

    # 4) Genlik Histogram
    ax4 = plt.subplot(3, 3, 4)
    ax4.hist(np.abs(sig), bins=50, alpha=0.9, density=True)
    ax4.set_xlabel('Genlik')
    ax4.set_ylabel('Olasılık Yoğunluğu')
    ax4.set_title('Genlik Histogramı')
    ax4.grid(True, alpha=0.3)

    # 5) Faz Histogram
    ax5 = plt.subplot(3, 3, 5)
    phase = np.angle(sig)
    ax5.hist(phase, bins=50, alpha=0.9, density=True)
    ax5.set_xlabel('Faz (radyan)')
    ax5.set_ylabel('Olasılık Yoğunluğu')
    ax5.set_title('Faz Histogramı')
    ax5.grid(True, alpha=0.3)

    # 6) Otokorelasyon (kompleks conjugate ile)
    ax6 = plt.subplot(3, 3, 6)
    corr = np.correlate(sig, np.conj(sig), mode='full')
    lags = np.arange(-len(sig)+1, len(sig))
    # Görsel rahatlık için merkez çevresini göster
    show = min(len(sig), 2000)
    mid = len(corr)//2
    view = corr[mid:mid+show]
    ax6.plot(np.arange(show), np.abs(view))
    ax6.set_xlabel('Gecikme (örnek)')
    ax6.set_ylabel('Korelasyon')
    ax6.set_title('Otokorelasyon (merkez pencere)')
    ax6.grid(True, alpha=0.3)

    # 7) İstatistikler
    ax7 = plt.subplot(3, 3, (7, 9))
    ax7.axis('off')

    mean_power = np.mean(np.abs(sig)**2)
    peak_power = np.max(np.abs(sig)**2)
    dr_db = 10 * np.log10((peak_power + 1e-20) / (mean_power + 1e-20))
    stats_text = [
        f"Sinyal İstatistikleri:",
        f"Toplam örnek: {len(sig):,}",
        f"Örnekleme hızı: {fs/1e6:.3f} MHz",
        "",
        f"Ortalama güç: {10*np.log10(mean_power+1e-20):.2f} dB",
        f"Tepe güç: {10*np.log10(peak_power+1e-20):.2f} dB",
        f"Güç dinamik aralığı: {dr_db:.2f} dB",
        "",
        f"I ort/std: {np.mean(np.real(sig)):.4f} / {np.std(np.real(sig)):.4f}",
        f"Q ort/std: {np.mean(np.imag(sig)):.4f} / {np.std(np.imag(sig)):.4f}",
        f"I/Q korelasyon: {np.corrcoef(np.real(sig), np.imag(sig))[0,1]:.4f}",
        "",
        f"Ortalama genlik: {np.mean(np.abs(sig)):.4f}",
        f"Genlik std: {np.std(np.abs(sig)):.4f}",
        f"Faz std: {np.std(np.angle(sig)):.4f} rad"
    ]
    ax7.text(0.1, 0.9, '\n'.join(stats_text), fontfamily='monospace',
             fontsize=10, verticalalignment='top')

    plt.tight_layout()
    plt.show()

    # Konsola da yaz
    print("\n" + "="*50)
    print(title)
    print("="*50)
    for line in stats_text:
        print(line)

def main():
    parser = argparse.ArgumentParser(description='Tek .dot dosyası için TX Sinyal Analizi')
    parser.add_argument('file', help='.dot dosya yolu')
    parser.add_argument('--fs', type=float, default=1e6, help='Örnekleme hızı (Hz)')
    args = parser.parse_args()

    if not os.path.exists(args.file):
        print(f"Hata: {args.file} dosyası bulunamadı!")
        return

    sig = read_dot_file(args.file)
    analyze_signal(sig, args.fs, f"Sinyal Analizi - {os.path.basename(args.file)}")

if __name__ == "__main__":
    main()
