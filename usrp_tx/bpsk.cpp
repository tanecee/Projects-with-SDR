#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/thread.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/types/time_spec.hpp>
#include <uhd/types/metadata.hpp>
#include <complex>
#include <vector>
#include <iostream>
#include <thread>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <random>
#include <algorithm>
#include <fstream>

void saveToFile(const std::vector<std::complex<float>>& signal, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    for (const auto& sample : signal) {
        float real_part = sample.real();
        float imag_part = sample.imag();
        file.write(reinterpret_cast<const char*>(&real_part), sizeof(float));
        file.write(reinterpret_cast<const char*>(&imag_part), sizeof(float));
    }
    file.close();
    std::cout << "Sinyal " << filename << " dosyasına kaydedildi. Toplam örnek: " << signal.size() << std::endl;
}

struct OFDMParameters {
    size_t Nfft = 256;
    size_t Nsym = 14;
    size_t actScs = Nfft / 2;
    size_t dataScs = actScs / 2;
    size_t pilotScs = actScs / 2;
    std::vector<size_t> dataInd;
    std::vector<size_t> pilotInd;
    std::vector<std::complex<float>> pilot1;
    std::vector<std::complex<float>> pilot2;
    std::vector<std::complex<float>> sync;
    size_t cpLength = Nfft / 4;
    size_t wformLength = (Nfft + cpLength) * Nsym;
    int M = 4;                 // <<< Varsayılan BPSK (2). QPSK için 4 yapın.
    double sample_rate = 20e6;

    OFDMParameters();
};

std::vector<std::complex<float>> generateSyncSequence(int u, int N) {
    std::vector<std::complex<float>> seq(N);
    for (int n = 0; n < N; n++) {
        double phase = -M_PI * u * n * (n + 1) / N;
        seq[n] = std::complex<float>(std::cos(phase), std::sin(phase));
    }
    return seq;
}

OFDMParameters::OFDMParameters() {
    // Basit bir orta bant ayrımı (örnek amaçlı): DC etrafında data/pilot indeksleri
    for (size_t i = 1; i <= 2 * dataScs; i += 2) {
        dataInd.push_back(i + Nfft / 4);
    }
    for (size_t i = 2; i <= 2 * pilotScs; i += 2) {
        pilotInd.push_back(i + Nfft / 4);
    }

    pilot1.resize(pilotScs);
    pilot2.resize(pilotScs);
    for (size_t i = 0; i < pilotScs / 2; i++) {
        pilot1[2*i]   = std::complex<float>(1.0f, 0.0f);
        pilot1[2*i+1] = std::complex<float>(-1.0f, 0.0f);
        pilot2[2*i]   = std::complex<float>(-1.0f, 0.0f);
        pilot2[2*i+1] = std::complex<float>(1.0f, 0.0f);
    }

    sync = generateSyncSequence(8, 255);
}

std::vector<int> generateRandomBits(size_t length) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    std::vector<int> bits(length);
    for (size_t i = 0; i < length; i++) {
        bits[i] = dis(gen);
    }
    return bits;
}

// --- Modülasyonlar ---
std::vector<std::complex<float>> bpskMod(const std::vector<int>& bits) {
    std::vector<std::complex<float>> symbols(bits.size());
    for (size_t i = 0; i < bits.size(); i++) {
        float val = (bits[i] == 0) ? 1.0f : -1.0f;
        symbols[i] = std::complex<float>(val, 0.0f); // reel eksen
    }
    return symbols;
}

std::vector<std::complex<float>> qpskMod(const std::vector<int>& bits) {
    // Her sembolde 2 bit
    if (bits.size() % 2 != 0) {
        throw std::runtime_error("QPSK: Bit sayısı 2'nin katı olmalı.");
    }
    std::vector<std::complex<float>> symbols(bits.size() / 2);
    const float scale = 1.0f / std::sqrt(2.0f);

    for (size_t i = 0; i < bits.size(); i += 2) {
        int b0 = bits[i];
        int b1 = bits[i+1];
        float i_part = (b0 == 0) ? 1.0f : -1.0f;
        float q_part = (b1 == 0) ? 1.0f : -1.0f;
        symbols[i/2] = std::complex<float>(i_part * scale, q_part * scale);
    }
    return symbols;
}

// Genel modülatör (BPSK/QPSK)
std::vector<std::complex<float>> modulate(const std::vector<int>& bits, int M) {
    if (M == 2) {
        return bpskMod(bits);
    } else if (M == 4) {
        return qpskMod(bits);
    } else {
        throw std::runtime_error("Desteklenmeyen M: Sadece M=2 (BPSK) veya M=4 (QPSK) desteklenir.");
    }
}

std::vector<std::complex<float>> IFFT(const std::vector<std::complex<float>>& input) {
    size_t N = input.size();
    std::vector<std::complex<float>> output(N);

    for (size_t k = 0; k < N; k++) {
        std::complex<float> sum(0.0f, 0.0f);
        for (size_t n = 0; n < N; n++) {
            float angle = 2.0f * M_PI * k * n / N;
            std::complex<float> twiddle(std::cos(angle), std::sin(angle));
            sum += input[n] * twiddle;
        }
        output[k] = sum;
    }
    return output;
}


std::vector<std::complex<float>> fftShift(const std::vector<std::complex<float>>& input) {
    size_t N = input.size();
    std::vector<std::complex<float>> output(N);
    size_t half = N / 2;

    for (size_t i = 0; i < half; i++) {
        output[i] = input[i + half];
        output[i + half] = input[i];
    }
    return output;
}

std::vector<std::complex<float>> generateOFDMSignal(const OFDMParameters& p) {
    // Sembol başına bit sayısı
    int bits_per_symbol = (p.M == 2) ? 1 : (p.M == 4 ? 2 : 0);
    if (bits_per_symbol == 0) throw std::runtime_error("M sadece 2 veya 4 olabilir.");

    size_t num_data_symbols = p.Nsym * p.dataScs;                 // toplam data RE adedi
    size_t bitStream_Length = num_data_symbols * bits_per_symbol; // BPSK/QPSK’ye göre

    // Bit üret
    auto tx_bits = generateRandomBits(bitStream_Length);

    // Modülasyon
    auto mod_syms = modulate(tx_bits, p.M); // boyut: num_data_symbols

    // OFDM grid (Nfft x Nsym)
    std::vector<std::complex<float>> ofdmGrid(p.Nfft * p.Nsym, std::complex<float>(0.0f, 0.0f));

    // Pilot ve data yerleştirme
    for (size_t symbol_idx = 0; symbol_idx < p.Nsym; symbol_idx++) {
        const auto& pilot = (symbol_idx % 2 == 0) ? p.pilot1 : p.pilot2;

        // Pilotlar
        for (size_t j = 0; j < p.pilotScs; j++) {
            size_t idx = symbol_idx * p.Nfft + p.pilotInd[j];
            if (p.pilotInd[j] < p.Nfft) ofdmGrid[idx] = pilot[j];
        }

        // Data
        for (size_t j = 0; j < p.dataScs; j++) {
            size_t grid_idx = symbol_idx * p.Nfft + p.dataInd[j];
            size_t sym_idx  = symbol_idx * p.dataScs + j;
            if (p.dataInd[j] < p.Nfft && sym_idx < mod_syms.size()) {
                ofdmGrid[grid_idx] = mod_syms[sym_idx];
            }
        }
    }

    // Zaman domenine geçiş + CP ekleme
    std::vector<std::complex<float>> tx_temp(p.wformLength);
    size_t sample_counter = 0;

    for (size_t symbol_idx = 0; symbol_idx < p.Nsym; symbol_idx++) {
        std::vector<std::complex<float>> symbol(p.Nfft);
        for (size_t i = 0; i < p.Nfft; i++) {
            symbol[i] = ofdmGrid[symbol_idx * p.Nfft + i];
        }

        auto shifted_symbol = fftShift(symbol);
        auto ifft_output = IFFT(shifted_symbol);

        // Güç ölçekleme (mevcut kod mantığını koruduk)
        for (auto& x : ifft_output) x *= std::sqrt(static_cast<float>(p.Nfft));
        // Alternatif: güç birimlemek isterseniz şu iki satırı kullanın:
        // for (auto& x : ifft_output) x /= static_cast<float>(p.Nfft);

        // CP ekle
        for (size_t i = 0; i < p.cpLength; i++) {
            tx_temp[sample_counter++] = ifft_output[p.Nfft - p.cpLength + i];
        }
        for (size_t i = 0; i < p.Nfft; i++) {
            tx_temp[sample_counter++] = ifft_output[i];
        }
    }

    // Normalizasyon
    float max_val = 0.0f;
    for (const auto& sample : tx_temp) max_val = std::max(max_val, std::abs(sample));
    if (max_val > 0) {
        for (auto& sample : tx_temp) sample /= max_val;
    }

    // Sync preamble + payload
    std::vector<std::complex<float>> txSignal(p.sync.size() + tx_temp.size());
    for (size_t i = 0; i < p.sync.size(); i++) txSignal[i] = p.sync[i];
    for (size_t i = 0; i < tx_temp.size(); i++) txSignal[i + p.sync.size()] = tx_temp[i];

    return txSignal;
}

int UHD_SAFE_MAIN(int argc, char* argv[]) {
    // USRP ayarları
    std::string device_args = "serial=327AB4A";
    std::string subdev = "A:A";
    std::string ant = "TX/RX";
    std::string ref = "internal";

    OFDMParameters p;
    p.M = 4; // <<< BPSK (2). QPSK için 4 yapın.

    double rate = p.sample_rate; // OFDM örnekleme hızı
    double freq = 3.1e9;
    double gain = 40;
    double period = 0.1; // 100 ms periyot

    // USRP başlatma
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(device_args);
    usrp->set_tx_rate(rate);
    usrp->set_tx_freq(freq);
    usrp->set_tx_gain(gain);

    // UHD versiyonuna göre gerekebilir:
    // usrp->set_tx_subdev_spec(uhd::usrp::subdev_spec_t(subdev));
    usrp->set_tx_subdev_spec(subdev);

    usrp->set_tx_antenna(ant);
    usrp->set_time_now(uhd::time_spec_t(0.0));

    uhd::stream_args_t stream_args("fc32", "sc16");
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);

    // OFDM sinyali oluştur
    std::vector<std::complex<float>> ofdmSignal = generateOFDMSignal(p);
    const size_t spb = ofdmSignal.size();

    std::cout << "OFDM sinyali oluşturuldu. Toplam örnek: " << spb << std::endl;
    std::cout << "Modülasyon M=" << p.M << (p.M==2?" (BPSK)":" (QPSK)") << std::endl;
    std::cout << "Örnek hızı: " << rate / 1e6 << " MHz" << std::endl;
    std::cout << "Merkez frekans: " << freq / 1e9 << " GHz" << std::endl;

    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst = true;
    md.has_time_spec = true;

    uhd::time_spec_t time_stamp = usrp->get_time_now() + uhd::time_spec_t(2.0);
    std::cout << "OFDM sinyali gönderimi başlıyor..." << std::endl;

    try {
        saveToFile(ofdmSignal, "tx_signal_ofdm.dot");

        while (true) {
            uhd::time_spec_t current_time = usrp->get_time_now();
            while (current_time.get_real_secs() < time_stamp.get_real_secs() - 0.01) {
                current_time = usrp->get_time_now();
                std::this_thread::sleep_for(std::chrono::microseconds(100));
            }

            md.time_spec = time_stamp;

            double full_seconds = time_stamp.get_full_secs();
            double frac_seconds = time_stamp.get_frac_secs();
            double total_seconds = full_seconds + frac_seconds;

            std::cout << "OFDM sinyali gönderildi: "
                      << std::fixed << std::setprecision(6) << total_seconds
                      << " saniye" << std::endl;

            size_t num_tx_samps = tx_stream->send(&ofdmSignal.front(), ofdmSignal.size(), md);

            (void)num_tx_samps; // kullanılmıyor uyarısını bastır
            time_stamp += uhd::time_spec_t(period);
        }

    } catch (const std::exception& e) {
        std::cerr << "Hata: " << e.what() << std::endl;
    }

    md.end_of_burst = true;
    md.has_time_spec = false;
    tx_stream->send("", 0, md);

    std::cout << "OFDM sinyal gönderimi tamamlandı." << std::endl;
    return 0;
}
