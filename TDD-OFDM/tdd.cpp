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
#include <cstddef>
#include <numeric>

// Fonksiyon: Sinyali ikili dosyaya kaydetme
void saveToFile(const std::vector<std::complex<float>>& signal, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Dosya acilamadi: " << filename << std::endl;
        return;
    }
    for (const auto& sample : signal) {
        float real_part = sample.real();
        float imag_part = sample.imag();
        file.write(reinterpret_cast<const char*>(&real_part), sizeof(float));
        file.write(reinterpret_cast<const char*>(&imag_part), sizeof(float));
    }
    file.close();
    std::cout << "Sinyal " << filename << " dosyasina kaydedildi. Toplam ornek: " << signal.size() << std::endl;
}

// Yapı: OFDM Parametreleri
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
    int M = 4; // QPSK
    double sample_rate = 20e6;

    OFDMParameters();
};

// Fonksiyon: Zadoff-Chu senkronizasyon dizisi oluşturma
std::vector<std::complex<float>> generateSyncSequence(int u, int N) {
    std::vector<std::complex<float>> seq(N);
    for (int n = 0; n < N; n++) {
        double phase = -M_PI * u * n * (n + 1) / N;
        seq[n] = std::complex<float>(std::cos(phase), std::sin(phase));
    }
    return seq;
}

// Fonksiyon: Data alt-taşıyıcı endekslerini oluşturma
std::vector<size_t> buildDataInd(size_t Nfft, size_t dataScs) {
    std::vector<size_t> idx;
    idx.reserve(dataScs);
    for (size_t i = 1; i <= 2 * dataScs; i += 2) {
        idx.push_back(i + Nfft / 4);
    }
    return idx;
}

// Fonksiyon: Pilot alt-taşıyıcı endekslerini oluşturma
std::vector<size_t> buildPilotInd(size_t Nfft, size_t pilotScs) {
    std::vector<size_t> idx;
    idx.reserve(pilotScs);
    for (size_t i = 2; i <= 2 * pilotScs; i += 2) {
        idx.push_back(i + Nfft / 4);
    }
    return idx;
}

// OFDMParameters yapı kurucusu
OFDMParameters::OFDMParameters() {
    dataInd = buildDataInd(Nfft, dataScs);
    pilotInd = buildPilotInd(Nfft, pilotScs);

    const size_t DC_1B = Nfft / 2 + 1;
    pilotInd.erase(std::remove(pilotInd.begin(), pilotInd.end(), DC_1B), pilotInd.end());
    dataInd.erase(std::remove(dataInd.begin(), dataInd.end(), DC_1B), dataInd.end());

    pilotScs = pilotInd.size();
    dataScs = dataInd.size();

    pilot1.resize(pilotScs);
    pilot2.resize(pilotScs);
    for (size_t k = 0; k < pilotScs; ++k) {
        pilot1[k] = (k % 2 == 0) ? std::complex<float>(+1.f, 0.f)
                                 : std::complex<float>(-1.f, 0.f);
        pilot2[k] = -pilot1[k];
    }

    sync = generateSyncSequence(8, 255);
}

// Fonksiyon: Rastgele bit dizisi oluşturma
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

// Fonksiyon: QPSK modülasyonu
std::vector<std::complex<float>> qpskMod(const std::vector<int>& bits) {
    std::vector<std::complex<float>> symbols(bits.size() / 2);
    const float scale = 1.0f / std::sqrt(2.0f);

    for (size_t i = 0; i < bits.size(); i += 2) {
        int bit1 = bits[i];
        int bit2 = bits[i + 1];

        float real_part = (bit1 == 0) ? 1.0f : -1.0f;
        float imag_part = (bit2 == 0) ? 1.0f : -1.0f;

        symbols[i / 2] = std::complex<float>(real_part * scale, imag_part * scale);
    }

    return symbols;
}

// Fonksiyon: IFFT (Kaba bir FFT uygulaması)
std::vector<std::complex<float>> IFFT(const std::vector<std::complex<float>>& X) {
    size_t N = X.size();
    std::vector<std::complex<float>> x(N);
    for (size_t n = 0; n < N; ++n) {
        std::complex<float> sum(0.f, 0.f);
        for (size_t k = 0; k < N; ++k) {
            float ang = 2.0f * M_PI * float(n) * float(k) / float(N);
            sum += X[k] * std::complex<float>(std::cos(ang), std::sin(ang));
        }
        x[n] = sum / static_cast<float>(N);
    }
    return x;
}

// Fonksiyon: FFT kaydırması
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

// Fonksiyon: OFDM sinyali oluşturma
std::vector<std::complex<float>> generateOFDMSignal(const OFDMParameters& p) {
    size_t num_data_symbols = p.Nsym * p.dataScs;
    size_t bitStream_Length = num_data_symbols * 2;
    auto tx_bits = generateRandomBits(bitStream_Length);

    auto qam_symbols = qpskMod(tx_bits);

    std::vector<std::complex<float>> ofdmGrid(p.Nfft * p.Nsym, 0.0f);

    for (size_t symbol_idx = 0; symbol_idx < p.Nsym; symbol_idx++) {
        const auto& pilot = (symbol_idx % 2 == 0) ? p.pilot1 : p.pilot2;

        for (size_t j = 0; j < p.pilotScs; j++) {
            ofdmGrid[symbol_idx * p.Nfft + (p.pilotInd[j] - 1)] = pilot[j];
        }

        for (size_t j = 0; j < p.dataScs; j++) {
            ofdmGrid[symbol_idx * p.Nfft + (p.dataInd[j] - 1)] =
                qam_symbols[symbol_idx * p.dataScs + j];
        }
    }

    std::vector<std::complex<float>> tx_temp(p.wformLength);
    size_t sample_counter = 0;

    for (size_t symbol_idx = 0; symbol_idx < p.Nsym; symbol_idx++) {
        std::vector<std::complex<float>> symbol(p.Nfft);
        for (size_t i = 0; i < p.Nfft; i++) {
            symbol[i] = ofdmGrid[symbol_idx * p.Nfft + i];
        }

        auto shifted_symbol = fftShift(symbol);
        auto ifft_output = IFFT(shifted_symbol);

        for (auto& x : ifft_output) {
            x *= std::sqrt(static_cast<float>(p.Nfft));
        }

        for (size_t i = 0; i < p.cpLength; i++) {
            tx_temp[sample_counter++] = ifft_output[p.Nfft - p.cpLength + i];
        }
        for (size_t i = 0; i < p.Nfft; i++) {
            tx_temp[sample_counter++] = ifft_output[i];
        }
    }

    float max_val = 0.0f;
    for (const auto& sample : tx_temp) {
        max_val = std::max(max_val, std::abs(sample));
    }
    if (max_val > 0.f) {
        for (auto& sample : tx_temp) {
            sample /= max_val;
        }
    }

    std::vector<std::complex<float>> txSignal(p.sync.size() + tx_temp.size());
    for (size_t i = 0; i < p.sync.size(); i++) {
        txSignal[i] = p.sync[i];
    }
    for (size_t i = 0; i < tx_temp.size(); i++) {
        txSignal[i + p.sync.size()] = tx_temp[i];
    }

    return txSignal;
}

std::pair<double, size_t> estimate_cfo_and_sync(const std::vector<std::complex<float>>& rx_buffer, const std::vector<std::complex<float>>& zc_seq, double sample_rate) {
    size_t zc_len = zc_seq.size();
    if (rx_buffer.size() < zc_len) {
        return {0.0, 0};
    }

    // ZC dizisinin konjugesini al
    std::vector<std::complex<float>> zc_conj(zc_len);
    for (size_t i = 0; i < zc_len; ++i) {
        zc_conj[i] = std::conj(zc_seq[i]);
    }
    
    // Alinan sinyal ve ZC dizisi arasindaki korelasyon
    std::vector<double> correlation_mag_sq(rx_buffer.size() - zc_len + 1);
    for (size_t i = 0; i <= rx_buffer.size() - zc_len; ++i) {
        std::complex<float> sum(0.f, 0.f);
        for (size_t j = 0; j < zc_len; ++j) {
            sum += rx_buffer[i + j] * zc_conj[j];
        }
        correlation_mag_sq[i] = std::norm(sum);
    }
    
    // En büyük korelasyon tepe noktasini bul
    size_t peak_index = 0;
    double max_mag_sq = 0.0;
    for (size_t i = 0; i < correlation_mag_sq.size(); ++i) {
        if (correlation_mag_sq[i] > max_mag_sq) {
            max_mag_sq = correlation_mag_sq[i];
            peak_index = i;
        }
    }
    
    // CFO tahmini için iki yarım ZC dizisinin faz farkını kullan
    std::complex<float> r1(0.f, 0.f);
    std::complex<float> r2(0.f, 0.f);
    size_t half_zc = zc_len / 2;

    for(size_t i = 0; i < half_zc; ++i) {
        r1 += rx_buffer[peak_index + i] * std::conj(zc_seq[i]);
    }

    for(size_t i = half_zc; i < zc_len; ++i) {
        r2 += rx_buffer[peak_index + i] * std::conj(zc_seq[i]);
    }

    double phase_diff = std::arg(r2) - std::arg(r1);
    double cfo_est = -phase_diff / (M_PI * zc_len) * sample_rate;
    
    std::cout << "Estimated CFO: " << cfo_est << " Hz" << std::endl;
    std::cout << "Estimated symbol start index: " << peak_index << std::endl;

    return {cfo_est, peak_index};
}

// UHD ana fonksiyonu
int UHD_SAFE_MAIN(int argc, char* argv[]) {
    // USRP ayarları
    std::string device_args = "serial=33A6D68";
    std::string subdev = "A:A";
    std::string ant = "TX/RX";
    std::string antr = "RX2";
    std::string ref = "external";

    OFDMParameters p;
    double rate = p.sample_rate;
    double freq = 3.1e9;
    double gain = 50;

    // USRP'yi başlatma
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(device_args);

    // TX ayarları
    usrp->set_tx_rate(rate);
    usrp->set_tx_freq(freq);
    usrp->set_tx_gain(gain);
    usrp->set_tx_subdev_spec(subdev);
    usrp->set_tx_antenna(ant);

    // RX ayarları (TX ile aynı)
    usrp->set_rx_rate(rate);
    usrp->set_rx_freq(freq);
    usrp->set_rx_gain(gain-10);
    usrp->set_rx_subdev_spec(subdev);
    usrp->set_rx_antenna(antr);

    usrp->set_time_now(uhd::time_spec_t(0.0));

    uhd::stream_args_t stream_args("fc32", "sc16");
    
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

    // OFDM sinyali oluştur
    std::vector<std::complex<float>> ofdmSignal = generateOFDMSignal(p);
    const size_t spb = ofdmSignal.size();

    std::cout << "OFDM sinyali olusturuldu. Toplam ornek: " << spb << std::endl;
    std::cout << "Ornek hizi: " << rate / 1e6 << " MHz" << std::endl;
    std::cout << "Merkez frekans: " << freq / 1e9 << " GHz" << std::endl;
    std::cout << "WF:" << p.wformLength <<std::endl;

    saveToFile(ofdmSignal, "tx.dot");
    
    uhd::time_spec_t next_burst_time = usrp->get_time_now() + uhd::time_spec_t(2.0);

    double burst_duration = static_cast<double>(ofdmSignal.size()) / rate;
    double guard_interval = 0.05;
    double period = 2 * burst_duration + 2 * guard_interval;

    std::cout << "Periyot: " << period << " saniye" << std::endl;

    bool first_rx_completed = false; 

    try {
        while (true) {
            // --- TX (İletim) İşlemi ---
            uhd::tx_metadata_t tx_md;
            tx_md.start_of_burst = true;
            tx_md.end_of_burst = false;
            tx_md.has_time_spec = true;
            tx_md.time_spec = next_burst_time;

            std::cout << "TX baslangic zamani: " << std::fixed << std::setprecision(6) << next_burst_time.get_real_secs() << " saniye" << std::endl;
            size_t num_tx_samps = tx_stream->send(&ofdmSignal.front(), ofdmSignal.size(), tx_md);
            
            uhd::time_spec_t tx_end_time = next_burst_time + uhd::time_spec_t(burst_duration);

            tx_md.start_of_burst = false;
            tx_md.end_of_burst = true;
            tx_md.has_time_spec = false;
            std::vector<std::complex<float>> empty_buffer;
            tx_stream->send(empty_buffer.data(), empty_buffer.size(), tx_md);
            std::cout << "Gonderilen ornek sayisi: " << num_tx_samps << std::endl;
            
            // --- RX (Alım) İşlemi ---
            uhd::time_spec_t rx_start_time = tx_end_time + uhd::time_spec_t(guard_interval);

            std::cout << "RX baslangic zamani: " << std::fixed << std::setprecision(6) << rx_start_time.get_real_secs() << " saniye" << std::endl;
            
            uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
            stream_cmd.num_samps = spb;
            stream_cmd.stream_now = false;
            stream_cmd.time_spec = rx_start_time;
            
            rx_stream->issue_stream_cmd(stream_cmd);
            
            std::vector<std::complex<float>> rx_buffer(spb);
            uhd::rx_metadata_t rx_md;
            size_t num_rx_samps = rx_stream->recv(&rx_buffer.front(), rx_buffer.size(), rx_md, 3.0);
            std::cout << "RX Metadata Time: " << rx_md.time_spec.get_real_secs() << std::endl;
            
            std::cout << "Alinan ornek sayisi: " << num_rx_samps << std::endl;
            
            if (rx_md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
                std::cerr << "RX hatasi: " << rx_md.strerror() << std::endl;
            } else if (!first_rx_completed) {
                // CFO ve sembol başlangıç indeksi tahmini
                auto results = estimate_cfo_and_sync(rx_buffer, p.sync, rate);
                double estimated_cfo = results.first;
                size_t start_index = results.second;

                // CFO düzeltmesini uygula ve sadece OFMD sembollerini al
                std::vector<std::complex<float>> corrected_signal(rx_buffer.size() - start_index);
                for (size_t i = 0; i < corrected_signal.size(); ++i) {
                    double phase_correction = 2.0 * M_PI * estimated_cfo * (static_cast<double>(i) / rate);
                    std::complex<float> rotation_factor(std::cos(phase_correction), -std::sin(phase_correction));
                    corrected_signal[i] = rx_buffer[start_index + i] * rotation_factor;
                }
                
                // Kaydedilecek sinyal boyutunu ayarlama
                if (corrected_signal.size() > p.wformLength) {
                    corrected_signal.resize(p.wformLength);
                } else {
                    std::cerr << "Uyari: Alinan sinyal beklenenden kisa!" << std::endl;
                }

                saveToFile(corrected_signal, "rx.dot");
                first_rx_completed = true;
                std::cout << "RX dosyası kaydedildi. Sonraki kayıtlar devre dışı." << std::endl;
            }
            
            // Bir sonraki döngü için zaman güncelleme
            uhd::time_spec_t rx_end_time = rx_start_time + uhd::time_spec_t(burst_duration);
            next_burst_time = rx_end_time + uhd::time_spec_t(guard_interval);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Hata: " << e.what() << std::endl;
    }

    return 0;
}
