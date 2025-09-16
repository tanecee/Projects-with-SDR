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

// Fonksiyon: CFO ve Senkronizasyon Tahmini (estimate_cfo_and_sync)
std::pair<double, size_t> estimate_cfo_and_sync(const std::vector<std::complex<float>>& rx_buffer, const std::vector<std::complex<float>>& zc_seq, double sample_rate) {
    size_t zc_len = zc_seq.size();
    if (rx_buffer.size() < zc_len) {
        std::cerr << "Hata: RX buffer Zadoff-Chu dizisinden kısa!" << std::endl;
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

    // --- BURAYA EKLENECEK KOD ---
    // correlation_mag_sq'yi dosyaya kaydet
    std::ofstream corr_file("correlation_output.txt");
    if (corr_file.is_open()) {
        for (size_t i = 0; i < correlation_mag_sq.size(); ++i) {
            corr_file << i << " " << correlation_mag_sq[i] << std::endl;
        }
        corr_file.close();
        std::cout << "Korelasyon sonuclari 'correlation_output.txt' dosyasina kaydedildi." << std::endl;
    } else {
        std::cerr << "Hata: 'correlation_output.txt' dosyasi acilamadi." << std::endl;
    }
    // --- EKLENECEK KOD BİTİŞİ ---
    
    // En büyük korelasyon tepe noktasini bul
    size_t peak_index = 0;
    double max_mag_sq = 0.0;
    for (size_t i = 0; i < correlation_mag_sq.size(); ++i) {
        if (correlation_mag_sq[i] > max_mag_sq) {
            max_mag_sq = correlation_mag_sq[i];
            peak_index = i;
        }
    }

    // Eğer peak_index çok büyükse, bu muhtemelen senkronizasyon hatasıdır.
    // Şimdilik sadece uyarı verelim.
    if (peak_index > rx_buffer.size() / 2 && rx_buffer.size() > 0) { // rx_buffer.size() > 0 kontrolü ekledim
        std::cerr << "Uyarı: Senkronizasyon dizisi beklenden çok geç algılandı. Peak index: " << peak_index << std::endl;
    }
    
    // CFO tahmini için iki yarım ZC dizisinin faz farkını kullan
    std::complex<float> r1(0.f, 0.f);
    std::complex<float> r2(0.f, 0.f);
    size_t half_zc = zc_len / 2;

    if (peak_index + zc_len > rx_buffer.size()) {
        std::cerr << "Hata: CFO tahmini için yeterli RX örneği yok. (peak_index + zc_len > rx_buffer.size())" << std::endl;
        return {0.0, peak_index};
    }

    for(size_t i = 0; i < half_zc; ++i) {
        r1 += rx_buffer[peak_index + i] * std::conj(zc_seq[i]);
    }

    for(size_t i = half_zc; i < zc_len; ++i) {
        r2 += rx_buffer[peak_index + i] * std::conj(zc_seq[i]);
    }

    double phase_diff = std::arg(r2) - std::arg(r1);
    // CFO formülü düzeltmesi: 2.0 * M_PI kullanılması önerilir.
    double cfo_est = -phase_diff * sample_rate / (2.0 * M_PI * static_cast<double>(half_zc)); 

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
            stream_cmd.num_samps = spb; // Gönderilen tüm sinyal uzunluğu (ZC + OFDM)
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

                // OFDM veri bölümünün başlangıç indeksi
                // ZC dizisi alındıktan sonra OFDM sembolleri başlıyor.
                // Eğer start_index ZC dizisinin başlangıcı ise, OFDM verisi start_index + p.sync.size()'dan başlar.
                size_t ofdm_data_start_idx = start_index + p.sync.size();

                // Eğer buffer'da yeterli örnek yoksa veya ofdm_data_start_idx çok büyükse
                if (ofdm_data_start_idx >= rx_buffer.size()) {
                    std::cerr << "Hata: OFDM veri başlangıç indeksi RX buffer boyutundan büyük veya eşit! "
                              << "ofdm_data_start_idx: " << ofdm_data_start_idx 
                              << ", rx_buffer.size(): " << rx_buffer.size() << std::endl;
                    first_rx_completed = true; // Hata durumunda tekrar kaydetmemek için
                    continue; // Bu iterasyonu atla
                }

                // İşlenecek maksimum örnek sayısı: ya buffer'da kalan örnekler ya da beklenen OFDM uzunluğu
                size_t actual_samples_to_process = std::min(p.wformLength, rx_buffer.size() - ofdm_data_start_idx);

                if (actual_samples_to_process == 0) {
                     std::cerr << "Hata: OFDM veri işlemek için hiç örnek yok. Senkronizasyon hatası çok ciddi." << std::endl;
                     first_rx_completed = true;
                     continue; 
                }

                std::vector<std::complex<float>> corrected_signal(actual_samples_to_process);

                // CFO düzeltmesini uygula
                for (size_t i = 0; i < actual_samples_to_process; ++i) {
                    double phase_correction = 2.0 * M_PI * estimated_cfo * (static_cast<double>(i) / rate);
                    std::complex<float> rotation_factor(std::cos(phase_correction), -std::sin(phase_correction));
                    corrected_signal[i] = rx_buffer[ofdm_data_start_idx + i] * rotation_factor;
                }
                
                // Kaydedilecek sinyal boyutunu ayarlama (gerekiyorsa)
                if (corrected_signal.size() < p.wformLength) {
                    std::cerr << "Uyari: Alinan ve duzeltilen OFDM sinyali beklenenden kisa! Beklenen: " << p.wformLength << ", Alinan: " << corrected_signal.size() << std::endl;
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

çıktı
OFDM sinyali olusturuldu. Toplam ornek: 4735
Ornek hizi: 20 MHz
Merkez frekans: 3.1 GHz
WF:4480
Sinyal tx.dot dosyasina kaydedildi. Toplam ornek: 4735
Periyot: 0.100474 saniye
TX baslangic zamani: 2.091254 saniye
Gonderilen ornek sayisi: 4735
RX baslangic zamani: 2.141491 saniye
RX Metadata Time: 2.141491
Alinan ornek sayisi: 4735
Korelasyon sonuclari 'correlation_output.txt' dosyasina kaydedildi.
Uyarı: Senkronizasyon dizisi beklenden çok geç algılandı. Peak index: 3592
Estimated CFO: 12894.064579 Hz
Estimated symbol start index: 3592
Uyari: Alinan ve duzeltilen OFDM sinyali beklenenden kisa! Beklenen: 4480, Alinan: 888
Sinyal rx.dot dosyasina kaydedildi. Toplam ornek: 888
RX dosyası kaydedildi. Sonraki kayıtlar devre dışı.
TX baslangic zamani: 2.191727 saniye
Gonderilen ornek sayisi: 4735
RX baslangic zamani: 2.241964 saniye
LLLRX Metadata Time: 2.241964
Alinan ornek sayisi: 4735
TX baslangic zamani: 2.292201 saniye
Gonderilen ornek sayisi: 4735

bu da correlation outputun birazının çıktısı:
0 1.46518e-06
1 5.66747e-05
2 4.46619e-05
3 7.32105e-05
4 0.000213096
5 0.000329132
6 0.000167028
7 0.00087012
8 0.000259158
9 0.000120539
10 1.5325e-05
11 0.00102425
12 4.03902e-05
13 0.000328149
14 9.74642e-05
15 0.000248931
16 0.000309294
17 8.64482e-05
18 0.000161164
19 7.97256e-05
20 0.000179677
21 0.000648373
22 0.000230209
23 2.26487e-05
24 6.68561e-05
25 0.00048367
26 0.000284618
27 4.11354e-05
28 0.000339072
29 3.6662e-05
30 5.52597e-05
31 0.000118718
32 8.12692e-06
33 0.000143533
34 0.000118938
35 0.00029095
36 7.81226e-05
37 0.00017507
38 7.60932e-05
39 0.000493602
40 5.89577e-05
41 0.000234793
42 4.35916e-05
43 0.000157792
44 9.2737e-05
45 0.000136394
46 0.000453182
47 6.57605e-05
48 0.000128372
49 0.000590715
50 0.000155466
51 0.000807593
52 0.000467576
53 8.21492e-05
54 0.000509053
55 0.00039269
56 0.000348644
57 0.000131798
58 0.000287294
59 6.43053e-05
60 0.000125685
61 0.000166112
62 9.44132e-05
63 2.62569e-05
64 6.5669e-05
65 0.000260489
66 5.66239e-05
67 0.00032206
68 0.000202139
69 0.000226552
70 7.72446e-05
71 8.41414e-05
72 0.000105724
73 8.84749e-05
74 0.000103567
75 1.30221e-05
76 0.000302944
77 3.73242e-06
78 2.25972e-05
79 0.000133242
80 0.000245073
81 2.73652e-05
82 0.000136201
83 0.000427267
84 0.000408021
85 8.75051e-05
86 0.00019862
87 0.00106538
88 0.000145389
89 0.000149752
90 6.19619e-05
91 5.16557e-05
92 2.34788e-06
93 9.37567e-05
94 0.000276363
95 0.000159078
96 0.000292866
97 0.000220056
98 0.00119792
99 0.00064399
100 0.000980618
101 0.000216949
102 0.000215549
103 0.000514858
104 0.000147997
105 0.000114536
106 8.43772e-05
107 0.000281822
108 0.000280229
109 2.86453e-06
110 8.0016e-05
111 0.000753775
112 0.000177695
113 3.16937e-05
114 0.00134309
115 0.000354523
116 0.000209897
117 1.35102e-05
118 0.000161171
119 0.000678294
120 0.000306155
121 0.000159704
122 0.000164405
123 0.000429333
124 0.000247531
125 4.01248e-05
126 5.94875e-05
127 0.00016118
128 0.000151416
129 0.000305838
130 0.000291446
131 0.00013534
132 0.000431103
133 1.2365e-05
134 0.000364986
135 0.000161631
136 0.000268218
137 0.000453919
138 0.000120458
139 0.000175454
140 8.20153e-05
141 2.75967e-06
142 5.83962e-05
143 8.0913e-05
144 0.000739026
145 0.000209047
146 0.000132651
147 0.000547067
148 0.000140951
149 0.000293388
150 9.69391e-05
151 0.00036327
152 0.000224249
153 0.000276098
154 0.000304073
155 0.00030704
156 2.45721e-05
157 1.33726e-05
158 0.000678299
159 0.000443065
160 5.78194e-05
161 0.000128304
162 0.00031768
163 1.55894e-05
164 0.000146393
165 0.000147408
166 0.000284172
167 5.4521e-05
168 2.56137e-05
169 0.000128007
170 4.41524e-05
171 1.84399e-05
172 8.15587e-05
