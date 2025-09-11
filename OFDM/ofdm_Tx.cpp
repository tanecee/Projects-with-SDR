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
    std::vector<size_t> dataInd;   // 1-tabanlı değerler tutulacak
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

std::vector<std::complex<float>> generateSyncSequence(int u, int N) {
    std::vector<std::complex<float>> seq(N);
    for (int n = 0; n < N; n++) {
        double phase = -M_PI * u * n * (n + 1) / N;
        seq[n] = std::complex<float>(std::cos(phase), std::sin(phase));
    }
    return seq;
}
std::vector<size_t> buildDataInd(size_t Nfft, size_t dataScs) {
    std::vector<size_t> idx;
    idx.reserve(dataScs); // belleği önceden ayır
    for (size_t i = 1; i <= 2 * dataScs; i += 2) {  // 1,3,5,...
        idx.push_back(i + Nfft / 4);
    }
    return idx;
}

std::vector<size_t> buildPilotInd(size_t Nfft, size_t pilotScs) {
    std::vector<size_t> idx;
    idx.reserve(pilotScs);
    for (size_t i = 2; i <= 2 * pilotScs; i += 2) {  // 2,4,6,...
        idx.push_back(i + Nfft / 4);
    }
    return idx;
}

OFDMParameters::OFDMParameters() {
    // 1-tabanlı indeks dizileri
    dataInd  = buildDataInd (Nfft, dataScs);
    pilotInd = buildPilotInd(Nfft, pilotScs);

    // DC (1-tabanlı 129) çıkar
    const size_t DC_1B = Nfft/2 + 1;
    pilotInd.erase(std::remove(pilotInd.begin(), pilotInd.end(), DC_1B), pilotInd.end());
    dataInd .erase(std::remove(dataInd .begin(), dataInd .end(), DC_1B), dataInd .end());

    // Güncel uzunluklar
    pilotScs = pilotInd.size();
    dataScs  = dataInd.size();

    // Pilot dizileri (+1,-1) ve (-1,+1) olarak sıra ile
    pilot1.resize(pilotScs);
    pilot2.resize(pilotScs);
    for (size_t k = 0; k < pilotScs; ++k) {
        pilot1[k] = (k % 2 == 0) ? std::complex<float>(+1.f, 0.f)
                                 : std::complex<float>(-1.f, 0.f);
        pilot2[k] = -pilot1[k];
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

std::vector<std::complex<float>> qpskMod(const std::vector<int>& bits) {
    std::vector<std::complex<float>> symbols(bits.size() / 2);
    const float scale = 1.0f / std::sqrt(2.0f);
    
    for (size_t i = 0; i < bits.size(); i += 2) {
        int bit1 = bits[i];
        int bit2 = bits[i+1];
        
        float real_part = (bit1 == 0) ? 1.0f : -1.0f;
        float imag_part = (bit2 == 0) ? 1.0f : -1.0f;
        
        symbols[i/2] = std::complex<float>(real_part * scale, imag_part * scale);
    }
    
    return symbols;
}

std::vector<std::complex<float>> IFFT(const std::vector<std::complex<float>>& X) {
    size_t N = X.size();
    std::vector<std::complex<float>> x(N);
    for (size_t n = 0; n < N; ++n) {
        std::complex<float> sum(0.f,0.f);
        for (size_t k = 0; k < N; ++k) {
            float ang = 2.0f * M_PI * float(n) * float(k) / float(N);
            sum += X[k] * std::complex<float>(std::cos(ang), std::sin(ang)); // e^{+j2πnk/N}
        }
        x[n] = sum / static_cast<float>(N); // <-- normalize (MATLAB ifft gibi)
    }
    return x;
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
    size_t num_data_symbols = p.Nsym * p.dataScs;
    size_t bitStream_Length = num_data_symbols * 2; // QPSK için 2 bit/sembol
    auto tx_bits = generateRandomBits(bitStream_Length);
    
    auto qam_symbols = qpskMod(tx_bits);
    
    std::vector<std::complex<float>> ofdmGrid(p.Nfft * p.Nsym, 0.0f);
    
    for (size_t symbol_idx = 0; symbol_idx < p.Nsym; symbol_idx++) {
        const auto& pilot = (symbol_idx % 2 == 0) ? p.pilot1 : p.pilot2;

        // [FIX] 1-tabanlı indeks → dizi erişiminde -1
        for (size_t j = 0; j < p.pilotScs; j++) {
            ofdmGrid[symbol_idx * p.Nfft + (p.pilotInd[j] - 1)] = pilot[j];
        }
        
        // [FIX] 1-tabanlı indeks → dizi erişiminde -1
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
        auto ifft_output    = IFFT(shifted_symbol);

        // [ADD] MATLAB ile aynı ölçek (ifft(...)*sqrt(Nfft))
        for (auto& x : ifft_output) {
            x *= std::sqrt(static_cast<float>(p.Nfft));
        }
        
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
    for (const auto& sample : tx_temp) {
        max_val = std::max(max_val, std::abs(sample));
    }
    if (max_val > 0.f) {
        for (auto& sample : tx_temp) {
            sample /= max_val;
        }
    }
    
    // Sync preamble ekleme
    std::vector<std::complex<float>> txSignal(p.sync.size() + tx_temp.size());
    for (size_t i = 0; i < p.sync.size(); i++) {
        txSignal[i] = p.sync[i];
    }
    for (size_t i = 0; i < tx_temp.size(); i++) {
        txSignal[i + p.sync.size()] = tx_temp[i];
    }
    
    return txSignal;
}

int UHD_SAFE_MAIN(int argc, char* argv[]) {
    // USRP ayarları
    std::string device_args = "serial=327AB4A"; 
    std::string subdev = "A:A";
    std::string ant = "TX/RX"; 
    std::string ref = "internal";
    
    OFDMParameters p;
    double rate = p.sample_rate; // OFDM sample rate'i kullanılır
    double freq = 3.3e9;
    double gain = 70;     
    double period = 0.1; // 100 ms periyot
    
    // USRP'yi başlatma
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(device_args);
    
    usrp->set_tx_rate(rate);
    usrp->set_tx_freq(freq);
    usrp->set_tx_gain(gain);
    usrp->set_tx_subdev_spec(subdev);
    usrp->set_tx_antenna(ant);
    
    usrp->set_time_now(uhd::time_spec_t(0.0));
    
    uhd::stream_args_t stream_args("fc32", "sc16");
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);
    
    // OFDM sinyali oluştur
    std::vector<std::complex<float>> ofdmSignal = generateOFDMSignal(p);
    const size_t spb = ofdmSignal.size(); // Tüm OFDM sinyalini buffer olarak kullanılır
    
    std::cout << "OFDM sinyali oluşturuldu. Toplam örnek: " << spb << std::endl;
    std::cout << "Örnek hızı: " << rate / 1e6 << " MHz" << std::endl;
    std::cout << "Merkez frekans: " << freq / 1e9 << " GHz" << std::endl;
    
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst = true;
    md.has_time_spec = true;

    
    uhd::time_spec_t time_stamp = usrp->get_time_now() + uhd::time_spec_t(2.0);

    std::cout << "WF:" << p.wformLength <<std::endl;
    std::cout << "Length in seconds: " << static_cast<double>(p.wformLength) / p.sample_rate << " s" << std::endl;

    std::cout << "OFDM sinyali gönderimi başlıyor..." << std::endl;
    std::cout << "Zaman damgaları:" << std::endl;
    
    try {
        saveToFile(ofdmSignal, "tx_ofdm.dot");

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
            
            time_stamp += uhd::time_spec_t(period);
        }
        

    }
    catch (const std::exception& e) {
        std::cerr << "Hata: " << e.what() << std::endl;
    }

    md.end_of_burst = true;
    md.has_time_spec = false;
    tx_stream->send("", 0, md);
    
    std::cout << "OFDM sinyal gönderimi tamamlandı." << std::endl;
    
    return 0;
}
