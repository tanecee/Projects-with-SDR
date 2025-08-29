#include <iio.h>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <thread>

int main(int argc, char* argv[]) {
    std::string uri = "usb:3.11.5";
    double freq = 3.1e9; // 3.1 GHz
    double gain = 20; // Kazanç, TX ile uyumlu
    size_t num_samples = static_cast<size_t>(rate); // 1 saniye için örnek sayısı
    size_t buf_size = 4096; // Buffer boyutu, verimli bir değer

    // IIO context oluştur
    iio_context* ctx = iio_create_context_from_uri(uri.c_str());
    if (!ctx) {
        std::cerr << "Context oluşturulamadı: " << uri << std::endl;
        return 1;
    }

    // AD9361 PHY cihazını bul
    iio_device* phy = iio_context_find_device(ctx, "ad9361-phy");
    if (!phy) {
        std::cerr << "ad9361-phy cihazı bulunamadı" << std::endl;
        iio_context_destroy(ctx);
        return 1;
    }

    // RX stream cihazını bul (cf-ad9361-lpc)
    iio_device* rx = iio_context_find_device(ctx, "cf-ad9361-lpc");
    if (!rx) {
        std::cerr << "cf-ad9361-lpc cihazı bulunamadı" << std::endl;
        iio_context_destroy(ctx);
        return 1;
    }

    // Frekans ayarla (RX LO)
    iio_channel* rx_lo = iio_device_find_channel(phy, "altvoltage0", true); // RX_LO
    if (!rx_lo) {
        std::cerr << "RX_LO kanalı bulunamadı" << std::endl;
        iio_context_destroy(ctx);
        return 1;
    }
    iio_channel_attr_write_longlong(rx_lo, "frequency", static_cast<long long>(freq));

    // Sample rate ayarla
    iio_channel* rx_samp = iio_device_find_channel(phy, "voltage0", false); // RX voltage0
    if (!rx_samp) {
        std::cerr << "RX voltage0 kanalı bulunamadı" << std::endl;
        iio_context_destroy(ctx);
        return 1;
    }
    iio_channel_attr_write_longlong(rx_samp, "sampling_frequency", static_cast<long long>(rate));

    // Gain control mode: manual
    iio_channel_attr_write(rx_samp, "gain_control_mode", "manual");

    // Hardware gain ayarla
    iio_channel_attr_write_longlong(rx_samp, "hardwaregain", static_cast<long long>(gain));

    // RX kanalları: voltage0 (I), voltage1 (Q)
    iio_channel* rx_i = iio_device_find_channel(rx, "voltage0", false);
    iio_channel* rx_q = iio_device_find_channel(rx, "voltage1", false);
    if (!rx_i || !rx_q) {
        std::cerr << "RX I/Q kanalları bulunamadı" << std::endl;
        iio_context_destroy(ctx);
        return 1;
    }
    iio_channel_enable(rx_i);
    iio_channel_enable(rx_q);

    // RX buffer oluştur
    iio_buffer* rx_buf = iio_device_create_buffer(rx, buf_size, false);
    if (!rx_buf) {
        std::cerr << "RX buffer oluşturulamadı" << std::endl;
        iio_context_destroy(ctx);
        return 1;
    }

    std::vector<std::complex<int16_t>> samples;
    samples.reserve(num_samples);

    std::cout << "Alım başlıyor... 1 saniye sürecek." << std::endl;

    size_t total_received = 0;
    while (total_received < num_samples) {
        int ret = iio_buffer_refill(rx_buf);
        if (ret < 0) {
            std::cerr << "Buffer refill hatası: " << ret << std::endl;
            break;
        }

        const int16_t* p_dat = static_cast<const int16_t*>(iio_buffer_start(rx_buf));
        size_t samples_in_buf = iio_buffer_end(rx_buf) - iio_buffer_start(rx_buf)) / (2 * sizeof(int16_t)); // I/Q per sample

        for (size_t j = 0; j < samples_in_buf && total_received < num_samples; ++j, ++total_received) {
            int16_t i = p_dat[2 * j];
            int16_t q = p_dat[2 * j + 1];
            samples.emplace_back(i, q);
        }

        iio_buffer_push(rx_buf);
    }

    std::cout << "Alım tamamlandı. Toplam örnek: " << samples.size() << std::endl;

    std::ofstream out("received.dat", std::ios::binary);
    if (!out) {
        std::cerr << "Dosya açılamadı: received.dat" << std::endl;
    } else {
        out.write(reinterpret_cast<const char*>(samples.data()), samples.size() * sizeof(std::complex<int16_t>));
        std::cout << "Veri received.dat dosyasına kaydedildi." << std::endl;
    }

    iio_buffer_destroy(rx_buf);
    iio_context_destroy(ctx);

    return 0;
}