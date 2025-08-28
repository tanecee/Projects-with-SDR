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

int UHD_SAFE_MAIN(int argc, char* argv[]) {
    std::string device_args = "serial=327AB4A"; 
    std::string subdev = "A:A";
    std::string ant = "TX/RX"; 
    std::string ref = "internal";
    
    double rate = 1e6; 
    double freq = 3.1e9;
    double gain = 20;     
    double period = 0.1; // 100 ms periyot
    
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(device_args);
    
    usrp->set_tx_rate(rate);
    usrp->set_tx_freq(freq);
    usrp->set_tx_gain(gain);
    usrp->set_tx_subdev_spec(subdev);
    usrp->set_tx_antenna(ant);
    
    usrp->set_time_now(uhd::time_spec_t(0.0));
    
    uhd::stream_args_t stream_args("fc32", "sc16");
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);
    
    const double wave_freq = 1000; // 1 kHz
    const size_t spb = 1024; // Buffer başına örnek sayısı
    const size_t bufs_per_burst = static_cast<size_t>(std::ceil(period * rate / spb));

    const double LEAD = 0.15; // 150 ms güvenlik payı

    std::vector<std::complex<float>> buff(spb);
    for (size_t i = 0; i < spb; i++) {
        double t = i / rate;
        buff[i] = std::complex<float>(
            std::cos(2 * M_PI * wave_freq * t),
            std::sin(2 * M_PI * wave_freq * t)
        );
    }
    
    uhd::tx_metadata_t md;
    md.start_of_burst = false;
    md.end_of_burst   = false;
    md.has_time_spec  = false;
    
    uhd::time_spec_t time_stamp = usrp->get_time_now() + uhd::time_spec_t(LEAD);
    
    std::cout << "Sinyal gönderimi başlıyor..." << std::endl;
    std::cout << "Zaman damgaları:" << std::endl;
    
    try {
        while (true) {

            while (time_stamp < usrp->get_time_now() + uhd::time_spec_t(LEAD)) {
                time_stamp += uhd::time_spec_t(period);
            }

            md.time_spec      = time_stamp;
            md.start_of_burst = true;
            md.end_of_burst   = false;
            md.has_time_spec  = true;

            double full_seconds = time_stamp.get_full_secs();
            double frac_seconds = time_stamp.get_frac_secs();
            double total_seconds = full_seconds + frac_seconds;
            std::cout << "Sinyal gönderildi: " 
                      << std::fixed << std::setprecision(6) << total_seconds 
                      << " saniye" << std::endl;

            for (size_t k = 0; k < bufs_per_burst; k++) {
                tx_stream->send(buff.data(), buff.size(), md);
                md.start_of_burst = false;  
                md.has_time_spec  = false;  
            }

            md.end_of_burst = true;
            tx_stream->send("", 0, md);
            md.end_of_burst = false;

            time_stamp += uhd::time_spec_t(period);

            while (usrp->get_time_now() + uhd::time_spec_t(0.050) < time_stamp) {
                std::this_thread::sleep_for(std::chrono::milliseconds(5));
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Hata: " << e.what() << std::endl;
    }
    
    md.end_of_burst = true;
    md.has_time_spec = false;
    tx_stream->send("", 0, md);
    
    std::cout << "Sinyal gönderimi tamamlandı." << std::endl;
    
    return 0;
}
