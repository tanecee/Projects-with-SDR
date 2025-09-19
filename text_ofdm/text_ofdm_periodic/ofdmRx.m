function [rx_bits, rx_data_symbols] = ofdmRx(rx_data, p)
    phase_offsets = zeros(p.Nsym, 1);
    for i = 1:p.Nsym
        sym_start_idx = (i-1)*(p.Nfft + p.cpLength) + 1;
        cp_samples = rx_data(sym_start_idx:sym_start_idx+p.cpLength-1);
        end_samples = rx_data(sym_start_idx+p.Nfft:sym_start_idx+p.Nfft+p.cpLength-1);
        phase_diff = median(angle(cp_samples .* conj(end_samples)));
        phase_offsets(i) = phase_diff;
    end
    residual_cfo = mean(phase_offsets) / (2 * pi * p.Nfft) * p.sample_rate;
    t = (0:length(rx_data)-1).' ./ p.sample_rate;
    rx_data = rx_data .* exp(-1j * 2 * pi * residual_cfo * t);
    
    rx_data = rx_data / max(abs(rx_data));
    
    rx_grid = zeros(p.Nfft, p.Nsym);
    for i = 1:p.Nsym
        sym_start_idx = (i-1)*(p.Nfft + p.cpLength) + p.cpLength + 1;
        symbol_data = rx_data(sym_start_idx:sym_start_idx+p.Nfft-1);
        fft_output = fftshift(fft(symbol_data) / sqrt(p.Nfft));
        rx_grid(:, i) = fft_output;
    end
    
    channel_grid = zeros(p.Nfft, p.Nsym);
    for i = 1:p.Nsym
        pilot_sequence = p.pilot1;
        if mod(i,2)==0
            pilot_sequence = p.pilot2;
        end
        rx_pilots = rx_grid(p.pilotInd, i);
        channel_est = rx_pilots ./ pilot_sequence;
        all_indices = min(p.dataInd):max(p.dataInd);
        channel_grid(all_indices,i) = interp1(p.pilotInd,channel_est,all_indices,'linear','extrap');
    end
    
    rx_data_symbols = zeros(p.dataScs, p.Nsym);
    for i = 1:p.Nsym
        rx_data_symbols(:, i) = rx_grid(p.dataInd, i) ./ channel_grid(p.dataInd, i);
    end
    rx_data_symbols = rx_data_symbols(:);
    
    rx_bits = zeros(2*length(rx_data_symbols),1);
    for k = 1:length(rx_data_symbols)
        bit1 = real(rx_data_symbols(k)) < 0; 
        bit2 = imag(rx_data_symbols(k)) < 0; 
        rx_bits(2*k-1) = double(bit1);
        rx_bits(2*k)   = double(bit2);
    end
end
