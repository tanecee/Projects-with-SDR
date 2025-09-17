function p = parametersOFDM()

    p.Nfft=256; 
    p.Nsym=14;
    p.actScs=p.Nfft/2; 
    p.dataScs=p.actScs/2;
    p.pilotScs=p.actScs/2;
    p.dataInd= (1:2:2*p.dataScs)+ p.Nfft/4;
    p.pilotInd = (2:2:2*p.pilotScs) + p.Nfft/4;
    p.pilot1=repmat([1 ; -1],p.pilotScs/2,1);
    p.pilot2=repmat([-1 ; 1],p.pilotScs/2,1);
    p.sync=zadoffChuSeq(8,255);
    p.cpLength=p.Nfft/4; 
    p.wformLength=(p.Nfft+p.cpLength)*p.Nsym;
    p.M=4; 
    p.sample_rate=20e6;
end

function [tx, tx_bits, ofdmGrid] = ofdmTx_1ms(p)
    [tx_normal, tx_bits, ofdmGrid] = ofdmTx(p);
    
    total_samples_1ms = round(1e-3 * p.sample_rate);
    
    current_length = length(tx_normal);
    
    if current_length < total_samples_1ms
        tx = [tx_normal; zeros(total_samples_1ms - current_length, 1)];
    elseif current_length > total_samples_1ms
        tx = tx_normal(1:total_samples_1ms);
        warning('Sinyal 1 ms''den uzun, kesildi');
    else
        % Tam 1 ms
        tx = tx_normal;
    end
end

function [tx, tx_bits, ofdmGrid] = ofdmTx(p)
p=parametersOFDM;

num_data_symbols = p.Nsym*p.dataScs; %Toplam veri sembölü
bitStream_Length = num_data_symbols * log2(p.M); %Üretilen bit sayısı
tx_bits = randi([0 1], bitStream_Length, 1); 
%bit dizisi QAM sembollerine çevrilir ve ortalama sembol gücü 1'e ölçeklenir.
qam_symbols = qammod(tx_bits, p.M, 'InputType', 'bit', 'UnitAveragePower', true); 
ofdmGrid = zeros(p.Nfft, p.Nsym);

for i = 1:p.Nsym %Her sütun bir OFDM sembolünün örnekelrini temsil eder.
    if mod(i, 2)
        ofdmGrid(p.pilotInd, i) = p.pilot1;
    else
        ofdmGrid(p.pilotInd, i) = p.pilot2;
    end
    ofdmGrid(p.dataInd, i) = qam_symbols((i-1)*p.dataScs+1:i*p.dataScs); %Veri sembollerini dataInd konumlarına yerleştirilir
end

tx = zeros(p.wformLength, 1);
sample_counter = 1;

for i = 1:p.Nsym
    ifft_output = ifft(fftshift(ofdmGrid(:, i)), p.Nfft) * sqrt(p.Nfft);
    cp_samples = ifft_output(end-p.cpLength+1:end);
    ofdm_symbol = [cp_samples; ifft_output];
    tx(sample_counter:sample_counter+length(ofdm_symbol)-1) = ofdm_symbol;
    sample_counter = sample_counter + length(ofdm_symbol);
end

tx = tx / max(abs(tx));
tx = [p.sync; tx];
%fftshift + ifft: grid'i merkez-DC'den IFFT'nin beklediği sıralamaya uygun şekilde kaydırıp zaman domenine geçilir.
%*sqrt(Nfft): enerji ölçekleme (IFFT'nin 1/√N normuna göre sabitleme).
%CP: sembol sonundaki cpLength örneği başa kopyalanır.
%Dalga formunu normalize edip başına sync (ZC) preamble eklenir.
end


function [ber, num_errors, rx_bits, rx_data_symbols] = ofdmRx(rx_data, tx_bits)
    p = parametersOFDM();

    % Preamble detection
    [corr, lag] = xcorr(rx_data, p.sync);
    [~, max_idx] = max(abs(corr));
    start_idx = lag(max_idx) + length(p.sync) + 1;

    if start_idx < 1 || start_idx + p.wformLength - 1 > length(rx_data)
        warning('Senkronizasyon indeksi geçersiz.');
        start_idx = 1;
        if length(rx_data) < p.wformLength
            rx_data = [rx_data; zeros(p.wformLength - length(rx_data), 1)];
        else
            rx_data = rx_data(1:p.wformLength);
        end
    else
        rx_data = rx_data(start_idx:start_idx + p.wformLength - 1);
    end

    % CFO estimation
    phase_offsets = zeros(p.Nsym, 1);
    for i = 1:p.Nsym
        sym_start_idx = (i-1)*(p.Nfft + p.cpLength) + 1;
        cp_samples = rx_data(sym_start_idx:sym_start_idx+p.cpLength-1);
        end_samples = rx_data(sym_start_idx+p.Nfft:sym_start_idx+p.Nfft+p.cpLength-1);
        phase_diff = median(angle(cp_samples .* conj(end_samples)));
        phase_offsets(i) = phase_diff;
    end
    residual_cfo = mean(phase_offsets) / (2 * pi * p.Nfft) * p.sample_rate;
    t = (0:length(rx_data)-1).'./ p.sample_rate;
    rx_data = rx_data .* exp(-1j * 2 * pi * residual_cfo * t);

    % Normalize
    rx_data = rx_data / max(abs(rx_data));

    % FFT processing
    rx_grid = zeros(p.Nfft, p.Nsym);
    for i = 1:p.Nsym
        sym_start_idx = (i-1)*(p.Nfft + p.cpLength) + p.cpLength + 1;
        symbol_data = rx_data(sym_start_idx:sym_start_idx+p.Nfft-1);
        fft_output = fftshift(fft(symbol_data) / sqrt(p.Nfft));
        rx_grid(:, i) = fft_output;
    end

    % Channel estimation
    channel_grid = zeros(p.Nfft, p.Nsym);
    for i = 1:p.Nsym
        pilot_sequence = p.pilot1;
        if mod(i, 2) == 0
            pilot_sequence = p.pilot2;
        end
        rx_pilots = rx_grid(p.pilotInd, i);
        channel_est = rx_pilots ./ pilot_sequence;
        all_indices = min(p.dataInd):max(p.dataInd);
        channel_grid(all_indices, i) = interp1(p.pilotInd, channel_est, all_indices, 'linear', 'extrap');
    end

    % Data extraction and equalization
    rx_data_symbols = zeros(p.dataScs, p.Nsym);
    for i = 1:p.Nsym
        rx_data_symbols(:, i) = rx_grid(p.dataInd, i);
        rx_data_symbols(:, i) = rx_data_symbols(:, i) ./ channel_grid(p.dataInd, i);
    end

    rx_data_symbols = rx_data_symbols(:); 
    
    % Demodulation
    if nargin >= 2 && ~isempty(tx_bits)
        rx_bits = qamdemod(rx_data_symbols, p.M, 'OutputType', 'bit', 'UnitAveragePower', true);
        num_errors = sum(rx_bits ~= tx_bits);
        ber = num_errors / length(tx_bits);
    else
        rx_bits = [];
        num_errors = [];
        ber = [];
    end

    % Visualization
    if ~isempty(rx_data_symbols)
        figure(1);
        scatter(real(rx_data_symbols), imag(rx_data_symbols), 'b.');
        title('Received Constellation');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        grid on;
    end

    figure(2);
    plot(20*log10(abs(fftshift(fft(rx_data, 1024)))), 'b');
    title('Received Signal Spectrum');
    xlabel('Frequency Bin');
    ylabel('Magnitude (dB)');
    grid on;
end

function main_USRP_TXRX()
    % Parametreleri yükle
    p = parametersOFDM();
    
    % TX sinyali üret
    [txFrame, txBits, ~] = ofdmTx_1ms(p);
    txFrame = single(txFrame);

    % --- TX Objesi ---
    txObj = comm.SDRuTransmitter( ...
        'Platform',             'B210', ...
        'SerialNum',            '33A6D68', ...
        'MasterClockRate',      20e6, ...
        'CenterFrequency',      3.1e9, ...
        'Gain',                 70, ...
        'InterpolationFactor',  1, ...
        'ChannelMapping',       1, ...
        'EnableBurstMode',      false, ...
        'ClockSource',          'Internal');

    % --- RX Objesi ---
    rxObj = comm.SDRuReceiver( ...
        'Platform',           'B210', ...
        'SerialNum',          '33A6D68', ...
        'MasterClockRate',    20e6, ...
        'CenterFrequency',    3.1e9, ...
        'Gain',               30, ...
        'DecimationFactor',    1, ...
        'SamplesPerFrame',    20000, ...
        'OutputDataType',     'double', ...
        'ChannelMapping',     1, ...
        'EnableBurstMode',    false); % Burst mode kapalı

    % --- TX/RX Döngüsü ---
    disp('TX/RX başlatılıyor. CTRL+C ile durdurabilirsiniz.');
    
    % Figure'ları hazırla
    fig1 = figure('Name', 'RX Spektrum');
    fig2 = figure('Name', 'RX Constellation');
    
    frame_count = 0;
    preamble_detected = false;
    
    try
        while true
            % --- 1. TX gönder ---
            underrun = txObj(txFrame);
            if underrun
                warning('TX underrun occurred');
            end
            
            % --- 2. RX ile gelen örnekleri oku ---
            [rxData, dataLen, overrun] = rxObj();
            
            if overrun
                warning('RX overrun occurred');
            end
            
            if dataLen > 0
                frame_count = frame_count + 1;
                
                % --- 3. OFDM demodülasyonu ---
                try
                    [ber, num_errors, rxBits, rx_symbols] = ofdmRx(rxData, txBits);
                    
                    if ~isempty(ber) && ber <= 1
                        preamble_detected = true;
                        fprintf('Frame %d: Preamble bulundu! BER: %.4f, Hatalar: %d/%d\n', ...
                            frame_count, ber, num_errors, length(txBits));
                        
                        % --- 4. Spektrumu göster ---
                        figure(fig1);
                        pwelch(rxData, 1024, 512, 1024, p.sample_rate, 'centered');
                        title(sprintf('RX Spektrumu - Frame %d (BER: %.4f)', frame_count, ber));
                        drawnow;
                        
                    else
                        if preamble_detected
                            fprintf('Frame %d: Preamble kayboldu\n', frame_count);
                            preamble_detected = false;
                        end
                    end
                    
                catch ME
                    fprintf('Frame %d: Hata - %s\n', frame_count, ME.message);
                end
            end
            
            % --- 5. Küçük bekleme ---
            pause(0.1);
        end
        
    catch ME
        fprintf('Program durduruldu: %s\n', ME.message);
    end
    
    % Temizlik
    release(txObj);
    release(rxObj);
    fprintf('USRP bağlantıları kapatıldı.\n');
end

% Script'i çalıştır
main_USRP_TXRX();

