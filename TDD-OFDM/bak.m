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

    num_data_symbols = p.Nsym*p.dataScs; 
    bitStream_Length = num_data_symbols * log2(p.M);
    tx_bits = randi([0 1], bitStream_Length, 1);

    qam_symbols = qammod(tx_bits, p.M, 'InputType', 'bit', 'UnitAveragePower', true);

    ofdmGrid = zeros(p.Nfft, p.Nsym);

    for i = 1:p.Nsym
        if mod(i,2)
            ofdmGrid(p.pilotInd,i) = p.pilot1;
        else
            ofdmGrid(p.pilotInd,i) = p.pilot2;
        end

        ofdmGrid(p.dataInd,i) = qam_symbols((i-1)*p.dataScs+1:i*p.dataScs);
    end

    % --- 3. DC taşıyıcısını boş bırak ---
    dc_idx = p.Nfft/2 + 1;  
    ofdmGrid(dc_idx, :) = 0;

    % --- 4. IFFT + CP ekle ---
    txPayload = zeros(p.wformLength,1);
    sample_counter = 1;
    for i = 1:p.Nsym
        ifft_out = ifft(fftshift(ofdmGrid(:,i)), p.Nfft) * sqrt(p.Nfft);
        cp_samples = ifft_out(end-p.cpLength+1:end);
        ofdm_symbol = [cp_samples; ifft_out];
        txPayload(sample_counter:sample_counter+length(ofdm_symbol)-1) = ofdm_symbol;
        sample_counter = sample_counter + length(ofdm_symbol);
    end
    txPayload = txPayload / max(abs(txPayload));

    txFrame = [p.sync; txPayload];

    % --- 6. 1 ms padding hesapla ---
    total_samples_1ms = round(1e-3 * p.sample_rate); 
    N_padding = total_samples_1ms - length(txFrame);
    if N_padding > 0
        txFrame = [txFrame; zeros(N_padding,1)];
    elseif N_padding < 0
        warning('Frame uzunluğu 1 ms');
    end

    % --- 7. Çıktı ---
    tx = txFrame;
end


function [ber, num_errors, rx_bits] = ofdmRx(rx_data, tx_bits, p)

[corr, lags] = xcorr(rx_data, p.sync);          
[~, peak_index] = max(abs(corr));
fprintf('Preamble pik indeksi: %d\n', peak_index);

% OFDM çerçeve başlangıcı (preamble sonu)
actual_ofdm_start_point = lags(peak_index) + length(p.sync);

% OFDM çerçeve uzunluğu
expected_full_signal_length = p.wformLength;
if actual_ofdm_start_point < 1 || actual_ofdm_start_point + expected_full_signal_length - 1 > length(rx_data)
    warning('OFDM sembol başlangıcı aralık dışında!');
    ber = 1; num_errors = length(tx_bits); rx_bits = zeros(size(tx_bits));
    return;
end

rx_data_synced = rx_data(actual_ofdm_start_point : actual_ofdm_start_point + expected_full_signal_length - 1);

%% --- 2) Kaba CFO düzeltme (CP üzerinden tahmin) ---
phase_offsets = zeros(p.Nsym,1);
for i = 1:p.Nsym
    idx = (i-1)*(p.Nfft+p.cpLength) + 1;
    cp_samples = rx_data_synced(idx:idx+p.cpLength-1);
    end_samples = rx_data_synced(idx+p.Nfft:idx+p.Nfft+p.cpLength-1);
    phase_offsets(i) = median(angle(cp_samples .* conj(end_samples)));
end
residual_cfo = mean(phase_offsets) / (2*pi*p.Nfft) * p.sample_rate;
n = (0:length(rx_data_synced)-1).';
rx_data_cfo_corrected = rx_data_synced .* exp(-1j*2*pi*residual_cfo*n/p.sample_rate);
fprintf('CFO düzeltildi: %.2f Hz\n', residual_cfo);

%% --- 3) OFDM FFT ---
rx_grid = zeros(p.Nfft, p.Nsym);
for i = 1:p.Nsym
    b = (i-1)*(p.Nfft+p.cpLength) + p.cpLength + 1;
    e = b + p.Nfft - 1;
    rx_grid(:,i) = fftshift(fft(rx_data_cfo_corrected(b:e))/sqrt(p.Nfft));
end

%% --- 4) CPE düzeltme ---
for i = 1:p.Nsym
    pilots_ref = p.pilot1; if mod(i,2)==0, pilots_ref = p.pilot2; end
    Yp = rx_grid(p.pilotInd,i);
    cpe = angle(mean(Yp ./ pilots_ref));
    rx_grid(:,i) = rx_grid(:,i) * exp(-1j*cpe);
end

%% --- 5) Kanal tahmini ve eşitleme ---
channel_grid = zeros(p.Nfft,p.Nsym);
for i = 1:p.Nsym
    pilots_ref = p.pilot1; if mod(i,2)==0, pilots_ref = p.pilot2; end
    Yp = rx_grid(p.pilotInd,i);
    Hp = Yp ./ pilots_ref;
    all_act = min(p.dataInd):max(p.dataInd);
    H = interp1(p.pilotInd, Hp, all_act, 'linear', 'extrap');
    H_full = ones(p.Nfft,1); H_full(all_act) = H;
    channel_grid(:,i) = H_full;
end

rx_data_symbols_eq = zeros(p.dataScs,p.Nsym);
for i = 1:p.Nsym
    rx_data_symbols_eq(:,i) = rx_grid(p.dataInd,i) ./ channel_grid(p.dataInd,i);
end
rx_data_symbols_eq = rx_data_symbols_eq(:);

%% --- 6) QAM demodülasyonu ve BER ---
rx_bits = qamdemod(rx_data_symbols_eq, p.M, 'gray', ...
                   'OutputType','bit', 'UnitAveragePower', true);
num_errors = sum(rx_bits ~= tx_bits(1:numel(rx_bits)));
ber = num_errors / numel(rx_bits);
end

%% --- Parametreleri al ---
p = parametersOFDM();
[txFrame, txBits, ~] = ofdmTx_1ms(p);   % 1ms OFDM çerçeve üret
txFrame = single(txFrame);               % USRP için single format

%% --- TX Objesi ---
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

%% --- RX Objesi ---
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
    'EnableBurstMode',    true, ...
    'NumFramesInBurst',   10);

%% --- TX/RX Döngüsü ---
disp('TX/RX başlatılıyor. CTRL+C ile durdurabilirsiniz.');
figure('Name','RX Spektrum');

while true
    % --- 1. TX gönder ---
    step(txObj, txFrame);

    % --- 2. RX ile gelen örnekleri oku ---
    [rxData, len] = step(rxObj);

    if len > 0
        % --- 3. ZC Preamble ile senkronizasyon ve OFDM çöz ---
        try
            [ber, num_errors, rxBits] = ofdmRx(rxData, txBits, p);

            % Eğer BER NaN veya 1 ise preamble bulunamamış demektir
            if ~isempty(ber) && ber <= 1
                fprintf('Preamble bulundu! BER: %.4f, Hata sayısı: %d\n', ber, num_errors);

                % --- 4. RX spektrumu pwelch ile çiz ---
                clf; % figure temizle
                pwelch(rxData, 1024, 512, 1024, p.sample_rate, 'centered');
                title('RX Spektrumu');
                drawnow;
            else
                % Preamble gelmedi, atla
                fprintf('Preamble gelmedi, bekleniyor...\n');
            end
        catch
            % Hata olursa (ör. preamble yok), atla
            fprintf('Preamble gelmedi veya analiz başarısız.\n');
            continue;
        end
    end

    % --- 5. 1 ms duraklama ---
    pause(1);
end
