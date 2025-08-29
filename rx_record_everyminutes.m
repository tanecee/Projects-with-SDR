% Basit Pluto RX -> 1 saniye IQ kaydet (.dot)
clear; clc;

% ---- PlutoSDR ayarları ----
rx = comm.SDRRxPluto( ...
    'RadioID','sn:1044734c9605001104000100d5934f698c', ...
    'CenterFrequency', 3.1e9, ...     % USRP ile aynı
    'BasebandSampleRate', 1e6, ...    % USRP ile aynı
    'SamplesPerFrame', 1024*100, ...  % ~0.1024 s / frame
    'OutputDataType','double', ...
    'GainSource','Manual');           % istersen 'AGC Fast Attack'

rx.ShowAdvancedProperties = true;
rx.Gain = 20;                         % Manual ise kazanç şart

% ---- Kayıt ayarları ----
recordDuration = 1;                   % 1 saniye
Fs = rx.BasebandSampleRate;
N  = Fs * recordDuration;             % hedef örnek sayısı
outDir = 'recorded_data';
if ~exist(outDir,'dir'); mkdir(outDir); end

fprintf('PlutoSDR Receiver hazir. Fc=%.1f MHz, Fs=%.1f MHz\n', ...
    rx.CenterFrequency/1e6, Fs/1e6);

packetCounter = 0;

try
    while true
        % 1 saniyelik tamponu doldur
        buf = complex(zeros(N,1));
        idx = 0;
        while idx < N
            x = rx();                      % bir frame oku (complex double)
            take = min(length(x), N-idx);  % hedefe taşmayalım
            buf(idx+1:idx+take) = x(1:take);
            idx = idx + take;
        end

        % Dosya adı (timestamp)
        ts = datestr(datetime('now'),'yyyymmdd_HHMMSS');
        outFile = fullfile(outDir, sprintf('usrp_rx_%s.dot', ts));

        % .dot'a yaz (float32 interleaved: I0 Q0 I1 Q1 ...)
        fid = fopen(outFile,'wb');
        iq  = single([real(buf) imag(buf)]');   % 2xN
        fwrite(fid, iq, 'float32');
        fclose(fid);

        packetCounter = packetCounter + 1;
        fprintf('[%d] Kaydedildi: %s (1 s, %d örnek)\n', packetCounter, outFile, N);

        % USRP'nin 100 ms periyodu ile çakışmaması için istersen küçük bekleme
        % pause(0.05);
    end

catch ME
    fprintf('Alım durduruldu: %s\n', ME.message);
end

release(rx);
fprintf('Toplam paket: %d\n', packetCounter);
