% PlutoSDR - Tek seferlik 1 saniye IQ kaydı (.dot) + Analiz

clear; clc;

% ---- PlutoSDR ayarları ----
rx = comm.SDRRxPluto( ...
    'RadioID','sn:1044734c9605001104000100d5934f698c', ...
    'CenterFrequency', 3.1e9, ...
    'BasebandSampleRate', 1e6, ...
    'SamplesPerFrame', 1024*100, ...
    'OutputDataType','double', ...
    'GainSource','Manual');

rx.ShowAdvancedProperties = true;
rx.Gain = 20;   % Manual modda kazanç

% ---- Kayıt ayarları ----
recordDuration = 1;            % saniye
Fs = rx.BasebandSampleRate;
N  = Fs * recordDuration;      % toplam örnek sayısı
outDir = 'recorded_data';
if ~exist(outDir,'dir'); mkdir(outDir); end

fprintf('Pluto RX basladi: Fc=%.1f MHz, Fs=%.1f MHz\n', ...
    rx.CenterFrequency/1e6, Fs/1e6);

% ---- 1 saniyelik veri topla ----
buf = complex(zeros(N,1));
idx = 0;
while idx < N
    x = rx();
    take = min(length(x), N-idx);
    buf(idx+1:idx+take) = x(1:take);
    idx = idx + take;
end

% ---- Dosyaya yaz ----
ts = datetime('now');
ts_str = datestr(ts,'yyyymmdd_HHMMSS');
outFile = fullfile(outDir, sprintf('usrp_rx_%s.dot', ts_str));

fid = fopen(outFile,'wb');
iq  = single([real(buf) imag(buf)]');   % 2xN -> interleaved
fwrite(fid, iq, 'float32');
fclose(fid);

fprintf('Kayit tamamlandi: %s (%d örnek, %.1f saniye)\n', ...
    outFile, N, recordDuration);

release(rx);

analyzeRecording(buf, Fs, ts);


function analyzeRecording(data, fs, timeStamp)
    % Detaylı sinyal analizi
    
    figure('Position', [100, 100, 1200, 800]);
    
    % 1. I/Q bileşenleri (zaman domeni)
    subplot(3,2,1);
    t = (0:length(data)-1)/fs;
    plot(t, real(data), 'b', t, imag(data), 'r');
    title('I ve Q Bileşenleri');
    xlabel('Zaman (s)'); ylabel('Genlik');
    legend('I','Q'); grid on;
    
    % 2. Spektrum (PSD)
    subplot(3,2,2);
    [Pxx,f] = pwelch(data,1024,512,1024,fs,'centered');
    plot(f/1e3,10*log10(Pxx));
    title('Güç Spektral Yoğunluğu');
    xlabel('Frekans (kHz)'); ylabel('dB/Hz'); grid on;
    xlim([-50 50]); % dar bant gösterim
    
    % 3. Constellation
    subplot(3,2,3);
    scatter(real(data), imag(data),1,'filled');
    title('Constellation Diagram');
    xlabel('I'); ylabel('Q'); axis equal; grid on;
    
    % 4. Histogram
    subplot(3,2,4);
    histogram(real(data),50,'Normalization','pdf'); hold on;
    histogram(imag(data),50,'Normalization','pdf');
    title('I/Q Histogram'); xlabel('Genlik'); ylabel('PDF');
    legend('I','Q'); grid on;
    
    % 5. Otokorelasyon
    subplot(3,2,5);
    [corrI,lags] = xcorr(real(data),100,'coeff');
    [corrQ,~]    = xcorr(imag(data),100,'coeff');
    plot(lags/fs*1000,corrI,'b',lags/fs*1000,corrQ,'r');
    title('Otokorelasyon');
    xlabel('Gecikme (ms)'); ylabel('Korelasyon'); grid on;
    legend('I','Q');
    
    % 6. İstatistikler
    subplot(3,2,6);
    stats = {
        sprintf('Ortalama Güç: %.2f dB', 10*log10(mean(abs(data).^2)))
        sprintf('Tepe Güç: %.2f dB', 10*log10(max(abs(data).^2)))
        sprintf('Dinamik Aralık: %.2f dB', 10*log10(max(abs(data).^2)/mean(abs(data).^2)))
        sprintf('I Std: %.4f', std(real(data)))
        sprintf('Q Std: %.4f', std(imag(data)))
        sprintf('SNR (tahmini): %.2f dB', estimateSNR(data))
    };
    text(0.1,0.8,stats,'FontSize',9);
    axis off; title('İstatistikler');
    
    sgtitle(sprintf('USRP/Pluto Sinyal Analizi - %s',datestr(timeStamp)),'FontSize',14);
end

function snr = estimateSNR(data)
    % Basit SNR tahmini
    signal_power = mean(abs(data).^2);
    noise_floor  = median(abs(data).^2); % Noise floor tahmini
    snr = 10*log10(signal_power/noise_floor);
end