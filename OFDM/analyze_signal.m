function analyze_signal
    [file, path] = uigetfile('*.dot', 'Kaydedilmiş IQ dosyasını seçin');
    if isequal(file, 0), error('Dosya seçilmedi'); end
    filename = fullfile(path, file);
    iqData = read_iq_file(filename);
    p = parametersOFDM();
    analyze_and_plot(iqData, datetime('now'), p);
end

function iqData = read_iq_file(filename)
    fprintf('Dosya okunuyor: %s\n', filename);
    fid = fopen(filename, 'rb');
    if fid == -1, error('Dosya açılamadı: %s', filename); end
    rawData = fread(fid, Inf, 'float32');
    fclose(fid);

    I = rawData(1:2:end);
    Q = rawData(2:2:end);
    iqData = complex(I, Q);

    fprintf('Okunan örnek sayısı: %d\n', length(iqData));
end

function p = parametersOFDM()
    p.Nfft   = 256;        % IFFT/FFT boyutu
    p.Nsym   = 14;         % OFDM sembol sayısı
    p.actScs = p.Nfft/2;   % aktif alt taşıyıcı (128)
    p.dataScs  = p.actScs/2;   % 64
    p.pilotScs = p.actScs/2;   % 64

    p.dataInd  = (1:2:(2*p.dataScs))  + p.Nfft/4; % tekler
    p.pilotInd = (2:2:(2*p.pilotScs)) + p.Nfft/4; % çiftler

    % Pilot desenleri (uzunluk = pilotScs)
    p.pilot1 = repmat([1; -1], p.pilotScs/2, 1);
    p.pilot2 = repmat([-1; 1], p.pilotScs/2, 1);
    p.altPilots = true;   % tek/çift sembolde alternans

    p.cpLength   = p.Nfft/4;                 % 64
    p.wformLength= (p.Nfft + p.cpLength)*p.Nsym;

    p.M  = 4;
    p.fs = 20e6;     % örnekleme hızı
    p.Fc = 3.3e9;    % taşıyıcı
    p.useFFTshift = true;
    p.constPointsToShow = 1e9;  % saçılımda gösterilecek nokta sayısı
end

function analyze_and_plot(data, timeStamp, p)
    Nfft    = p.Nfft;
    cpLen   = p.cpLength;
    symLen  = Nfft + cpLen;
    Nsym    = p.Nsym;
    useShift= p.useFFTshift;
    Fc      = p.Fc;
    fs      = p.fs;

    dataInd  = p.dataInd(:);
    pilotInd = p.pilotInd(:);
    dcIdx    = Nfft/2 + 1;
    dataInd  = dataInd( dataInd>=1 & dataInd<=Nfft & dataInd~=dcIdx );
    pilotInd = pilotInd( pilotInd>=1 & pilotInd<=Nfft & pilotInd~=dcIdx );

    % CP-korelasyon
    searchSpan = symLen;
    NsymProbe  = Nsym;
    tauMetric  = zeros(searchSpan+1,1);
    for tau = 0:searchSpan
        acc = 0;
        for s = 1:NsymProbe
            i0 = tau + (s-1)*symLen + 1;
            i1 = i0 + symLen - 1;
            if i1 > numel(data), break; end
            sym_td = data(i0:i1);
            pre  = sym_td(1:cpLen);
            tail = sym_td(Nfft+1:end);
            acc  = acc + abs(sum(conj(pre).*tail));
        end
        tauMetric(tau+1) = acc;
    end
    [~,bestTauIdx] = max(tauMetric);
    fineOff = bestTauIdx-1;

    need = Nsym*symLen;
    if fineOff + need <= length(data)
        xWin = data(fineOff+1 : fineOff+need);
    else
        Nsym = floor((length(data)-fineOff)/symLen);
        need = Nsym*symLen;
        xWin = data(fineOff+1 : fineOff+need);
        warning('Kayıt kısa geldi. Nsym=%d olarak ayarlandı.', Nsym);
        if Nsym < 1, warning('Veri yetersiz.'); return; end
    end

    % CFO tahmini
    eps_hat = estimate_cfo_cp(xWin, Nfft, cpLen, Nsym);
    f_off   = eps_hat * fs / Nfft;
    fprintf('CFO: eps=%.6f (Δf), f_off=%.2f Hz\n', eps_hat, f_off);
    if isfinite(Fc)
        fprintf('CFO ~ %.2f ppm\n', 1e6 * f_off / Fc);
    end

    % Zaman alanında CFO düzelt
    n = (0:length(xWin)-1).';
    xCorr = xWin .* exp(-1j*2*pi*eps_hat * n / Nfft);

    % Grafikler 
    figure('Position',[100,100,1400,900],'Name','Kayıtlı Sinyal Analizi');

    % (1) Ham I/Q
    subplot(3,3,1);
    t = (0:length(xWin)-1)/fs;
    plot(t, real(xWin), 'b', t, imag(xWin), 'r'); grid on;
    title('I/Q (ham)'); xlabel('s');
    xlim([0, min(0.1, t(end))]); legend('I','Q');

    % (2) Spektrum Analyze
    subplot(3,3,2);
    [Pxx,f] = pwelch(xCorr,1024,512,1024,fs,'centered');
    plot(f/1e3,10*log10(Pxx)); grid on;
    title('Spektrum Analizi'); xlabel('kHz'); ylabel('dB/Hz');
    xlim([-fs/2/1e3, fs/2/1e3]);


    % (3) Histogram 
    subplot(3,3,4);
    histogram(real(xCorr),50,'Normalization','pdf'); hold on;
    histogram(imag(xCorr),50,'Normalization','pdf'); grid on;
    title('I/Q Histogram (CFO düzeltmeli)'); xlabel('Genlik'); ylabel('PDF'); legend('I','Q');

    % (4) Autocorrelation
    subplot(3,3,5);
    [corrI,lags] = xcorr(real(xCorr),100,'coeff');
    [corrQ,~]    = xcorr(imag(xCorr),100,'coeff');
    plot(lags/fs*1000,corrI,'b',lags/fs*1000,corrQ,'r'); grid on;
    title('Otokorelasyon (CFO düzeltmeli)'); xlabel('ms'); ylabel('\rho'); legend('I','Q');

    % (5) Statistics
    subplot(3,3,6); axis off;
    snr_est = estimateSNR(xCorr);
    lines = {
        sprintf('Örnek sayısı: %d', length(xCorr))
        sprintf('Fs: %.1f MHz', fs/1e6)
        sprintf('CFO (MHz): %s', ternary(isfinite(Fc), sprintf('%.6f MHz', f_off/1e6), 'N/A'))
        sprintf('Güç ort.: %.2f dB', 10*log10(mean(abs(xCorr).^2)))
        sprintf('SNR tahmini: %.2f dB', snr_est)
    };
    text(0.1,0.7,lines,'FontSize',10);
    title('İstatistikler');

    % OFDM: Discarding CP, Applying FFT
    S = ofdm_fft(xCorr, Nfft, cpLen, Nsym, useShift);  % boyut: (Nfft × Nsym)

    % H[k] estimate
    Npil = numel(pilotInd);
    knownPil = ones(Npil, Nsym);
    if p.altPilots
        knownPil(:,1:2:end) = repmat(p.pilot1, 1, numel(1:2:Nsym));
        if Nsym >= 2
            knownPil(:,2:2:end) = repmat(p.pilot2, 1, numel(2:2:Nsym));
        end
    else
        knownPil = repmat(p.pilot1, 1, Nsym);
    end

    Yp  = S(pilotInd, :);       % (Npil × Nsym)
    Hp  = Yp ./ knownPil;       % kaba kanal (CPE içerir)
    CPE = angle(mean(Hp, 1));   % sembol başına ortak faz
    Hp_c = Hp .* exp(-1j*CPE);  % CPE çıkar
    S_c  = S  .* exp(-1j*CPE);  % tüm taşıyıcılara uygula

    % Frekans içi interpolasyon (her sembol)
    xq   = (1:Nfft).';
    H_est = zeros(Nfft, Nsym);
    for s = 1:Nsym
        H_est(:,s) = interp1(pilotInd, Hp_c(:,s), xq, 'linear', 'extrap');
    end

    % Gürültü varyansı (aktif bant dışından, kaba)
    excl = true(Nfft,1);
    excl(min([dataInd; pilotInd]) : max([dataInd; pilotInd])) = false; % aktif bandı hariç tut
    excl(pilotInd) = true;  % pilotlar hariç kalsın
    if any(~excl,'all')
        noiseVar = median(abs(S_c(excl,:)).^2,'all');
    else
        noiseVar = 0;
    end

    % Eşitleme (MMSE)
    He2 = abs(H_est).^2;
    W   = conj(H_est) ./ (He2 + noiseVar);
    Xeq   = S_c .* W;
    Xd_eq = Xeq(dataInd, :);

    % Constellation
    subplot(3,3,8);
    Xd_flat = Xd_eq(:);
    NshowF  = min(numel(Xd_flat), p.constPointsToShow);
    idxF    = round(linspace(1, numel(Xd_flat), NshowF));
    scatter(real(Xd_flat(idxF)), imag(Xd_flat(idxF)), 3, 'filled'); grid on; axis equal;
    title('Constellation'); xlabel('I'); ylabel('Q');

    sgtitle(sprintf('Kayıtlı Sinyal Analizi - %s', datestr(timeStamp)), 'FontSize', 14);
end

function eps_hat = estimate_cfo_cp(x, Nfft, cpLen, NsymUse)
    symLen = Nfft + cpLen;
    eps_list = zeros(NsymUse,1);
    for k = 1:NsymUse
        i0 = (k-1)*symLen + 1;
        i1 = i0 + symLen - 1;
        if i1 > numel(x), eps_list = eps_list(1:k-1); break; end
        sym_td = x(i0:i1);
        pre  = sym_td(1:cpLen);
        tail = sym_td(Nfft+1:end);
        M = sum(conj(pre).*tail);
        eps_list(k) = angle(M)/(2*pi);
    end
    eps_med = median(eps_list);
    dev  = abs(eps_list - eps_med);
    mask = dev <= 3*median(dev+eps);
    eps_hat = mean(eps_list(mask));
end

function snr = estimateSNR(data)
    signal_power = mean(abs(data).^2);
    noise_floor  = median(abs(data).^2);
    snr = 10*log10(signal_power/noise_floor);
end

function S = ofdm_fft(x, Nfft, cpLen, NsymUse, doShift)
    symLen = Nfft + cpLen;
    S = complex(zeros(Nfft, NsymUse));
    for k = 1:NsymUse
        i0 = (k-1)*symLen + 1;
        i1 = i0 + symLen - 1;
        if i1 > numel(x), S = S(:,1:k-1); break; end
        sym_td = x(i0 + cpLen : i1); % CP'yi at
        Xf = fft(sym_td, Nfft);
        if doShift, Xf = fftshift(Xf); end
        S(:,k) = Xf;
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
