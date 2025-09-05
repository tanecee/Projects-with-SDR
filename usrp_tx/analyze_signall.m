function analyze_signall()
    % --- Dosya seç ---
    [filename, pathname] = uigetfile({'*.dot;*.dat','IQ Data Files'}, 'Analiz için dosya seçin');
    if isequal(filename,0), fprintf('Dosya seçilmedi.\n'); return; end
    filepath = fullfile(pathname, filename);
    fprintf('Analiz ediliyor: %s\n', filepath);

    % --- Örnekleme hızı ---
    answer = inputdlg({'Örnekleme Hızı (Hz):'}, 'Sample Rate', [1 40], {'20e6'});
    Fs = 20e6; if ~isempty(answer), Fs = str2double(answer{1}); end

    % --- IQ oku (float32 interleaved) ---
    sig = read_dot_file(filepath);

    % --- Parametreler (TX kodunla uyumlu) ---
    Nfft  = 256;
    cpLen = Nfft/4;             % 64
    Nsym  = 14;
    Lsync = 255;  zc_u = 8;

    activeBins = (Nfft/4+1):(3*Nfft/4); % 65:192 (128 aktif)
    dataBins   = activeBins(1:2:end);   % tekler: data

    % ---- ZC ile payload başlangıcını bul (gösterme, sadece hizalama) ----
    zc = zadoffChuSeq(zc_u, Lsync);
    h  = conj(flipud(zc));
    c  = conv(sig, h, 'valid');
    [~, loc] = max(abs(c));
    k_payload = min(loc + Lsync, numel(sig)); % ZC biter bitmez
    Lsym = Nfft + cpLen;

    % ---- Zaman bölgesi (ilk 4000 örnek) ----
    figure('Position',[80 80 1200 800],'Name','OFDM TX Analysis','NumberTitle','off');

    subplot(2,2,1);
    Ntime = min(4000, numel(sig));
    t_ms  = (0:Ntime-1)/Fs*1e3;
    plot(t_ms, real(sig(1:Ntime)), 'LineWidth',1); hold on;
    plot(t_ms, imag(sig(1:Ntime)), 'LineWidth',1);
    xlabel('Time (ms)'); ylabel('Amplitude');
    title('Time Domain (first 4000 samples)'); grid on; legend('I','Q');

    % ---- PSD (Welch) two-sided, eksen −10..+10 MHz ----
    subplot(2,2,2);
    [Pxx, f] = pwelch(sig, hann(1024), 512, 1024, Fs, 'centered'); % power/Hz
    plot(f/1e6, 10*log10(Pxx+1e-20), 'LineWidth',1.3);
    xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)'); title('Power Spectrum (Welch)'); grid on;
    xlim([-Fs/2 Fs/2]/1e6);   % Fs=20e6 için tam −10..+10 MHz

    % ---- CP'yi at, FFT al, sadece data SC'lerini çiz (eşitlemesiz) ----
    subplot(2,2,[3 4]);  % altta geniş tek eksen
    constel = [];

    % payload içinde mümkün olduğu kadar sembol işle
    remain = numel(sig) - k_payload + 1;
    nSymAvail = max(0, floor(remain / Lsym));
    nUse = min(Nsym, nSymAvail);  % dosya kısa ise kısmi

    for m = 1:nUse
        s0  = k_payload + (m-1)*Lsym + cpLen;           % CP'den sonraki ilk örnek
        td  = sig(s0 : s0+Nfft-1);                      % Nfft örnek
        fd  = fftshift(fft(td, Nfft));                  % frekans alanı
        constel = [constel; fd(dataBins).'];            % data taşıyıcıları
    end

    % Eğer sembol bulunamadıysa (ör. hizalama tutmadı) fallback: baştan bir pencere dene
    if isempty(constel) && numel(sig) >= (cpLen+Nfft)
        td  = sig(1+cpLen : cpLen+Nfft);
        fd  = fftshift(fft(td, Nfft));
        constel = fd(dataBins).';
    end

    scatter(real(constel), imag(constel), 8, 'filled'); axis equal; grid on;
    xlabel('I'); ylabel('Q'); title('Constellation (Data SCs) [CP removed]');
end

% ---------- yardımcı ----------
function sig = read_dot_file(fname)
    fid = fopen(fname,'rb'); assert(fid>0,'Dosya açılamadı.');
    raw = fread(fid,'float32'); fclose(fid);
    assert(mod(numel(raw),2)==0,'IQ sayısı çift değil.');
    I = raw(1:2:end); Q = raw(2:2:end);
    sig = double(I) + 1j*double(Q);
end
