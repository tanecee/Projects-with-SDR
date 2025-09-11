function [rx_file, info] = pluto_record_signal()

outDir = 'recorded_data'; 
Fs = 20e6;
Fc = 3.3e9;
gain = 70;
durSec = 1.0;
zc_u = 8;
Lsync = 255;
timeout_s = 60;
pwr_threshold = 1e-3; 
peak_to_avg_ratio = 5; %Sahte pikleri yok eder.
k_med = 12;
baseThr = 0.15;

 
    Nneed   = round(Fs * durSec); %Yazılacak toplam örnek sayısı
    rx_file = "";
    info    = struct('found',false,'peak',NaN,'avg',NaN,'thr',NaN,'pwr',NaN,'loc',NaN);

    if ~exist(outDir, 'dir'); mkdir(outDir); end
    sync = zadoffChuSeq(zc_u, Lsync);

    mf   = conj(flipud(sync));

    rx = comm.SDRRxPluto(...
        'RadioID',              'sn:1044734c9605001104000100d5934f698c', ...
        'CenterFrequency',      Fc, ...
        'BasebandSampleRate',   Fs, ...
        'SamplesPerFrame',      1024*100, ...
        'OutputDataType',       'double', ...
        'GainSource',           'Manual');
    rx.ShowAdvancedProperties = true;
    rx.Gain = gain;

    fprintf('=== Pluto RX (%.1f MHz @ %.1f MHz), gain=%d dB ===\n', Fs/1e6, Fc/1e6, gain);
    fprintf('Senkron bekleniyor (timeout=%ds)...\n', timeout_s);

    % Önceki frame'den kalan örnekleri saklamak için 'prev' değişkeni.
    prev  = complex(zeros(Lsync-1,1)); %Frame sınırında kaçabilecek piki yakalamak için, bir önceki frame'den Lsync-1 örneği saklanıyor
    found = false;
    wrote = 0;
    t0    = tic;
    last_msg_time = tic; % Mesaj yazdırma için zamanlayıcı

    try
        while true
            % Zaman aşımı kontrolü
            if toc(t0) > timeout_s
                fprintf('!!! Zaman aşımı: Senkron bulunamadı. Kayıt yapılmadı.\n');
                break;
            end
            
            % Sinyal aranırken periyodik mesaj yazdırma
            if toc(last_msg_time) > 2 % Her 2 saniyede bir mesaj yazdır
                fprintf('Senkron sinyali aranıyor...\n');
                last_msg_time = tic;
            end

            % PlutoSDR'dan bir frame sinyal okuma
            x = rx();

            % DC bastırma
            %x = x - mean(x);

            % Güç kapısı (Power Gate) kontrolü: Sinyal gücü eşiğin altındaysa, işlem yapma.
            pwr = mean(abs(x).^2);
            if pwr < pwr_threshold
                prev = x(end-(Lsync-2):end);
                continue;
            end

            % Önceki frame ile birleştirme (senkronizasyonun frame sınırında kaçmaması için)
            xcat = [prev; x];

            % Eşleşmiş filtre (korelasyon)
            c    = conv(xcat, mf, 'valid');
            magc = abs(c);
            
            % Medyan ve ortalama hesaplama
            med  = median(magc);
            avg  = mean(magc);

            % Uyarlanabilir eşik
            thr  = max(k_med*med, baseThr);

            % Zirve bulma
            [pks, locs] = findpeaks(magc, 'MinPeakHeight', thr, ...
                                        'MinPeakDistance', round(0.5*Lsync), ...
                                        'MinPeakProminence', 0.1);

            if ~found && ~isempty(pks)
                % Ek güvenlik kontrolü: Tepe/Ortalama oranı kontroü
                if pks(1)/avg < peak_to_avg_ratio
                    % Zayıf veya sahte pik; bir sonraki frame'e geç
                    prev = x(end-(Lsync-2):end);
                    continue;
                end

                % Senkron sinyalinin başlangıç indeksi
                k = locs(1);
                start_in_xcat = k + Lsync - 1;
                start_in_x    = start_in_xcat - length(prev);
                start_in_x    = max(start_in_x, 1);

                % İlk kaydedilecek sinyal parçasını hazırla
                seg = x(start_in_x:end);

                % Dosyayı sadece senkron bulunduğunda AÇ
                ts      = datestr(datetime('now'), 'yyyymmdd_HHMMSS');
                rx_file = fullfile(outDir, sprintf('pluto_rx_%s_%.0fMHz.dot', ts, Fs/1e6));
                fid     = fopen(rx_file, 'wb');
                assert(fid > 0, 'Dosya açılamadı: %s', rx_file);

                % İlk parçayı diske yaz
                fwrite(fid, [real(seg(:)) imag(seg(:))].', 'float32');
                wrote = wrote + numel(seg);

                % Bilgilendirme
                info.found = true; info.peak = pks(1); info.avg = avg; info.thr = thr; info.pwr = pwr; info.loc = k;
                fprintf('Senkron bulundu (peak=%.2f, avg=%.2f, p/avg=%.2f, pwr=%.2e). Kayda geçildi.\n', ...
                        pks(1), avg, pks(1)/avg, pwr);

                % Kalan 1 saniyelik veriyi tamamla
                while wrote < Nneed
                    xi = rx();
                    xi = xi - mean(xi); % DC bastırma
                    nwrite = min(length(xi), Nneed - wrote);
                    if nwrite > 0
                        fwrite(fid, [real(xi(1:nwrite)) imag(xi(1:nwrite))].', 'float32');
                        wrote = wrote + nwrite;
                    end
                end

                fclose(fid);
                fprintf('Kayıt tamam: %s (%.2f s, %d örnek)\n', rx_file, durSec, Nneed);
                found = true;
                break; % İşlem bitti, döngüden çık
            end

            % Sonraki frame için kuyruk
            prev = x(end-(Lsync-2):end);
        end

    catch ME
        % Hata durumunda açık dosyayı kapat ve kısmi dosyayı sil
        if exist('fid','var') && fid > 0
            fclose(fid);
            if strlength(rx_file) > 0 && isfile(rx_file)
                delete(rx_file);
            end
        end
        release(rx);
        rethrow(ME);
    end

    release(rx);
    
    % Senkron bulunamadıysa, dosya yolu boş döner
    if ~found
        rx_file = '';
    end
end
