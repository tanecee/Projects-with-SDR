function main_RX()
    % Bu p yapısı, sadece Nsym başlığını okumak için kullanılır.
    p_base = parametersOFDM(1);

    % USRP alıcı nesnesini bir kere oluştur
    rx = comm.SDRuReceiver( ...
        'Platform','B210', ...
        'SerialNum','33A6D57', ...
        'MasterClockRate',p_base.sample_rate, ...
        'CenterFrequency',2.9e9, ...
        'Gain',75, ...
        'DecimationFactor',1, ...
        'SamplesPerFrame', 20000, ...
        'OutputDataType','single');
    
    sync = zadoffChuSeq(8,255);
    disp('--- Dinleme başlıyor ---');
    
    % OFDM parametrelerini ekrana bir kere yazdır
    fprintf('OFDM Parametreleri:\n');
    fprintf('Nfft: %d\n', p_base.Nfft);
    fprintf('Nsym: %d\n', p_base.Nsym);
    fprintf('Data Alt-taşıyıcı Sayısı: %d\n', p_base.dataScs);
    fprintf('Pilot Alt-taşıyıcı Sayısı: %d\n', p_base.pilotScs);
    fprintf('CP Uzunluğu: %d\n', p_base.cpLength);
    fprintf('Örnek Hızı: %.2f Msps\n', p_base.sample_rate / 1e6);
    disp('-------------------------');

    % Log dosyasını oluştur ve başlığı yaz
    % Dosya adını o anki tarih ve saate göre otomatik oluşturur.
    log_filename = sprintf('OFDM_Log_%s.txt', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
    fileID = fopen(log_filename, 'a');
    if fileID == -1
        error('Hata: Log dosyası oluşturulamadı.');
    end
    
    % Log dosyasına başlangıç bilgilerini yaz
    fprintf(fileID, '--- OFDM Alıcı Log Dosyası ---\n');
    fprintf(fileID, 'Oluşturulma Zamanı: %s\n', datestr(now));
    fprintf(fileID, '---------------------------------\n\n');
    
    % Grafik penceresini bir kere oluştur
    h_fig = figure;
    set(h_fig, 'Name', 'Gerçek Zamanlı OFDM Analizi');
    
    % Alt-grafikleri oluştur
    h_constellation = subplot(2, 1, 1);
    title(h_constellation, 'Alınan Takımyıldız (Constellation)');
    xlabel(h_constellation, 'Gerçek Kısım');
    ylabel(h_constellation, 'Sanal Kısım');
    grid(h_constellation, 'on');
    
    h_spectrum = subplot(2, 1, 2);
    title(h_spectrum, 'Alınan Sinyal Spektrumu');
    xlabel(h_spectrum, 'Frekans Bini');
    ylabel(h_spectrum, 'Genlik (dB)');
    grid(h_spectrum, 'on');
    
    while true
        [rxSig, len] = rx();
        if len > 0
            corr = abs(conv(rxSig, conj(flipud(sync))));
            [peak, idx] = max(corr);
            
            if peak > 1   % Senkron eşiği
                % Komut penceresinde senkronizasyon bilgisini göstermeye gerek yok
                %fprintf('Senkron bulundu! peak=%.2f, idx=%d\n', peak, idx);
                
                header_length = (p_base.Nfft + p_base.cpLength) * 1; 
                start_idx_header = idx + 1;
                
                if start_idx_header + header_length - 1 <= length(rxSig)
                    rxFrame_header = rxSig(start_idx_header:start_idx_header+header_length-1);
                    
                    [rx_bits_header, ~] = ofdmRx(rxFrame_header, p_base);
                    
                    if length(rx_bits_header) >= 16
                        nsym_bits = rx_bits_header(1:16);
                        actual_Nsym = bi2de(nsym_bits', 'left-msb');
                        
                        p_full = parametersOFDM(actual_Nsym);
                        
                        if start_idx_header + p_full.wformLength - 1 <= length(rxSig)
                            rxFrame_full = rxSig(start_idx_header:start_idx_header+p_full.wformLength-1);
                            
                            [rx_bits_full, rx_data_symbols] = ofdmRx(rxFrame_full, p_full);
                            
                            rx_message_bits = rx_bits_full(17:end);
                            rx_message = bitsToText(double(rx_message_bits));
                            
                            % Mesajı hem komut penceresine hem de dosyaya yaz
                            fprintf('(%s) Alınan Mesaj (Nsym: %d): %s\n', datestr(now, 'HH:MM:SS.FFF'), actual_Nsym, rx_message);
                            fprintf(fileID, '[%s] Nsym: %d, Mesaj: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), actual_Nsym, rx_message);
                            
                            % Grafiklerin güncellenmesi
                            if ishandle(h_fig)
                                % Konstelasyon grafiğini güncelle
                                cla(h_constellation);
                                scatter(h_constellation, real(rx_data_symbols), imag(rx_data_symbols), 'b.');
                                title(h_constellation, 'Alınan Takımyıldız (Constellation)');
                                
                                % Spektrum grafiğini güncelle
                                cla(h_spectrum);
                                plot(h_spectrum, 20*log10(abs(fftshift(fft(rxFrame_full)))), 'b');
                                title(h_spectrum, 'Alınan Sinyal Spektrumu');
                                
                                drawnow;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Script durduğunda log dosyasını kapat
    fclose(fileID);
end
