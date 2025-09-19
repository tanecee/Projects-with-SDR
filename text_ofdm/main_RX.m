function main_RX()

    p_base = parametersOFDM(1); % Nsym=1 ile temel parametreler oluşturulur

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
    
    while true
        [rxSig, len] = rx();
        if len > 0
            corr = abs(conv(rxSig, conj(flipud(sync))));
            [peak, idx] = max(corr);
            
            if peak > 1   % Senkron eşiği
                fprintf('Senkron bulundu! peak=%.2f, idx=%d\n', peak, idx);
                
                header_length = (p_base.Nfft + p_base.cpLength) * 1; % Sadece 1 sembol alalım
                start_idx_header = idx + 1;
                
                if start_idx_header + header_length - 1 <= length(rxSig)
                    rxFrame_header = rxSig(start_idx_header:start_idx_header+header_length-1);
                    
                    [rx_bits_header, ~] = ofdmRx(rxFrame_header, p_base);
                    
                    % İlk 16 biti Nsym bilgisi olarak alınır
                    % Burada bir kontrol eklemek iyi olur, eğer rx_bits_header 16 bitten kısaysa hata almamak için
                    if length(rx_bits_header) >= 16
                        nsym_bits = rx_bits_header(1:16);
                        actual_Nsym = bi2de(nsym_bits', 'left-msb');
                        fprintf('Alınan Nsym degeri: %d\n', actual_Nsym);
                    else
                        disp('Nsym başlığı okunamadı. Yeterli bit yok.');
                        continue; % Döngüyü tekrarla
                    end
                    
                    % Yeni p yapısını actual_Nsym'e göre oluşturuldu
                    p_full = parametersOFDM(actual_Nsym);
                    
                    if start_idx_header + p_full.wformLength - 1 <= length(rxSig)
                        rxFrame_full = rxSig(start_idx_header:start_idx_header+p_full.wformLength-1);
                        
                        [rx_bits_full, ~] = ofdmRx(rxFrame_full, p_full);
                        
                        % Nsym bilgisini içeren ilk 16 bit'i atla
                        rx_message_bits = rx_bits_full(17:end);
                        
                        rx_message = bitsToText(double(rx_message_bits));
                        fprintf('Alınan Mesaj: %s\n', rx_message);
                        
                        save('rxFrame.mat','rxFrame_full','rx_bits_full');
                        disp('Frame decode edildi ve kaydedildi.');
                        break;
                    else
                        disp('Frame tamamı arabellekte değil, sonraki frame bekleniyor.');
                    end
                else
                    disp('Başlık arabellekte değil, sonraki frame bekleniyor.');
                end
            end
        end
    end
end
