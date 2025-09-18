function main_RX()
    p = parametersOFDM();

    rx = comm.SDRuReceiver( ...
        'Platform','B210', ...
        'SerialNum','33A6D57', ...
        'MasterClockRate',p.sample_rate, ...
        'CenterFrequency',3.1e9, ...
        'Gain',55, ...
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

            if peak > 1   % eşik, deneme ile ayarlanmalı
                fprintf('Senkron bulundu! peak=%.2f, idx=%d\n', peak, idx);

                % Senkron sonrası frame uzunluğunu çıkar
                start_idx = idx + 1;
                if start_idx + p.wformLength - 1 <= length(rxSig)
                    rxFrame = rxSig(start_idx:start_idx+p.wformLength-1);

                    % OFDM çözümü
                    [~, rx_bits, rx_symbols] = ofdmRx(rxFrame);

                    % Kaydet
                    save('rxFrame.mat','rxFrame','rx_bits','rx_symbols');
                    disp('Frame decode edildi ve kaydedildi.');
                    
                    % Tek seferlik çalıştırmak istersen break koy
                    % break;
                else
                    disp('Frame tamamı bufferda değil, sonraki frame bekleniyor.');
                end
            end
        end
    end
end
