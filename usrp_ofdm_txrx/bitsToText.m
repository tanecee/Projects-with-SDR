function message = bitsToText(bits)
    message = '';
    
    % 8 bit gruplarına ayır
    for i = 1:8:length(bits)
        if i+7 > length(bits)
            break;
        end
        
        % 8 biti al
        byte_bits = bits(i:i+7);
        
        % ASCII karaktere dönüştür
        char_code = 0;
        for j = 1:8
            char_code = char_code + byte_bits(j) * 2^(8-j);
        end
        
        % Null terminator bulunca dur
        if char_code == 0
            break;
        end
        
        % Karakteri string'e ekle
        message = [message, char(char_code)];
    end
end
