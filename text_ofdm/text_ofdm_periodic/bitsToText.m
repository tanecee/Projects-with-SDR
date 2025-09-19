function message = bitsToText(bits)
    message = '';
    
    for i = 1:8:length(bits)
        if i+7 > length(bits)
            break;
        end
        
        byte_bits = bits(i:i+7);
        
        char_code = 0;
        for j = 1:8
            char_code = char_code + byte_bits(j) * 2^(8-j);
        end
        
        if char_code == 0
            break;
        end
        
        message = [message, char(char_code)];
    end
end
