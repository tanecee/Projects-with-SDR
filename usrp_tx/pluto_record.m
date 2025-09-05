
function rx_file = pluto_record()

Fs = 20e6;            
Fc = 3.1e9;        
Lsync = 255;
zc_u = 8;
recordDuration = 1.0;
Nneed = round(Fs * recordDuration);
timeout = 30;
pwr_threshold = 1e-3;
peak_to_avg_ratio = 5;


sync = zadoffChuSeq(zc_u, Lsync);
mf = conj(flipud(sync));

rx = comm.SDRRxPluto( ...
    'RadioID', 'sn:1044734c9605001104000100d5934f698c', ...
    'CenterFrequency', Fc, ...
    'BasebandSampleRate', Fs, ...
    'SamplesPerFrame', 1024*100, ...
    'OutputDataType', 'double', ...
    'GainSource', 'Manual');
rx.ShowAdvancedProperties = true;
rx.Gain = 40;


outDir = 'recorded_data';
if ~exist(outDir, 'dir'); mkdir(outDir); end
ts = datestr(datetime('now'), 'yyyymmdd_HHMMSS');
rx_file = fullfile(outDir, ['pluto_rx_' ts '.dot']);
[fid, msg] = fopen(rx_file, 'wb');
assert(fid > 0, msg);

fprintf('Senkron bekleniyor...\n');
buf = complex([], []);
prev = complex(zeros(Lsync-1, 1));
found = false;
start_time = tic;

try
    while true
        if toc(start_time) > timeout
            error('Timeout: Synchronization sequence not found within %d seconds.', timeout);
        end
        x = rx();
        pwr = mean(abs(x).^2);
        if pwr < pwr_threshold
            fprintf('Low signal power (%.2e), waiting for signal...\n', pwr);
            prev = x(end-(Lsync-2):end);
            continue;
        end
        xcat = [prev; x];
        c = conv(xcat, mf, 'valid');
        magc = abs(c);
        med = median(magc);
        avg = mean(magc);
        thr = max(15*med, 0.2); % Stricter threshold
        [pks, locs] = findpeaks(magc, 'MinPeakHeight', thr, 'MinPeakDistance', round(0.5*Lsync), 'MinPeakProminence', 0.1);

        if ~found && ~isempty(pks)
            if pks(1)/avg < peak_to_avg_ratio
                fprintf('Peak detected but too weak (peak=%.2f, avg=%.2f, ratio=%.2f), continuing...\n', ...
                    pks(1), avg, pks(1)/avg);
                prev = x(end-(Lsync-2):end);
                continue;
            end
            k = locs(1);
            start_in_xcat = k + Lsync - 1;
            start_in_x = start_in_xcat - length(prev);
            start_in_x = max(start_in_x, 1);
            seg = x(start_in_x:end);
            fwrite(fid, [real(seg(:)) imag(seg(:))].', 'float32');
            buf = seg;
            wrote = length(seg);
            fprintf('Senkron bulundu (peak=%.2f, power=%.2e, peak/avg=%.2f). Kayda geçildi.\n', ...
                pks(1), pwr, pks(1)/avg);
            found = true;

            while wrote < Nneed
                xi = rx();
                nwrite = min(length(xi), Nneed - wrote);
                if nwrite > 0
                    fwrite(fid, [real(xi(1:nwrite)) imag(xi(1:nwrite))].', 'float32');
                    buf = [buf; xi(1:nwrite)];
                end
                wrote = wrote + nwrite;
            end
            break;
        end
        prev = x(end-(Lsync-2):end);
    end
catch ME
    fclose(fid);
    release(rx);
    rethrow(ME);
end
fclose(fid);
release(rx);
fprintf('Kayıt tamam: %s (1.0 s)\n', rx_file);

end
