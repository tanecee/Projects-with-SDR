function record_ofdm_signal(record_duration, center_freq, sample_rate)
% RECORD_OFDM_SIGNAL - Records OFDM signal with Zadoff-Chu synchronization
%
% Parameters:
%   record_duration: Recording duration in seconds (default: 1.0)
%   center_freq: Center frequency in Hz (default: 3.3e9)
%   sample_rate: Sample rate in Hz (default: 20e6)

% Default parameters
if nargin < 1, record_duration = 1.0; end
if nargin < 2, center_freq = 3.3e9; end
if nargin < 3, sample_rate = 20e6; end

% OFDM parameters (from TX code)
Nfft = 256;
Nsym = 14;
cpLength = Nfft / 4;
wformLength = (Nfft + cpLength) * Nsym;
sync_seq_length = 255;
zadoff_chu_root = 8;

% Generate Zadoff-Chu sequence for synchronization
sync_sequence = generate_zadoff_chu_sequence(zadoff_chu_root, sync_seq_length);

% Setup SDR
rx = sdrrx('Pluto', ...
    'RadioID', 'ip:pluto.local', ...
    'CenterFrequency', center_freq, ...
    'BasebandSampleRate', sample_rate, ...
    'SamplesPerFrame', 2^15, ...
    'GainSource', 'Manual', ...
    'Gain', 30, ...
    'OutputDataType', 'double');

% Create record_data directory if it doesn't exist
if ~exist('record_data', 'dir')
    mkdir('record_data');
end

% Generate filename with timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('record_data/ofdm_signal_%s.mat', timestamp);

fprintf('Recording OFDM signal for %.1f seconds...\n', record_duration);

% Record signal
num_samples = round(record_duration * sample_rate);
rx_signal = zeros(num_samples, 1);
num_frames = ceil(num_samples / rx.SamplesPerFrame);

for i = 1:num_frames
    samples_to_read = min(rx.SamplesPerFrame, num_samples - (i-1)*rx.SamplesPerFrame);
    rx_signal((i-1)*rx.SamplesPerFrame + 1 : (i-1)*rx.SamplesPerFrame + samples_to_read) = rx();
end

% Find synchronization position
[sync_pos, corr_magnitude] = find_sync_position(rx_signal, sync_sequence);

fprintf('Synchronization found at position: %d\n', sync_pos);

% Extract complete frame (sync + OFDM data)
frame_start = sync_pos;
frame_end = frame_start + sync_seq_length + wformLength;

if frame_end > length(rx_signal)
    warning('Not enough samples for a complete frame');
    frame_end = length(rx_signal);
end

complete_frame = rx_signal(frame_start:frame_end);

% Save signal and metadata
save(filename, 'complete_frame', 'sample_rate', 'center_freq', ...
    'sync_pos', 'timestamp', 'corr_magnitude', 'sync_sequence', ...
    'Nfft', 'Nsym', 'cpLength', 'wformLength');

fprintf('Signal saved to %s\n', filename);

% Plot synchronization results
plot_synchronization(rx_signal, sync_pos, sync_seq_length, wformLength, corr_magnitude);

release(rx);
end

function sync_sequence = generate_zadoff_chu_sequence(u, N)
% Generate Zadoff-Chu sequence for synchronization
n = 0:N-1;
phase = -pi * u * n .* (n + 1) / N;
sync_sequence = exp(1j * phase).';
end

function [sync_pos, corr_magnitude] = find_sync_position(rx_signal, sync_sequence)
% Find synchronization position using cross-correlation
corr = xcorr(rx_signal, sync_sequence);
corr = corr(length(sync_sequence):end); % Keep valid part
corr_magnitude = abs(corr);

% Find the peak
[~, sync_pos] = max(corr_magnitude);
sync_pos = sync_pos(1); % In case of multiple peaks

% Calculate threshold
threshold = 5.0 * median(corr_magnitude);

if corr_magnitude(sync_pos) < threshold
    fprintf('Warning: Weak synchronization peak. Peak: %.2f, Threshold: %.2f\n', ...
        corr_magnitude(sync_pos), threshold);
end
end

function plot_synchronization(rx_signal, sync_pos, sync_seq_length, wform_length, corr_magnitude)
% Plot synchronization results
figure('Position', [100, 100, 1200, 800]);

% Plot received signal
subplot(2, 1, 1);
plot(abs(rx_signal));
hold on;
xline(sync_pos, 'r--', 'Sync Position', 'LineWidth', 1.5);
xline(sync_pos + sync_seq_length, 'g--', 'Sync End', 'LineWidth', 1.5);
xline(sync_pos + sync_seq_length + wform_length, 'm--', 'OFDM End', 'LineWidth', 1.5);
title('Received Signal with Sync Position');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Signal', 'Location', 'best');
grid on;

% Plot correlation
subplot(2, 1, 2);
plot(corr_magnitude);
hold on;
xline(sync_pos, 'r--', 'Peak Position', 'LineWidth', 1.5);
title('Cross-Correlation with Zadoff-Chu Sequence');
xlabel('Sample Index');
ylabel('Correlation Magnitude');
legend('Correlation', 'Location', 'best');
grid on;

% Save plot
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
saveas(gcf, sprintf('record_data/sync_plot_%s.png', timestamp));
end
