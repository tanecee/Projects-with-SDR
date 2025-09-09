function analyze_ofdm_signal(filename)
% ANALYZE_OFDM_SIGNAL - Analyzes recorded OFDM signal
%
% Parameters:
%   filename: Path to the recorded signal file

% Load recorded data
data = load(filename);
signal_data = data.complete_frame;
sample_rate = data.sample_rate;
sync_pos = data.sync_pos;
sync_sequence = data.sync_sequence;
Nfft = data.Nfft;
Nsym = data.Nsym;
cpLength = data.cpLength;
wformLength = data.wformLength;

fprintf('Analyzing signal from %s\n', filename);
fprintf('Sample rate: %.1f MHz\n', sample_rate/1e6);
fprintf('Signal length: %d samples\n', length(signal_data));
fprintf('Duration: %.3f seconds\n', length(signal_data)/sample_rate);
fprintf('Sync position: %d\n', sync_pos);
fprintf('OFDM parameters: Nfft=%d, Nsym=%d, cpLength=%d\n', Nfft, Nsym, cpLength);

% Separate sync and OFDM data
sync_part = signal_data(1:length(sync_sequence));
ofdm_part = signal_data(length(sync_sequence)+1:end);

% Create analysis plots
figure('Position', [100, 100, 1400, 1000]);

% PSD plot using pwelch
subplot(3, 2, 1);
nfft_psd = 1024;
window = hamming(nfft_psd);
noverlap = nfft_psd/2;
[pxx, f] = pwelch(signal_data, window, noverlap, nfft_psd, sample_rate, 'centered');
plot(f/1e6, 10*log10(pxx));
title('Power Spectral Density (PSD)');
xlabel('Frequency (MHz)');
ylabel('Power Spectral Density (dB/Hz)');
grid on;

% Time domain plot
subplot(3, 2, 2);
time_axis = (0:length(signal_data)-1) / sample_rate;
plot(time_axis, real(signal_data), 'b', 'DisplayName', 'Real');
hold on;
plot(time_axis, imag(signal_data), 'r', 'DisplayName', 'Imag');
title('Time Domain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('show');
grid on;

% Constellation plot of OFDM part
subplot(3, 2, 3);
scatter(real(ofdm_part), imag(ofdm_part), 1, 'filled', 'MarkerAlpha', 0.6);
title('OFDM Data Constellation Diagram');
xlabel('In-phase');
ylabel('Quadrature');
axis equal;
grid on;

% Constellation plot of sync part
subplot(3, 2, 4);
scatter(real(sync_part), imag(sync_part), 1, 'filled', 'MarkerAlpha', 0.6);
title('Sync Sequence Constellation Diagram');
xlabel('In-phase');
ylabel('Quadrature');
axis equal;
grid on;

% Histogram of amplitudes
subplot(3, 2, 5);
amplitudes = abs(signal_data);
histogram(amplitudes, 50, 'FaceAlpha', 0.7);
title('Amplitude Distribution');
xlabel('Amplitude');
ylabel('Count');
grid on;

% Phase histogram
subplot(3, 2, 6);
phases = angle(signal_data);
histogram(phases, 50, 'FaceAlpha', 0.7);
title('Phase Distribution');
xlabel('Phase (radians)');
ylabel('Count');
grid on;

% Save plot
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
saveas(gcf, sprintf('record_data/analysis_%s.png', timestamp));

% Additional analysis: OFDM symbols
analyze_ofdm_symbols(ofdm_part, Nfft, cpLength, Nsym, sample_rate);
end

function analyze_ofdm_symbols(ofdm_signal, Nfft, cpLength, Nsym, sample_rate)
% Analyze individual OFDM symbols

symbol_length = Nfft + cpLength;
num_symbols = floor(length(ofdm_signal) / symbol_length);

fprintf('\nOFDM signal analysis:\n');
fprintf('  Expected symbols: %d\n', Nsym);
fprintf('  Detected symbols: %d\n', num_symbols);
fprintf('  Symbol length: %d samples\n', symbol_length);
fprintf('  CP length: %d samples\n', cpLength);

if num_symbols == 0
    fprintf('No complete OFDM symbols found\n');
    return;
end

% Analyze first few symbols
for i = 1:min(3, num_symbols)
    start_idx = (i-1) * symbol_length + 1;
    end_idx = start_idx + symbol_length - 1;
    symbol = ofdm_signal(start_idx:end_idx);
    
    % Remove CP
    symbol_no_cp = symbol(cpLength+1:end);
    
    % Apply FFT
    symbol_freq = fft(symbol_no_cp);
    symbol_freq = fftshift(symbol_freq);
    
    % Plot symbol in frequency domain
    figure('Position', [100, 100, 1200, 500]);
    
    subplot(1, 2, 1);
    plot(abs(symbol_freq));
    title(sprintf('OFDM Symbol %d - Frequency Domain', i));
    xlabel('Subcarrier Index');
    ylabel('Magnitude');
    grid on;
    
    subplot(1, 2, 2);
    scatter(real(symbol_freq), imag(symbol_freq), 'filled', 'MarkerAlpha', 0.6);
    title(sprintf('OFDM Symbol %d - Constellation', i));
    xlabel('In-phase');
    ylabel('Quadrature');
    axis equal;
    grid on;
    
    % Save plot
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    saveas(gcf, sprintf('record_data/symbol_%d_%s.png', i, timestamp));
end
end

% Function to analyze a specific recorded file
function analyze_specific_file()
% Example: Analyze a specific recorded file
filename = 'record_data/ofdm_signal_20231201_123045.mat';
analyze_ofdm_signal(filename);
end
