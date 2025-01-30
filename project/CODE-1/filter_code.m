clc;
clear;
close all;

%% Part 1: Read the audio file
[audio_signal, sampling_rate] = audioread("projectDSP.wav");  % Read the audio file
% Display audio properties
fprintf('Sampling Rate: %d Hz\n', sampling_rate);
fprintf('Number of Channels: %d\n', size(audio_signal, 2));
fprintf('Audio Duration: %.2f seconds\n', length(audio_signal)/sampling_rate);

% Play the original audio file
disp('Playing the original audio file...');
%sound(audio_signal, sampling_rate);
%pause(length(audio_signal)/sampling_rate); % Wait for playback to complete

%% Part 2: Frequency analysis
% Length of the audio signal
num_samples = length(audio_signal);

% Apply FFT to the first channel of the audio signal (if stereo)
frequency_data = fft(audio_signal(:, 1), num_samples);

% Calculate magnitude spectrum and frequency axis
magnitude_spectrum = abs(frequency_data(1:num_samples/2));  % Only positive frequencies
frequency_axis = (0:num_samples/2-1) * (sampling_rate / num_samples);  % Frequency axis

% Plot the magnitude spectrum
figure;
plot(frequency_axis, magnitude_spectrum);
title('Magnitude Spectrum of the Audio Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Detect prominent peaks in the frequency spectrum
[peaks, locations] = findpeaks(magnitude_spectrum, frequency_axis, 'MinPeakHeight', max(magnitude_spectrum)/10);
fprintf('Detected Noise Frequencies:\n');
disp(locations);

%% Part 3: Low-pass filtering and fading
% Sampling frequency
fs = sampling_rate;  % Define the sampling frequency

% Low-pass filter parameters
filter_order = 300;  % Filter order
cutoff_frequency = 1300;  % Low-pass cutoff frequency (Hz)
normalized_cutoff = cutoff_frequency / (fs / 2);  % Normalize cutoff frequency

% Design the low-pass filter using Kaiser window
low_pass_filter = fir1(filter_order, normalized_cutoff, 'low', kaiser(filter_order + 1, 5));

% Apply the low-pass filter to the audio signal
filtered_signal = filter(low_pass_filter, 1, audio_signal);

% Generate a linear fading signal
signal_length = length(filtered_signal);  % Length of the filtered signal
fade_duration = 1.5;  % Fade duration in seconds
fade_samples = round(fade_duration * fs);  % Number of samples for fading

% Create the fade-in signal
fading_signal = ones(1, signal_length);  % Initialize with ones
fading_signal(1:fade_samples) = linspace(0, 1, fade_samples);  % Linear fade-in

% Apply fading to the filtered signal
filtered_signal = filtered_signal .* fading_signal(:);

% Compute the magnitude spectrum of the filtered signal
num_samples_filtered = length(filtered_signal);
fft_filtered = fft(filtered_signal, num_samples_filtered);
magnitude_filtered = abs(fft_filtered(1:num_samples_filtered/2));
frequency_axis_filtered = (0:num_samples_filtered/2-1) * (fs / num_samples_filtered);

% Cap the magnitude for better visualization
magnitude_filtered = min(magnitude_filtered, 5);

% Plot the magnitude spectrum of the filtered signal
figure;
plot(frequency_axis_filtered, magnitude_filtered, 'b');
title('Magnitude Spectrum After Low-pass Filtering and Fading');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Plot the fading signal
time_vector = (0:length(fading_signal)-1) / fs;  % Time vector in seconds
figure;
plot(time_vector, fading_signal, 'r', 'LineWidth', 1.5);
title('Linear Fading Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;
xlim([0 max(time_vector)]);

% Play and save the processed audio
disp('Playing the filtered and faded audio...');
sound(filtered_signal, fs);
audiowrite('filtered_sound_1.wav', filtered_signal, fs);
