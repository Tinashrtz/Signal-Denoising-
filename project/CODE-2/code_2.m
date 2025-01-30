clc
clear all;

%% Data Acquisition and Initial Listening

[x, Fs] = audioread('projectDSP.wav');
[samples,channels] = size(x);

info = audioinfo('projectDSP.wav');
Nbits = info.BitsPerSample; % Bit depth

sound(x, Fs);
%%%%%%%%% uncomment later if needed
t = (0:length(x)-1) / Fs; % Time vector
plot(t, x);
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform of the Audio Signal');
hold;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Spectrum Analysis
N = length(x); % Length of the signal
X = fft(x); % Compute the DFT
magnitude = abs(X) / N; % Magnitude spectrum
f = (0:N-1) * (Fs / N); % Frequency vector for one-sided spectrum

magnitude_one_sided = magnitude(1:N/2+1);
f_one_sided = f(1:N/2+1);

%%%%%%%%%%%%%uncomment later
figure;
plot(f_one_sided, magnitude_one_sided);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum of the Audio Signal');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Highest Magnitudes that is the sound of the siren noise
[peaks, locs] = findpeaks(abs(X), f, 'MinPeakHeight', max(abs(X))*0.1);
peaks = peaks(1:2);
locs = locs(1:2);
disp(['first peak: ', num2str(locs(1)), 'Hz, second peak: ', num2str(locs(2)), 'Hz']);

% Plot the spectrum and highlight the peaks
%%%%%%%%%%%% uncomment later
% hold on;
% plot(f_one_sided, magnitude_one_sided, locs, peaks, 'r*');
% plot(locs, peaks, 'r*'); % Mark peaks
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('Magnitude Spectrum with Identified Peaks');
% legend('Magnitude Spectrum', 'Identified Peaks');
% grid on;
% hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filter Design and Frecuency Resoponse Accuracy

% Parameters
Fs = 8000; % Sampling rate
N = 102; % Filter order (adjust as needed)

% Normalize frequencies (relative to Nyquist frequency)
f1 = [2050 2150] / (Fs/2); % Stopband around 2100 Hz
f2 = [3250 3350] / (Fs/2); % Stopband around 3300 Hz

% Combine stopbands into one specification
bands = [0 f1(1) f1(2) f2(1) f2(2) 1];
desired = [1 1 0 0 1 1]; % Desired response (1 = pass, 0 = stop)

% Design with different windows
hamming_filter = fir2(N, bands, desired, hamming(N+1));
blackman_filter = fir2(N, bands, desired, blackman(N+1));
kaiser_filter = fir2(N, bands, desired, kaiser(N+1, 3));

%%%%%%%%%%uncomment later
% Frequency response accuracy plots
fvtool(hamming_filter, 1, blackman_filter, 1, kaiser_filter, 1, ...
    'Fs', Fs);
legend('Hamming', 'Blackman', 'Kaiser');


[H1, W1] = freqz(hamming_filter, 1, 1024, Fs); % Frequency response
[H2, W2] = freqz(blackman_filter, 1, 1024, Fs); % Frequency response
[H3, W3] = freqz(kaiser_filter, 1, 1024, Fs); % Frequency response
% Phase Response (in radians)
phase_response_hamming = angle(H1);
phase_response_blackman = angle(H2);
phase_response_kaiser = angle(H3);


figure;
subplot(3,1,1)
plot(W1, phase_response_hamming, 'r');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase Response of hamming');
grid on;
subplot(3,1,2)
plot(W2, phase_response_blackman, 'g');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase Response of blackman');
grid on;
subplot(3,1,3)
plot(W3,phase_response_kaiser,'b');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase Response of kaiser');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% applying filters
filtered_signal_hamming = filter(hamming_filter, 1, x);
filtered_signal_blackman = filter(blackman_filter, 1, x);
filtered_signal_kaiser = filter(kaiser_filter, 1, x);


%%%%%%%% uncomment later
t = (0:length(x)-1) / Fs;
subplot(4, 1, 1);
plot(t, x);
title('Original Signal');
subplot(4, 1, 2);
plot(t, filtered_signal_hamming);
title('Filtered with Hamming Window');
subplot(4, 1, 3);
plot(t, filtered_signal_blackman);
title('Filtered with Blackman Window');
subplot(4, 1, 4);
plot(t, filtered_signal_kaiser);
title('Filtered with Kaiser Window');
xlabel('Time (s)');
%%%%%%%%%%%%%%%%%%%%%%%%%




% Magnitude spectrum of filtered signals
figure;
[original_fft, f] = compute_fft(x, Fs);
[hamming_fft, ~] = compute_fft(filtered_signal_hamming, Fs);
[blackman_fft, ~] = compute_fft(filtered_signal_blackman, Fs);
[kaiser_fft, ~] = compute_fft(filtered_signal_kaiser, Fs);


%%%%%%%%%% uncommect later
subplot(1,3,1)
plot(f, original_fft, 'k',f, hamming_fft, 'r');
legend('Original', 'Hamming');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain Comparison');
subplot(1,3,2)
plot(f, original_fft, 'k',f, blackman_fft, 'g');
legend('Original', 'blackman');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain Comparison');
subplot(1,3,3)
plot(f, original_fft, 'k',f, kaiser_fft, 'b');
legend('Original', 'kaiser');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain Comparison');
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% kaiser is our choice %%%


%% Signal Processing
% Filtered signals after using kaiser - Equiripple

% Equiripple multiband filter (Order: 100)
% freq = [0 1940 2040 2160 2260 3150 3250 3350 3450 4000], mag = [10 10 0 0 10 10 0 0 10 10]
load('equiripple_multiband.mat');
signal_multiband = filter(multibandO100, filtered_signal_kaiser);
dft_multiband = fft(signal_multiband);
% sound(signal_multiband, Fs);

% High-order Equiripple multiband filter (Order: 350) - Narrow stopband (Not practical)
% freq = [0 1800 1900 2300 2400 3000 3100 3500 3600 4000], mag = [-20 -20 0 0 -20 -20 0 0 -20 -20]
load('equiripple_high_order.mat');
signal_multiband_high = filter(multibandO350, filtered_signal_kaiser);
dft_multiband_high = fft(signal_multiband_high);
% sound(signal_multiband_high, Fs);

% Plot the magnitude spectrum of the filtered signals
%%%%%%%%%%% uncomment later if needed
subplot(2, 1, 1);
plot(f(1:floor(N/2)+1), abs(dft_multiband(1:floor(N/2)+1)));
title('Equiripple (Order: 100)');
xlabel('Frequency (Hz)');
ylabel('|DFT|');

subplot(2, 1, 2);
plot(f(1:floor(N/2)+1), abs(dft_multiband_high(1:floor(N/2)+1)));
title('Equiripple High (Order: 350)');
xlabel('Frequency (Hz)');
ylabel('|DFT|');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Spectrum Analysis Results:', newline, ...
    '1. Equiripple Filter (Multiband: Order 100)', newline, ...
    '2. High-Order Equiripple Filter (Order 350)']);

% comparing the two filtered audio signal
audiowrite('filtered_sound_2.wav', signal_multiband_high / max(abs(signal_multiband_high)), Fs);
audiowrite('filtered_sound_1.wav', signal_multiband / max(abs(signal_multiband)), Fs);

%%%%% filtered_sound_2.wav removed the siren noise ...
%%%%% but filtered_sound_1.wav couldnt do it ...
%%%%% so using multiband_high after kaiser is the best result


%% Function Parts

% Helper Function
function [magnitude, f] = compute_fft(signal, Fs)
    N = length(signal);
    X = fft(signal);
    magnitude = abs(X(1:N/2+1)) / N;
    f = (0:N/2) * Fs / N;
end

