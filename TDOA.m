

%% trying to find the TDOA on gong audio dataset of matlab (Fs = 8192 in gong dataset) 

clear;
load gong;
refsig = y;
delay1 = 5;
delay2 = 25;

delay_sig1 = delayseq(refsig,delay1);
delay_sig2 = delayseq(refsig,delay2);
sig1= awgn(delay_sig1, -5,"measured");
sig2= awgn(delay_sig2, -5,"measured");
refsig= awgn(refsig, -5,"measured");
tau_est = gccphat([sig1,sig2],refsig);
disp(tau_est);

%%
Y = 20*log10(abs(fft(y)));
L = length(y);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

% Plot spectrum
figure;

subplot(4,1,1);
plot(y);
xlabel('samples')
ylabel('amplitude')
title('original signal');
subplot(4,1,2);
plot(delay_sig1);
xlabel('samples')
ylabel('amplitude')
title('delayed signal1 signal');
subplot(4,1,3);
plot(delay_sig2);
xlabel('samples')
ylabel('amplitude')
title('delayed signal1 signal');
subplot(4,1,4);
plot(sig1);
xlabel('samples')
ylabel('amplitude')
title('noisy delayed signal');

figure;
plot(f,P1)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum ');

%% finding the TDOA of two signals

% noise generation of distribution in frequency domain within 1000Hz 

% Parameters
fs = 512000;  
N = 512000;  
f_high = 1000;

% Frequency resolution
df = fs / N;

% Create frequency bins
freq_bins = (0:N-1)' * df;

% Create a frequency domain signal with random phases and amplitudes
% only within the desired frequency range
freq_noise = zeros(N, 1);
valid_bins = (freq_bins <= f_high) & (freq_bins >= 0); % Frequency range 0 to 1000 Hz

% Random amplitude and phase for the valid frequency bins
amplitude = randn(sum(valid_bins), 1);
phase = rand(sum(valid_bins), 1) * 2 * pi;
freq_noise(valid_bins) = amplitude .* exp(1i * phase);

% Ensure the signal is real by mirroring the frequencies
freq_noise = [freq_noise; conj(freq_noise(end-1:-1:2))];

% Convert to time domain using IFFT
time_noise = ifft(freq_noise, 'symmetric');

% Take only the first N samples
time_noise = time_noise(1:N);


L = length(time_noise);
Y = fft(time_noise);
P2 = abs(Y / L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2 * P1(2:end-1);
f = fs * (0:(L/2)) / L;
P1_dB = 20 * log10(P1);

figure;
plot(f, P1_dB);
title('Single-Sided Amplitude Spectrum of Band-Limited Noise in dB');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;


%%  noise generation of distribution in frequency domain within 1000Hz 
% Given signal y
y = randn(512000, 1); % Example signal, replace with your actual signal

% Parameters
fs = 512000; % Sampling frequency in Hz
fmax_noise = 1000; % Maximum frequency of noise in Hz
SNR_dB = 10; % Desired SNR in dB

% Generate noise in frequency domain within the desired frequency range
N = length(y);
f = (0:N-1)' * (fs/N); % Frequency vector

% Create a frequency mask for the noise within 1000 Hz
noise_mask = (f <= fmax_noise) ;

% Generate random noise in frequency domain
noise_freq_domain = (randn(N, 1) + 1i*randn(N, 1)) .* noise_mask;

% Convert noise back to time domain
noise = ifft(noise_freq_domain, 'symmetric');

% Normalize noise power to achieve desired SNR
signal_power = mean(y.^2);
noise_power = mean(noise.^2);
K = sqrt(signal_power / noise_power / (10^(SNR_dB / 10)));
noise = noise * K;

% Add noise to the original signal
y_noisy = y + noise;

% Plot frequency distribution of the noise in dB
noise_freq_domain_magnitude = abs(noise_freq_domain);
noise_freq_domain_dB = 20 * log10(noise_freq_domain_magnitude + eps); % Add eps to avoid log(0)

figure;
subplot(2, 1, 1);
plot(f, noise_freq_domain_dB);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Distribution of Noise (in dB)');
xlim([0, fs/2]); % Plot up to Nyquist frequency

% Plot the original and noisy signals in time domain
subplot(2, 1, 2);
plot(y, 'b');
hold on;
plot(y_noisy, 'r');
xlabel('Sample');
ylabel('Amplitude');
legend('Original Signal', 'Noisy Signal');
title('Original and Noisy Signals in Time Domain');

% Optional: Verify SNR
actual_SNR_dB = 10 * log10(mean(y.^2) / mean(noise.^2));
disp(['Actual SNR: ', num2str(actual_SNR_dB), ' dB']);



%% noise adding in signal and finding the TDOA (trying on drone recorded signal)
%
clear;
% reading the microphone recording file
file= {"E:/IITD/sem6/Project TDOA/main/240412_Drone recordings from Kapil/Channel_0.tdms"};
datashell = cell(numel(file), 1);

datashell{1}  = tdmsread(file{1}, 'ChannelGroupName', 'Untitled', 'ChannelNames', 'Filtered');


%%

start_sample = 20 * 64000+1; 
end_sample = 25* 64000;   
             
Raw = datashell{1}{1}.Filtered(start_sample:end_sample);
upsample_factor = 8;

% Upsample received signals
received_signal1_upsampled = resample(Raw, upsample_factor, 1);


% Input signal
fs = 64000*8;
dt = 1/fs; 
t_end = 1;
t = 0 : dt : t_end-dt;
length_t = length(t);

% Setting of delay
% del_time =0.001;

% del_time =0.000014576895;

del_time = 0.0000291516 ;  % for testing of estimated TDOA ( set original equal to del_time)
%the value of delay between signals is the TDOA when two microphones are placed 1m far and
%a source signal is at 50 m from the microphones (microphones and source are in perpendicular fashion, not linear). 
x0_num = round(del_time/dt);               



% Setting the reference signal
rs1=received_signal1_upsampled;


y = rs1(1:length_t);

%x = rs1(end-512000-x0_num+1:end-x0_num);
x = [zeros(x0_num,1); y(1:end-x0_num)]; % Delay the signal by td

%% passing signal and refrance signal to low pass filter( if we don't pass signal to LPF, it don't affect result because signal has frequency spactrum within 1200Hz) and ploting the spectrum of y

% Sampling frequency (Hz)
fs = 64000*8; % Example sampling frequency, you should replace this with your actual sampling frequency

% Define the cutoff frequency
Fc = 1000; % Cutoff frequency of 2000 Hz

% Design a low-pass filter using designfilt
lpFilt = designfilt('lowpassfir', 'PassbandFrequency', Fc, 'StopbandFrequency', Fc+500, ...
                    'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', fs);

% Apply the filter to the signal
y_filtered = filter(lpFilt, y);

x_filtered = filter(lpFilt, x);


L = length(y_filtered);
Y = fft(y_filtered);

% Compute the two-sided spectrum P2
P2 = abs(Y / L);

% Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define the frequency domain f
f = fs*(0:(L/2))/L;

% Plot the single-sided amplitude spectrum
figure;
plot(f, P1)
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
grid on;

%%
X = fft(x_filtered);
% Compute the two-sided spectrum P2
P2 = abs(X / L);

% Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(f,P1);
title('Single-Sided Amplitude Spectrum of x(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
grid on;

%%
function noisy_signal = add_noise(signal, snr_db,fs)

    signal_power = rms(signal)^2;
    snr_linear = 10^(snr_db / 10); % Correct conversion
    noise_power = signal_power / snr_linear;
    N = length(signal);  
    

    fmax_noise = 1000; % Maximum frequency of noise in Hz
    f = (0:N-1)' * (fs/N); % Frequency vector
    
    % Create a frequency mask for the noise within 1000 Hz
    noise_mask = (f <= fmax_noise) ;
    
    %disp(noise_mask);
    %fprintf('dimension of noise_mask: %f x %f \n',size(noise_mask));
    
    % Generate random noise in frequency domain
    noise_freq_domain = (randn(N, 1) + 1i*randn(N, 1)) .* noise_mask;
    
    % Convert noise back to time domain
    time_noise = ifft(noise_freq_domain, 'symmetric');
    

    % Plot frequency distribution of the noise in dB
    noise_freq_domain_magnitude = abs(noise_freq_domain);
    noise_freq_domain_dB = 20 * log10(noise_freq_domain_magnitude + eps); % Add eps to avoid log(0)
    
    figure;
    plot(f, noise_freq_domain_dB);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Frequency Distribution of Noise (in dB)');
    xlim([0, fs/2]); % Plot up to Nyquist frequency

    
    generated_noise_power = rms(time_noise)^2;

    % Scale the filtered noise to achieve the desired noise power
    noise = time_noise * sqrt(noise_power / generated_noise_power);
    noisy_signal = signal + noise;

    % Calculate the actual SNR to verify the result
    actual_snr_linear = rms(signal)^2 / rms(noise)^2;
    actual_snr_db = 10 * log10(actual_snr_linear);
    
    % Display the actual SNR
    fprintf('Desired SNR (dB): %f, Actual SNR (dB): %f\n', snr_db, actual_snr_db);

    


end


%%
% Define the desired SNR in dB
desired_snr_db = -10; 
% Add noise to the signals
noisy_x = add_noise(x_filtered, desired_snr_db,fs);
noisy_y = add_noise(y_filtered, desired_snr_db,fs);



% Plot the signals for visual verification
figure;
subplot(2,1,1);
plot(t,x_filtered);
hold on;

plot(t,noisy_x);
title('Signal 1 and Noisy Signal 1');
legend('Original Signal 1', 'Noisy Signal 1');

subplot(2,1,2);
plot(t,y_filtered);
hold on;
plot(t,noisy_y);
title('Signal 2 and Noisy Signal 2');
legend('Original Signal 2', 'Noisy Signal 2');

%
figure
subplot(4,1,1);  
plot(t,x_filtered);
xlabel('time(s)')
ylabel('amplitude')
title('Delayed signal without noise')

subplot(4,1,2);  
plot(t,noisy_x);
xlabel('time(s)')
ylabel('amplitude')
title('Delayed signal')

subplot(4,1,3);  
plot(t,y_filtered);
xlabel('time(s)')
ylabel('amplitude')
title('Reference signal without noise')

subplot(4,1,4);  
plot(t,noisy_y);
xlabel('time(s)')
ylabel('amplitude')
title('Reference signal')

% finding the TDOA between noisy signal x and noisy signal y. 
% The signals noisy_x and noisy_y are multiplied by a Hanning window (hann) to taper the signals before further processing.
% The Hanning window reduces spectral leakage by smoothly tapering the signal to zero at the boundaries, which is especially useful in frequency domain analyses.

sig = noisy_x.*hann(length(noisy_x));    
refsig = noisy_y.*hann(length(noisy_y)); 

% GCC-PHAT Time Delay Estimation:
% 
% [tau_est, R, lags] = gccphat(sig, refsig, fs);
% The gccphat function is used to estimate the time delay (tau_est) between the two signals sig and refsig.
% R is the cross-correlation function, and lags represents the time lags corresponding to the cross-correlation values.
% The function returns tau_est, the estimated time delay between the signals, which indicates the difference in arrival times of a signal at two sensors or microphones.

[tau_est,R,lags] = gccphat(sig,refsig,fs);  


figure;
plot(lags,10*log10(real(R(:,1))))

xlabel('Lag Times (ms)')
ylabel('Cross-correlation in dB scale')

disp(tau_est);
