
filePaths = {
    "E:/IITD/sem6/Project TDOA/main/240412_Drone recordings from Kapil/Channel_0.tdms", ...
    "E:/IITD/sem6/Project TDOA/main/240412_Drone recordings from Kapil/Channel_1.tdms", ...
    "E:/IITD/sem6/Project TDOA/main/240412_Drone recordings from Kapil/Channel_2.tdms", ...
    "E:/IITD/sem6/Project TDOA/main/240412_Drone recordings from Kapil/Channel_3.tdms", ...
    % Add more file paths as needed
};

% Preallocate cell array to store data from each file
dataCell = cell(numel(filePaths), 1);

% Read data from each file
for i = 1:numel(filePaths)
    dataCell{i} = tdmsread(filePaths{i}, 'ChannelGroupName', 'Untitled', 'ChannelNames', 'Filtered');
end
%%
for i = 1:numel(filePaths)    
    ttData = dataCell{i};
    start_sample = 1 * 64000; % 1 is the starting time in seconds
    end_sample = 60 * 64000;   % 60 is the end time in seconds
    fs = 64000;                % Sampling frequency

    RawData = ttData{1}.Filtered(start_sample:end_sample);
    time_segment = (start_sample:end_sample) / fs;

    % Plot original signal
    figure;
    
    subplot(4, 1, 1);
    plot(time_segment, RawData);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Original Signal from File ' num2str(i)]);
    
    newarr=RawData(20*64000+1:40*64000);
    f=64000;
    [pxx,f2] = pwelch(newarr,64000,32000,f,fs);
    subplot(4, 1, 2);
    plot(f2,10*log10(pxx));
    title('pwelch plot');

    % Compute spectrogram
    window = 64000; % Specify the window length
    overlap = window/2; % Specify the overlap between adjacent segments
    nfft = 64000; % Specify the number of FFT points
    [S,F,T] = spectrogram(RawData, window, overlap, nfft, fs);

    % Plot spectrogram
    subplot(4, 1, 3);
    imagesc(T, F, 10*log10(abs(S)));
    axis xy;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title(['Spectrogram from File ' num2str(i)]);
    colorbar;

    % Extract line frequencies (tonals)
    [f0, timeInSec] = pitch(RawData, fs);
    
    % plot(20*log10(abs(fft(RawData(20*64000+1:40*64000)))));
    % Compute FFT
    Y = fft(RawData);
    L = length(RawData);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;

    % Plot spectrum
    figure;
    plot(f,P1)
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title(['Spectrum from File ' num2str(i)]);

    
  

end 
