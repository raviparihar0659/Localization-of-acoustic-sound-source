% This MATLAB code simulates the Time Difference of Arrival (TDOA) quantization error as the position of a sound source varies. 
% The code calculates and compares the original TDOA values with quantized TDOA values for a range of source distances.



function tdoa_quantization_error_varying_source_position
    % Constants
    speed_of_sound = 343; % Speed of sound in m/s
    mic1_position = [0, 0, 0]; % Microphone 1 coordinates in m
    mic2_position_y = 1; % Position of mic2 on y-axis in m (example fixed position)

    % Sampling frequency
    fs = 64000*8; % Sampling frequency in Hz
    Ts = 1 / fs; % Sampling period in seconds

    % The code simulates source distances ranging from 5 meters to 200 meters with a step size of 1 meter.
    % Source position range (in m)
    source_min = 5; % Minimum source distance in m (5 meters)
    source_max = 200; % Maximum source distance in m (200 meters)
    source_step = 1; % Step size in m
    source_distances = source_min:source_step:source_max; % Vector of source distances

    % Preallocate result arrays
    original_tdoas = zeros(size(source_distances));
    quantized_tdoas = zeros(size(source_distances));

    % Fixed position for microphone 2
    mic2_position = [0, mic2_position_y, 0];

    for i = 1:length(source_distances)
        % Calculate source position
        source_position = [source_distances(i), 0, 0];

        % Calculate distances from source to microphones
        d1 = norm(source_position - mic1_position);
        d2 = norm(source_position - mic2_position);

        % Calculate times to reach microphones
        t1 = d1 / speed_of_sound;
        t2 = d2 / speed_of_sound;

        % Original TDOA
        tau = t2 - t1;
        original_tdoas(i) = tau;

        % Quantized TDOA
        tau_m = round(tau / Ts) * Ts;
        quantized_tdoas(i) = tau_m;
    end

    % Plot results
    figure;
    plot(source_distances , original_tdoas * 1e6, 'b', 'LineWidth', 1.5); hold on;
    plot(source_distances , quantized_tdoas * 1e6, 'r--', 'LineWidth', 1.5);
    xlabel('Source Distance (m)');
    ylabel('TDOA (\mus)');
    legend('Original TDOA', 'Quantized TDOA');
    title('original TDOA and estimated TDOA due Quantization Error vs Source Distance ');
    grid on;
end

tdoa_quantization_error_varying_source_position;


%%

% This MATLAB code examines the impact of quantization error on Time Difference of Arrival (TDOA) measurements by varying both
%  the source distance and the distance between two microphones. 
% The code calculates and compares the original and quantized TDOA values for a range of source distances and microphone spacings.

% TDOA Calculation:
% 
% For each combination of source distance and microphone spacing, the code calculates the distance from the source to both microphones (d1 and d2).
% It computes the original TDOA (tau) based on the difference in arrival times of the sound at each microphone.
% The TDOA is then quantized using the sampling period (Ts).
function tdoa_quantization_error_varying_source_and_distance_from_mic
    % Constants
    speed_of_sound = 34300; % Speed of sound in cm/s
    mic1_position = [0, 0, 0]; % Microphone 1 coordinates in cm
    
    % Sampling frequency
    fs = 64000*8; % Sampling frequency in Hz
    Ts = 1 / fs; % Sampling period in seconds
    
    % Source position range (in cm)
    source_min = 500; % Minimum source distance in cm (5 meters)
    source_max = 20000; % Maximum source distance in cm (200 meters)
    source_step = 10; % Step size in cm
    source_distances = source_min:source_step:source_max; % Vector of source distances
    
    % Distance range for mic2 (in cm)
    d_min = 20; % Minimum distance in cm
    d_max = 100; % Maximum distance in cm
    d_step = 10; % Step size in cm
    mic2_distances = d_min:d_step:d_max; % Vector of distances from 20 cm to 150 cm
    
    % Preallocate result matrices
    original_tdoas = zeros(length(mic2_distances), length(source_distances));
    quantized_tdoas = zeros(length(mic2_distances), length(source_distances));
    
    for j = 1:length(mic2_distances)
        % Calculate position of mic2
        mic2_position = [0, mic2_distances(j), 0];
        
        for i = 1:length(source_distances)
            % Calculate source position
            source_position = [source_distances(i), 0, 0];
            
            % Calculate distances from source to microphones
            d1 = norm(source_position - mic1_position);
            d2 = norm(source_position - mic2_position);
            
            % Calculate times to reach microphones
            t1 = d1 / speed_of_sound;
            t2 = d2 / speed_of_sound;
            
            % Original TDOA
            tau = t2 - t1;
            original_tdoas(j, i) = tau;
            
            % Quantized TDOA
            tau_m = round(tau / Ts) * Ts;
            quantized_tdoas(j, i) = tau_m;
        end
    end
    
    % Plot results
    figure;
    subplot(2, 1, 1);
    for j = 1:length(mic2_distances)
        plot(source_distances / 100, original_tdoas(j, :) * 1e6, 'LineWidth', 1.5); hold on;
    end
    xlabel('Source Distance (m)');
    ylabel('Original TDOA (\mus)');
    title('Original TDOA vs Source Distance for different values of distance btw mics');
    legend(arrayfun(@(x) sprintf('d = %d cm', x), mic2_distances, 'UniformOutput', false));
    grid on;
    
    subplot(2, 1, 2);
    for j = 1:length(mic2_distances)
        plot(source_distances / 100, quantized_tdoas(j, :) * 1e6, 'LineWidth', 1.5); hold on;
    end
    xlabel('Source Distance (m)');
    ylabel('Quantized TDOA (\mus)');
    title('Quantized TDOA vs Source Distance for different values of distance btw mics');
    legend(arrayfun(@(x) sprintf('d = %d cm', x), mic2_distances, 'UniformOutput', false));
    grid on;
end
tdoa_quantization_error_varying_source_and_distance_from_mic;
