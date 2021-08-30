%% This file runs a frequency tagging analysis on EEG data files defined in a JSON file.
% The pipeline was created by Jessica Figueira, then I adapted it to this dataset.
% Benjamin Velie, 7/29/2021
% veliebm@ufl.edu


%% Load our parameters.
parameters_file = 'processed/freqtageeg/parameters.json';
do_all(parameters_file)
delete_lock_file(mfilename('fullpath'))


%% Driver functions.
function do_all(parameters_file)
    % Read a JSON file to get parameters, then use them process all subjects.
    all_parameters = read_json(parameters_file);
    for i = 1:numel(all_parameters)
        parameters = all_parameters(i);
        do_one(parameters.in_eeg_name, parameters.in_eeg_dir, parameters.out_fft_path, parameters.out_hilbert_path, parameters.out_sliding_window_prefix)
    end
end
function do_one(in_eeg_name, in_eeg_dir, out_fft_path, out_hilbert_path, out_sliding_window_prefix)
    % Process a subject.
    %% Load our data.
    EEG = load_dataset(in_eeg_name, in_eeg_dir);
    
    % Convert array to a double rather than single because freqtag_HILB requires that.
    dataset = double(EEG.data);
    
    % Oz sensor number
    oz_id = 20;
    
    num_time_points = EEG.pnts;

    
    %% Define all frequencies contained in the data, see section 3.1 for a guide on how to delimit the frequency axis.
    epoch_duration = EEG.xmax - EEG.xmin;
    frequency_resolution = 1/epoch_duration;
    sampling_rate = EEG.srate;
    num_frequencies_all = sampling_rate / 2;
    faxisall = 0:frequency_resolution:num_frequencies_all;

    
    %% Eliminate unnecessary frequencies 
    faxis = faxisall(12:2:265);
    num_freqs = length(faxis);
    
    
    %% Plot the data, averaging the trials 
    %plot_average_of_trials(dataset)
    %plot_sensor(dataset, oz_id)
    

    %% Run FFT on the data
    baseline = 4/5*500;
    stimulus_start = baseline + 500;
    stimulus_end = 2483;
    
    [pow, phase, freqs] = freqtag_FFT(mean(dataset(:, stimulus_start:stimulus_end, :),3), 500); %The fft function takes a 2-D matrix, it is necessary to average over the trials (third dimension) 
    to_tsv(pow, out_fft_path)
    to_tsv(freqs, 'processed/freqtageeg/freqs.tsv')
    
    
    %% Plot the spectrum. For visualization purposes, only plot the frequencies of interest.
    %plot_fft_all_sensors(pow, freqs, num_freqs)
    %plot_fft_sensor(pow, freqs, num_freqs, oz_id)
    
    
    %% Run FFT on single-trials.
    [spec] = freqtag3D_FFT(dataset, stimulus_start:stimulus_end, sampling_rate); %No need to average over the third dimension here.
    
    %% Run sliding window analysis.
    [trialpow,winmat3d,phasestabmat,trialSNR] = freqtag_slidewin(dataset, 0, stimulus_start:stimulus_end, stimulus_start:stimulus_end, 12, 600, 500, out_sliding_window_prefix);


    %% Averaging across sliding windows.
    meanwinmat = mean(winmat3d, 3);
    %plot_sliding_window_average(meanwinmat, 'Mean of moving windows at 12 Hz', 'Sample points', 'Voltage')
    

    %% Project sliding window average into frequency domain.
    % Note the sample rate is 600 Hz at 6 Hz
    [meanwinmat_pow, meanwinmat_phase, meanwinmat_freqs] = freqtag_FFT(meanwinmat, 600); 
    %plot_sliding_window_average_FFT(meanwinmat_freqs, meanwinmat_pow, 'Power spectrum of the mean window shifted at 12 Hz', 'Frequency (Hz)')


    %% Plot spectrum
    %plot_single_fft_all_sensors(faxisall, spec)
    %plot_single_fft_one_sensor(spec, faxisall, oz_id)
    
    
    %% Run Hilbert Transform. Check the 12Hz frequency.
    [powermat, phasemat, complexmat] = freqtag_HILB(mean(dataset(:, stimulus_start:stimulus_end, :), 3), 12, 8, oz_id, 0, sampling_rate);
    to_tsv(powermat, out_hilbert_path)
    
    
    %% Plot frequency over time
    %plot_hilbert(powermat, oz_id, 'Hilbert Transform (12Hz Stimuli)')
    
    
    %% Close all figures and clean the worskpace
    %close all
    %clear all
    
    
end


%% Plot functions.
function plot_sliding_window_average(meanwinmat, fig_title, fig_x_label, fig_y_label)
    % Plot that fancy sliding window average you calculated.
    figure
    plot(meanwinmat')
    title(fig_title)
    xlabel(fig_x_label)
    ylabel(fig_y_label)
end
function plot_sliding_window_average_FFT(meanwinmat_freqs, meanwinmat_pow, fig_title, fig_x_label)
    % Plot that fancy sliding window average you calculated - but now after it's been FFT'd!
    figure,
    bar(meanwinmat_freqs(2:20), meanwinmat_pow(32,2:20)');
    title(fig_title)
    xlabel(fig_x_label)
end


function plot_hilbert(powermat, sensor, fig_title)
    % Plot frequency over time
    figure(), plot(powermat(sensor, :)')
    ax = gca;
    ax.FontSize = 18;
    ax.Box = 'off';
    xlabel('Time (ms)'), ylabel('Frequency')
    title([fig_title])
    legend ('Oz sensor'), legend boxoff
end

function plot_single_fft_one_sensor(spec, faxisall, sensor)
    figure(), plot(faxisall(:, 1:100), spec(sensor, 1:100))
    ax = gca;
    ax.FontSize = 18;
    ax.Box = 'off';
    xlabel('Frequency (Hz)'), ylabel('Power (?V�)')
    title([ 'Average Frequency Spectrum'])
    legend ('Oz sensor'), legend boxoff
end
function plot_single_fft_all_sensors(faxisall, spec)
    % Plot our single FFT for every sensor.
    figure(), plot(faxisall(:, 1:100), spec(:, 1:100))
    ax = gca;%%editing the plot
    ax.FontSize = 18;
    ax.Box = 'off'
    xlabel('Frequency (Hz)'), ylabel('Power (?V�)')
    title([ 'Average Frequency Spectrum' ])
end

function plot_fft_sensor(pow, freqs, num_freqs, sensor)
    % Plot the spectrum of a sensor of interest. Only plot the frequencies of interest
    figure(), plot(freqs(1:num_freqs), pow(sensor, 1:num_freqs)) %select Oz sensor
    ax = gca;
    ax.FontSize = 18;
    ax.Box = 'off'
    xlabel('Frequency (Hz)'), ylabel('Power (?V�)')
    title([ 'Frequency Spectrum' ])
    legend ('Oz sensor'), legend boxoff
end
function plot_fft_all_sensors(pow, freqs, num_freqs)
    % Plot the FFT of all sensors
    figure(), plot(freqs(1:num_freqs), pow(:, 1:num_freqs))
    ax = gca;%%editing the plot
    ax.FontSize = 18;
    ax.Box = 'off'
    xlabel('Frequency (Hz)'), ylabel('Power (?V�)')
    title(['Frequency Spectrum'])
end

function plot_sensor(dataset, sensor)
    % Plot the data of a specific sensor.
    figure(), plot(mean(dataset(sensor,:,:), 3)')
    ax = gca;
    ax.FontSize = 18;
    ax.Box = 'off';
    xlabel('Time (ms)'), ylabel('Amplitude (?V)');
    title([ 'Time domain plot']);
    legend ('Oz sensor'), legend boxoff;
end
function plot_average_of_trials(dataset)
    % Plot the data, averaging all trials together.
    figure(), plot(mean(dataset,3)') 
    ax = gca;%%editing the plot
    ax.FontSize = 18;
    ax.Box = 'off';
    xlabel('Time (ms)'), ylabel('Amplitude (?V)');
    title([ 'Time domain plot'])
end


%% Other functions.
function to_tsv(mat, out_filename)
    % Write a matrix to tsv.
    dlmwrite(out_filename, mat, '\t');
end
function [EEG] = load_dataset(file_name, directory)
    % Load a dataset.
    eeglab;
    EEG = pop_loadset('filename', file_name, 'filepath', directory);
    EEG = eeg_checkset(EEG);
end
function [data] = read_json(in_path)
    % Read a JSON file.
    fname = in_path; 
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    data = jsondecode(str);
end
