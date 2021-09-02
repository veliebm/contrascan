%% This file runs a frequency tagging analysis on EEG data files defined in a JSON file.
% The pipeline was created by Jessica Figueira and Andreas Keil, then I adapted it to this dataset.
% Benjamin Velie, 7/29/2021
% veliebm@ufl.edu


%% Load our parameters.
parameters_file = 'processed/freqtageeg/parameters.json';
do_all(parameters_file)
delete_lock_file(mfilename('fullpath'))


%% Driver functions.
function do_all(parameters_file)
    % Read a JSON file to get parameters, then use them to process all subjects.
    eeglab;
    all_parameters = read_json(parameters_file);
    for i = 1:numel(all_parameters)
        parameters = all_parameters(i);
        do_one(parameters.in_eeg_name, parameters.in_eeg_dir, parameters.out_fft_path, parameters.out_hilbert_path, parameters.out_sliding_window_prefix, str2num(parameters.frequency), parameters.sliding_window_average_plot, parameters.sliding_window_average_fft_plot, parameters.out_faxisall_path, parameters.out_spec_path, parameters.out_meanwinmat_pow_path, parameters.out_meanwinmat_freqs_path)
    end
end
function do_one(in_eeg_name, in_eeg_dir, out_fft_path, out_hilbert_path, out_sliding_window_prefix, frequency, sliding_window_average_plot, sliding_window_average_fft_plot, out_faxisall_path, out_spec_path, out_meanwinmat_pow_path, out_meanwinmat_freqs_path)
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
    to_csv(pow, out_fft_path)
    to_csv(freqs, 'processed/freqtageeg/freqs.csv')
    
    
    %% Plot the spectrum. For visualization purposes, only plot the frequencies of interest.
    %plot_fft_all_sensors(pow, freqs, num_freqs)
    %plot_fft_sensor(pow, freqs, num_freqs, oz_id)
    
    
    %% Run FFT on single-trials.
    [spec] = freqtag3D_FFT(dataset, stimulus_start:stimulus_end, sampling_rate); %No need to average over the third dimension here.
    
    %% Run sliding window analysis.
    [trialpow,winmat3d,phasestabmat,trialSNR] = freqtag_slidewin(dataset, 0, stimulus_start:stimulus_end, stimulus_start:stimulus_end, frequency, 600, 500, out_sliding_window_prefix);


    %% Averaging across sliding windows.
    meanwinmat = mean(winmat3d, 3);
    plot_sliding_window_average(meanwinmat, sprintf('Mean of moving windows at %i Hz', frequency), 'Sample points', 'Voltage', sliding_window_average_plot)
    

    %% Project sliding window average into frequency domain.
    % Note the sample rate is 600 Hz at 6 Hz
    [meanwinmat_pow, meanwinmat_phase, meanwinmat_freqs] = freqtag_FFT(meanwinmat, 600); 
    plot_sliding_window_average_FFT(meanwinmat_freqs, meanwinmat_pow, sprintf('Power spectrum of the mean window shifted at %i Hz', frequency), 'Frequency (Hz)', sliding_window_average_fft_plot)
    to_csv(meanwinmat_pow, out_meanwinmat_pow_path)
    to_csv(meanwinmat_freqs, out_meanwinmat_freqs_path)

    %% Plot spectrum
    %plot_single_fft_all_sensors(faxisall, spec)
    %plot_single_fft_one_sensor(spec, faxisall, oz_id)
    to_csv(faxisall, out_faxisall_path)
    to_csv(spec, out_spec_path)
    
    
    %% Run Hilbert Transform. Check the frequency of interest.
    [powermat, phasemat, complexmat] = freqtag_HILB(mean(dataset(:, stimulus_start:stimulus_end, :), 3), frequency, 8, oz_id, 0, sampling_rate);
    to_csv(powermat, out_hilbert_path)
    
    
    %% Plot frequency over time
    %plot_hilbert(powermat, oz_id, 'Hilbert Transform ({FREQUENCY} Stimuli)')
    
    
    %% Close all figures and clean the worskpace
    %close all
    %clear all
end


%% Plot functions.
function plot_sliding_window_average(meanwinmat, fig_title, fig_x_label, fig_y_label, out_path)
    % Plot that fancy sliding window average you calculated.
    fig = figure('visible', 'off');
    
    plot(meanwinmat')
    title(fig_title)
    xlabel(fig_x_label)
    ylabel(fig_y_label)

    exportgraphics(fig, out_path, 'Resolution', 300)
    close(fig)
end
function plot_sliding_window_average_FFT(meanwinmat_freqs, meanwinmat_pow, fig_title, fig_x_label, out_path)
    % Plot that fancy sliding window average you calculated - but now after it's been FFT'd!
    fig = figure('visible', 'off');

    bar(meanwinmat_freqs(2:20), meanwinmat_pow(32,2:20)');
    title(fig_title)
    xlabel(fig_x_label)

    exportgraphics(fig, out_path, 'Resolution', 300)
    close(fig)
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
    ax = gca;
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
function to_csv(mat, out_filename)
    % Write a matrix to a csv file.
    writematrix(mat, out_filename);

    % Pause to hopefully allow MatLab to actually write more than one file :(
    pause(2)
end
function [EEG] = load_dataset(file_name, directory)
    % Load a dataset.
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
