%% This script plots the spectrums we used for our alpha analyses.

freqs = [0,1.25000000000000,2.50000000000000,3.75000000000000,5,6.25000000000000,7.50000000000000,8.75000000000000,10,11.2500000000000,12.5000000000000,13.7500000000000,15,16.2500000000000,17.5000000000000,18.7500000000000,20,21.2500000000000,22.5000000000000,23.7500000000000,25,26.2500000000000,27.5000000000000,28.7500000000000,30,31.2500000000000,32.5000000000000,33.7500000000000,35,36.2500000000000,37.5000000000000,38.7500000000000,40,41.2500000000000,42.5000000000000,43.7500000000000,45,46.2500000000000,47.5000000000000,48.7500000000000,50,51.2500000000000,52.5000000000000,53.7500000000000,55,56.2500000000000,57.5000000000000,58.7500000000000,60,61.2500000000000,62.5000000000000,63.7500000000000,65,66.2500000000000,67.5000000000000,68.7500000000000,70,71.2500000000000,72.5000000000000,73.7500000000000,75,76.2500000000000,77.5000000000000,78.7500000000000,80,81.2500000000000,82.5000000000000,83.7500000000000,85,86.2500000000000,87.5000000000000,88.7500000000000,90,91.2500000000000,92.5000000000000,93.7500000000000,95,96.2500000000000,97.5000000000000,98.7500000000000,100,101.250000000000,102.500000000000,103.750000000000,105,106.250000000000,107.500000000000,108.750000000000,110,111.250000000000,112.500000000000,113.750000000000,115,116.250000000000,117.500000000000,118.750000000000,120,121.250000000000,122.500000000000,123.750000000000,125,126.250000000000,127.500000000000,128.750000000000,130,131.250000000000,132.500000000000,133.750000000000,135,136.250000000000,137.500000000000,138.750000000000,140,141.250000000000,142.500000000000,143.750000000000,145,146.250000000000,147.500000000000,148.750000000000,150,151.250000000000,152.500000000000,153.750000000000,155,156.250000000000,157.500000000000,158.750000000000,160,161.250000000000,162.500000000000,163.750000000000,165,166.250000000000,167.500000000000,168.750000000000,170,171.250000000000,172.500000000000,173.750000000000,175,176.250000000000,177.500000000000,178.750000000000,180,181.250000000000,182.500000000000,183.750000000000,185,186.250000000000,187.500000000000,188.750000000000,190,191.250000000000,192.500000000000,193.750000000000,195,196.250000000000,197.500000000000,198.750000000000,200,201.250000000000,202.500000000000,203.750000000000,205,206.250000000000,207.500000000000,208.750000000000,210,211.250000000000,212.500000000000,213.750000000000,215,216.250000000000,217.500000000000,218.750000000000,220,221.250000000000,222.500000000000,223.750000000000,225,226.250000000000,227.500000000000,228.750000000000,230,231.250000000000,232.500000000000,233.750000000000,235,236.250000000000,237.500000000000,238.750000000000,240,241.250000000000,242.500000000000,243.750000000000,245,246.250000000000,247.500000000000,248.750000000000];

pretrial = get_pretrial_spectrums();
posttrial = get_posttrial_spectrums();

plot_pretrial(freqs(1:30), pretrial(1:30,:));
plot_posttrial(freqs(1:30), posttrial(1:30,:));


function [pretrial] = get_pretrial_spectrums()
    %% Get average occipital posttrial spectrum of each subject.
    pow_paths = dir('processed/eeg_trial_alphas/sub-*_data-averagebaselinepower_trial_alpha.mat');
    cells_of_data = cell(size(pow_paths));
    for ii=1:length(pow_paths)
        file = pow_paths(ii);
        aFileData = load(fullfile(file.folder, file.name));
        spectrum = aFileData.averagebaselinepower;
        occipital_spectrum = spectrum([20 31 19 7 8 9 10], :)';
        average_occipital_spectrum = mean(occipital_spectrum, 2);
        cells_of_data{ii} = average_occipital_spectrum;
    end
    pretrial = cat(2, cells_of_data{:});
end


function [posttrial] = get_posttrial_spectrums()
    %% Get average occipital pretrial spectrum of each subject.
    pow_paths = dir('processed/eeg_trial_alphas/sub-*_data-averagerawpower_trial_alpha.mat');
    cells_of_data = cell(size(pow_paths));
    for ii=1:length(pow_paths)
        file = pow_paths(ii);
        aFileData = load(fullfile(file.folder, file.name));
        spectrum = aFileData.averagerawpower;
        occipital_spectrum = spectrum([20 31 19 7 8 9 10], :)';
        average_occipital_spectrum = mean(occipital_spectrum, 2);
        cells_of_data{ii} = average_occipital_spectrum;
    end
    posttrial = cat(2, cells_of_data{:});
end


function plot_pretrial(x_values, y_values)
    %% Plot post trial data.
    figure;
    plot(x_values, y_values, 'LineWidth', 3);
    fig = gcf;
    ax = gca;
    axis(gca, 'tight')
    box off;

    patch([8.75 8.75 12.5 12.5], [0 20 20 0], [.7 .7 .7])
    set(gca,'children',flipud(get(gca,'children')))
    set(gca, 'Layer', 'top')

    ylabel('Power (µV^2)');
    xlabel('Frequency (Hz)');

    yticks([10 20]);

    ax.FontSize = 30;
    ax.LineWidth = 3;

    text(7, 23, ['Analysis Window'], 'FontSize', 30)
    text(8, 21, ['8.75-12.5Hz'], 'FontSize', 30)

    % Make window size enormous.
    set(gcf, 'Position',  [0, 0, 2000, 2000])
    p = get(gca,'Position');
    set(gca,'Position',[p(1) p(2)+0.1 p(3)-0.1 p(4)-0.2]);

    mkdir 'processed/plots/'
    saveas(fig, 'processed/plots/alpha_pretrial.png')
end


function plot_posttrial(x_values, y_values)
    %% Plot post trial data.
    figure;
    plot(x_values, y_values, 'LineWidth', 3);
    fig = gcf;
    ax = gca;
    axis(gca, 'tight')
    box off;

    patch([8.75 8.75 12.5 12.5], [0 20 20 0], [.7 .7 .7])
    set(gca,'children',flipud(get(gca,'children')))
    set(gca, 'Layer', 'top')

    ylabel('Power (µV^2)');
    xlabel('Frequency (Hz)');

    yticks([10 20]);

    ax.FontSize = 30;
    ax.LineWidth = 3;

    text(7, 22.5, ['Analysis Window'], 'FontSize', 30)
    text(8, 21, ['8.75-12.5Hz'], 'FontSize', 30)

    % Make window size enormous.
    set(gcf, 'Position',  [0, 0, 2000, 2000])
    p = get(gca,'Position');
    set(gca,'Position',[p(1) p(2)+0.1 p(3)-0.1 p(4)-0.2]);

    mkdir 'processed/plots/'
    saveas(fig, 'processed/plots/alpha_posttrial.png')
end