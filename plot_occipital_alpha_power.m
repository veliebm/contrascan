%% This script plots the spectrums we used for our alpha analyses.

freqs = [0,1.25,2.50,3.75,5,6.25,7.50,8.75,10,11.25,12.5,13.75,15,16.25,17.5,18.75,20,21.25,22.5,23.75,25,26.25,27.5,28.75,30,31.25,32.5,33.75,35,36.25,37.5,38.75,40,41.25,42.5,43.75,45,46.25,47.5,48.75,50,51.25,52.5,53.75,55,56.25,57.5,58.75,60,61.25,62.5,63.75,65,66.25,67.5,68.75,70,71.25,72.5,73.75,75,76.25,77.5,78.75,80,81.25,82.5,83.75,85,86.25,87.5,88.75,90,91.25,92.5,93.75,95,96.25,97.5,98.75,100,101.25,102.5,103.75,105,106.25,107.5,108.75,110,111.25,112.5,113.75,115,116.25,117.5,118.75,120,121.25,122.5,123.75,125,126.25,127.5,128.75,130,131.25,132.5,133.75,135,136.25,137.5,138.75,140,141.25,142.5,143.75,145,146.25,147.5,148.75,150,151.25,152.5,153.75,155,156.25,157.5,158.75,160,161.25,162.5,163.75,165,166.25,167.5,168.75,170,171.25,172.5,173.75,175,176.25,177.5,178.75,180,181.25,182.5,183.75,185,186.25,187.5,188.75,190,191.25,192.5,193.75,195,196.25,197.5,198.75,200,201.25,202.5,203.75,205,206.25,207.5,208.75,210,211.25,212.5,213.75,215,216.25,217.5,218.75,220,221.25,222.5,223.75,225,226.25,227.5,228.75,230,231.25,232.5,233.75,235,236.25,237.5,238.75,240,241.25,242.5,243.75,245,246.25,247.5,248.75];

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

    ylabel('Power (??V^2)');
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

    ylabel('Power (??V^2)');
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
