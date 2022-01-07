%% This script plots the spectrums we used for our continuous alpha analyses.

freqs = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28, 28.5, 29, 29.5, 30, 30.5, 31, 31.5, 32, 32.5, 33, 33.5, 34, 34.5, 35, 35.5, 36, 36.5, 37, 37.5, 38, 38.5, 39, 39.5, 40, 40.5, 41, 41.5, 42, 42.5, 43, 43.5, 44, 44.5, 45, 45.5, 46, 46.5, 47, 47.5, 48, 48.5, 49, 49.5, 50, 50.5, 51, 51.5, 52, 52.5, 53, 53.5, 54, 54.5, 55, 55.5, 56, 56.5, 57, 57.5, 58, 58.5, 59, 59.5, 60, 60.5, 61, 61.5, 62, 62.5, 63, 63.5, 64, 64.5, 65, 65.5, 66, 66.5, 67, 67.5, 68, 68.5, 69, 69.5, 70, 70.5, 71, 71.5, 72, 72.5, 73, 73.5, 74, 74.5, 75, 75.5, 76, 76.5, 77, 77.5, 78, 78.5, 79, 79.5, 80, 80.5, 81, 81.5, 82, 82.5, 83, 83.5, 84, 84.5, 85, 85.5, 86, 86.5, 87, 87.5, 88, 88.5, 89, 89.5, 90, 90.5, 91, 91.5, 92, 92.5, 93, 93.5, 94, 94.5, 95, 95.5, 96, 96.5, 97, 97.5, 98, 98.5, 99, 99.5, 100, 100.5, 101, 101.5, 102, 102.5, 103, 103.5, 104, 104.5, 105, 105.5, 106, 106.5, 107, 107.5, 108, 108.5, 109, 109.5, 110, 110.5, 111, 111.5, 112, 112.5, 113, 113.5, 114, 114.5, 115, 115.5, 116, 116.5, 117, 117.5, 118, 118.5, 119, 119.5, 120, 120.5, 121, 121.5, 122, 122.5, 123, 123.5, 124, 124.5, 125, 125.5, 126, 126.5, 127, 127.5, 128, 128.5, 129, 129.5, 130, 130.5, 131, 131.5, 132, 132.5, 133, 133.5, 134, 134.5, 135, 135.5, 136, 136.5, 137, 137.5, 138, 138.5, 139, 139.5, 140, 140.5, 141, 141.5, 142, 142.5, 143, 143.5, 144, 144.5, 145, 145.5, 146, 146.5, 147, 147.5, 148, 148.5, 149, 149.5, 150, 150.5, 151, 151.5, 152, 152.5, 153, 153.5, 154, 154.5, 155, 155.5, 156, 156.5, 157, 157.5, 158, 158.5, 159, 159.5, 160, 160.5, 161, 161.5, 162, 162.5, 163, 163.5, 164, 164.5, 165, 165.5, 166, 166.5, 167, 167.5, 168, 168.5, 169, 169.5, 170, 170.5, 171, 171.5, 172, 172.5, 173, 173.5, 174, 174.5, 175, 175.5, 176, 176.5, 177, 177.5, 178, 178.5, 179, 179.5, 180, 180.5, 181, 181.5, 182, 182.5, 183, 183.5, 184, 184.5, 185, 185.5, 186, 186.5, 187, 187.5, 188, 188.5, 189, 189.5, 190, 190.5, 191, 191.5, 192, 192.5, 193, 193.5, 194, 194.5, 195, 195.5, 196, 196.5, 197, 197.5, 198, 198.5, 199, 199.5, 200, 200.5, 201, 201.5, 202, 202.5, 203, 203.5, 204, 204.5, 205, 205.5, 206, 206.5, 207, 207.5, 208, 208.5, 209, 209.5, 210, 210.5, 211, 211.5, 212, 212.5, 213, 213.5, 214, 214.5, 215, 215.5, 216, 216.5, 217, 217.5, 218, 218.5, 219, 219.5, 220, 220.5, 221, 221.5, 222, 222.5, 223, 223.5, 224, 224.5, 225, 225.5, 226, 226.5, 227, 227.5, 228, 228.5, 229, 229.5, 230, 230.5, 231, 231.5, 232, 232.5, 233, 233.5, 234, 234.5, 235, 235.5, 236, 236.5, 237, 237.5, 238, 238.5, 239, 239.5, 240, 240.5, 241, 241.5, 242, 242.5, 243, 243.5, 244, 244.5, 245, 245.5, 246, 246.5, 247, 247.5, 248, 248.5, 249, 249.5];

spectrums = get_spectrums();
average_spectrum = get_average_spectrum(spectrums)

plot_spectrums(freqs(1:80), spectrums(1:80,:));
plot_average_spectrum(freqs(1:80), average_spectrum(1:80));


function [spectrums] = get_spectrums()
    %% Get average occipital spectrum of each subject for all 2s intervals.
    pow_paths = dir('processed/eeg_alphas/sub-*_average_power.mat');
    cells_of_data = cell(size(pow_paths));
    for ii=1:length(pow_paths)
        file = pow_paths(ii);
        aFileData = load(fullfile(file.folder, file.name));
        spectrum = aFileData.average_power;
        occipital_spectrum = spectrum([20 31 19 7 8 9 10], :)';
        average_occipital_spectrum = mean(occipital_spectrum, 2);
        cells_of_data{ii} = average_occipital_spectrum;
    end
    spectrums = cat(2, cells_of_data{:});
end


function [average_spectrum] = get_average_spectrum(spectrums)
    %% Get the average of our spectrums.
    average_spectrum = mean(spectrums, 2);
end


function plot_spectrums(x_values, y_values)
    %% Plot spectrums.
    figure;
    plot(x_values, y_values, 'LineWidth', 3);
    fig = gcf;
    ax = gca;
    axis(gca, 'tight')
    box off;

    patch([8.5 8.5 12.5 12.5], [0 6 6 0], [.7 .7 .7])
    set(gca,'children',flipud(get(gca,'children')))
    set(gca, 'Layer', 'top')

    ylabel('Power (µV^2)');
    xlabel('Frequency (Hz)');

    yticks([10 20]);

    ax.FontSize = 30;
    ax.LineWidth = 3;

    text(6.5, 8.5, ['Analysis Window'], 'FontSize', 30)
    text(8, 7, ['8.5-12.5Hz'], 'FontSize', 30)

    % Make window size enormous.
    set(gcf, 'Position',  [0, 0, 2000, 2000])
    p = get(gca,'Position');
    set(gca,'Position',[p(1) p(2)+0.1 p(3)-0.1 p(4)-0.2]);

    mkdir 'processed/plots/'
    saveas(fig, 'processed/plots/alpha_spectrum_continuous.png')
end


function plot_average_spectrum(x_values, y_values)
    %% Plot average spectrum.
    figure;
    plot(x_values, y_values, 'LineWidth', 3);
    fig = gcf;
    ax = gca;
    axis(gca, 'tight')
    box off;

    patch([8.5 8.5 12.5 12.5], [0 4 4 0], [.7 .7 .7])
    set(gca,'children',flipud(get(gca,'children')))
    set(gca, 'Layer', 'top')

    ylabel('Power (µV^2)');
    xlabel('Frequency (Hz)');

    yticks([5]);

    ax.FontSize = 30;
    ax.LineWidth = 3;

    text(6.5, 5, ['Analysis Window'], 'FontSize', 30)
    text(8, 4.5, ['8.5-12.5Hz'], 'FontSize', 30)

    % Make window size enormous.
    set(gcf, 'Position',  [0, 0, 2000, 2000])
    p = get(gca,'Position');
    set(gca,'Position',[p(1) p(2)+0.1 p(3)-0.1 p(4)-0.2]);

    mkdir 'processed/plots/'
    saveas(fig, 'processed/plots/alpha_average_spectrum_continuous.png')
end
