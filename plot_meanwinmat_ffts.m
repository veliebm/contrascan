% Get the mean of the FFTs of meanwinmat we calculated.

%% Load all pows into a single array.
FFT_dir = 'processed/freqtag_sliding_window_dir/FFT';
pow_paths = dir('processed/freqtag_sliding_window_dir/FFT/*pow.mat');
allData = cell(size(pow_paths));
for ii=1:length(pow_paths)
    file = pow_paths(ii);
    aFileData = load(fullfile(file.folder, file.name));
    allData{ii} = aFileData.pow;
end
data_array = cat(3, allData{:});


%% Load frequency data.
load('processed/freqtag_sliding_window_dir/FFT/sub-104_freqs.mat')


%% Mean data together.
mean_pow = mean(data_array, 3);
num_subjects = length(data_array(1,1,:));
stderror= std(data_array, 0, 3) / sqrt(num_subjects);


%% Plot data.
x_values = freqs(2:11);
y_values = mean_pow(20, 2:11);
error = stderror(20, 2:11);

fig = bar(x_values, y_values)
ax = gca;

% Edit the ticks.
yticks([0, .1]);
xticks([6:6:36])
ax.XAxis.TickLength = [0, 0];
ax.LineWidth = 2;

% Make figure grayscale.
fig.FaceColor = '#808080'

% Make window size enormous.
set(gcf, 'Position',  [0, 0, 1000, 800])

% Remove top and right side of box.
box off;

% Label the axes.
ylabel('Power (µV^2)', 'FontSize', 40);
xlabel('Frequency (Hz)', 'FontSize', 40);
ax.FontSize = 40;

% Plot error bars.
hold on
err = errorbar(x_values, y_values, error, 'LineWidth', 2);
err.Color = [0 0 0];
err.LineStyle = 'none';
hold off

saveas(fig, 'processed/plots/oz_spectrum.png')
