% Plot our EEG only analyses.

%% Load variables.
eeg_analyses = load('./data/misc/eeg_alone_analyses.mat');
alpha = eeg_analyses.alpha(1:end-0);
ssvep = eeg_analyses.ssVEP(1:end-0);
taxis = eeg_analyses.taxis(1:end-0)/1000;

%% Make plot.
set(0, 'DefaultLineLineWidth', 8);
fig = plot(taxis, alpha);

hold on;
plot(taxis, ssvep);
plot(taxis, zeros(length(taxis),1), 'k--', 'LineWidth',3)
hold off;

ylabel('Voltage (µV)')
xlabel('Time (s)')
ax = gca;
set(gca,'fontname','calibri')
axis(ax, 'tight')
ax.FontSize = 40;
ax.ColorOrder = [0 0 0; .7 .7 .7; 0 0 0; 0 0 0];
box off;
ax.LineWidth = 3;
xticks([0 1 2 3 4 5])

% Make window size enormous.
set(gcf, 'Position',  [0, 0, 1000, 800])

ax.Children = ax.Children([3 2 1]);
legend('', 'ssVEP Amplitude', 'Alpha Amplitude', 'location','northwest', 'Box','off')

saveas(fig, './processed/plots/eeg_alone_analyses.svg')
