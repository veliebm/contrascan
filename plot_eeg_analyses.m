% Plot our EEG only analyses.

%% Load variables.
eeg_analyses = load('./data/misc/eeg_alone_analyses.mat');
alpha = eeg_analyses.alpha(1:end-30);
ssvep = eeg_analyses.ssVEP(1:end-30);
taxis = eeg_analyses.taxis(1:end-30)/1000;

%% Make plot.
set(0, 'DefaultLineLineWidth', 12);
fig = plot(taxis, alpha);

hold on;
plot (taxis, ssvep);
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

% Plot a dashed line at y=0.
yline(0,'k--', 'LineWidth', 3)

legend('Alpha Amplitude', 'ssVEP Amplitude', '', 'location','northwest')

saveas(fig, './processed/plots/eeg_alone_analyses.png')
