% Plot our EEG only analyses.

%% Load variables.
eeg_analyses = load('./data/misc/eeg_alone_analyses.mat');
alpha = eeg_analyses.alpha;
ssvep = eeg_analyses.ssVEP;
taxis = eeg_analyses.taxis;

%% Make plot.
set(0, 'DefaultLineLineWidth', 4);
fig = plot(taxis, alpha);

hold on;
plot (taxis, ssvep);
hold off;

ylabel('Electrical Activity (µV)')
xlabel('Time (ms)')
ax = gca;
axis(ax, 'tight')
ax.FontSize = 40;
ax.ColorOrder = [0 0 0; .7 .7 .7; 0 0 0; 0 0 0];
box off;
ax.LineWidth = 3;

% Make window size enormous.
set(gcf, 'Position',  [0, 0, 2000, 2000])

% Plot a dashed line at y=0.
yline(0,'k--', 'LineWidth', 3)

legend('Alpha', 'ssVEP', '')

saveas(fig, './processed/plots/eeg_alone_analyses.png')
