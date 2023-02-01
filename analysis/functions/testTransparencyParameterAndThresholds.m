close all; clearvars;

% Define thresholds
thresholds = 1:0.1:2.5;

% Prepare data
fits = cell(numel(thresholds), 1);

% Loop over data sets
for k = 1:numel(thresholds)

    % Load the fit
    path = sprintf('../fits/TransparencyParameter/Params_%s.mat', strrep(sprintf('%.1f', thresholds(k)), '.', '_'));
    load(path, 'y');
    fits{k} = y;

end

if ~exist('../../figures/Figure S2', 'dir')
    mkdir('../../figures/Figure S2')
end

% Prepare figure
fh1 = figure(1); clf; hold on; box on;
ax1 = gca;
ax1.FontSize = 20;
ax1.LineWidth = 1;
% ax1.TickLength = [0.02 0.05];
% ax1.XLim = [0 25];
% ax1.YLim = [0 400];

plot(ax1, thresholds, arrayfun(@(i)fits{i}(3), 1:numel(thresholds)), 'LineWidth', 2)
plot(ax1, [thresholds(1) thresholds(end)], [1 1], '--k', 'Linewidth', 2)

% ax1.XTick = 1:0.2:2;

ytickformat('%.1f');

xlabel(ax1, 'Threshold value')
ylabel(ax1, '\nu')

saveas(fh1, '../../figures/Figure S2/FigS2.png')