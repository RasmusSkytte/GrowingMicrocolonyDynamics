close all; clearvars; clc;

% Load the fitted values
load('../fits/GrowthParams.mat');

% Get configration struct
p = getConfiguration(x);

% Load the fitting data
p.T_end = 24;

% Solve model with current parameters
[t, C, ~, N, p] = solveModel(p);

% Plot the best result
fh = figure(); clf; hold on;
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1;
ax.Box = 'on';

% Compute the velocity of expantion
v = C;
tt = t;
v(diff(C) == 0)  = [];
tt(diff(C) == 0) = [];

v = diff(v) ./ diff(tt);
tt = tt(1:(end-1)) + diff(tt) / 2;

% Determine the surface nutrient level
f = @(C)min(1,max(0, (C-p.r)./([p.r(2:end); p.rNp]-p.r)));

nC = nan(size(t));
for i = 1:numel(nC)
    fC = f(C(i));
    I = find(fC, 1, 'last');
    nC(i) = N(i, I-1) * (1-fC(I)) + N(i, I) * fC(I);
end

% Load the fitted values
load('../fits/PhageAttackParams_onlyDelta.mat');
p = getConfiguration(x, y, 3);

yyaxis left
plot(ax, tt, v, 'LineWidth', 2);
ylabel('Colony growth rate ({\mu}m / h)', FontSize = 18)
ax.XLim = [0 24];
ax.YLim = [0 80];

yyaxis right
plot(ax, t, p.epsilon * p.dR * p.g_max * nC ./ (nC + p.K), 'LineWidth', 2);
ylabel('Phage killing rate ({\mu}m / h)', FontSize = 20)
ax.YLim = [0 80];

xlabel('Time (h)')

if ~exist('../../figures/Figure 6/' , 'dir')
    mkdir('../../figures/Figure 6')
end

saveas(fh, '../../figures/Figure 6/Fig6.png')