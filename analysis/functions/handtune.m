%close all; clearvars; clc;

% Load the fitting data
load('../../experiments/GrowthData.mat')

T = dataset{1, 1};
BF = [dataset{:, 2}];

% First dataset is only zeros
BF = BF(2:end, :);
T  = T(2:end);

% Get mean and error
BF(BF == 0) = nan;
mBF = nanmean(BF, 2);
dBF = nanstd(BF, [], 2); % ./ sqrt(sum(~isnan(BF), 2));

% Load the fitted values
%load('../fits/GrowthParams.mat');

%     n_0     log10(p.D)  g_max    K       mu        R0
D = 0.8 * 1e-9 * 3600 * 1e12; % in micrometer^2 / h
D_max = 0.9 * 1e-9 * 3600 * 1e12; % in micrometer^2 / h


tau = 20 / 60; % 20 minute doubling time
% N(t) = 2^t
% N(0.333) = 2
% N(t) = exp(g*t)
% N(0.3333) = 2 = exp(g*0.333)
% l2(2) = g * 0.333
g_max = log(2) / tau

% Set startpoint and bounds
%     n_0    K    mu    R0  log10(p.D)   g_max
x = [0.005   0.05  0.75  100   log10(D)];

% Get configration struct
p = getConfiguration(x)

% Load the fitting data
p.T_end = T;

% Solve model with current parameters
[t, C, ~, N, p] = solveModel(p);

% Plot the best result
fh = figure(1); clf; hold on;
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1;
ax.Box = 'on';

plot(ax, t, C, 'b', 'LineWidth', 2);
errorbar(T, mBF, dBF, '.k', 'LineWidth', 2)

xlabel('Time (h)')
ylabel('Colony Radius ({\mu}m)')
