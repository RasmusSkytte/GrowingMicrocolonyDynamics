%close all; clearvars;
rng(0)

% This script run a minimizer to find optimal fit values
D = 0.8 * 1e-9 * 3600 * 1e12; % in micrometer^2 / h

tau = 20 / 60; % 20 minute doubling time
g_max = log(2) / tau;

% Set startpoint and bounds
%     n_0    log10(p.D)   mu*g_max  K    
x0 = [0.01   log10(D)     g_max/0.64  1e-1 ];

% Load the fitting data
load('../../experiments/GrowthData.mat')
T  = dataset{1, 1};
BF = [dataset{:, 2}];

% First dataset is only zeros
BF = BF(2:end, :);
T  = T(2:end);

% Get mean and error
BF(BF==0) = nan;
mBF = nanmean(BF, 2);
dBF = nanstd(BF, [], 2); 

% Truncate reference data
I = T > 16;
mBF(I) = [];
dBF(I) = [];
T(I)   = [];

% Check for previous fit
path = sprintf('../fits/GrowthParams.mat');
if exist(path, 'file')
    load(path, 'x');
    x0 = x
end

% Call minimizier
x = fminsearch(@(x)fitGrowthCurve(x, T, mBF, dBF, true), x0, optimset('MaxIter', 10*numel(x0), 'UseParallel', true));

% Create output folder
if ~exist('../fits', 'dir')
    mkdir('../fits')
end

save('../fits/GrowthParams.mat', 'x');
