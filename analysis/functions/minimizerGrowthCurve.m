close all; clearvars;

% This script run a minimizer to find optimal fit values

% Running the fit using all parameters
% n = 1;
%       n_0 log10(p.D) g_max   K     delta   mu     R0
% lb = [0.001    4     2.00    0.01   0      0.2    0.1  ];
% x0 = [0.026    5.25  3.20    0.7    1e-9   0.80   1.7  ];
% ub = [0.500    7     4.00    2.0    0.100  1.00   10   ];

% Running the fit without intrinsic decay of bacteria
n = 2;
%     n_0     log10(p.D)  g_max     K      delta    mu        R0
x0 = [0.007   6.07       4         0.57     0      0.75      1.4  ];
lb = [0.001    4         1.0000    0.01     0      0.2       0.1  ];
ub = [0.050    7         5.0000    2.0      0      1.00      10   ];


% Load the fitting data
load('../../experiments/Datasets/DataSet_1')
BF = [dataset{1}{2} dataset{2}{2} dataset{3}{2} dataset{4}{2} dataset{5}{2}];
T = dataset{1}{1};

% First dataset is only zeros
BF = BF(2:end, :);
T  = T(2:end);

% Get mean and error
t = T;
BF(BF==0) = nan;
mBF = nanmean(BF, 2);
dBF = nanstd(BF, [], 2) ./ sqrt(sum(~isnan(BF), 2));

% Check for previous fit
path = sprintf('../fits/GrowthParams_Model_%d.mat', n);
if exist(path, 'file')
    load(path, 'x');
    x0 = x;
end

% Call minimizier
% x = fmincon(@(x)fitGrowthCurve(x, T, mBF, dBF, true), x0, [], [], [], [], lb, ub, [], optimoptions('fmincon', 'UseParallel', true));
x = fminsearch(@(x)fitGrowthCurve(x, T, mBF, dBF, true), x0);

% Create output folder
if ~exist('../fits', 'dir')
    mkdir('../fits')
end

save(sprintf('../fits/GrowthParams_Model_%d.mat', n), 'x');