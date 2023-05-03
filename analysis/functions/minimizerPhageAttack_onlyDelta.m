%close all; clearvars;
rng(0)

if isempty(gcp('nocreate'))
    parpool('local', 12);
end


% Set starting parameters

%    dR        epsilon  nu
y0 = [30   1    0  ];

% Load the best fitting parameters (bacteria parameters)
load('../fits/GrowthParams.mat')

% Load the fitting data
load('../../experiments/PhageData.mat')

% Allocate cells for fitting data
tmp1 = [dataset{:, 1}]; % Time Data
tmp2 = [dataset{:, 2}]; % BF data
tmp3 = [dataset{:, 3}]; % GDP Data

T_i  = tmp1(1, :) + 10^0.5;

nd = length(dataset);
Time = cell(nd, 1);
BF   = cell(nd, 1);
GFP  = cell(nd, 1);
for d = 1:nd
    Time{d} = tmp1(:, d);
    BF{d}   = tmp2(:, d);
    GFP{d}  = tmp3(:, d);
end

% Check for previous fit
path = '../fits/PhageAttackParams_onlyDelta.mat';
if exist(path, 'file')
    load(path, 'y', 'fitImprovement');
    y0 = y(1:3);
    T_i = y(4:end);
    if (fitImprovement < 5)
        error('Fit has  converged');
    end
end

% Add T_i to fitting parameters
y0 = [y0 T_i];

% Determine start fit value
y_old = y0;
d_old = fitPhageAttack(x, y_old, 3, Time, BF, GFP, false);

% Call minimizier
y = fminsearch(@(y)fitPhageAttack(x, y, 3, Time, BF, GFP, true), y0, optimset('MaxIter', 10*numel(y0)));

% Report change in conditions
fprintf('Total parameter change: %.3f\n', norm(y-y_old))

d_new = fitPhageAttack(x, y, 3, Time, BF, GFP, false);
fitImprovement = 100*(d_old-d_new)/d_old;
fprintf('Fit improvement: %.3f %%\n', fitImprovement)

% Store fit
d = d_new;
save(path, 'y', 'd', 'fitImprovement');