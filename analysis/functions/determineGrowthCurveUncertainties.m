close all; clearvars;

% parpool('local', 16);

% After minimizing, the script checks if the fit is converged and estimates
% the uncertainties on the parameters.
% The script uses the previous (if any) run as the starting parameters for
% the min search. If the values are not converged, try repeating the script
% a few times.


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
end

% Reimpose bounds
x(x < lb) = lb(x < lb);
x(x > ub) = ub(x > ub);

% Determine fit uncertainties
dx = nan(size(x));

% Get current chi2 value
chi2 = fitGrowthCurve(x, T, mBF, dBF, true);
chi2_p = nan(size(x));
chi2_m = nan(size(x));

% Perturb all parameters
for i = 1:numel(x)
    if x(i) == 0
        continue
    end
    
    % Create perubation vector
    p = ones(size(x));
    p(i) = 0;
    
    % Vary parameters a little
    pertubation = 10.^(-4:0.5:-0.5);

    % Determine the size of the positve pertubation
    positivePertubations = arrayfun(@(f)fitGrowthCurve(x.*p + f*x.*(1-p), T, mBF, dBF, false), 1 + pertubation);
    
    % Determine the size of the positve pertubation
    negativePertubations = arrayfun(@(f)fitGrowthCurve(x.*p + f*x.*(1-p), T, mBF, dBF, false), 1 - pertubation);
    
    % Assert chi2 are higher
    if any([positivePertubations negativePertubations] < chi2)
        warning('Fit has not (fully) converged (chi2 lower at x_%d%+.1g)', i, [pertubation(find(positivePertubations < chi2, 1)) -pertubation(find(negativePertubations < chi2, 1))]);
    end
    
    % Assert chi2 are high enough
    if ~any(positivePertubations > chi2 + 1)
        error('Positive pertubation not large enough');
    end
    if ~any(negativePertubations > chi2 + 1)
        error('Positive pertubation not large enough');
    end
    
    % Set the initial pertubations
    fp = 1 + pertubation(find((positivePertubations - chi2) > 1, 1));
    fm = 1 - pertubation(find((negativePertubations - chi2) > 1, 1));
    
    % Run a quick search for minimal value
    fp = bisection(@(f)fitGrowthCurve(x.*p + f*x.*(1-p), T, mBF, dBF, false)-(chi2+1), 1, fp);
    fm = bisection(@(f)fitGrowthCurve(x.*p + f*x.*(1-p), T, mBF, dBF, false)-(chi2+1), fm, 1);
    
    % Check minimal values are "close"
    chi2_p(i) = fitGrowthCurve(x.*p + fp*x.*(1-p), T, mBF, dBF, false);
    chi2_m(i) = fitGrowthCurve(x.*p + fm*x.*(1-p), T, mBF, dBF, false);
    
    % Store the value
    dx(i) = mean([fp*x(i) - x(i) x(i) - fm*x(i)]);
    
end

% Report warnings
for i = 1:numel(x) 
    if x(i) == 0
        continue
    end
    fprintf('chi2_p - chi2 = %.3f\n', chi2_p(i) - chi2);
    fprintf('chi2_m - chi2 = %.3f\n', chi2_m(i) - chi2);
    
    if any(([chi2_p(i) chi2_m(i)] - 1 - chi2) / chi2 > 1e-2)
        warning('Uncertianty value %d deviates more than 1 percent', i);  
    elseif any(([chi2_p(i) chi2_m(i)] - 1 - chi2) / chi2 > 1e-3)
        warning('Uncertianty value %d deviates more than 1 permille', i);  
    end
end

save(sprintf('../fits/GrowthParams_Model_%d.mat', n), 'x', 'dx');