close all; clearvars;

% After minimizing, the script checks if the fit is converged and estimates
% the uncertainties on the parameters.
% The script uses the previous (if any) run as the starting parameters for
% the min search. If the values are not converged, try repeating the
% minimizer script a few times.

% Load the fitting data
load('../../experiments/GrowthData.mat')

T = dataset{1, 1};
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

% Load the growth rate fit parameters
load('../fits/GrowthParams.mat', 'x');

% Determine fit uncertainties
dx = nan(size(x));

% Get current chi2 value
chi2 = fitGrowthCurve(x, T, mBF, dBF, true)
chi2_p = nan(size(x));
chi2_m = nan(size(x));

% Perturb all parameters
for i = 1:numel(x)

    % Create perubation vector
    p = ones(size(x));
    p(i) = 0;

    % Vary parameters a little
    pertubation = 10.^(-4:0.5:0);

    % Determine the size of the positve pertubation
    positivePertubations = nan(size(pertubation));
    for j = 1:numel(pertubation)
        positivePertubations(j) = fitGrowthCurve(x.*p + (1+pertubation(j))*x.*(1-p), T, mBF, dBF, false);
        if positivePertubations(j) > chi2 + 1
            break
        end
    end
    positivePertubations(j) - chi2;

    % Determine the size of the positve pertubation
    negativePertubations = nan(size(pertubation));
    for j = 1:numel(pertubation)
        negativePertubations(j) = fitGrowthCurve(x.*p + (1-pertubation(j))*x.*(1-p), T, mBF, dBF, false);
        if negativePertubations(j) > chi2 + 1
            break
        end
    end
    negativePertubations(j) - chi2;

    % Assert chi2 are higher
    if any(positivePertubations < chi2)
        error('Fit has not (fully) converged (chi2 lower at x_%d%+.1g)', i,  pertubation(find(positivePertubations < chi2, 1)));
    end
    if any(negativePertubations < chi2)
        error('Fit has not (fully) converged (chi2 lower at x_%d%+.1g)', i, -pertubation(find(negativePertubations < chi2, 1)));
    end

    % Assert chi2 are high enough
    if ~any(positivePertubations > (chi2 + 1))
        error('Positive pertubation for x_%d not large enough', i);
    end
    if ~any(negativePertubations > (chi2 + 1))
        error('Positive pertubation for x_%d not large enough', i);
    end

    % Set the initial pertubations
    fp = 1 + pertubation(find(positivePertubations > (chi2 + 1), 1));
    fm = 1 - pertubation(find(negativePertubations > (chi2 + 1), 1));

    % Run a quick search for minimal value
    fp = bisection(@(f)fitGrowthCurve(x.*p + f*x.*(1-p), T, mBF, dBF, false)-(chi2+1), 1, fp);
    fm = bisection(@(f)fitGrowthCurve(x.*p + f*x.*(1-p), T, mBF, dBF, false)-(chi2+1), fm, 1);

    % Check minimal values are "close"
    chi2_p(i) = fitGrowthCurve(x.*p + fp*x.*(1-p), T, mBF, dBF, false);
    chi2_m(i) = fitGrowthCurve(x.*p + fm*x.*(1-p), T, mBF, dBF, false);

    % Store the value
    dx(i) = mean([fp*x(i) - x(i) x(i) - fm*x(i)]);

    fprintf('chi2_p - chi2 = %.3f\n', chi2_p(i) - chi2);
    fprintf('chi2_m - chi2 = %.3f\n', chi2_m(i) - chi2);

    if any(([chi2_p(i) chi2_m(i)] - 1 - chi2) / chi2 > 1e-2)
        warning('Uncertianty value %d deviates more than 1 percent', i);
    end
end

save('../fits/GrowthParams.mat', 'x', 'dx');
