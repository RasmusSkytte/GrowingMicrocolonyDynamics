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

chi2 = fitGrowthCurve(x, T, mBF, dBF, false);
reduced_chi2 = chi2 / (length(mBF) - 4);

print(reduced_chi2)
