function chi2 = fitGrowthCurve(x, T, mBF, dBF, reportChi2)

% Construct the confugration
p = getConfiguration(x);

% Load the fitting data
p.T_end = T;

% Solve model with current parameters
[t, C] = solveModel(p);

% Locate the points to compare:
I = find(ismember(t, T));
I = I(1:2:end);

% Determine ndf
ndf = numel(T)-sum(x > 0);

% Compute loss
y = C(I);
chi2 = sum(  (mBF - y).^2./dBF.^2 );
if reportChi2
    fprintf('\t\\chi^2 / ndf. = %.2f, x = [%s]\n', chi2 / ndf,  sprintf([repmat('%.3f, ', 1, numel(x)-1) '%.3f'], x) )
end

end
