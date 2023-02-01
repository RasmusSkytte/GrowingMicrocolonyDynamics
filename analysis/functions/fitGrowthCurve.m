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
chi2 = (mBF - y).^2./dBF.^2;
chi2 = sum(chi2);% ./ (1:length(chi2))');
if reportChi2
    fprintf('\t\\chi^2 / ndf. = %.5f, x = [%s]\n', chi2 / ndf,  sprintf([repmat('%.4f, ', 1, numel(x)-1) '%.3f'], x) )
    fprintf('\tmin_doubling_time = %.1f\n', 60*log(2) / (x(2) / (1 + x(3))))
end

end
