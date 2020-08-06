close all; clearvars; clc;

% Select growth model
n = 2;

% Load the fitting data
load('../../experiments/Datasets/DataSet_1')
BF = [dataset{1}{2} dataset{2}{2} dataset{3}{2} dataset{4}{2} dataset{5}{2}];
T = dataset{1}{1};

% First dataset is only zeros
BF = BF(2:end, :);
T  = T(2:end);

% Get mean and error
BF(BF == 0) = nan;
mBF = nanmean(BF, 2);
dBF = nanstd(BF, [], 2) ./ sqrt(sum(~isnan(BF), 2));

% Load the fitted values
load(sprintf('../fits/GrowthParams_Model_%d.mat', n));

% Get configration struct
p = getConfiguration(x);

% Load the fitting data
p.T_end = T;

% Solve model with current parameters
[t, C] = solveModel(p);

% Plot the best result
fh = figure(n); clf; hold on;
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1;
ax.Box = 'on';

plot(ax, t, C, 'b', 'LineWidth', 2);
errorbar(T, mBF, dBF, '.k', 'LineWidth', 2)

xlabel('Time (h)')
ylabel('Colony Radius ({\mu}m)')

saveas(fh, sprintf('../fits/GrowthParams_Model_%d.png', n))
saveas(fh, '../../figures/fig1.png')

% Get current chi2 value
% Locate the points to compare:
I = find(ismember(t, T));
I = I(1:2:end);

% Determine ndf
ndf = numel(T)-sum(x > 0);

% Compute loss
y = C(I);
chi2 = sum(  (mBF - y).^2./dBF.^2 );
fprintf('--- Fit quality ---\n')
fprintf('\\chi^2 / ndf. = %.3f\n\n', chi2 / ndf)


% Report fit values
fprintf('--- Fit values ---\n')

fprintf('Diffusion constant: \n')
v = 10^x(2) / 10^floor(x(2)); % Value
d = (10^x(2) * log(10) * dx(2)) / 10^floor(x(2)); % Error
e = floor(x(2)); % Exponent

fprintf('(%.3f +- %.3f) * 10^%d {\\mu}m^2 / h\n', v, d, e)

fprintf('Maximal growth rate: \n')
v = x(3) * 1 / (1 + x(4)); % Value
d = sqrt(dx(3) / (1 + x(4))^2 + (dx(4) * x(3) * 1 / (1 + x(4))^2)^2); % Error
fprintf('(%.3f +- %.3f) 1 / h\n', v, d)

fprintf('Initial nutrients: \n')
v = x(1); % Value
d = dx(1); % Error
e = floor(log10(x(1))); % Exponent
fprintf('(%.3f +- %.3f) * 10^%d\n', v/10^e, d/10^e, e)

fprintf('Monod constant: \n')
v = x(4)*x(1); % Value
d = sqrt((x(1)*dx(4)).^2 + (x(4)*dx(1)).^2); % Error
e = floor(log10(v)); % Exponent
fprintf('(%.3f +- %.3f)* 10^%d\n', v/10^e, d/10^e, e)

fprintf('Packing density: \n')
v = x(6); % Value
d = dx(6); % Error
fprintf('(%.4f +- %.4f)\n', v, d)

fprintf('Initial radius: \n')
v = (3/(4*pi)*x(7)*1.3)^(1/3) / x(6); % Value
d = sqrt(((3/(4*pi)*x(7)*1.3)^(1/3) / x(6)^2 * dx(6))^2 + (dx(7)*1/3*x(7)^(-2/3)*(3/(4*pi)*1.3)^(1/3) / x(6))^2);
fprintf('(%.4f +- %.4f)\n', v, d)