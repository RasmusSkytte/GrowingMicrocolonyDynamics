close all; clearvars; clc;

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
dBF = nanstd(BF, [], 2);

% Load the fitted values
load('../fits/GrowthParams.mat');



% Get configration struct
p = getConfiguration(x);

% Load the fitting data
p.T_end = T;

% Solve model with current parameters
[t, C, ~, N, p] = solveModel(p);

% Plot the best result
fh = figure(); clf; hold on;
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1;
ax.Box = 'on';

%I = mBF > 300;
I = T > 16;
plot(ax, t, C, 'b', 'LineWidth', 2);
errorbar(T(I),  mBF(I),  dBF(I),  '.', 'LineWidth', 2)
errorbar(T(~I), mBF(~I), dBF(~I), '.k', 'LineWidth', 2)

xlabel('Time (h)')
ylabel('Colony Radius ({\mu}m)')

saveas(fh, '../fits/GrowthParams.png')

sdir = '../../figures/Figure 3';
if ~exist(sdir, 'dir')
    mkdir(sdir)
end
saveas(fh, sprintf('%s/Fig3.png', sdir));

% Get current chi2 value
% Locate the points to compare:
II = find(ismember(t, T(~I)));
II = II(1:2:end);

% Determine ndf
ndf = numel(T)-sum(x > 0);

% Compute loss
y = C(II);
chi2 = sum(  (mBF(~I) - y).^2./dBF(~I).^2 );
fprintf('--- Fit quality ---\n')
fprintf('\\chi^2 / ndf. = %.3f\n\n', chi2 / ndf)


% Report fit values
fprintf('--- Fit values ---\n')

fprintf('Diffusion constant: \n')
v = 10^x(2) / 10^floor(x(2)); % Value
d = (10^x(2) * log(10) * dx(2)) / 10^floor(x(2)); % Error
e = floor(x(2)); % Exponent
fprintf('(%.3f +- %.3f) * 10^%d {\\mu}m^2 / h\n', v, d, e)

fprintf('mu g_max: \n')
v = x(3); % Value
d = dx(3); % Error
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