close all; clearvars;

% load the data
load ../../experiments/data

GFP_radius      = GFP_radius_1_80;
GFP_radius_alt1 = GFP_radius_1_60;
GFP_radius_alt2 = GFP_radius_2_20;


% Replace NaN values
GFP_radius(isnan(GFP_radius)) = 0;

% Get the number of colonies and number of time points
nT = size(radius, 1);
nC = size(radius, 2);

% B4X2, series 22
%    Problem: First 2 frames includes bubbles from phage spray
%    Solution: Ommit first 2 frames
GFP_radius(1:2, 22) = nan;
radius(1:2, 22) = nan;

% B6X1, series 26
%    Problem: Frame 14 failed to capture properly
%    Solution: Ommit frame 14
GFP_radius(14, 26) = nan;
radius(14, 26) = nan;

% B6X_, series 26 - 30
%    Problem: Smallest GFP radius behaves as outlier
%    Solution: Ommit frame 1
GFP_radius(1, 26:30) = nan;
radius(1, 26:30) = nan;


if ~exist('../../experiments/Radius Curves', 'dir')
    mkdir('../../experiments/Radius Curves')
end

% Loop over all colonies and plot radius curves
B = [1 1 1 1 1 3 3 3 3 3 5 5 5 5 5 2 2 2 2 2 4 4 4 4 4 6 6 6 6 6];
X = arrayfun(@(k)sum(B(1:k)==B(k)), 1:numel(B));

for n = 1:25

    fh = figure(1); clf; hold on;
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 1;
    ax.Box = 'on';

    % Check for filtered data:
    f = ~isnan(radius(:, n));

    % Plot curves
    plot(ax, TT(f)+dT(n), radius(f, n), 'k', 'Linewidth', 2)
    plot(ax, TT(f)+dT(n), GFP_radius(f, n), 'g', 'Linewidth', 2)
    %     plot([T_i(n) T_i(n)], [0 400], ':k')

    % Label graph
    xlabel('Time (hours)')
    ylabel('Radius ({\mu}m)')

    xlim([0 25])
    ylim([0 400])

    % Save graph
    saveas(fh, sprintf('../../experiments/Radius Curves/B%dX%d.png', B(n), X(n)))

end

% Manuel grouping
for n = 1:5

    switch n
        case 1 % Growth rate runs
            d = [61 62 63 64 65];
            spath = 'GrowthData.mat';

        case 2 % Phage attack runs
            d = [12 13 14 21 22 23 24 25 31 32 33 35 41 42 43 44 51 52 53 54];
            spath = 'PhageData.mat';

        case 3 % "Outliers"
            d = [11 15 34 45 55];
            spath = 'OutlierData.mat';

        case 4 % Phage attack runs with higher threshold
            d = [12 13 14 21 22 23 24 25 31 32 33 35 41 42 43 44 51 52 53 54];
            spath = 'PhageData_higherTreshold_1_60.mat';
            %GFP_radius = GFP_radius_alt2;

        case 5 % Phage attack runs with higher threshold
            d = [12 13 14 21 22 23 24 25 31 32 33 35 41 42 43 44 51 52 53 54];
            spath = 'PhageData_higherTreshold_2_20.mat';
           % GFP_radius = GFP_radius_alt4;
    end

    % Map the odd indicies
    l = mod(floor(d/10), 2) == 1;
    I = nan(size(d));
    I(l) = floor(floor(d(l) / 10 ) / 2) * 5 + mod(d(l), 10);
    % Map the even indicies
    I(~l) = (5-(6-floor(d(~l)/10)) / 2)*5 + mod(d(~l), 10);

    % Prepare data set
    dataset = cell(numel(I), 3);
    cc = lines(numel(I));
    for i = 1:numel(I)

        % Store data
        dataset{i, 1} = TT' + dT(I(i));
        dataset{i, 2} = radius(:, I(i));
        dataset{i, 3} = GFP_radius(:, I(i));

    end

    % Save fitting data
    save(sprintf('../../experiments/%s', spath), 'dataset')

end