clearvars;

% load the data
load ../../experiments/data

% Get the number of colonies and number of time points
nT = size(radius, 1);
nC = size(radius, 2);

% Store the invasion times
T_i = nan(nC, 1);

% Some runs need some immediate filtering:
% 42, series 22,
%    Problem: GFP is not being detected by threshold level
%    Solution: Ommit data set

% 45, series 25
%    Problem: Image field of view contain two colonies
%    Solution: Ommit data set

% 61, series 26
%    Problem: Frame 14 failed to capture properly
%    Solution: Ommit frame 14
GFP_radius(14, 26) = nan;
radius(14, 26) = nan;

if ~exist('../../experiments/Radius Curves', 'dir')
    mkdir('../../experiments/Radius Curves')
end

% Loop over all colonies and plot radius curves
for n = 1:30

    fh = figure(1); clf; hold on; box on;

    % Check for filtered data:
    f = ~isnan(radius(:, n));

    % Plot curves
    plot(TT(f)+dT(n), radius(f, n), 'k')
    plot(TT(f)+dT(n), GFP_radius(f, n), 'g')
    %     plot([T_i(n) T_i(n)], [0 400], ':k')

    % Label graph
    xlabel('Time (hours)')
    ylabel('Radius ({\mu}m)')

    xlim([0 25])
    ylim([0 400])

    % Save graph
    saveas(fh, sprintf('../../experiments/Radius Curves/B%dX%d.png', ceil(n/5), n-(ceil(n/5)-1)*5))

end


if ~exist('../../experiments/Datasets', 'dir')
    mkdir('../../experiments/Datasets')
end

% Manuel grouping
for n = 1:7

    switch n
        case 1 % Growth rate runs
            d = [61 62 63 64 65];
        case 2 % Decaying runs (Untill atleast T = 20)
            d = [14 31 32 35 21 22 23 24 25 44];
        case 3 % Increasing runs
            d = [33 51 53 54 41 43];
        case 4 % Decaying and increasing
            d = [14 31 32 35 21 22 23 24 25 44 33 51 53 54 41 43];
        case 5 % Resistance ? (Inncreases again before T = 20)
            d = [11 12 13 15 34 52 55];
        case 6 % Data from set for with low z
            d = [14 31 32 35 21 23 25 41 43 44];
        case 7 % Data from set for with intermediate - high z
            d = [33 22 24 51 53 54];

        % Included runs:
        % All of W1
        % All of W2
        % All of W3
        % All of W4 - 42, 45
        % All of W5
        % All of W6

        % Omitted runs :
        % 42, series 22
        % 45, series 25
    end
    % Map the odd indicies
    l = mod(floor(d/10), 2) == 1;
    I = nan(size(d));
    I(l) = floor(floor(d(l) / 10 ) / 2) * 5 + mod(d(l), 10);
    % Map the even indicies
    I(~l) = (5-(6-floor(d(~l)/10)) / 2)*5 + mod(d(~l), 10);

    % Prepare figure
    fh = figure(1); clf; hold on; box on;

    dataset = cell(numel(I), 1);
    cc = lines(numel(I));
    for i = 1:numel(I)

        % Check for filtered data:
        f = ~isnan(radius(:, I(i)));

        % Find the first time the GFP signal decreases
        J = find(diff(GFP_radius(f, I(i)))<0);

        % Detmine time of deviation
        if  ~isempty(J)
            T_i(I(i)) = TT(J(1)) + dT(I(i));
        end

        % Plot curves
        plot(TT(f) + dT(I(i)), radius(f, I(i)),     '-', 'color', cc(i, :), 'LineWidth', 2, 'DisplayName', num2str(I(i)));
        plot(TT(f) + dT(I(i)), GFP_radius(f, I(i)), ':', 'color', cc(i, :), 'LineWidth', 2, 'DisplayName', num2str(I(i)));

        % Label graph
        xlabel('Time (hours)')
        ylabel('Radius ({\mu}m)')

        xlim([0 25])
        ylim([0 400])

        % Store data
        dataset{i} = {TT' + dT(I(i)), radius(:, I(i)), GFP_radius(:, I(i)), T_i(I(i))};
    end

    % Save graph
    saveas(fh, sprintf('../../experiments/Datasets/DataSet_%d.png', n))


    % Save fitting data
    save(sprintf('../../experiments/Datasets/DataSet_%d.mat', n), 'dataset')

end