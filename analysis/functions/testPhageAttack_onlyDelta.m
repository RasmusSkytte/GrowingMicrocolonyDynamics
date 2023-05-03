close all; clearvars;

% Load the fitting data
load('../../experiments/PhageData.mat')
nDataSets = numel(dataset);
labels = {'a', 'b', 'c', 'd', 'e'};

% Load growth curve fit
load('../fits/GrowthParams.mat')

% Load the phage attack fit
load('../fits/PhageAttackParams_onlyDelta.mat')

% Prepare folders
if ~exist('../../figures/Figure 5', 'dir')
    mkdir('../../figures/Figure 5')
end

if ~exist('../../figures/Figure S4', 'dir')
    mkdir('../../figures/Figure S4')
end

% Get configration struct
p = getConfiguration(x, y, 3);

% Loop over data sets (sorted by T_i)
[T_i, I] = sort(round(y(4:end), 2));

[~, I7_5]  = min(abs(y(4:end) - 7.5));
[~, I8_0]  = min(abs(y(4:end) - 8.0));
[~, I8_5]  = min(abs(y(4:end) - 8.5));
[~, I9_0]  = min(abs(y(4:end) - 9.0));
[~, I9_5]  = min(abs(y(4:end) - 9.5));

for i = 1:numel(I)
    n = I(i);

    % Prepare figure
    fh1 = figure(1); clf; hold on; box on;
    ax1 = gca;
    ax1.FontSize = 18;
    ax1.Units = 'pixels';
    ax1.LineWidth = 1;
    ax1.TickLength = [0.02 0.05];
    ax1.XLim = [0 25];
    ax1.YLim = [0 400];

    % Determine size
    if mod(i-1, 4) == 0 % Left edge
        dx = 95;
        ylabel(ax1, 'Colony Radius ({\mu}m)')
        ax1.YLabel.Position(1) = -5.5;
    else
        dx = 0;
    end

    if ceil(i/4) == 5 % Bottom edge
        dy = 75;
        xlabel(ax1, 'Time (h)')
    else
        dy = 0;
    end

    % Set size
    ax1.Position = [10+dx 15+dy 290 315];
    fh1.Position(3:4) = [325+dx 345+dy];

    % Create aliases
    Time = dataset{n, 1}; % Time Data
    BF   = dataset{n, 2}; % BF data
    GFP  = dataset{n, 3}; % GDP Data

    % Set the invasion time
    p.T_i = y(3+n);

    % Simulation sample times
    p.T_end = Time;

    % Solve model with current parameters
    [t, C, R] = solveModel(p);

    % Plot individual data and fit
    plot(ax1, Time, BF,       'Color', 'k', 'LineWidth', 2);
    plot(ax1, t,    C,  '--', 'Color', 'k', 'LineWidth', 2);

    plot(ax1, Time, GFP,      'Color', 'b', 'LineWidth', 2);
    plot(ax1, t,    R,  '--', 'Color', 'b', 'LineWidth', 2);

    text(ax1, 2.5, 350, sprintf('T_i = %.2f h', p.T_i), 'FontSize', 20);

    % Save figure
    saveas(fh1, sprintf('../../figures/Figure S4/Run_%d.png', i))


    % Add plots to the main figure
    if any(I(i) == [I7_5 I8_0 I8_5 I9_0 I9_5 I7_5])

        fh3 = figure(3); clf; hold on; box on;
        ax3 = gca;
        ax3.FontSize = 18;
        ax3.Units = 'pixels';
        ax3.LineWidth = 1;
        ax3.TickLength = [0.02 0.05];
        ax3.XLim = [0 25];
        ax3.YLim = [0 400];

        % Determine size
        if I(i) == I7_5 % Left edge
            dx = 95;
            ylabel(ax3, 'Colony Radius ({\mu}m)')
            ax3.YLabel.Position(1) = -7.5;
        else
            dx = 0;
        end

        dy = 75;
        xlabel(ax3, 'Time (h)')

        % Set size
        ax3.Position = [10+dx 15+dy 200 315];
        fh3.Position(3:4) = [235+dx 345+dy];

        % Plot together
        plot(ax3, Time, BF,       'Color', 'k', 'LineWidth', 2);
        plot(ax3, t,    C,  '--', 'Color', 'k', 'LineWidth', 2);

        plot(ax3, Time, GFP,      'Color', 'b', 'LineWidth', 2);
        plot(ax3, t,    R,  '--', 'Color', 'b', 'LineWidth', 2);

        text(ax3, 2.5, 350, sprintf('T_i = %.2f h', p.T_i), 'FontSize', 20);

        saveas(fh3, sprintf('../../figures/Figure 5/Fig5%s.png', labels{I(i) == [I7_5 I8_0 I8_5 I9_0 I9_5]}))
    end

    drawnow;

end
