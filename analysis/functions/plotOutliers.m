close all; clearvars;

% Load the fitting data
load('../../experiments/OutlierData.mat')
nDataSets = numel(dataset);
    
% Prepare folders
if ~exist('../../figures/Figure S3', 'dir')
    mkdir('../../figures/Figure S3')
end

for n = 1:nDataSets
    
    % Prepare figure
    fh1 = figure(1); clf; hold on; box on;
    ax1 = gca;
    ax1.FontSize = 20;
    ax1.Units = 'pixels';
    ax1.LineWidth = 1;
    ax1.TickLength = [0.02 0.05];
    ax1.XLim = [0 25];
    ax1.YLim = [0 400];
    
    % Determine size
    if mod(n-1, 3) == 0 % Left edge
        dx = 95;
        ylabel(ax1, 'Colony Radius ({\mu}m)')
        ax1.YLabel.Position(1) = -5.5;
    else
        dx = 0;
    end
    
    if ceil(n/3) == 2 % Bottom edge
        dy = 75;
        xlabel(ax1, 'Time')
    else
        dy = 0;
    end
    
    % Set size
    ax1.Position = [10+dx 15+dy 290 315];
    fh1.Position(3:4) = [325+dx 345+dy];
    
    % Create aliases
    Time = dataset{n}{1};
    BF   = dataset{n}{2};
    GFP  = dataset{n}{3};
    
    % Plot individual data
    plot(ax1, Time, BF,       'Color', 'k', 'LineWidth', 2);
    plot(ax1, Time, GFP,      'Color', 'b', 'LineWidth', 2);
    
    drawnow;
    
    % Save figure
    saveas(fh1, sprintf('../../figures/Figure S3/Run_%d.png', n))
    
end
