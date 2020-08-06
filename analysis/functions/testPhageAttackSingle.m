close all; clearvars;

% Choose data and model
DataSetID = 4;
model = 2;

detectionThreshold = 0;

% Load fits
load('../fits/GrowthParams_Model_2.mat')
sdir = sprintf('../fits/PhageAttack_Model_%d_DataSetID_%d', model, DataSetID);
load(sprintf('%s/Params_single.mat', sdir));

% Load the fitting data
load(sprintf('../../experiments/Datasets/DataSet_%d', DataSetID))

% Loop over data sets
nDataSets = numel(dataset);

if ~exist(sprintf('%s/Single', sdir), 'dir')
    mkdir(sprintf('%s/Single', sdir))
end

fh1 = figure(1); clf; hold on; box on;
ax1 = gca;
ax1.FontSize = 20;

fh2 = figure(2); clf; hold on; box on;
ax2 = gca;
ax2.FontSize = 20;

fh3 = figure(3); clf; hold on; box on;
ax3 = gca;
ax3.FontSize = 20;
s = scatter(ax3, [], [], 'filled');

cc = lines(nDataSets);

% Add labels to figure 2
xlabel(ax2, 'Time')
ylabel(ax2, 'Colony Radius')


dist = 0;
for n = 1:numel(dataset)
    
    % Get configration struct
    p = getConfiguration(x, yy(:, n), model);

    % Create aliases
    Time = dataset{n}{1};
    BF   = dataset{n}{2};
    GFP  = dataset{n}{3};
    
    % Adjust GFP level
    
    % Scale turning point value
    ind = find(diff(GFP) < 0, 1);
    k = min(BF(1:ind)./GFP(1:ind));
    
    if isempty(ind)
        ind = 1;
        k = 1;
    end

    % Scale inital value
    % k = dataset{1}{2}(find(I, 1)) / dataset{1}{3}(find(I, 1));
    
    % No scale
    % k = 1;
    
    % Apply scale
    GFP = GFP * k;
    
    % Set the invasion time
    p.T_i = y(5+n);
    
    % Simulation sample times
    p.T_end = Time;
    
    % Solve model with current parameters
    [t, C, R] = solveModel(p);
    
    % Filter out the sample points
    I = ismember(t, Time);
    I = find(I);
    I(find(diff(I)==1)+1) = [];
    
    % Plot individual data and fit
    delete(ax2.Children)
    plot(ax2, Time, BF,        'Color', 'k', 'LineWidth', 2);
    plot(ax2, t,    C,   '--', 'Color', 'k', 'LineWidth', 2);
    
    plot(ax2, Time, GFP,       'Color', 'b', 'LineWidth', 2);
    plot(ax2, t,    R,   '--', 'Color', 'b', 'LineWidth', 2);
    
    ax2.XLim = [0 25];
    ax2.YLim = [0 400];
    
    saveas(fh2, sprintf('%s/Single/Run_%d.png', sdir, n))
    
    % Plot a few together
    if any( n == (1:floor(nDataSets/4):nDataSets))
        plot(ax1, Time, BF,        'Color', cc(n, :), 'LineWidth', 2);
        plot(ax1, t,    C,   '--', 'Color', cc(n, :), 'LineWidth', 2);
    
        plot(ax1, Time, GFP, '-.', 'Color', cc(n, :), 'LineWidth', 2);
        plot(ax1, t,    R,   '--', 'Color', cc(n, :), 'LineWidth', 2);
    end
    
    % Plot the scaling parameter
    s.XData = [s.XData k];
    s.YData = [s.YData ind];
    
    drawnow;
    
    % Keep track of the loss function
    R = R(I);
    C = C(I);
 
    % Compute loss stemming from BF Radius
    dist = dist + sum( ( C - BF ).^2 );

    % Compute loss stemming from GFP Radius
    dist = dist + sum( ( R(GFP > detectionThreshold) - GFP(GFP > detectionThreshold) ).^2 );


end


%saveas(fh, sprintf('data/FitPhageAttack/Model_%d.fig', n))
%saveas(fh, sprintf('data/FitPhageAttack/Model_%d.png', n))

fprintf('\tfitting distance = %.3f, y = [%s]\n', dist,  sprintf([repmat('%.5f, ', 1, numel(y)-1) '%.5f'], y) )