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
sdir='../../figures/Figure 6'
if ~exist(sdir, 'dir')
    mkdir(sdir)
end


% Get configration struct
p = getConfiguration(x, y, 1);

p.T_end = linspace(0.5, 24, 50) % X-axis should linearly cover between 0 and 24 hours
T_i = linspace(1, 15,20) % Z-axis should linearly cover attack timings from 0 to 15(?)
                         % hours. I use 10 steps so it should not be too slow. Can be increased.

% Loop over attack timings and solve model
X = p.T_end;
for i = 1: 20;
    % Set the invasion time
    p.T_i = T_i(i)

    % Solve model with current parameters
    [t, C, R] = solveModel(p);

    % Extract the alive radius at the sample times
    rr = [];

    j = 1;
    for ii = 1:length(p.T_end);
        iflag = 0;
        while (iflag == 0);
            if (t(j) < p.T_end(ii));
                j = j + 1;
            else
                if(t(j)==p.T_end(ii));
                    rr = [rr, R(j)];
                    iflag = 1;
                else
                    print('problem')
                end
            end
        end
    end

    if(i==1)
        Z=rr
        Y=p.T_i
    else
        Z = [Z; rr]
        Y =[Y p.T_i]
    end
end

% Construct the grid to solve on
%[X, Y] = meshgrid(p.T_end, T_i)

% Plot surface
fh = figure(); clf;
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1;
ax.Box = 'on'

surf(X, Y, Z, "Linestyle", "none")
xlabel('Time (h)')
ylabel('T_i (h)')
zlabel('Alive radius ({\mu}m)')

saveas(fh, sprintf('%s/Fig6.png', sdir));
