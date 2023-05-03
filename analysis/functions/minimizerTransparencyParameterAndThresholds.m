close all; clearvars;

if isempty(gcp('nocreate'))
    parpool('local', 32);
end

% Load the best fitting parameters (bacteria parameters)
load('../fits/GrowthParams.mat')

% load the data
load('../../experiments/data.mat')

% Get the number of colonies and number of time points
nT = size(radius, 1);
nC = size(radius, 2);

% Optimize over all runs
nd = 25;

% Store the invasion times
T_i = nan(nd, 1);

% Define thresholds
thresholds = unique([1:0.1:2.5 1.8:0.025:2.0]); % Adding extra fidelity around where \nu seems to be

% Prepare data sets
datasets  = cell(numel(thresholds), 1);
threshold = zeros(numel(thresholds), 1);

% Remove outliers
d = [11 13 31 45 55];

% Map the odd indicies
l = mod(floor(d/10), 2) == 1;
I = nan(size(d));
I(l) = floor(floor(d(l) / 10 ) / 2) * 5 + mod(d(l), 10);
% Map the even indicies
I(~l) = (5-(6-floor(d(~l)/10)) / 2)*5 + mod(d(~l), 10);


% Loop over data sets
for k = 1:numel(thresholds)

    % Loop over all thresholds
    eval(['GFP_radius = GFP_radius_' strrep(sprintf('%.2f', thresholds(k)), '.', '_')])

    % Replace NaN values
    GFP_radius(isnan(GFP_radius)) = 0;

    % Set threshold value
    threshold = thresholds(k);

    % Prepare data set
    datasets{k} = cell(nd, 1);
    for i = 1:nd

        % Check if point is outlier
        if any(i == I)
            continue
        end

        % Check if any GFP signal is detected
        if ~any(GFP_radius(:, i) > 0)
            continue
        end

        % Check for filtered data:
        f = ~isnan(radius(:, i));

        % Find the first time the GFP signal decreases
        J = find(diff(GFP_radius(f, i))<0);

        % Detmine time of deviation
        if  ~isempty(J)
            T_i(i) = TT(J(1)) + dT(i);
        else
            T_i(i) = TT(1) + dT(i) + 5;
        end

        % Store data
        datasets{k}{i} = {TT' + dT(i), radius(:, i), GFP_radius(:, i), T_i(i), threshold};

    end
end


for k = 1:numel(thresholds)

    % Create output folder
    sdir = '../fits/TransparencyParameter';
    if ~exist(sdir, 'dir')
        mkdir(sdir)
    end

    % Check for previous fit
    path = sprintf('%s/Params_%s.mat', sdir, strrep(sprintf('%.1f', thresholds(k)), '.', '_'));
    if exist(path, 'file')

        load(path, 'y', 'fitImprovement');
        y0  = y;

    else
        fitImprovement = 1;
        if k < numel(thresholds)
            tpath = sprintf('%s/Params_%s.mat', sdir, strrep(sprintf('%.1f', thresholds(k+1)), '.', '_'));
            if exist(tpath, 'file')
                load(tpath, 'y');
                y0  = y;
            end
        
        elseif k > 1
            tpath = sprintf('%s/Params_%s.mat', sdir, strrep(sprintf('%.1f', thresholds(k-1)), '.', '_'));
            if exist(tpath, 'file')
                load(tpath, 'y');
                y0  = y;
            end
        else
            y0  = [10^0.5  10  1];
        end
    end


    % Perform minimization
    if fitImprovement < 0.1
        continue
    end

    dataset = datasets{k};

    % Allocate cells for fitting data
    I = find(cellfun(@(x)~isempty(x), dataset));
    n = numel(I);
    Time = cell(1, n);
    BF   = cell(1, n);
    GFP  = cell(1, n);

    % Store the data into the cells
    for i = 1:n
        Time{i} = dataset{I(i)}{1};
        BF{i}   = dataset{I(i)}{2};
        GFP{i}  = dataset{I(i)}{3};
    end

    % Add the values of T_i to y0 if not already
    if length(y0) == 3
        y0 = [y0 T_i(I)'];
    end

    % Determine start fit value
    d_old = fitPhageAttack(x, y0, 2, Time, BF, GFP, false);

    % Call minimizier
    %options = optimset('MaxIter', 100*numel(y0));
    options = optimset('MaxIter', 400*numel(y0),'MaxFunEvals', 400*length(y0));
    y = fminsearch(@(y)fitPhageAttack(x, y, 2, Time, BF, GFP, true), y0, options);

    % Report change in conditions
    fprintf('Total parameter change: %.3f\n', norm(y-y0))

    d_new = fitPhageAttack(x, y, 2, Time, BF, GFP, false);
    fitImprovement = 100*(d_old-d_new)/d_old;
    fprintf('Fit improvement: %.3f %%\n', fitImprovement)

    path = sprintf('%s/Params_%s.mat', sdir, strrep(sprintf('%.1f', thresholds(k)), '.', '_'));

    save(path, 'y', 'fitImprovement');
end

