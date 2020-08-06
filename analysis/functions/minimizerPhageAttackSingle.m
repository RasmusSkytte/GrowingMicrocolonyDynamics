close all; clearvars;
pause(60)
parpool('local', 16);

% Minimum detectable GFP size
detectionThreshold = 0;

% Choose dataset to fit
for DataSetID = 4
    
    fprintf('Analyzing DataSetID %d\n', DataSetID)
    
    % Select model
    for model = 3
        
        switch model
            case 1
                %    dR        epsilon  tau      gamma    nu
                lb = [1        0.05     0        0        0  ];
                y0 = [4        12.5     0        0        0  ];
                ub = [50       100.0    0        0        0  ];
                
                %T_i = [8.02, 8.68, 8.62, 9.23, 10.45, 10.50, 11.36];
                
            case 2
                %    dR        epsilon  tau      gamma    nu
                lb = [1        0.05     0        0        0.01  ];
                y0 = [4        12.5     0        0        1     ];
                ub = [20       50.0     0        0        10    ];
                
                %T_i = [7.69, 8.15, 8.28, 9.27, 9.66, 10.11, 10.99];
                
            case 3
                %    dR        epsilon  tau      gamma    nu
                lb = [1        0.05     0.1      0.01     0.1 ];
                y0 = [5.33     1.18     1        5e4      1    ];
                ub = [20       50.0     1        10       10.0 ];
                
                %T_i = [7.69, 8.15, 8.28, 9.27, 9.66, 10.11, 10.99];
        end
        
        % Load the best fitting parameters (bacteria parameters)
        load('../fits/GrowthParams_Model_2.mat')
        
        % Load the fitting data
        load(sprintf('../../experiments/Datasets/DataSet_%d', DataSetID))
        
        % Allocate cells for fitting data
        nd = numel(dataset);
        Time = cell(1, nd);
        BF   = cell(1, nd);
        GFP  = cell(1, nd);
        
        T_i  = nan(1, nd);
        
        % Store the data into the cells
        parfor i = 1:nd
            
            % Adjust GFP level
            
            % Scale turning point value
            k = find(diff(dataset{i}{3}) < 0, 1);
            k = min(dataset{i}{2}(1:k)./dataset{i}{3}(1:k));
            
            if isempty(k)
                k = 1;
            end
            
            % Scale inital value
            % k = dataset{1}{2}(find(I, 1)) / dataset{1}{3}(find(I, 1));
            
            % No scale
            if model == 3
                k = 1;
            end
            
            % Get estimate for T_i
            if ~isnan(dataset{i}{4})
                T_i(i) = dataset{i}{4};
            else
                T_i(i) = dataset{i}{1}(1)+5;
            end
            
            Time{i} = dataset{i}{1};
            BF{i}   = dataset{i}{2};
            GFP{i}  = dataset{i}{3} * k;
        end
        
        % Create output folder
        sdir = sprintf('../fits/PhageAttack_Model_%d_DataSetID_%d', model, DataSetID);
        if ~exist(sdir, 'dir')
            mkdir(sdir)
        end

        % Check for previous fit
        path = sprintf('%s/Params_single.mat', sdir);
        if exist(path, 'file')
            load(path, 'yy');
        else
            y = [y0 T_i];
            
            % Allocate fitting values
            yy = [repmat(y(1:5)', 1, nd); T_i];
        
        end
        
        % Determine start fit value
        yy_old = yy;
        
        d_old = 0;
        parfor i = 1:nd
            d_old = d_old + fitPhageAttackSingle(x,  yy(:, i), model, Time{i}, BF{i}, GFP{i}, detectionThreshold, false);
        end
        fprintf('Initial fit distance: %.3f\n', d_old)
        
        % FMINCON OR GENETIC ALGORITHM
        lb = [lb 0]';
        ub = [ub 24]';
        
        parfor i = 1:nd
            y0 = yy(:, i);
        
            % Assert bounds are set
            y0(y0 < lb) = lb(y0 < lb);
            y0(y0 > ub) = ub(y0 > ub);

            % Call minimizier
            % FMINCON
%             options = optimoptions('fmincon', 'UseParallel', true);
%             yy(:, i) = fmincon(@(y)fitPhageAttackSingle(x, y, model, Time{i}, BF{i}, GFP{i}, detectionThreshold, true), y0, [], [], [], [], lb, ub, [], options);
            options = optimset('MaxIter', 200*numel(y0));
            yy(:, i) = fminsearch(@(y)fitPhageAttackSingle(x, y, model, Time{i}, BF{i}, GFP{i}, detectionThreshold, true), y0, options);
%               yy(:, i) = ga(@(y)fitPhageAttackSingle(x, y, model, Time{i}, BF{i}, GFP{i}, detectionThreshold, true), numel(y0), [], [], [], [], lb, ub);
        end

        % Convert to the right format
        y = [mean(yy(1:5, :), 2)' yy(end, :)];

        % Report change in conditions
        fprintf('Total parameter change: %.3f\n', norm(yy(:)-yy_old(:)))
        
        d_new = 0;
        parfor i = 1:nd
            d_new = d_new + fitPhageAttackSingle(x,  yy(:, i), model, Time{i}, BF{i}, GFP{i}, detectionThreshold, false);
        end
        
        fprintf('Fit distance: %.3f\n', d_new)
        fprintf('Fit improvement: %.3f %%\n', 100*(d_old-d_new)/d_old)
        
        % Store fit
        d = d_new;
        save(path, 'y', 'yy', 'd');
        
    end
end