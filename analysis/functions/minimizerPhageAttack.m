close all; clearvars;

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
                ub = [20       50.0     0        0        0  ];
                
                %T_i = [8.02, 8.68, 8.62, 9.23, 10.45, 10.50, 11.36];
                
            case 2
                %    dR        epsilon  tau      gamma    nu
                lb = [1        0.05     0        0        0.01  ];
                y0 = [4        12.5     0        0        1     ];
                ub = [20       50.0     0        0        10    ];
                
                %T_i = [8.02, 8.68, 8.62, 9.23, 10.45, 10.50, 11.36];
                
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
        path = sprintf('%s/Params.mat', sdir);
        if exist(path, 'file')
            load(path, 'y');
            y0 = y(1:5);
            T_i = y(6:end);
        end
        
        % Add T_i to fitting parameters
        y0 = [y0 T_i];
        
        % Determine start fit value
        y_old = y0;
        d_old = fitPhageAttack(x, y_old, model, Time, BF, GFP, detectionThreshold, false);
         
        % FMINCON OR GENETIC ALGORITHM
        lb = [lb zeros(size(T_i))];
        ub = [ub 24*ones(size(T_i))];
        
        % Assert bounds are set
        y0(y0 < lb) = lb(y0 < lb);
        y0(y0 > ub) = ub(y0 > ub);
                
        % Call minimizier
        % FMINCON
%         options = optimoptions('fmincon', 'UseParallel', true);
%         y = fmincon(@(y)fitPhageAttack(x, y, model, Time, BF, GFP, detectionThreshold, true), y0, [], [], [], [], lb, ub, [], options);

        % GENETIC ALGORITHM
%         options = optimoptions('ga', 'UseParallel', true);
%         y = ga(@(y)fitPhageAttack(x, y, model, Time, BF, GFP, true), numel(y0), [], [], [], [], lb, ub, [], options);
        
        % FMINSEARCH
        options = optimset('MaxIter', 200*numel(y0));
        y = fminsearch(@(y)fitPhageAttack(x, y, model, Time, BF, GFP, detectionThreshold, true), y0, options);
        
        % Report change in conditions
        fprintf('Total parameter change: %.3f\n', norm(y-y_old))
        
        d_new = fitPhageAttack(x, y, model, Time, BF, GFP, detectionThreshold, false);
        fprintf('Fit improvement: %.3f %%\n', 100*(d_old-d_new)/d_old)
        
        % Store fit
        d = d_new;
        save(path, 'y', 'd');
        
    end
end