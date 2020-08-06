function dist = fitPhageAttackCoupled(x, y, model, Time, BF, GFP, detectionThreshold, reportDistance)

% Define dist
dist = 0;

% Loop over data
parfor n = 1:size(Time, 2)
    
    % Construct the configration
    p = getConfiguration(x, [y(4+(2*n-1)) y(1:4)], model);

    % Set the invasion time
    p.T_i = y(4+2*n);

    % Simulation sample times
    p.T_end = Time{n};

    % Solve model with current parameters
    [t, C, R] = solveModel(p);

    % Filter out the sample points
    I = ismember(t, Time{n});
    I = find(I);
    I(find(diff(I)==1)+1) = [];

    if model ~= 3
        R = R(I);
        C = C(I);
    end
    
    if model == 2
        % Volume of dead matter
        % Vd = p.nu * 4 * pi / 3 * (C^3 - R^3)
        
        % Volume of alive matter
        % Va = 4 * pi / 3 R^3
        
        % Colony volume
        % Vc = Va + Vd = 4 * pi / 3 * (p.nu * (C^3 - R^3) + R^3)
        
        % Colony radius
        C = (p.nu * (C.^3 - R.^3) + R.^3).^(1/3);
    end
    
    if model == 3
        
        % Choose C as normally                
        C = C(I);
        
        % GFP has maturation time
        % Choose R closest to maturation time
        I = arrayfun(@(T)find(min(abs(t-(T-p.nu))) == abs(t-(T-p.nu)), 1), Time{n});
        R = R(I);
    end
 
    % Compute loss stemming from BF Radius
    dist = dist + sum( ( C - BF{n} ).^2 );

    % Compute loss stemming from GFP Radius
    dist = dist + sum( ( R(GFP{n} > detectionThreshold) - GFP{n}(GFP{n} > detectionThreshold) ).^2 );

end

if reportDistance
    fprintf('\tfitting distance = %.1f, y = [%s]\n', dist,  sprintf([repmat('%.1f, ', 1, numel(y)-1) '%.1f'], y) )
end

end
