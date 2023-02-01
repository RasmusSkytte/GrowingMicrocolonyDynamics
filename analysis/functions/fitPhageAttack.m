function distance = fitPhageAttack(x, y, model, Time, BF, GFP, reportDistance)

% Define dist
distances = nan(1, numel(Time));

% Loop over data
parfor n = 1:numel(Time)

    % Construct the configration
    p = getConfiguration(x, y, model);

    % Set the invasion time
    p.T_i = y(3+n);

    % Simulation sample times
    p.T_end = Time{n};

    % Solve model with current parameters
    [t, C, R] = solveModel(p);

    % Filter out the sample points
    I = ismember(t, Time{n});
    I = find(I);
    I(find(diff(I)==1)+1) = [];

    R = R(I);
    C = C(I);

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

    % Compute loss stemming from BF Radius
    I = ~isnan(BF{n}) & BF{n} > 0;
    dist = sum( ( C(I) - BF{n}(I) ).^2 )

    % Compute loss stemming from GFP Radius
    I = ~isnan(GFP{n}) & GFP{n} > 0;
    dist = dist + sum( ( R(I) - GFP{n}(I) ).^2 );

    if isnan(dist)
        error('Nans produced in fitPhageAttack')
    end

    distances(n) = dist

end


distance = sum(distances);

if reportDistance
    fprintf('\tfitting distance = %10.1f, y = [%s]\n', distance,  sprintf([repmat('%.3f, ', 1, numel(y)-1) '%.3f'], y) )
end

end
