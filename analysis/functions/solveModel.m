function [T, C, R, N, p] = solveModel(p)

% Compute grid properties: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the number of grid points
if p.k > 1
    % The grid uses dynamics length scales
    N = ceil(log(1-p.r_max*(1-p.k) / p.dr) / log(p.k))+1;
else
    % The grid uses linear length scales
    N = p.r_max/p.dr+1;
end

p.r = cumsum([0 p.k.^(0:N-2)*p.dr])';     % Location of gridpoints
p.rNp = p.r(end) + p.k^(N-2)*p.dr;          % Location of phantom point at N+1

% Define laplacian (on a potentially non-uniform grid)
dn0 = @(n, R) 6 * ( n(2) - n(1) ) / p.r(2).^2;
dni = @(n, R) 2 ./ p.r(2:end-1) .* ( ( 2*p.r(2:end-1)-p.r(3:end) ) .* n(1:end-2) ./ ( (p.r(1:end-2)-p.r(2:end-1)) .* (p.r(1:end-2)-p.r(3:end)) ) + ( 3*p.r(2:end-1)-p.r(1:end-2)-p.r(3:end)) .* n(2:end-1) ./ ( (p.r(2:end-1)-p.r(1:end-2)) .* (p.r(2:end-1)-p.r(3:end)) ) + ( 2*p.r(2:end-1)-p.r(1:end-2)) .* n(3:end) ./ ( (p.r(3:end)-p.r(2:end-1)) .* (p.r(3:end)-p.r(1:end-2)) ));
dnN = @(n, R) 2 / (p.r(end)-p.r(end-1)).^2 * ( n(end-1) - n(end) );

lap = @(t, y) [ dn0(y(1:N), y(N+1));
                dni(y(1:N), y(N+1));
                dnN(y(1:N), y(N+1))];

% Initialize nutrient density on the grid
n = ones(N, 1)*p.n_0;


% Define helper equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = @(y) y(N+1);                    % Colony Radius (Alive)

% Dynamics of the penetration depth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are a few models for the penetration depth dR
heaviside = @(t) (t >= p.T_i);                    % Heaviside function

% The order of definitions is important.
switch p.model
    case 0
        % Model 0; no phage pressure
        dR = @(t, y) 0;

    otherwise
        % Model 1: Instant phage pressure:
        dR = @(t, y) heaviside(t) * p.dR;
end


% Define more helper equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
monod    = @(t, y) p.g_max * y(1:N)./(y(1:N)+p.K);
core     = @(t, y) min(1, max(0, (R(y)-dR(t, y)-p.r)./([p.r(2:end); p.rNp]-p.r)));
% core     = @(t, y)     p.r<(B(y)-dR(t, R(y), y));
shell    = @(t, y) min(1, max(0, (R(y)-p.r)./([p.r(2:end); p.rNp]-p.r))) - core(t, y);
% shell    = @(t, y) and(p.r>(B(y)-dR(t, R(y), y)), p.r<B(y));
p.dV =  p.r.^2.*p.k.^(0:N-1)'*p.dr;

% Define problem equations
dndt    = @(t, y) p.D.*lap(t, y) - p.mu * monod(t, y) .* core(t, y);

growth  = @(t, y)             sum( monod(t, y) .* core(t, y)  .* p.dV );
phage   = @(t, y) p.epsilon * sum( monod(t, y) .* shell(t, y) .* p.dV );

dRdt    = @(t, y) p.mu / R(y).^2 * (growth(t, y) - phage(t, y));
dVdt    = @(t, y) 4 * pi * p.mu * growth(t, y);

% Assemble problem
dfdt    = @(t, y) [dndt(t, y); ...
                   dRdt(t, y);
                   dVdt(t, y)];

% Define starting conditions
y0 = [n; p.R0; 4 * pi / 3 * p.R0^3];

% options = odeset('NonNegative', ones(size(y0)), 'RelTol', 1e-2, 'AbsTol', 1e-4);
options = odeset('NonNegative', ones(size(y0)), 'RelTol', 1e-5, 'AbsTol', 1e-6);
% options = odeset('Refine', 4);
Y = [];
T = [];
T_prev = 0;
for t = 1:length(p.T_end)
    [tt, yy] = ode15s(dfdt, [T_prev p.T_end(t)], y0, options);
    T_prev = p.T_end(t);
    Y = [Y; yy];
    T = [T; tt];
    y0 = Y(end, :);
end
% [T, Y] = ode15s(dfdt, [0 p.T_end], y0, options);
V = Y(:, N+2);
C = (3 * V / (4 * pi)).^(1/3);
R = Y(:, N+1);
N = Y(:, 1:N);
