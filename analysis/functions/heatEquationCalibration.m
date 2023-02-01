close all; 
clearvars;

% Define parameters
T0 = 200;
T_end = 10;
D = 1/100;

% Define the grid
p.dr = 0.01;
p.r = 0:p.dr:1;

N = length(p.r);

q = 0.25;

% Define exact solution
l_n = (1:1000)'*pi;
U = @(r, t) arrayfun(@(r)2*T0./r .* sum((sin(l_n)./l_n.^2 - cos(l_n)./l_n).*sin(l_n.*r).*exp(-l_n.^2*t)), r);

% Plot exact solution
plot(p.r, U(p.r, T_end/100), '-', 'LineWidth', 6, 'DisplayName', 'Exact')

% Run IEEE solution
n_ieee  = ones(N, 1)*T0; n_ieee(end) = 0;

i = 1:N-2;
dt = 0.0025;
for t = 0:dt:T_end

    % IEEE
    n_p = n_ieee(3:end);
    n_m = n_ieee(1:end-2);
    
    n0 = n_ieee(1) + 6*q.*(n_ieee(2)   - n_ieee(1));
    
    n_ieee  =(1-2*q)*n_ieee(2:end-1) + q*((1 + 1./i').*n_p + (1 - 1./i').*n_m);
    
    n_ieee = [n0; n_ieee; 0];
    
end

figure(1); hold on;
plot(p.r, n_ieee, '-', 'LineWidth', 4, 'DisplayName', 'IEEE')


% Run our solution
p.r_max = 1;
p.dr    = 0.01;
k = 1.05;                                      % Grid scaling factor
if k > 1
    N = ceil(log(1-p.r_max*(1-k)/p.dr)/log(k))+1; % Number of grid points
else
    N = p.r_max/p.dr+1;
end
p.r = cumsum([0 k.^(0:N-2)*p.dr])';     % Location of gridpoints
p.r(end) = p.r_max;

n_0     = ones(N, 1)*T0; n_0(end) = 0;


% Define laplacian (on a potentially non-uniform grid)
dn0 = @(n) 6 * ( n(2) - n(1) ) / p.r(2).^2;
dni = @(n) 2 ./ p.r(2:end-1) .* ( ( 2*p.r(2:end-1)-p.r(3:end) ) .* n(1:end-2) ./ ( (p.r(1:end-2)-p.r(2:end-1)) .* (p.r(1:end-2)-p.r(3:end)) ) + ( 3*p.r(2:end-1)-p.r(1:end-2)-p.r(3:end)) .* n(2:end-1) ./ ( (p.r(2:end-1)-p.r(1:end-2)) .* (p.r(2:end-1)-p.r(3:end)) ) + ( 2*p.r(2:end-1)-p.r(1:end-2)) .* n(3:end) ./ ( (p.r(3:end)-p.r(2:end-1)) .* (p.r(3:end)-p.r(1:end-2)) ));
dnN = @(n) 0;%2 / (p.r(end)-p.r(end-1)).^2 * ( n(end-1) - n(end) );

lap = @(t, n) [  dn0(n);
                dni(n);
                dnN(n)];
            
options = odeset('NonNegative', 1, 'Refine', 1, 'RelTol', 1e-8, 'AbsTol', 1e-10);
[T, Y] = ode15s(@(t, y)D.*lap(t, y), [0 T_end], n_0, options);
n_1 = Y(end, :)';

plot(p.r, n_1, '-', 'LineWidth', 2, 'DisplayName', 'Our Method')
legend('Location', 'Best')