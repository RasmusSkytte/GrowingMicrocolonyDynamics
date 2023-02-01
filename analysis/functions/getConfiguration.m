function p = getConfiguration(x,varargin)

% Nutrient parameters
p.n_0     = x(1);       % Nutrient desity
p.D       = 10^x(2);    % Diffusion constant

% Cell parameters
p.K       = x(1)*min(1, x(4));  % Monod growth
p.mu      = 1;          % Packing density
p.R0      = (3/(4*pi)*1.4*1.3)^(1/3)/p.mu; % Initial colony size
p.g_max   = x(3);       % Maximal growth rate


% Phage parameters (Disabled when fitting growth curve)
if numel(varargin) > 0
    p.dR      = varargin{1}(1);  % Shell layer thickness
    p.epsilon = varargin{1}(2);  % Killing rate
    p.nu      = varargin{1}(3);  % dead density effect

    p.model   = varargin{2};
else
    p.dR      = 0;
    p.epsilon = 0;
    p.nu      = 0;

    p.model   = 0;
end

% Grid parameters
p.r_max = 2000;
p.dr = 0.05;
p.k = 1.02;
end
