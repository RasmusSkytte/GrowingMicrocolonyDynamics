function p = getConfiguration(x,varargin)

% Nutrient parameters
p.n_0     = x(1);       % Nutrient desity
p.D       = 10^x(2);    % Diffusion constant

% Cell parameters
p.g_max   = x(3);       % Maximal growth rate
p.K       = x(1)*x(4);  % Monod growth
p.delta   = x(5);       % Decay rate
p.mu      = x(6);       % Packing density
p.R0      = (3/(4*pi)*x(7)*1.3)^(1/3)/p.mu; % Initial colony size

% Phage parameters (Disabled when fitting growth curve)
if  numel(varargin) > 0
    p.dR      = varargin{1}(1);  % Shell layer thickness
    p.epsilon = varargin{1}(2);  % Killing rate
    p.tau     = varargin{1}(3);  % Establishing time
    p.gamma   = varargin{1}(4);  % Adsorption effect
    p.nu      = varargin{1}(5);  % dead density effect
    
    p.model   = varargin{2};
else
    p.dR      = 0;
    p.epsilon = 0;
    p.tau     = 0;
    p.gamma   = 0;
    p.nu      = 0;
    
    p.model   = 0;
end

% Grid parameters
p.r_max = 2000;
% p.dr = 0.1; p.k = 1.05; % chi_2 / ndf ~ 8.7
% p.dr = 0.1; p.k = 1.03; % chi_2 / ndf ~ 10
% p.dr = 0.1; p.k = 1.02; % chi_2 / ndf ~ 2.5
% p.dr = 0.05; p.k = 1.05; % chi_2 / ndf ~ 7.9
% p.dr = 0.05; p.k = 1.03; % chi_2 / ndf ~ 1.7
p.dr = 0.05; p.k = 1.02; % chi_2 / ndf ~ 0.4
% p.dr = 0.01; p.k = 1.05; % chi_2 / ndf ~ 0.74
% p.dr = 0.01; p.k = 1.03; % chi_2 / ndf ~ 1.1
end