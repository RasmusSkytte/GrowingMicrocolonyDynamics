close all; clearvars; clc;

% Load the fitted values
load('../fits/GrowthParams.mat');
load('../fits/PhageAttackParams.mat');

% Plot the reference and the pertubation
fh = figure(); clf;
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1;
ax.Box = 'on';


[t_ref, tt_ref, v_ref, nC_ref] = getVelocityOfExpansion(x);


x_param_labels = {'n_0', 'D', '{\mu}g_{max}', 'K'};
for x_ind = 1:numel(x)
    perm = zeros(size(x));
    perm(x_ind) = dx(x_ind);

    [t_p,   tt_p,   v_p,   nC_p]   = getVelocityOfExpansion(x + perm);
    [t_m,   tt_m,   v_m,   nC_m]   = getVelocityOfExpansion(x - perm);

    yyaxis left
    plot(ax, tt_ref, v_ref, 'LineWidth', 2, 'color', [0, 0.4470, 0.7410], Marker='none');
    hold on;
    plot(ax, tt_p,   v_p,   'LineWidth', 2, 'color', [0, 0.4470, 0.7410], 'LineStyle', ':', Marker='none');
    plot(ax, tt_m,   v_m,   'LineWidth', 2, 'color', [0, 0.4470, 0.7410], 'LineStyle', ':', Marker='none');
    ylabel('Colony growth rate ({\mu}m / h)', fontsize = 18)
    ax.XLim = [0 24];
    ax.YLim = [0 80];

    % Get configuration with phage
    p_ref = getConfiguration(x, y, 1);
    p_p   = getConfiguration(x + perm, y, 1);
    p_m   = getConfiguration(x - perm, y, 1);

    yyaxis right
    plot(ax, t_ref, p_ref.epsilon * p_ref.dR * p_ref.g_max * nC_ref ./ (nC_ref + p_ref.K), 'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980], Marker='none');
    plot(ax, t_p,   p_p.epsilon   * p_p.dR   * p_p.g_max   * nC_p   ./ (nC_p   + p_p.K),   'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980], Marker='none');
    plot(ax, t_m,   p_m.epsilon   * p_m.dR   * p_m.g_max   * nC_m   ./ (nC_m   + p_m.K),   'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980], Marker='none');
    ylabel('Phage killing rate ({\mu}m / h)', fontsize = 20)
    ax.YLim = [0 80];

    xlabel('Time (h)')

    title(x_param_labels{x_ind})

    if ~exist('../../figures/Figure S5', 'dir')
        mkdir('../../figures/Figure S5')
    end

    saveas(fh, sprintf('../../figures/Figure S5/FigS5_%d_.png', x_ind))
    hold off;
end


y_param_labels = {'{\Delta}R', '\epsilon'};
for y_ind = 1:2

    y_p = y; y_p(y_ind) = y(y_ind)*1.1;
    y_m = y; y_m(y_ind) = y(y_ind)*0.9;

    yyaxis left
    plot(ax, tt_ref, v_ref, 'LineWidth', 2, 'color', [0, 0.4470, 0.7410]);
    hold on;
    ylabel('Colony growth rate ({\mu}m / h)', fontsize = 18)
    ax.XLim = [0 24];
    ax.YLim = [0 80];

    % Get configuration with phage
    p_ref = getConfiguration(x, y, 1);
    p_p   = getConfiguration(x, y_p, 1);
    p_m   = getConfiguration(x, y_m, 1);

    yyaxis right
    plot(ax, t_ref, p_ref.epsilon * p_ref.dR * p_ref.g_max * nC_ref ./ (nC_ref + p_ref.K), 'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980], Marker='none');
    plot(ax, t_ref, p_p.epsilon   * p_p.dR   * p_ref.g_max * nC_ref ./ (nC_ref + p_ref.K), 'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980], Marker='none', LineStyle=':');
    plot(ax, t_ref, p_m.epsilon   * p_m.dR   * p_ref.g_max * nC_ref ./ (nC_ref + p_ref.K), 'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980], Marker='none', LineStyle=':');
    ylabel('Phage killing rate ({\mu}m / h)', fontsize = 20)
    ax.YLim = [0 80];

    xlabel('Time (h)')

    title(y_param_labels{y_ind})

    if ~exist('../../figures/Figure S5', 'dir')
        mkdir('../../figures/Figure S5')
    end

    saveas(fh, sprintf('../../figures/Figure S5/FigS5_%d_.png', y_ind+numel(x)))
    hold off;
end




% Define function that computes velocity of expansion
function [t, tt, v, nC] = getVelocityOfExpansion(x)

    % Get configration struct
    p = getConfiguration(x);

    % Load the fitting data
    p.T_end = 24;

    % Solve model with current parameters
    [t, C, ~, N, p] = solveModel(p);

    % Compute the velocity of expantion
    v = C;
    tt = t;
    v(diff(C) == 0)  = [];
    tt(diff(C) == 0) = [];

    v = diff(v) ./ diff(tt);
    tt = tt(1:(end-1)) + diff(tt) / 2;

    % Determine the surface nutrient level
    f = @(C)min(1,max(0, (C-p.r)./([p.r(2:end); p.rNp]-p.r)));
    nC = nan(size(t));
    for i = 1:numel(nC)
        fC = f(C(i));
        I = find(fC, 1, 'last');
        nC(i) = N(i, I-1) * (1-fC(I)) + N(i, I) * fC(I);
    end
end

