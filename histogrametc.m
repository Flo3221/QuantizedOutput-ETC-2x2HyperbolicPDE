% Histogram of ETC inter-event times (tau_{k+1} - tau_k)
% Figure 1 : mu0   = 2  (red) vs mu0   = 10 (blue), Omega fixed = 0.6
% Figure 2 : Omega = 0.1(red) vs Omega = 0.6(blue), mu0   fixed = 8
% Computed over 100 initial conditions, T_sim = 15 s
n_IC  = 100;
T_sim = 15;
alphas = linspace(0.2, 20, n_IC);

iet_mu2_all     = [];
iet_mu10_all    = [];
iet_Omega01_all = [];
iet_Omega06_all = [];

% Model workspace handle — needed to override Omega (lives there, not in base)
mws = get_param('q2x2sim', 'ModelWorkspace');

for ic = 1:n_IC
    fprintf('IC %d/%d\n', ic, n_IC);

    v0_ic = @(x) alphas(ic) * (1 - x);
    z0    = [p * v0_ic(x1); v0_ic(x0)];
    z0hat = 1.5 * z0;

    assignin('base', 'z0',    z0);
    assignin('base', 'z0hat', z0hat);

    % ── Figure 1: vary mu0, fix Omega = 0.6 ──────────────────────────────
    mws.assignin('Omega', 0.6);          % Omega -> model workspace

    mu0 = 2;
    assignin('base', 'mu0', mu0);
    out  = sim('q2x2sim.slx', 'StopTime', num2str(T_sim));
    tetc = out.Uetc.Time(:);  Uetc = out.Uetc.Data(:);
    Utg  = find(diff(Uetc) ~= 0);
    if ~isempty(Utg)
        t_ev = tetc(Utg+1);
        iet_mu2_all = [iet_mu2_all; diff([0; t_ev])];
    end

    mu0 = 10;
    assignin('base', 'mu0', mu0);
    out  = sim('q2x2sim.slx', 'StopTime', num2str(T_sim));
    tetc = out.Uetc.Time(:);  Uetc = out.Uetc.Data(:);
    Utg  = find(diff(Uetc) ~= 0);
    if ~isempty(Utg)
        t_ev = tetc(Utg+1);
        iet_mu10_all = [iet_mu10_all; diff([0; t_ev])];
    end

    % ── Figure 2: vary Omega, fix mu0 = 8 ────────────────────────────────
    assignin('base', 'mu0', 8);          % mu0 -> base workspace

    mws.assignin('Omega', 0.1);
    out  = sim('q2x2sim.slx', 'StopTime', num2str(T_sim));
    tetc = out.Uetc.Time(:);  Uetc = out.Uetc.Data(:);
    Utg  = find(diff(Uetc) ~= 0);
    if ~isempty(Utg)
        t_ev = tetc(Utg+1);
        iet_Omega01_all = [iet_Omega01_all; diff([0; t_ev])];
    end

    mws.assignin('Omega', 0.6);
    out  = sim('q2x2sim.slx', 'StopTime', num2str(T_sim));
    tetc = out.Uetc.Time(:);  Uetc = out.Uetc.Data(:);
    Utg  = find(diff(Uetc) ~= 0);
    if ~isempty(Utg)
        t_ev = tetc(Utg+1);
        iet_Omega06_all = [iet_Omega06_all; diff([0; t_ev])];
    end
end

% Restore default values
mws.assignin('Omega', 0.6);
assignin('base', 'mu0', 8);

fprintf('Done.\n  mu0=2    : %d IETs\n  mu0=10   : %d IETs\n  Omega=0.1: %d IETs\n  Omega=0.6: %d IETs\n', ...
    numel(iet_mu2_all), numel(iet_mu10_all), ...
    numel(iet_Omega01_all), numel(iet_Omega06_all));

%% ── Clean: remove zero/negative artefacts ────────────────────────────────
iet_mu2_all     = iet_mu2_all(iet_mu2_all       > 0);
iet_mu10_all    = iet_mu10_all(iet_mu10_all     > 0);
iet_Omega01_all = iet_Omega01_all(iet_Omega01_all > 0);
iet_Omega06_all = iet_Omega06_all(iet_Omega06_all > 0);

%% ── Shared plot settings ─────────────────────────────────────────────────
red  = [1,   0,   0  ];
blue = [0.1, 0.1, 0.8];
nbins = 100;

% Helper: compute log edges and ticks from a pair of IET vectors
make_axes = @(v1, v2) deal( ...
    floor(log10(min([v1;v2]))*10)/10, ...
    ceil( log10(max([v1;v2]))*10)/10 );

%% ── Figure 1 : mu0 comparison ────────────────────────────────────────────
[e_min, e_max] = make_axes(iet_mu2_all, iet_mu10_all);
edges      = logspace(e_min, e_max, nbins);
xtick_exp  = e_min:0.2:e_max;
xtick_vals = 10.^xtick_exp;
xtick_lbls = arrayfun(@(v) sprintf('$10^{%.1f}$',v), xtick_exp, 'UniformOutput', false);

figure('Color', 'w');
ax1 = gca;  hold(ax1, 'on');

h2 = histogram(iet_mu10_all, edges, ...   % blue — background
    'Normalization', 'probability', 'FaceColor', blue, 'FaceAlpha', 0.95, 'EdgeColor', 'none');
h1 = histogram(iet_mu2_all,  edges, ...   % red  — foreground
    'Normalization', 'probability', 'FaceColor', red,  'FaceAlpha', 0.95, 'EdgeColor', 'none');

ax1.XScale = 'log';  ax1.FontSize = 30;  ax1.LineWidth = 2;  ax1.Box = 'on';
ax1.XTick = xtick_vals;  ax1.TickLabelInterpreter = 'latex';
xticklabels(ax1, xtick_lbls);
xlim(ax1, [10^e_min, 10^e_max]);
ylabel(ax1, 'Density','Interpreter', 'latex', 'FontSize', 50);
xlabel(ax1, 'Inter-execution times', 'Interpreter', 'latex', 'FontSize', 50);
legend(ax1, [h1,h2], '$\nu_0 = 2$', '$\nu_0 = 10$', ...
    'Interpreter', 'latex', 'FontSize', 50, 'Location', 'northwest');
grid(ax1, 'off');

%% ── Figure 2 : Omega comparison ──────────────────────────────────────────
[e_min, e_max] = make_axes(iet_Omega01_all, iet_Omega06_all);
edges      = logspace(e_min, e_max, nbins);
xtick_exp  = e_min:0.2:e_max;
xtick_vals = 10.^xtick_exp;
xtick_lbls = arrayfun(@(v) sprintf('$10^{%.1f}$',v), xtick_exp, 'UniformOutput', false);

figure('Color', 'w');
ax2 = gca;  hold(ax2, 'on');

h4 = histogram(iet_Omega06_all, edges, ...  % blue — background
    'Normalization', 'probability', 'FaceColor', blue, 'FaceAlpha', 0.95, 'EdgeColor', 'none');
h3 = histogram(iet_Omega01_all, edges, ...  % red  — foreground
    'Normalization', 'probability', 'FaceColor', red,  'FaceAlpha', 0.95, 'EdgeColor', 'none');

ax2.XScale = 'log';  ax2.FontSize = 30;  ax2.LineWidth = 2;  ax2.Box = 'on';
ax2.XTick = xtick_vals;  ax2.TickLabelInterpreter = 'latex';
xticklabels(ax2, xtick_lbls);
xlim(ax2, [10^e_min, 10^e_max]);
ylabel(ax2, 'Density',               'Interpreter', 'latex', 'FontSize', 50);
xlabel(ax2, 'Inter-execution times', 'Interpreter', 'latex', 'FontSize', 50);
legend(ax2, [h3,h4], '$\Omega = 0.1$', '$\Omega = 0.6$', ...
    'Interpreter', 'latex', 'FontSize', 50, 'Location', 'northwest');
grid(ax2, 'off');
