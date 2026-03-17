%% Profile Likelihood Analysis v3  (parfor-accelerated)
%  Why a_rs and a_sr Cannot Be Shared Across MOI Treatments
%
%  RUNTIME ESTIMATE
%    Serial equivalent : ~10-40 min  (320 lsqnonlin x ~100 ode15s each)
%    2 workers         : ~5-20 min
%    4 workers         : ~3-10 min    (typical laptop)
%    8+ workers        : ~2-5 min
%
%    Each sweep (a_rs, a_sr) runs as its own parfor over 40 grid points.
%    The pool stays alive between sweeps so there is no re-launch overhead.
%    The inner loop over 4 treatments stays serial inside each worker
%    (only 4 iterations -- overhead of nested parfor would exceed benefit).
%    NOTE: both sweeps cannot be merged into one 80-task parfor because
%    MATLAB's parfor classifier requires sliced outputs to be indexed with
%    the loop variable DIRECTLY -- conditional indexing (gi = task - n_grid)
%    violates this rule and causes a variable-classification error.
%
%  OUTPUTS (all written automatically):
%    profile_likelihood_results.mat            <- all computed profile data
%    profile_likelihood_interspecific.pdf      <- Figure 1: R^2 profile curves
%    profile_optimal_values.pdf               <- Figure 2: bar chart
%    profile_likelihood_interpretation.pdf    <- written interpretation + ms text
%
%  PREREQUISITES:
%    fitting_results_4treatments_Y25747F.mat           (M1 results)
%    fitting_results_4treatments_Y25747F_revised.mat   (M2 results)
%    dataset_Y25747F.xlsx
%    generate_interpretation_pdf.py                    (companion script)
% -------------------------------------------------------------------------

clear; clc; close all;

% =========================================================================
% 0.  Parallel pool setup
% =========================================================================
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting parallel pool...\n');
    pool = parpool('local');   % uses all available cores by default
end
fprintf('Parallel pool active: %d workers\n\n', pool.NumWorkers);

% =========================================================================
% 1.  Load data and fixed parameters
% =========================================================================
excelFileName = 'dataset_25747F.xlsx';
sheetName     = '25747_Fern_Grouped';
data_table    = readtable(excelFileName, 'Sheet', sheetName);
time          = data_table.Time / 60;          % minutes -> hours
OD_to_CFU     = 1e9;

all_data_mean = {data_table.Low_Group_Mean        * OD_to_CFU, ...
                 data_table.Medium_Group_Mean      * OD_to_CFU, ...
                 data_table.Medium_High_Group_Mean * OD_to_CFU, ...
                 data_table.High_Group_Mean        * OD_to_CFU};

P0_values       = [0.203, 2.03, 20.3, 203];
treatment_names = {'Low', 'Medium', 'Med-High', 'High'};
n_treatments    = 4;
n_grid          = 40;

K          = 1e9;
r_fixed    = 0.3730;
a_ss_fixed = 1.0963;
S0_fixed   = 7e6;

% --- M1 shared reference values ---
m1             = load('fitting_results_4treatments_25747F_5shared.mat');
a_rs_M1_shared = m1.results.shared.a_rs;
a_sr_M1_shared = m1.results.shared.a_sr;
M1_R2          = m1.results.R_squared(:);

% --- M2 starting points ---
m2          = load('fitting_results_4treatments_25747F_3shared.mat');
r2          = m2.results;
a_rr_shared = r2.shared.a_rr;
beta_shared = r2.shared.beta;
m_shared    = r2.shared.m;
a_rs_M2     = r2.treatment_specific.a_rs(:);
a_sr_M2     = r2.treatment_specific.a_sr(:);
eps_init    = r2.treatment_specific.epsilon(:);
gamma_init  = r2.treatment_specific.gamma(:);
delta_init  = r2.treatment_specific.delta(:);
phi_s_init  = r2.treatment_specific.phi_s(:);
R0_init     = r2.treatment_specific.R0(:);
M2_R2       = r2.R_squared(:);

% =========================================================================
% 2.  Profile grids  (log-spaced, centred on M2 optimal range)
% =========================================================================
a_rs_grid = logspace(log10(max(0.01, min(a_rs_M2)*0.1)), ...
                     log10(max(a_rs_M2)*5), n_grid);
a_sr_grid = logspace(log10(max(0.001, min(a_sr_M2)*0.1)), ...
                     log10(max(a_sr_M2)*5), n_grid);

fprintf('a_rs grid : [%.3f , %.1f]  (%d log-spaced points)\n', ...
    a_rs_grid(1), a_rs_grid(end), n_grid);
fprintf('a_sr grid : [%.4f , %.2f]  (%d log-spaced points)\n\n', ...
    a_sr_grid(1), a_sr_grid(end), n_grid);

% =========================================================================
% 3.  Optimisation options
%     Use 'Display','off' — essential for parfor (parallel workers cannot
%     write to shared console without garbling output).
% =========================================================================
opts = optimoptions('lsqnonlin', ...
    'Display','off', 'MaxIterations',300, ...
    'MaxFunctionEvaluations',6e4, ...
    'FunctionTolerance',1e-10, 'StepTolerance',1e-10);

lb_ts = [1e-5; 1e-5; 1e-6; 1e-9; 1.0];
ub_ts = [0.50; 1.00; 1e-3; 1e-6; 1e7];

% =========================================================================
% 4.  Broadcast copies of variables used inside parfor
%     (MATLAB requires these to be sliced or broadcast; listing them
%      explicitly as local copies clarifies intent and avoids warnings)
% =========================================================================
% These are read-only inside the parfor loops:
time_par       = time;
data_par       = all_data_mean;       % cell array — broadcast
P0_par         = P0_values;
r_par          = r_fixed;
a_ss_par       = a_ss_fixed;
S0_par         = S0_fixed;
K_par          = K;
a_rr_par       = a_rr_shared;
beta_par       = beta_shared;
m_par          = m_shared;
a_rs_M2_par    = a_rs_M2;            % treatment-specific M2 values
a_sr_M2_par    = a_sr_M2;
eps_par        = eps_init;
gamma_par      = gamma_init;
delta_par      = delta_init;
phi_s_par      = phi_s_init;
R0_par         = R0_init;
nt_par         = n_treatments;
opts_par       = opts;
lb_par         = lb_ts;
ub_par         = ub_ts;

% =========================================================================
% 5.  Compute a_rs profile  (parfor over 40 grid points)
%
%     MATLAB parfor requires that sliced output variables are indexed with
%     the loop variable DIRECTLY (e.g. R2_profile_ars(gi,:)).
%     That is why both sweeps run as separate parfor loops rather than a
%     single combined loop: a combined loop would need conditional indexing
%     (gi = task - n_grid), which the parfor classifier rejects.
%     The pool stays alive between the two loops so there is no re-launch
%     overhead.
% =========================================================================
fprintf('Computing a_rs profile  (%d grid pts x %d treatments, %d workers)...\n', ...
    n_grid, n_treatments, pool.NumWorkers);

R2_profile_ars = nan(n_grid, n_treatments);   % pre-allocate (sliced output)
t_start = tic;

parfor gi = 1:n_grid                          % loop variable indexes output directly
    a_rs_val = a_rs_grid(gi);                 % scalar broadcast -- fine in parfor
    R2_row   = nan(1, nt_par);               % local accumulator for this gi

    for ti = 1:nt_par                         % serial inner loop (only 4 iters)
        p0  = local_clip([eps_par(ti); gamma_par(ti); delta_par(ti); ...
                          phi_s_par(ti); R0_par(ti)], lb_par, ub_par);
        obj = @(p) residuals_ts(p, time_par, data_par{ti}, P0_par(ti), ...
                       r_par, a_ss_par, S0_par, K_par, ...
                       a_rs_val, a_sr_M2_par(ti), ...
                       a_rr_par, beta_par, m_par);
        try
            [~, rn] = lsqnonlin(obj, p0, lb_par, ub_par, opts_par);
            sst = sum((data_par{ti} - mean(data_par{ti})).^2);
            R2_row(ti) = 1 - rn / sst;
        catch
        end
    end
    R2_profile_ars(gi, :) = R2_row;          % sliced write: gi indexes row directly
end

t_ars = toc(t_start);
fprintf('  a_rs complete in %.1f s  (~%.1f s/task on %d workers)\n\n', ...
    t_ars, t_ars * pool.NumWorkers / n_grid, pool.NumWorkers);

% =========================================================================
% 6.  Compute a_sr profile  (parfor over 40 grid points)
% =========================================================================
fprintf('Computing a_sr profile  (%d grid pts x %d treatments, %d workers)...\n', ...
    n_grid, n_treatments, pool.NumWorkers);

R2_profile_asr = nan(n_grid, n_treatments);
t_start = tic;

parfor gi = 1:n_grid
    a_sr_val = a_sr_grid(gi);
    R2_row   = nan(1, nt_par);

    for ti = 1:nt_par
        p0  = local_clip([eps_par(ti); gamma_par(ti); delta_par(ti); ...
                          phi_s_par(ti); R0_par(ti)], lb_par, ub_par);
        obj = @(p) residuals_ts(p, time_par, data_par{ti}, P0_par(ti), ...
                       r_par, a_ss_par, S0_par, K_par, ...
                       a_rs_M2_par(ti), a_sr_val, ...
                       a_rr_par, beta_par, m_par);
        try
            [~, rn] = lsqnonlin(obj, p0, lb_par, ub_par, opts_par);
            sst = sum((data_par{ti} - mean(data_par{ti})).^2);
            R2_row(ti) = 1 - rn / sst;
        catch
        end
    end
    R2_profile_asr(gi, :) = R2_row;
end

t_asr = toc(t_start);
fprintf('  a_sr complete in %.1f s  (~%.1f s/task on %d workers)\n\n', ...
    t_asr, t_asr * pool.NumWorkers / n_grid, pool.NumWorkers);

% =========================================================================
% 7.  Derive summary statistics
% =========================================================================
optimal_ars   = zeros(1, n_treatments);
optimal_asr   = zeros(1, n_treatments);
R2_at_opt_ars = zeros(1, n_treatments);
R2_at_opt_asr = zeros(1, n_treatments);
R2_at_M1_ars  = zeros(1, n_treatments);
R2_at_M1_asr  = zeros(1, n_treatments);

for ti = 1:n_treatments
    [R2_at_opt_ars(ti), ia] = max(R2_profile_ars(:,ti));
    [R2_at_opt_asr(ti), is] = max(R2_profile_asr(:,ti));
    optimal_ars(ti) = a_rs_grid(ia);
    optimal_asr(ti) = a_sr_grid(is);
    [~, im1a] = min(abs(a_rs_grid - a_rs_M1_shared));
    [~, im1s] = min(abs(a_sr_grid - a_sr_M1_shared));
    R2_at_M1_ars(ti) = R2_profile_ars(im1a, ti);
    R2_at_M1_asr(ti) = R2_profile_asr(im1s, ti);
end

fold_range_ars = max(optimal_ars) / min(optimal_ars);
fold_range_asr = max(optimal_asr) / min(optimal_asr);
R2_cost_ars    = R2_at_opt_ars - R2_at_M1_ars;
R2_cost_asr    = R2_at_opt_asr - R2_at_M1_asr;

% =========================================================================
% 8.  Console summary
% =========================================================================
fprintf('%s\n', repmat('=',1,74));
fprintf('PROFILE LIKELIHOOD SUMMARY\n');
fprintf('%s\n', repmat('=',1,74));
fprintf('%-10s  %10s  %10s  %10s  %10s  %10s\n', ...
    'Treatment','Opt a_rs','M1 a_rs','R2@opt','R2@M1','R2 cost');
fprintf('%s\n', repmat('-',1,74));
for ti = 1:n_treatments
    fprintf('%-10s  %10.2f  %10.2f  %10.4f  %10.4f  %10.4f\n', ...
        treatment_names{ti}, optimal_ars(ti), a_rs_M1_shared, ...
        R2_at_opt_ars(ti), R2_at_M1_ars(ti), R2_cost_ars(ti));
end
fprintf('  Fold-range a_rs : %.1fx\n\n', fold_range_ars);

fprintf('%-10s  %10s  %10s  %10s  %10s  %10s\n', ...
    'Treatment','Opt a_sr','M1 a_sr','R2@opt','R2@M1','R2 cost');
fprintf('%s\n', repmat('-',1,74));
for ti = 1:n_treatments
    fprintf('%-10s  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
        treatment_names{ti}, optimal_asr(ti), a_sr_M1_shared, ...
        R2_at_opt_asr(ti), R2_at_M1_asr(ti), R2_cost_asr(ti));
end
fprintf('  Fold-range a_sr : %.1fx\n\n', fold_range_asr);

% =========================================================================
% 9.  Save .mat
% =========================================================================
profile_results = struct( ...
    'a_rs_grid',          a_rs_grid, ...
    'a_sr_grid',          a_sr_grid, ...
    'R2_profile_ars',     R2_profile_ars, ...
    'R2_profile_asr',     R2_profile_asr, ...
    'n_grid',             n_grid, ...
    'n_treatments',       n_treatments, ...
    'treatment_names',    {treatment_names}, ...
    'P0_values',          P0_values, ...
    'optimal_ars',        optimal_ars, ...
    'optimal_asr',        optimal_asr, ...
    'R2_at_opt_ars',      R2_at_opt_ars, ...
    'R2_at_opt_asr',      R2_at_opt_asr, ...
    'fold_range_ars',     fold_range_ars, ...
    'fold_range_asr',     fold_range_asr, ...
    'a_rs_M1_shared',     a_rs_M1_shared, ...
    'a_sr_M1_shared',     a_sr_M1_shared, ...
    'a_rs_M2',            a_rs_M2, ...
    'a_sr_M2',            a_sr_M2, ...
    'R2_at_M1_ars',       R2_at_M1_ars, ...
    'R2_at_M1_asr',       R2_at_M1_asr, ...
    'R2_cost_ars',        R2_cost_ars, ...
    'R2_cost_asr',        R2_cost_asr, ...
    'M1_R2',              M1_R2, ...
    'M2_R2',              M2_R2, ...
    'n_workers',          pool.NumWorkers, ...
    'compute_time_ars_s', t_ars, ...
    'compute_time_asr_s', t_asr ...
);
save('ATCC25747_profile_likelihood_results.mat', 'profile_results');
fprintf('Saved : ATCC25747_profile_likelihood_results.mat\n\n');

% =========================================================================
% 10. Figure 1 — Profile R^2 curves
% =========================================================================
colors = [0.20 0.60 1.00;
          0.10 0.65 0.10;
          1.00 0.50 0.00;
          0.80 0.20 0.20];

fig1 = figure('Units','inches','Position',[0.5 0.5 14 6]);
tl1  = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl1,'Profile Likelihood: Per-Treatment Optimal a_{rs} and a_{sr}', ...
    'FontSize',13,'FontWeight','bold');

%  --- Panel A : a_rs ---
nexttile; hold on;
patch([min(optimal_ars) max(optimal_ars) max(optimal_ars) min(optimal_ars)], ...
      [-0.4 -0.4 1.1 1.1],[1 0.8 0.8],'EdgeColor','none', ...
      'FaceAlpha',0.30,'DisplayName','Conflict zone');
for ti=1:n_treatments
    plot(a_rs_grid, R2_profile_ars(:,ti),'-', ...
        'Color',colors(ti,:),'LineWidth',2.5,'DisplayName',treatment_names{ti});
end
for ti=1:n_treatments
    xline(optimal_ars(ti),'--','Color',colors(ti,:),'LineWidth',1.5,'Alpha',0.8);
    text(optimal_ars(ti), 0.04+0.07*(ti-1), sprintf('  %.1f',optimal_ars(ti)), ...
        'Color',colors(ti,:),'FontSize',8.5,'FontWeight','bold');
end
xline(a_rs_M1_shared,'k-','LineWidth',2.5, ...
    'DisplayName',sprintf('M1 shared = %.2f',a_rs_M1_shared));
text(a_rs_M1_shared, 0.37, ...
    sprintf('  M1 shared\n  (%.1f)',a_rs_M1_shared), ...
    'FontSize',8.5,'FontWeight','bold','Color','k');
set(gca,'XScale','log','FontSize',10);
xlabel('a_{rs} (fixed value)','FontSize',11);
ylabel('Best achievable R^2','FontSize',11);
title({sprintf('(A) Profile R^2 vs fixed a_{rs}'), ...
       sprintf('Optimal values span %.1f-fold range',fold_range_ars)},'FontSize',10);
ylim([-0.35 1.08]); xlim([a_rs_grid(1)*0.8 a_rs_grid(end)*1.2]);
legend('Conflict zone',treatment_names{:},'M1 shared', ...
    'Location','southwest','FontSize',8.5);
grid on;

%  --- Panel B : a_sr ---
nexttile; hold on;
patch([min(optimal_asr) max(optimal_asr) max(optimal_asr) min(optimal_asr)], ...
      [-0.4 -0.4 1.1 1.1],[1 0.8 0.8],'EdgeColor','none', ...
      'FaceAlpha',0.30,'DisplayName','Conflict zone');
for ti=1:n_treatments
    plot(a_sr_grid, R2_profile_asr(:,ti),'-', ...
        'Color',colors(ti,:),'LineWidth',2.5,'DisplayName',treatment_names{ti});
end
for ti=1:n_treatments
    xline(optimal_asr(ti),'--','Color',colors(ti,:),'LineWidth',1.5,'Alpha',0.8);
    text(optimal_asr(ti), 0.04+0.07*(ti-1), sprintf('  %.3f',optimal_asr(ti)), ...
        'Color',colors(ti,:),'FontSize',8.5,'FontWeight','bold');
end
xline(a_sr_M1_shared,'k-','LineWidth',2.5, ...
    'DisplayName',sprintf('M1 shared = %.3f',a_sr_M1_shared));
text(a_sr_M1_shared, 0.37, ...
    sprintf('  M1 shared\n  (%.2f)',a_sr_M1_shared), ...
    'FontSize',8.5,'FontWeight','bold','Color','k');
set(gca,'XScale','log','FontSize',10);
xlabel('a_{sr} (fixed value)','FontSize',11);
ylabel('Best achievable R^2','FontSize',11);
title({sprintf('(B) Profile R^2 vs fixed a_{sr}'), ...
       sprintf('Optimal values span %.1f-fold range',fold_range_asr)},'FontSize',10);
ylim([-0.35 1.08]); xlim([a_sr_grid(1)*0.8 a_sr_grid(end)*1.2]);
legend('Conflict zone',treatment_names{:},'M1 shared', ...
    'Location','southwest','FontSize',8.5);
grid on;

exportgraphics(fig1,'ATCC25747_profile_likelihood_interspecific.pdf', ...
    'ContentType','vector','Resolution',300);
fprintf('Saved : ATCC25747_profile_likelihood_interspecific.pdf\n');

% =========================================================================
% 11. Figure 2 — Optimal-values bar chart
% =========================================================================
fig2 = figure('Units','inches','Position',[0.5 0.5 12 5]);
tl2  = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl2,['Treatment-Specific Optimal a_{rs} and a_{sr}: ' ...
           'Incompatible with a Single Shared Value'], ...
    'FontSize',12,'FontWeight','bold');

nexttile; hold on;
bh1 = bar(1:n_treatments, optimal_ars, 0.55,'FaceColor','flat');
for ti=1:n_treatments; bh1.CData(ti,:)=colors(ti,:); end
plot(1:n_treatments, a_rs_M2,'ko','MarkerSize',9,'MarkerFaceColor','w', ...
    'LineWidth',2,'DisplayName','M2 fitted');
yline(a_rs_M1_shared,'r-','LineWidth',2.5, ...
    'DisplayName',sprintf('M1 shared = %.2f',a_rs_M1_shared));
for ti=1:n_treatments
    text(ti, optimal_ars(ti)*1.10, sprintf('%.1f',optimal_ars(ti)), ...
        'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end
set(gca,'YScale','log','XTick',1:n_treatments, ...
    'XTickLabel',treatment_names,'FontSize',10);
ylabel('Optimal a_{rs}  (log scale)','FontSize',11);
title('(A)  Optimal a_{rs} per treatment','FontSize',11);
legend('Profile optimum','M2 fitted','M1 shared', ...
    'Location','northeast','FontSize',9);
text(0.03,0.07,sprintf('Fold-range: %.0fx',fold_range_ars), ...
    'Units','normalized','FontSize',9,'Color',[0.7 0 0],'FontWeight','bold');
grid on;

nexttile; hold on;
bh2 = bar(1:n_treatments, optimal_asr, 0.55,'FaceColor','flat');
for ti=1:n_treatments; bh2.CData(ti,:)=colors(ti,:); end
plot(1:n_treatments, a_sr_M2,'ko','MarkerSize',9,'MarkerFaceColor','w', ...
    'LineWidth',2,'DisplayName','M2 fitted');
yline(a_sr_M1_shared,'r-','LineWidth',2.5, ...
    'DisplayName',sprintf('M1 shared = %.3f',a_sr_M1_shared));
for ti=1:n_treatments
    text(ti, optimal_asr(ti)*1.25, sprintf('%.3f',optimal_asr(ti)), ...
        'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end
set(gca,'YScale','log','XTick',1:n_treatments, ...
    'XTickLabel',treatment_names,'FontSize',10);
ylabel('Optimal a_{sr}  (log scale)','FontSize',11);
title('(B)  Optimal a_{sr} per treatment','FontSize',11);
legend('Profile optimum','M2 fitted','M1 shared', ...
    'Location','northeast','FontSize',9);
text(0.03,0.07,sprintf('Fold-range: %.0fx',fold_range_asr), ...
    'Units','normalized','FontSize',9,'Color',[0.7 0 0],'FontWeight','bold');
grid on;

exportgraphics(fig2,'25747_profile_optimal_values.pdf', ...
    'ContentType','vector','Resolution',300);
fprintf('Saved : 25747_profile_optimal_values.pdf\n');

% =========================================================================
% 12. Call Python to generate interpretation PDF
% =========================================================================
fprintf('\nGenerating interpretation PDF...\n');
[status, msg] = system('python3 25747_generate_interpretation_pdf.py');
if status == 0
    fprintf('Saved : 25747_profile_likelihood_interpretation.pdf\n');
else
    fprintf('Python call returned status %d.\n  %s\n', status, msg);
    fprintf('Run manually: python3 generate_interpretation_pdf.py\n');
end

fprintf('\n=== All done. Total wall time: a_rs %.1fs  a_sr %.1fs ===\n', ...
    t_ars, t_asr);

% =========================================================================
%% Local helper functions
%  NOTE: local functions are visible to parfor workers when they are
%  defined in the same file as the parfor loop.
% =========================================================================

function v = local_clip(v, lo, hi)
    v = max(min(v, hi), lo);
end

function res = residuals_ts(p_ts, time, data_mean, P0, ...
        r_fixed, a_ss_fixed, S0_fixed, K, ...
        a_rs_val, a_sr_val, a_rr, beta, m)
    epsilon = p_ts(1);  gamma_v = p_ts(2);
    delta   = p_ts(3);  phi_s   = p_ts(4);  R0 = p_ts(5);
    y0   = [S0_fixed; R0; P0];
    ode_opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'NonNegative',[1 2 3]);
    f = @(t,y) [ ...
        r_fixed*y(1)*(1-(a_ss_fixed*y(1)+a_rs_val*y(2))/K) ...
            - delta*y(1) - phi_s*y(3)*y(1); ...
        gamma_v*r_fixed*y(2)*(1-(a_sr_val*y(1)+a_rr*y(2))/K) ...
            + delta*y(1) - epsilon*phi_s*y(3)*y(2); ...
        beta*phi_s*(y(1)+epsilon*y(2))*y(3) - m*y(3) ];
    try
        [~,Y] = ode15s(f, time, y0, ode_opts);
        res   = Y(:,1) + Y(:,2) - data_mean;
    catch
        res = 1e15*ones(size(time));
    end
end