%% Model Comparison Evidence for Reviewer 2 Response
%
% This script generates four pieces of evidence to justify that
% a_rs and a_sr must be treatment-specific rather than shared:
%
%   FIGURE 1 – Side-by-side fit comparison (M1 vs M2), 2x4 panels
%   FIGURE 2 – Residual analysis (M1 vs M2)
%   TABLE    – AIC / BIC / R² formal model selection
%
% MODEL DEFINITIONS
% -----------------
%   M1 (shared a_rs, a_sr):  4 fixed + 5 shared + 20 treatment-specific = 25 params
%   M2 (free   a_rs, a_sr):  4 fixed + 3 shared + 28 treatment-specific = 31 params
%
% IMPORTANT: Run your two fitting scripts first so the .mat result files exist:
%   fitting_results_4treatments_25747F_5shared.mat          <- M1 results
%   fitting_results_4treatments_25747F_3shared.mat          <- M2 results
% -------------------------------------------------------------------------

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

excelFileName = 'dataset_25747F.xlsx';
sheetName     = '25747_Fern_Grouped';
data_table    = readtable(excelFileName, 'Sheet', sheetName);
time          = data_table.Time / 60;          % minutes → hours
OD_to_CFU     = 1e9;

all_data_mean = {data_table.Low_Group_Mean         * OD_to_CFU, ...
                 data_table.Medium_Group_Mean       * OD_to_CFU, ...
                 data_table.Medium_High_Group_Mean  * OD_to_CFU, ...
                 data_table.High_Group_Mean         * OD_to_CFU};

all_data_std  = {data_table.Low_Group_Std          * OD_to_CFU, ...
                 data_table.Medium_Group_Std        * OD_to_CFU, ...
                 data_table.Medium_High_Group_Std   * OD_to_CFU, ...
                 data_table.High_Group_Std          * OD_to_CFU};

P0_values       = [0.203, 2.03, 20.3, 203];
treatment_names = {'Low', 'Medium', 'Med-High', 'High'};
n_treatments    = 4;

% Fixed parameters (same for both models)
K          = 1e9;
r_fixed    = 0.3730;
a_ss_fixed = 1.0963;
S0_fixed   = 7e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load Fitted Parameters from Both Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- M1: shared a_rs, a_sr ---
m1 = load('fitting_results_4treatments_25747F_5shared.mat');
r1 = m1.results;

% Reconstruct M1 parameter vector (25 params)
% Layout: [a_rs, a_sr, a_rr, beta, m | eps x4 | gamma x4 | delta x4 | phi_s x4 | R0 x4]
params_M1 = [r1.shared.a_rs;  r1.shared.a_sr;  r1.shared.a_rr; ...
             r1.shared.beta;  r1.shared.m; ...
             r1.treatment_specific.epsilon(:); ...
             r1.treatment_specific.gamma(:);   ...
             r1.treatment_specific.delta(:);   ...
             r1.treatment_specific.phi_s(:);   ...
             r1.treatment_specific.R0(:)];

% --- M2: treatment-specific a_rs, a_sr ---
m2 = load('fitting_results_4treatments_25747F_3shared.mat');
r2 = m2.results;

% Reconstruct M2 parameter vector (31 params)
% Layout: [a_rr, beta, m | a_rs x4 | a_sr x4 | eps x4 | gamma x4 | delta x4 | phi_s x4 | R0 x4]
params_M2 = [r2.shared.a_rr; r2.shared.beta; r2.shared.m; ...
             r2.treatment_specific.a_rs(:);   ...
             r2.treatment_specific.a_sr(:);   ...
             r2.treatment_specific.epsilon(:); ...
             r2.treatment_specific.gamma(:);   ...
             r2.treatment_specific.delta(:);   ...
             r2.treatment_specific.phi_s(:);   ...
             r2.treatment_specific.R0(:)];

n_params_M1 = 25;   % fitted (excluding 4 fixed)
n_params_M2 = 31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Compute Fit Quality + AIC/BIC for Both Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stats_M1, resid_M1, sim_M1] = compute_all_stats(params_M1, 'M1', ...
    time, all_data_mean, P0_values, r_fixed, a_ss_fixed, S0_fixed, K, n_params_M1);

[stats_M2, resid_M2, sim_M2] = compute_all_stats(params_M2, 'M2', ...
    time, all_data_mean, P0_values, r_fixed, a_ss_fixed, S0_fixed, K, n_params_M2);

% Print formal model selection table
fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('FORMAL MODEL SELECTION TABLE\n');
fprintf('%s\n', repmat('=', 1, 78));
fprintf('%-12s | %6s | %10s | %10s | %10s | %10s | %10s\n', ...
    'Model', 'Params', 'R²_Low', 'R²_Med', 'R²_MedHi', 'R²_High', 'R²_Overall');
fprintf('%s\n', repmat('-', 1, 78));
fprintf('%-12s | %6d | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
    'M1 (shared)', n_params_M1, ...
    stats_M1.R2(1), stats_M1.R2(2), stats_M1.R2(3), stats_M1.R2(4), ...
    stats_M1.R2_overall);
fprintf('%-12s | %6d | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
    'M2 (free)', n_params_M2, ...
    stats_M2.R2(1), stats_M2.R2(2), stats_M2.R2(3), stats_M2.R2(4), ...
    stats_M2.R2_overall);
fprintf('%s\n', repmat('=', 1, 78));
fprintf('%-12s | %6s | %10s | %10s | %10s\n', ...
    'Model', 'Params', 'AIC', 'BIC', 'ΔAIC vs M1');
fprintf('%s\n', repmat('-', 1, 78));
fprintf('%-12s | %6d | %10.2f | %10.2f | %10s\n', ...
    'M1 (shared)', n_params_M1, stats_M1.AIC, stats_M1.BIC, 'reference');
fprintf('%-12s | %6d | %10.2f | %10.2f | %10.2f\n', ...
    'M2 (free)', n_params_M2, stats_M2.AIC, stats_M2.BIC, ...
    stats_M2.AIC - stats_M1.AIC);
fprintf('%s\n', repmat('=', 1, 78));
fprintf('Interpretation: ΔAIC < −10 = strong evidence for M2; ΔBIC < −10 = very strong.\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1 – Side-by-Side Fit Comparison (2 × 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors_treat = [0.20, 0.60, 1.00;   % Blue   Low
                0.10, 0.65, 0.10;   % Green  Medium
                1.00, 0.50, 0.00;   % Orange Med-High
                0.80, 0.20, 0.20];  % Red    High

fig1 = figure('Units', 'inches', 'Position', [0.5, 0.5, 16, 9]);
tl1  = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, 'Model Fit Comparison for 25747 Fern: M1 (shared a_{rs}, a_{sr}) vs M2 (treatment-specific a_{rs}, a_{sr})', ...
    'FontSize', 18, 'FontWeight', 'bold');

row_labels = {'M1 (shared)', 'M2 (free)'};
all_params = {params_M1, params_M2};
all_sims   = {sim_M1, sim_M2};
model_ids  = {'M1', 'M2'};

tspan_plot = linspace(0, max(time), 500)';

for row = 1:2
    for col = 1:n_treatments
        nexttile((row-1)*4 + col);

        % Data
        errorbar(time, all_data_mean{col}, all_data_std{col}, 'o', ...
            'MarkerSize', 4, 'Color', colors_treat(col,:), ...
            'MarkerFaceColor', colors_treat(col,:), ...
            'LineWidth', 1, 'CapSize', 2.5, 'DisplayName', 'Data');
        hold on;

        % Simulated total (S+R)
        [~, S_p, R_p, ~] = simulate_model_dispatch(all_params{row}, ...
            tspan_plot, P0_values(col), col, r_fixed, a_ss_fixed, S0_fixed, K, model_ids{row});
        plot(tspan_plot, S_p + R_p, 'k-', 'LineWidth', 3, 'DisplayName', 'S+R');
        plot(tspan_plot, S_p,       'b--','LineWidth', 1.5, 'DisplayName', 'S');
        plot(tspan_plot, R_p,       'r:', 'LineWidth', 1.5, 'DisplayName', 'R');

        % Labels
        R2_val = stats_M1.R2(col);
        if row == 2; R2_val = stats_M2.R2(col); end

        % Colour R² red when poor
        r2_color = 'k';
        if R2_val < 0.7; r2_color = [0.8 0 0]; end

        title(sprintf('%s | %s ', treatment_names{col}, row_labels{row}), ...
            'FontSize', 14);
        % Overlay R² as text with colour
        ax = gca;
        text(0.97, 0.97, sprintf('R^2 = %.3f', R2_val), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold', ...
            'Color', r2_color);

        xlabel('Time (h)', 'FontSize', 12);
        ylabel('CFU/mL',   'FontSize', 12);
        xlim([0 max(time)]);
        grid on;
        set(gca, 'FontSize', 10);

        if row == 1 && col == 1
            legend('Data','S+R','S','R', 'Location','northwest','FontSize',7);
        end
    end
end

exportgraphics(fig1, '25747F_fit_comparison_M1vsM2.pdf', ...
    'ContentType','vector','Resolution',300);
fprintf('Figure 1 saved: 25747F_fit_comparison_M1vsM2.pdf\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2 – Residual Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig2 = figure('Units', 'inches', 'Position', [0.5, 0.5, 16, 10]);
tl2  = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Residual Analysis for 25747 Fern: M1 vs M2', ...
    'FontSize', 18, 'FontWeight', 'bold');

% Row 1: M1 residuals vs time
% Row 2: M2 residuals vs time
% Row 3: Absolute residual comparison bar chart per treatment

max_abs_resid = 0;
for col = 1:n_treatments
    max_abs_resid = max(max_abs_resid, ...
        max(abs([resid_M1{col}; resid_M2{col}])));
end
y_lim_resid = max_abs_resid * 1.15;

for row = 1:2
    all_r = {resid_M1, resid_M2};
    mdl_label = row_labels{row};

    for col = 1:n_treatments
        nexttile((row-1)*4 + col);

        res = all_r{row}{col};
        % Colour bars: blue positive, red negative
        pos_idx = res >= 0;
        bar(time(pos_idx),  res(pos_idx),  'FaceColor', [0.25 0.55 0.95], ...
            'EdgeColor','none','BarWidth',0.8);
        hold on;
        bar(time(~pos_idx), res(~pos_idx), 'FaceColor', [0.95 0.30 0.25], ...
            'EdgeColor','none','BarWidth',0.8);
        yline(0, 'k-', 'LineWidth', 1);

        % Overlay smoothed trend to reveal systematic pattern
        if length(time) > 5
            smooth_r = movmean(res, max(3, round(length(time)/10)));
            plot(time, smooth_r, 'k-', 'LineWidth', 2.5);
        end

        ylim([-y_lim_resid, y_lim_resid]);
        xlabel('Time (h)', 'FontSize', 12);
        ylabel('Residual (CFU/mL)', 'FontSize', 12);
        title(sprintf('%s  |  %s', mdl_label, treatment_names{col}), 'FontSize', 14);

        % Annotate with RMSE and sign-run test p-value
        rmse_val = sqrt(mean(res.^2));
        text(0.97, 0.97, sprintf('RMSE=%.2e', rmse_val), ...
            'Units','normalized','HorizontalAlignment','right', ...
            'VerticalAlignment','top','FontSize',10);

        grid on; set(gca,'FontSize',12);
    end
end

% Row 3: RMSE bar comparison per treatment
nexttile(9);  % position (3,1) in 3×4 layout - but we'll use all 4 tiles for bar charts
for col = 1:n_treatments
    nexttile(8 + col);

    rmse_m1 = sqrt(mean(resid_M1{col}.^2));
    rmse_m2 = sqrt(mean(resid_M2{col}.^2));

    bar_data = [rmse_m1, rmse_m2];
    b = bar(bar_data, 0.6);
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.70 0.85 1.00];   % light blue  M1
    b.CData(2,:) = [0.30 0.70 0.30];   % green       M2

    set(gca, 'XTickLabel', {'M1 (shared)','M2 (free)'}, 'FontSize', 12);
    ylabel('RMSE (CFU/mL)', 'FontSize', 8);
    title(sprintf('RMSE: %s', treatment_names{col}), 'FontSize', 14);
    text(1, rmse_m1*1.02, sprintf('%.2e', rmse_m1), ...
        'HorizontalAlignment','center','FontSize',10);
    text(2, rmse_m2*1.02, sprintf('%.2e', rmse_m2), ...
        'HorizontalAlignment','center','FontSize',10);
    grid on; set(gca,'FontSize',12);
end

exportgraphics(fig2, '25747F_residual_comparison_M1vsM2.pdf', ...
    'ContentType','vector','Resolution',300);
fprintf('Figure 2 saved: 25747F_residual_comparison_M1vsM2.pdf\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print Final Summary for Paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('SUMMARY FOR MANUSCRIPT / RESPONSE LETTER\n');
fprintf('%s\n', repmat('=', 1, 78));
fprintf('\nModel M1 (a_rs, a_sr shared):  %d free parameters\n', n_params_M1);
fprintf('  Per-treatment R²: %.3f | %.3f | %.3f | %.3f\n', ...
    stats_M1.R2(1), stats_M1.R2(2), stats_M1.R2(3), stats_M1.R2(4));
fprintf('  Overall R² = %.4f,  AIC = %.2f,  BIC = %.2f\n', ...
    stats_M1.R2_overall, stats_M1.AIC, stats_M1.BIC);

fprintf('\nModel M2 (a_rs, a_sr free):    %d free parameters\n', n_params_M2);
fprintf('  Per-treatment R²: %.3f | %.3f | %.3f | %.3f\n', ...
    stats_M2.R2(1), stats_M2.R2(2), stats_M2.R2(3), stats_M2.R2(4));
fprintf('  Overall R² = %.4f,  AIC = %.2f,  BIC = %.2f\n', ...
    stats_M2.R2_overall, stats_M2.AIC, stats_M2.BIC);

fprintf('\nΔAIC (M2 – M1) = %.2f\n', stats_M2.AIC - stats_M1.AIC);
fprintf('ΔBIC (M2 – M1) = %.2f\n', stats_M2.BIC - stats_M1.BIC);
fprintf('\nNote: ΔAIC < −10 indicates very strong evidence favouring M2.\n');
fprintf('%s\n\n', repmat('=', 1, 78));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stats, residuals_cell, sims_cell] = compute_all_stats(params, model_id, ...
    time, all_data_mean, P0_values, r_fixed, a_ss_fixed, S0_fixed, K, n_params)
    n_treatments   = length(all_data_mean);
    n_points       = length(time);
    R2             = zeros(1, n_treatments);
    RMSE           = zeros(1, n_treatments);
    SS_res_total   = 0;
    SS_tot_total   = 0;
    n_total        = 0;
    residuals_cell = cell(1, n_treatments);
    sims_cell      = cell(1, n_treatments);

    for i = 1:n_treatments
        [~, S_s, R_s, ~] = simulate_model_dispatch(params, time, P0_values(i), i, ...
            r_fixed, a_ss_fixed, S0_fixed, K, model_id);
        total_s = S_s + R_s;
        sims_cell{i} = total_s;

        d    = all_data_mean{i};
        res  = total_s - d;
        residuals_cell{i} = res;

        ss_r = sum(res.^2);
        ss_t = sum((d - mean(d)).^2);
        R2(i)   = 1 - ss_r / ss_t;
        RMSE(i) = sqrt(mean(res.^2));

        SS_res_total = SS_res_total + ss_r;
        SS_tot_total = SS_tot_total + ss_t;
        n_total      = n_total + n_points;
    end

    R2_overall  = 1 - SS_res_total / SS_tot_total;
    sigma2      = SS_res_total / n_total;
    AIC_val     = n_total * log(sigma2) + 2 * n_params;
    BIC_val     = n_total * log(sigma2) + n_params * log(n_total);

    stats.R2         = R2;
    stats.RMSE       = RMSE;
    stats.R2_overall = R2_overall;
    stats.AIC        = AIC_val;
    stats.BIC        = BIC_val;
end

function [t_out, S_out, R_out, P_out] = simulate_model_dispatch(params, time, P0, ...
    tidx, r_fixed, a_ss_fixed, S0_fixed, K, model_id)
    % Routes to the correct parameter-extraction scheme depending on model.

    if strcmp(model_id, 'M1')
        % M1 layout: [a_rs, a_sr, a_rr, beta, m | eps x4 | gamma x4 | delta x4 | phi_s x4 | R0 x4]
        a_rs    = params(1);
        a_sr    = params(2);
        a_rr    = params(3);
        beta    = params(4);
        m       = params(5);
        epsilon = params(5  + tidx);
        gamma   = params(9  + tidx);
        delta   = params(13 + tidx);
        phi_s   = params(17 + tidx);
        R0      = params(21 + tidx);

    else
        % M2 layout: [a_rr, beta, m | a_rs x4 | a_sr x4 | eps x4 | gamma x4 | delta x4 | phi_s x4 | R0 x4]
        a_rr    = params(1);
        beta    = params(2);
        m       = params(3);
        a_rs    = params(3  + tidx);
        a_sr    = params(7  + tidx);
        epsilon = params(11 + tidx);
        gamma   = params(15 + tidx);
        delta   = params(19 + tidx);
        phi_s   = params(23 + tidx);
        R0      = params(27 + tidx);
    end

    y0       = [S0_fixed; R0; P0];
    ode_opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'NonNegative',[1,2,3]);

    odefun = @(t, y) [ ...
        r_fixed * y(1) * (1 - (a_ss_fixed*y(1) + a_rs*y(2)) / K) ...
            - delta*y(1) - phi_s*y(3)*y(1); ...
        gamma * r_fixed * y(2) * (1 - (a_sr*y(1) + a_rr*y(2)) / K) ...
            + delta*y(1) - epsilon*phi_s*y(3)*y(2); ...
        beta * phi_s * (y(1) + epsilon*y(2)) * y(3) - m*y(3) ...
    ];

    try
        [t_out, Y] = ode15s(odefun, time, y0, ode_opts);
        S_out = Y(:,1);  R_out = Y(:,2);  P_out = Y(:,3);
    catch
        t_out = time;
        S_out = 1e15*ones(size(time));
        R_out = S_out;  P_out = S_out;
    end
end