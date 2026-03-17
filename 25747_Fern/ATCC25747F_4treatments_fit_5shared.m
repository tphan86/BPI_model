%% Global Fitting Including All 4 Treatments for 25747 + Fern
%
% This script fits the REVISED model to ALL 4 phage treatments simultaneously
% (Low, Medium, Med-High, High) following the same structure as Y3650F.
%
% REVISED MODEL (from BPI_revision.pdf):
%   dS/dt = r*S*(1 - (a_ss*S + a_rs*R)/K) - delta*S - phi_s*P*S
%   dR/dt = gamma*r*R*(1 - (a_sr*S + a_rr*R)/K) + delta*S - epsilon*phi_s*P*R
%   dP/dt = beta*phi_s*(S + epsilon*R)*P - m*P
%
% where epsilon is the resistance efficiency (epsilon << 1 means strong resistance)
%
% FIXED (from control): r, a_ss, S0, K
% SHARED: a_rs, a_sr, a_rr, beta, m
% TREATMENT-SPECIFIC: epsilon, gamma, delta, phi_s, R0 (5 × 4 = 20)
%
% Total: 5 shared + 20 treatment-specific = 25 parameters

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

excelFileName = 'dataset_25747F.xlsx';
sheetName = '25747_Fern_Grouped';

data_table = readtable(excelFileName, 'Sheet', sheetName);
time = data_table.Time / 60;  % Convert minutes to hours

OD_to_CFU = 1e9;

% Load ALL 4 phage treatments (excluding Control which is used for fixed params)
data_low_mean = data_table.Low_Group_Mean * OD_to_CFU;
data_low_std = data_table.Low_Group_Std * OD_to_CFU;

data_med_mean = data_table.Medium_Group_Mean * OD_to_CFU;
data_med_std = data_table.Medium_Group_Std * OD_to_CFU;

data_medhigh_mean = data_table.Medium_High_Group_Mean * OD_to_CFU;
data_medhigh_std = data_table.Medium_High_Group_Std * OD_to_CFU;

data_high_mean = data_table.High_Group_Mean * OD_to_CFU;
data_high_std = data_table.High_Group_Std * OD_to_CFU;

% P0 values for all 4 treatments (PFU/mL)
% From manuscript: "Phages Fern IDv1 were introduced at a maximum concentration of 2.0E2 PFU/mL"
% Adjust these based on your actual dilution series
P0_values = [0.203, 2.03, 20.3, 203];  % Low, Medium, Med-High, High

all_data_mean = {data_low_mean, data_med_mean, data_medhigh_mean, data_high_mean};
all_data_std = {data_low_std, data_med_std, data_medhigh_std, data_high_std};
treatment_names = {'Low', 'Medium', 'Med-High', 'High'};
n_treatments = 4;

fprintf('=== Fitting ALL 4 Treatments for Y25747 + Fern ===\n');
fprintf('Treatments: Low (P0=0.2), Medium (P0=2), Med-High (P0=20), High (P0=200)\n');
fprintf('Data points per treatment: %d\n', length(time));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed Parameters (from control fitting)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These should be determined from fitting the control (no phage) data
% You may need to adjust these values based on your control fitting results

K = 1e9;              % Carrying capacity (CFU/mL)
r_fixed = 0.3730;     % Intrinsic growth rate from control (h^-1)
a_ss_fixed = 1.0963;  % Intraspecific competition of S from control
S0_fixed = 7e6;       % Initial susceptible density from control (CFU/mL)

fprintf('\n=== Fixed Parameters (from control) ===\n');
fprintf('  r = %.4f h^-1, a_ss = %.4f, S0 = %.2e CFU/mL, K = %.2e CFU/mL\n', ...
    r_fixed, a_ss_fixed, S0_fixed, K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params(1:5) = SHARED [a_rs, a_sr, a_rr, beta, m]
% params(6:9) = epsilon for [Low, Med, MedHigh, High]
% params(10:13) = gamma for [Low, Med, MedHigh, High]
% params(14:17) = delta for [Low, Med, MedHigh, High]
% params(18:21) = phi_s for [Low, Med, MedHigh, High]
% params(22:25) = R0 for [Low, Med, MedHigh, High]

% SHARED parameter initial guesses (more moderate values)
a_rs_init = 66.2876;       % Effect of R on S - moderate
a_sr_init = 2.3421;       % Effect of S on R - moderate  
a_rr_init = 1.003;       % R intraspecific - similar to a_ss
beta_init = 38.3358;        % Burst size - moderate
m_init = 0.2733;          % Phage decay - moderate

% TREATMENT-SPECIFIC initial guesses (4 treatments: Low, Med, MedHigh, High)
eps_init = [0.0276, 0.0433, 0.0525, 0.0682];     % Resistance efficiency
gamma_init = [0.9756, 1, 0.402, 0.9547];           % Fitness cost (higher MOI = lower fitness?)
delta_init = [9.98e-4, 1e-3, 1e-3, 1e-3];       % Resistance acquisition
phi_s_init = [1e-9, 1.26e-8, 1e-9, 1.29e-8];       % Adsorption rate (increase with MOI effect)
R0_init = [7.2e6, 4.6e6, 1.4e6, 6e6];              % Initial resistant (may decrease with MOI)
% Assemble initial parameter vector
params_init = [a_rs_init; a_sr_init; a_rr_init; beta_init; m_init; ...
               eps_init(:); gamma_init(:); delta_init(:); phi_s_init(:); R0_init(:)];
% BOUNDS - More flexible to allow optimizer to find good solutions
lb_shared = [0.01;    % a_rs min
             0.01;   % a_sr min
             0.1;    % a_rr min
             5;      % beta min
             0.01];  % m min

ub_shared = [500;    % a_rs max (allow higher)
             500;     % a_sr max
             100;    % a_rr max
             100;    % beta max
             0.5];   % m max (allow higher decay)

% Treatment-specific bounds
lb_eps = 1e-5 * ones(4, 1);    ub_eps = 0.5 * ones(4, 1);
lb_gamma = 1e-5 * ones(4, 1);   ub_gamma = 1.0 * ones(4, 1);
lb_delta = 1e-6 * ones(4, 1);   ub_delta = 1e-3 * ones(4, 1);
lb_phi_s = 1e-9 * ones(4, 1);  ub_phi_s = 1e-6 * ones(4, 1);
lb_R0 = 1 * ones(4, 1);         ub_R0 = 1e7 * ones(4, 1);

% Assemble bounds
lower_bound = [lb_shared; lb_eps; lb_gamma; lb_delta; lb_phi_s; lb_R0];
upper_bound = [ub_shared; ub_eps; ub_gamma; ub_delta; ub_phi_s; ub_R0];

% Ensure initial params are within bounds
params_init = max(params_init, lower_bound);
params_init = min(params_init, upper_bound);

fprintf('\n=== Parameter Structure ===\n');
fprintf('  SHARED (5): a_rs, a_sr, a_rr, beta, m\n');
fprintf('  TREATMENT-SPECIFIC (20): epsilon, gamma, delta, phi_s, R0 × 4\n');
fprintf('  TOTAL: 25 parameters\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-Start Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_restarts = 1;  % Number of random restarts
best_params = params_init;
best_resnorm = Inf;

opts = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'MaxIterations', 100, ...
    'MaxFunctionEvaluations', 1e5, ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-10);

fprintf('\n=== Starting Multi-Start Optimization (%d restarts) ===\n', n_restarts);

tic;

for restart = 1:n_restarts
    fprintf('\n--- Restart %d/%d ---\n', restart, n_restarts);
    
    if restart == 1
        params_start = params_init;
    else
        % Random perturbation for multi-start
        perturbation = 0.5 + rand(size(best_params));
        params_start = best_params .* perturbation;
        params_start = max(params_start, lower_bound);
        params_start = min(params_start, upper_bound);
    end
    
    % Define objective function
    objective = @(params) compute_residuals(params, time, all_data_mean, ...
        P0_values, r_fixed, a_ss_fixed, S0_fixed, K);
    
    try
        [params_opt, resnorm, ~, exitflag] = ...
            lsqnonlin(objective, params_start, lower_bound, upper_bound, opts);
        
        fprintf('Restart %d: resnorm = %.4e, exitflag = %d\n', restart, resnorm, exitflag);
        
        if resnorm < best_resnorm
            best_resnorm = resnorm;
            best_params = params_opt;
            fprintf('  *** New best found! ***\n');
        end
    catch ME
        fprintf('Restart %d failed: %s\n', restart, ME.message);
    end
end

optimization_time = toc;
params_opt = best_params;

fprintf('\n=== Optimization Complete ===\n');
fprintf('Total time: %.2f seconds\n', optimization_time);
fprintf('Best resnorm: %.4e\n', best_resnorm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shared parameters
a_rs_fit = params_opt(1);
a_sr_fit = params_opt(2);
a_rr_fit = params_opt(3);
beta_fit = params_opt(4);
m_fit = params_opt(5);

% Treatment-specific parameters
eps_fit = params_opt(6:9);
gamma_fit = params_opt(10:13);
delta_fit = params_opt(14:17);
phi_s_fit = params_opt(18:21);
R0_fit = params_opt(22:25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Quality Assessment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== Model Quality ===\n');

R_squared_all = zeros(1, n_treatments);
RMSE_all = zeros(1, n_treatments);
SS_res_total = 0;
SS_tot_total = 0;
n_points_total = 0;

for i = 1:n_treatments
    [~, S_sim, R_sim, ~] = simulate_model(params_opt, time, P0_values(i), ...
        i, r_fixed, a_ss_fixed, S0_fixed, K);
    total_sim = S_sim + R_sim;
    
    data_i = all_data_mean{i};
    ss_res = sum((data_i - total_sim).^2);
    ss_tot = sum((data_i - mean(data_i)).^2);
    R_squared_all(i) = 1 - ss_res / ss_tot;
    RMSE_all(i) = sqrt(mean((data_i - total_sim).^2));
    
    SS_res_total = SS_res_total + ss_res;
    SS_tot_total = SS_tot_total + ss_tot;
    n_points_total = n_points_total + length(data_i);
    
    fprintf('  %s: R² = %.4f, RMSE = %.2e CFU/mL\n', ...
        treatment_names{i}, R_squared_all(i), RMSE_all(i));
end

% Overall R-squared
R_squared_overall = 1 - SS_res_total / SS_tot_total;
fprintf('\n  Overall R² = %.4f\n', R_squared_overall);
fprintf('  Mean R² = %.4f\n', mean(R_squared_all));

% AIC and BIC
n_params_fit = length(params_opt);
sigma2 = SS_res_total / n_points_total;
AIC = n_points_total * log(sigma2) + 2 * n_params_fit;
BIC = n_points_total * log(sigma2) + n_params_fit * log(n_points_total);
fprintf('\n  AIC = %.2f, BIC = %.2f\n', AIC, BIC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('FITTING RESULTS: Y25747 + Fern (All 4 Treatments)\n');
fprintf('%s\n', repmat('=', 1, 80));

fprintf('\n--- FIXED Parameters (from control) ---\n');
fprintf('  r = %.4f h^-1\n', r_fixed);
fprintf('  a_ss = %.4f\n', a_ss_fixed);
fprintf('  S0 = %.2e CFU/mL\n', S0_fixed);
fprintf('  K = %.2e CFU/mL\n', K);

fprintf('\n--- SHARED Parameters ---\n');
fprintf('  a_rs (effect of R on S) = %.4f\n', a_rs_fit);
fprintf('  a_sr (effect of S on R) = %.4f\n', a_sr_fit);
fprintf('  a_rr (R intraspecific)  = %.4f\n', a_rr_fit);
fprintf('  beta (burst size)       = %.2f PFU/CFU\n', beta_fit);
fprintf('  m (phage decay)         = %.4f h^-1\n', m_fit);

fprintf('\n--- TREATMENT-SPECIFIC Parameters ---\n');
fprintf('%-10s %8s %10s %10s %12s %12s %12s %8s\n', ...
    'Treatment', 'P0', 'epsilon', 'gamma', 'delta', 'phi_s', 'R0', 'R²');
fprintf('%s\n', repmat('-', 1, 90));
for i = 1:n_treatments
    fprintf('%-10s %8.1f %10.4f %10.4f %12.2e %12.2e %12.2e %8.4f\n', ...
        treatment_names{i}, P0_values(i), eps_fit(i), gamma_fit(i), ...
        delta_fit(i), phi_s_fit(i), R0_fit(i), R_squared_all(i));
end
fprintf('%s\n', repmat('-', 1, 90));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = [0.2, 0.6, 1.0;    % Blue (Low)
          0.2, 0.8, 0.2;    % Green (Medium)
          1.0, 0.5, 0.0;    % Orange (Med-High)
          0.8, 0.2, 0.2];   % Red (High)

fig1 = figure('Units', 'inches', 'Position', [1, 1, 14, 10]);

t_layout = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_layout, sprintf('25747 + Fern: Global Fit (4 Treatments)\nOverall R² = %.3f', R_squared_overall), ...
    'FontSize', 14, 'FontWeight', 'bold');

% Plot each treatment
for i = 1:n_treatments
    nexttile;
    
    % Simulate for plotting (use finer time grid)
    tspan_plot = linspace(0, max(time), 500)';
    [~, S_sim, R_sim, P_sim] = simulate_model(params_opt, tspan_plot, P0_values(i), ...
        i, r_fixed, a_ss_fixed, S0_fixed, K);
    total_sim = S_sim + R_sim;
    
    % Plot data with error bars
    errorbar(time, all_data_mean{i}, all_data_std{i}, 'o', ...
        'MarkerSize', 3, 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
        'LineWidth', 0.5, 'CapSize', 2);
    hold on;
    
    % Plot model
    plot(tspan_plot, total_sim, 'k-', 'LineWidth', 2, 'DisplayName', 'S+R');
    plot(tspan_plot, S_sim, 'b--', 'LineWidth', 1.5, 'DisplayName', 'S');
    plot(tspan_plot, R_sim, 'r:', 'LineWidth', 1.5, 'DisplayName', 'R');
    
    xlabel('Time (hours)', 'FontSize', 10);
    ylabel('CFU/mL', 'FontSize', 10);
    title(sprintf('%s (P_0=%.3f)\nR²=%.3f, \\epsilon=%.4f, \\gamma=%.2f', ...
        treatment_names{i}, P0_values(i), R_squared_all(i), eps_fit(i), gamma_fit(i)), ...
        'FontSize', 10);
    
    if i == 1
        legend('Data', 'S+R', 'S', 'R', 'Location', 'best', 'FontSize', 8);
    end
    
    grid on;
    set(gca, 'FontSize', 9);
    xlim([0, max(time)]);
end

% Panel 5: All treatments overlaid
nexttile;
plot(time, data_low_mean, '-', 'LineWidth', 1.5, 'Color', colors(1,:), 'DisplayName', 'Low');
hold on;
plot(time, data_med_mean, '-', 'LineWidth', 1.5, 'Color', colors(2,:), 'DisplayName', 'Medium');
plot(time, data_medhigh_mean, '-', 'LineWidth', 1.5, 'Color', colors(3,:), 'DisplayName', 'Med-High');
plot(time, data_high_mean, '-', 'LineWidth', 1.5, 'Color', colors(4,:), 'DisplayName', 'High');
xlabel('Time (hours)', 'FontSize', 10);
ylabel('CFU/mL', 'FontSize', 10);
title('All Treatments (Data)', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 9);

% Panel 6: Summary
nexttile;
axis off;

summary_text = {
    '\bf{Y25747 + Fern Global Fit}', ...
    '', ...
    '\bf{Fixed (from control):}', ...
    sprintf('  r = %.4f h^{-1}', r_fixed), ...
    sprintf('  a_{ss} = %.4f', a_ss_fixed), ...
    sprintf('  S_0 = %.2e', S0_fixed), ...
    '', ...
    '\bf{Shared:}', ...
    sprintf('  a_{rs} = %.4f, a_{sr} = %.4f', a_rs_fit, a_sr_fit), ...
    sprintf('  a_{rr} = %.4f', a_rr_fit), ...
    sprintf('  \\beta = %.1f, m = %.4f', beta_fit, m_fit), ...
    '', ...
    '\bf{R² by Treatment:}', ...
    sprintf('  Low: %.3f', R_squared_all(1)), ...
    sprintf('  Medium: %.3f', R_squared_all(2)), ...
    sprintf('  Med-High: %.3f', R_squared_all(3)), ...
    sprintf('  High: %.3f', R_squared_all(4)), ...
    '', ...
    sprintf('\\bf{Overall R² = %.3f}', R_squared_overall), ...
};

text(0.02, 0.98, summary_text, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 9, 'Interpreter', 'tex');
box on;
title('Summary', 'FontSize', 11);

% Save figure
exportgraphics(fig1, '25747_Fern_4treatment_fit_5shared.pdf', 'ContentType', 'vector', 'Resolution', 300);
% exportgraphics(fig1, 'Y25747_Fern_4treatment_fit.png', 'Resolution', 300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Trends Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig2 = figure('Units', 'inches', 'Position', [1, 1, 12, 8]);

t_layout2 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_layout2, 'Y25747 + Fern: Treatment-Specific Parameters vs MOI', ...
    'FontSize', 14, 'FontWeight', 'bold');

params_to_plot = {eps_fit, gamma_fit, delta_fit, phi_s_fit, R0_fit};
param_labels = {'\epsilon', '\gamma', '\delta', '\phi_s', 'R_0'};
param_titles = {'Resistance Efficiency', 'Fitness Cost', 'Resistance Acquisition', ...
    'Adsorption Rate', 'Initial Resistant'};

for i = 1:5
    nexttile;
    loglog(P0_values, params_to_plot{i}, 'ko-', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'b', 'LineWidth', 2);
    xlabel('P_0 (PFU/mL)', 'FontSize', 10);
    ylabel(param_labels{i}, 'FontSize', 10);
    title(param_titles{i}, 'FontSize', 11);
    grid on;
    set(gca, 'FontSize', 9);
    
    % Add trend line
    hold on;
    x_log = log10(P0_values);
    y_log = log10(params_to_plot{i});
    p = polyfit(x_log, y_log, 1);
    x_fit = linspace(min(x_log), max(x_log), 100);
    y_fit = polyval(p, x_fit);
    plot(10.^x_fit, 10.^y_fit, 'r--', 'LineWidth', 1);
    text(0.05, 0.95, sprintf('slope = %.2f', p(1)), 'Units', 'normalized', ...
        'FontSize', 8, 'VerticalAlignment', 'top');
end

% R-squared (linear y-axis)
nexttile;
semilogx(P0_values, R_squared_all, 'ko-', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'g', 'LineWidth', 2);
xlabel('P_0 (PFU/mL)', 'FontSize', 10);
ylabel('R^2', 'FontSize', 10);
title('Fit Quality', 'FontSize', 11);
ylim([0, 1]);
yline(0.9, 'r--', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 9);

% Save figure
exportgraphics(fig2, '25747_Fern_4treatment_params_5shared.pdf', 'ContentType', 'vector', 'Resolution', 300);
% exportgraphics(fig2, 'Y25747_Fern_4treatment_params.png', 'Resolution', 300);

fprintf('\nFigures saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results.strain = '25747';
results.phage = 'Fern_IDv1';

results.fixed.r = r_fixed;
results.fixed.a_ss = a_ss_fixed;
results.fixed.S0 = S0_fixed;
results.fixed.K = K;

results.shared.a_rs = a_rs_fit;
results.shared.a_sr = a_sr_fit;
results.shared.a_rr = a_rr_fit;
results.shared.beta = beta_fit;
results.shared.m = m_fit;

results.treatment_specific.epsilon = eps_fit;
results.treatment_specific.gamma = gamma_fit;
results.treatment_specific.delta = delta_fit;
results.treatment_specific.phi_s = phi_s_fit;
results.treatment_specific.R0 = R0_fit;

results.treatment_names = treatment_names;
results.P0_values = P0_values;
results.R_squared = R_squared_all;
results.R_squared_overall = R_squared_overall;
results.RMSE = RMSE_all;
results.AIC = AIC;
results.BIC = BIC;
results.note = 'Global fit with revised model (epsilon for partial resistance)';

save('fitting_results_4treatments_25747F_5shared.mat', 'results');
fprintf('Results saved to: fitting_results_4treatments_25747F_5shared.mat\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manuscript Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s\n', repmat('=', 1, 90));
fprintf('TABLE FOR MANUSCRIPT: 25747 + Fern IDv1 Global Fitting Results\n');
fprintf('%s\n', repmat('=', 1, 90));

fprintf('\nTable: Fixed parameters (from control)\n');
fprintf('%-15s %15s %20s\n', 'Parameter', 'Value', 'Units');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%-15s %15.4f %20s\n', 'r', r_fixed, 'h^-1');
fprintf('%-15s %15.4f %20s\n', 'a_ss', a_ss_fixed, 'dimensionless');
fprintf('%-15s %15.2e %20s\n', 'S0', S0_fixed, 'CFU/mL');
fprintf('%-15s %15.2e %20s\n', 'K', K, 'CFU/mL');

fprintf('\nTable: Shared parameters (across all MOI treatments)\n');
fprintf('%-15s %15s %20s\n', 'Parameter', 'Value', 'Units');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%-15s %15.4f %20s\n', 'a_rs', a_rs_fit, 'dimensionless');
fprintf('%-15s %15.4f %20s\n', 'a_sr', a_sr_fit, 'dimensionless');
fprintf('%-15s %15.4f %20s\n', 'a_rr', a_rr_fit, 'dimensionless');
fprintf('%-15s %15.2f %20s\n', 'beta', beta_fit, 'PFU/CFU');
fprintf('%-15s %15.4f %20s\n', 'm', m_fit, 'h^-1');

fprintf('\nTable: Treatment-specific parameters\n');
fprintf('%-10s | %8s | %10s | %10s | %12s | %12s | %12s | %6s\n', ...
    'Treatment', 'P0', 'epsilon', 'gamma', 'delta', 'phi_s', 'R0', 'R^2');
fprintf('%s\n', repmat('-', 1, 90));
for i = 1:n_treatments
    fprintf('%-10s | %8.1f | %10.4f | %10.4f | %12.2e | %12.2e | %12.2e | %6.3f\n', ...
        treatment_names{i}, P0_values(i), eps_fit(i), gamma_fit(i), ...
        delta_fit(i), phi_s_fit(i), R0_fit(i), R_squared_all(i));
end
fprintf('%s\n', repmat('-', 1, 90));
fprintf('Overall R²: %.4f\n', R_squared_overall);
fprintf('AIC: %.2f, BIC: %.2f\n', AIC, BIC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function residuals = compute_residuals(params, time, all_data_mean, ...
    P0_values, r_fixed, a_ss_fixed, S0_fixed, K)
    %
    % Compute residuals for ALL treatments simultaneously
    % This ensures shared parameters are truly shared across MOI treatments
    %
    
    n_treatments = length(all_data_mean);
    n_timepoints = length(time);
    residuals = zeros(n_treatments * n_timepoints, 1);
    
    for i = 1:n_treatments
        [~, S_sim, R_sim, ~] = simulate_model(params, time, P0_values(i), ...
            i, r_fixed, a_ss_fixed, S0_fixed, K);
        total_sim = S_sim + R_sim;
        
        idx_start = (i-1) * n_timepoints + 1;
        idx_end = i * n_timepoints;
        residuals(idx_start:idx_end) = total_sim - all_data_mean{i};
    end
end

function [t_out, S_out, R_out, P_out] = simulate_model(params, time, P0, ...
    treatment_idx, r_fixed, a_ss_fixed, S0_fixed, K)
    %
    % Simulate the REVISED phage-bacteria model
    %
    % REVISED MODEL (from BPI_revision.pdf):
    %   dS/dt = r*S*(1 - (a_ss*S + a_rs*R)/K) - delta*S - phi_s*P*S
    %   dR/dt = gamma*r*R*(1 - (a_sr*S + a_rr*R)/K) + delta*S - epsilon*phi_s*P*R
    %   dP/dt = beta*phi_s*(S + epsilon*R)*P - m*P
    %
    % Parameter indexing:
    % params(1:5) = [a_rs, a_sr, a_rr, beta, m]  (SHARED)
    % params(6:9) = epsilon for 4 treatments
    % params(10:13) = gamma for 4 treatments
    % params(14:17) = delta for 4 treatments
    % params(18:21) = phi_s for 4 treatments
    % params(22:25) = R0 for 4 treatments
    
    % Extract SHARED parameters
    a_rs = params(1);
    a_sr = params(2);
    a_rr = params(3);
    beta = params(4);
    m = params(5);
    
    % Extract TREATMENT-SPECIFIC parameters
    epsilon = params(5 + treatment_idx);      % params(6:9)
    gamma = params(9 + treatment_idx);        % params(10:13)
    delta = params(13 + treatment_idx);       % params(14:17)
    phi_s = params(17 + treatment_idx);       % params(18:21)
    R0 = params(21 + treatment_idx);          % params(22:25)
    
    % Compute adsorption rate to resistant bacteria
    phi_r = epsilon * phi_s;
    
    % Initial conditions
    y0 = [S0_fixed; R0; P0];
    
    % ODE options
    ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NonNegative', [1, 2, 3]);
    
    % Define ODE system (REVISED MODEL)
    odefun = @(t, y) [
        % dS/dt = r*S*(1 - (a_ss*S + a_rs*R)/K) - delta*S - phi_s*P*S
        r_fixed * y(1) * (1 - (a_ss_fixed * y(1) + a_rs * y(2)) / K) ...
            - delta * y(1) - phi_s * y(3) * y(1);
        
        % dR/dt = gamma*r*R*(1 - (a_sr*S + a_rr*R)/K) + delta*S - epsilon*phi_s*P*R
        gamma * r_fixed * y(2) * (1 - (a_sr * y(1) + a_rr * y(2)) / K) ...
            + delta * y(1) - phi_r * y(3) * y(2);
        
        % dP/dt = beta*phi_s*(S + epsilon*R)*P - m*P
        beta * phi_s * (y(1) + epsilon * y(2)) * y(3) - m * y(3)
    ];
    
    % Solve ODE
    try
        [t_out, Y] = ode15s(odefun, time, y0, ode_opts);
        S_out = Y(:, 1);
        R_out = Y(:, 2);
        P_out = Y(:, 3);
    catch
        % If ODE solver fails, return large values to penalize
        t_out = time;
        S_out = 1e15 * ones(size(time));
        R_out = 1e15 * ones(size(time));
        P_out = 1e15 * ones(size(time));
    end
end