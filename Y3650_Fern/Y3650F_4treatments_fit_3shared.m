%% Global Fitting Including All 4 Treatments
%
% This script fits the model to ALL 4 phage treatments simultaneously
% (Low, Medium, Med-High, High).
%
% FIXED (from control): r, a_ss, S0, K
% SHARED (3): a_rr, beta, m
% TREATMENT-SPECIFIC (7 × 4 = 28): a_rs, a_sr, epsilon, gamma, delta, phi_s, R0
%
% Total: 3 shared + 28 treatment-specific = 31 parameters

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

excelFileName = 'dataset_Y3650F.xlsx';
sheetName = 'Y3650_Fern_Grouped';

data_table = readtable(excelFileName, 'Sheet', sheetName);
time = data_table.Time / 60;  % hours

OD_to_CFU = 1e9;

data_low_mean    = data_table.Low_Group_Mean         * OD_to_CFU;
data_low_std     = data_table.Low_Group_Std          * OD_to_CFU;

data_med_mean    = data_table.Medium_Group_Mean      * OD_to_CFU;
data_med_std     = data_table.Medium_Group_Std       * OD_to_CFU;

data_medhigh_mean = data_table.Medium_High_Group_Mean * OD_to_CFU;
data_medhigh_std  = data_table.Medium_High_Group_Std  * OD_to_CFU;

data_high_mean   = data_table.High_Group_Mean        * OD_to_CFU;
data_high_std    = data_table.High_Group_Std         * OD_to_CFU;

% P0 values for all 4 treatments: Low, Medium, Med-High, High
P0_values = [0.347, 3.47, 34.7, 347];

all_data_mean  = {data_low_mean, data_med_mean, data_medhigh_mean, data_high_mean};
all_data_std   = {data_low_std,  data_med_std,  data_medhigh_std,  data_high_std};
treatment_names = {'Low', 'Medium', 'Med-High', 'High'};
n_treatments = 4;

fprintf('=== Fitting ALL 4 Treatments (3 Shared Parameters) ===\n');
fprintf('Treatments: Low (P0=0.347), Medium (P0=3.47), Med-High (P0=34.7), High (P0=347)\n');
fprintf('Data points per treatment: %d\n', length(time));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed Parameters (from control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K          = 1e9;
r_fixed    = 0.4699;
a_ss_fixed = 0.9411;
S0_fixed   = 9e6;

fprintf('\n=== Fixed Parameters ===\n');
fprintf('  r = %.4f, a_ss = %.4f, S0 = %.2e, K = %.2e\n', ...
    r_fixed, a_ss_fixed, S0_fixed, K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params(1)   = a_rr   (SHARED)
% params(2)   = beta   (SHARED)
% params(3)   = m      (SHARED)
%
% params(4:7)   = a_rs  for [Low, Med, MedHigh, High]
% params(8:11)  = a_sr  for [Low, Med, MedHigh, High]
% params(12:15) = epsilon for [Low, Med, MedHigh, High]
% params(16:19) = gamma for [Low, Med, MedHigh, High]
% params(20:23) = delta for [Low, Med, MedHigh, High]
% params(24:27) = phi_s for [Low, Med, MedHigh, High]
% params(28:31) = R0    for [Low, Med, MedHigh, High]

% --- Initial values: SHARED ---
a_rr_init  = 0.9231;
beta_init  = 18.77;
m_init     = 0.48;

% --- Initial values: TREATMENT-SPECIFIC ---
a_rs_init  = [4.4942, 162.534, 271.6803, 8.8929];   % effect of R on S
a_sr_init  = [0.0338,  0.3356, 232.7256, 5.1228];    % effect of S on R
eps_init   = [0.0370,  0.0129, 0.1292,   0.0088];    % resistance efficiency
gamma_init = [0.9988,  0.6206,  0.7700,  0.3378];    % fitness cost
delta_init = [1e-3,    6.5683e-4, 1.3426e-5,  2.4489e-4];      % resistance acquisition
phi_s_init = [1e-9,    1.7951e-9, 2.676e-9,  3.7685e-9];   % adsorption rate
R0_init    = [3.4569e6, 1.8634e4, 1.4508e6,  6.1508e6];       % initial resistant density

params_init = [a_rr_init; beta_init; m_init; ...
               a_rs_init(:); a_sr_init(:); eps_init(:); gamma_init(:); ...
               delta_init(:); phi_s_init(:); R0_init(:)];

% --- Bounds ---
lb_shared = [1e-2; 10;   1.5e-3];   % a_rr, beta, m
ub_shared = [100;  100;  0.5];

lb_a_rs  = 1e-2 * ones(4,1);   ub_a_rs  = 500  * ones(4,1);
lb_a_sr  = 1e-2 * ones(4,1);   ub_a_sr  = 500  * ones(4,1);
lb_eps   = 1e-4 * ones(4,1);   ub_eps   = 0.5  * ones(4,1);
lb_gamma = 1e-4 * ones(4,1);   ub_gamma = 1.0  * ones(4,1);
lb_delta = 1e-6 * ones(4,1);   ub_delta = 1e-3 * ones(4,1);
lb_phi_s = 1e-9 * ones(4,1);   ub_phi_s = 1e-6 * ones(4,1);
lb_R0    = 1    * ones(4,1);   ub_R0    = 1e7  * ones(4,1);

lower_bound = [lb_shared; lb_a_rs; lb_a_sr; lb_eps; lb_gamma; lb_delta; lb_phi_s; lb_R0];
upper_bound = [ub_shared; ub_a_rs; ub_a_sr; ub_eps; ub_gamma; ub_delta; ub_phi_s; ub_R0];

params_init = max(params_init, lower_bound);
params_init = min(params_init, upper_bound);

fprintf('\n=== Parameter Structure ===\n');
fprintf('  SHARED (3):              a_rr, beta, m\n');
fprintf('  TREATMENT-SPECIFIC (28): a_rs, a_sr, epsilon, gamma, delta, phi_s, R0 x4\n');
fprintf('  TOTAL: 31 parameters\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-Start Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_restarts = 1;
best_params  = params_init;
best_resnorm = Inf;

opts = optimoptions('lsqnonlin', ...
    'Display',               'iter', ...
    'MaxIterations',         1, ...
    'MaxFunctionEvaluations', 1e5, ...
    'FunctionTolerance',     1e-10, ...
    'StepTolerance',         1e-10);

fprintf('\n=== Starting Optimization ===\n');

tic;

for restart = 1:n_restarts
    fprintf('\n--- Restart %d/%d ---\n', restart, n_restarts);
    
    if restart == 1
        params_start = params_init;
    else
        perturbation  = 0.5 + rand(size(best_params));
        params_start  = best_params .* perturbation;
        params_start  = max(params_start, lower_bound);
        params_start  = min(params_start, upper_bound);
    end
    
    objective = @(params) compute_residuals(params, time, all_data_mean, ...
        P0_values, r_fixed, a_ss_fixed, S0_fixed, K);
    
    try
        [params_opt, resnorm, ~, exitflag] = ...
            lsqnonlin(objective, params_start, lower_bound, upper_bound, opts);
        
        fprintf('Restart %d: resnorm = %.4e, exitflag = %d\n', restart, resnorm, exitflag);
        
        if resnorm < best_resnorm
            best_resnorm = resnorm;
            best_params  = params_opt;
            fprintf('  *** New best! ***\n');
        end
    catch ME
        fprintf('Restart %d failed: %s\n', restart, ME.message);
    end
end

optimization_time = toc;
params_opt = best_params;

fprintf('\n=== Optimization Complete ===\n');
fprintf('Time: %.2f seconds\n', optimization_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_rr_fit  = params_opt(1);
beta_fit  = params_opt(2);
m_fit     = params_opt(3);

a_rs_fit  = params_opt(4:7);
a_sr_fit  = params_opt(8:11);
eps_fit   = params_opt(12:15);
gamma_fit = params_opt(16:19);
delta_fit = params_opt(20:23);
phi_s_fit = params_opt(24:27);
R0_fit    = params_opt(28:31);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Quality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== Model Quality ===\n');

R_squared_all = zeros(1, n_treatments);
RMSE_all      = zeros(1, n_treatments);

for i = 1:n_treatments
    [~, S_sim, R_sim, ~] = simulate_model(params_opt, time, P0_values(i), ...
        i, r_fixed, a_ss_fixed, S0_fixed, K);
    total_sim = S_sim + R_sim;
    
    data_i   = all_data_mean{i};
    ss_res   = sum((data_i - total_sim).^2);
    ss_tot   = sum((data_i - mean(data_i)).^2);
    R_squared_all(i) = 1 - ss_res / ss_tot;
    RMSE_all(i)      = sqrt(mean((data_i - total_sim).^2));
    
    fprintf('  %s: R² = %.4f, RMSE = %.2e\n', ...
        treatment_names{i}, R_squared_all(i), RMSE_all(i));
end

fprintf('\n  Mean R² = %.4f\n', mean(R_squared_all));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information Criteria (AIC and BIC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate total residual sum of squares across all treatments
total_rss = 0;
total_n   = 0;
for i = 1:n_treatments
    [~, S_sim, R_sim, ~] = simulate_model(params_opt, time, P0_values(i), ...
        i, r_fixed, a_ss_fixed, S0_fixed, K);
    total_sim = S_sim + R_sim;
    data_i = all_data_mean{i};
    rss_i = sum((data_i - total_sim).^2);
    total_rss = total_rss + rss_i;
    total_n = total_n + length(data_i);
end

% Number of fitted parameters
k_params = 31;  % 3 shared + 28 treatment-specific

% Calculate AIC and BIC
% AIC = 2k + n*ln(RSS/n)
% BIC = k*ln(n) + n*ln(RSS/n)
log_likelihood = -0.5 * total_n * log(total_rss / total_n);
AIC = 2 * k_params - 2 * log_likelihood;
BIC = k_params * log(total_n) - 2 * log_likelihood;

fprintf('\n=== Information Criteria ===\n');
fprintf('  Total data points: %d\n', total_n);
fprintf('  Number of parameters: %d\n', k_params);
fprintf('  Total RSS: %.6e\n', total_rss);
fprintf('  Log-likelihood: %.6f\n', log_likelihood);
fprintf('  AIC: %.6f\n', AIC);
fprintf('  BIC: %.6f\n', BIC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('FITTING RESULTS (All 4 Treatments, 3 Shared)\n');
fprintf('%s\n', repmat('=', 1, 70));

fprintf('\n--- SHARED Parameters ---\n');
fprintf('  a_rr = %.4f, beta = %.2f, m = %.4f\n', a_rr_fit, beta_fit, m_fit);

fprintf('\n--- TREATMENT-SPECIFIC Parameters ---\n');
fprintf('%-10s %8s %8s %8s %8s %10s %12s %12s %10s\n', ...
    'Treatment', 'P0', 'a_rs', 'a_sr', 'epsilon', 'gamma', 'delta', 'phi_s', 'R0');
fprintf('%s\n', repmat('-', 1, 100));
for i = 1:n_treatments
    fprintf('%-10s %8.2f %8.4f %8.4f %8.4f %8.4f %10.2e %12.2e %12.2e\n', ...
        treatment_names{i}, P0_values(i), a_rs_fit(i), a_sr_fit(i), ...
        eps_fit(i), gamma_fit(i), delta_fit(i), phi_s_fit(i), R0_fit(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting: Fitted Trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure('Units', 'inches', 'Position', [1, 1, 14, 10]);

t_layout1 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_layout1, 'Y3650 Fern: Global Fitting Results (All 4 Treatments, 3 Shared)', ...
    'FontSize', 14, 'FontWeight', 'bold');

colors = lines(n_treatments);
tspan_plot = linspace(min(time), max(time), 500)';

for i = 1:n_treatments
    nexttile;
    
    [~, S_sim, R_sim, ~] = simulate_model(params_opt, tspan_plot, P0_values(i), ...
        i, r_fixed, a_ss_fixed, S0_fixed, K);
    total_sim = S_sim + R_sim;
    
    % Recompute R² on original time grid for label
    [~, S_q, R_q, ~] = simulate_model(params_opt, time, P0_values(i), ...
        i, r_fixed, a_ss_fixed, S0_fixed, K);
    data_i  = all_data_mean{i};
    ss_res  = sum((data_i - (S_q + R_q)).^2);
    ss_tot  = sum((data_i - mean(data_i)).^2);
    r2_label = 1 - ss_res / ss_tot;
    
    errorbar(time, all_data_mean{i}, all_data_std{i}, 'o', ...
        'MarkerSize', 3, 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
        'LineWidth', 0.5, 'CapSize', 2);
    hold on;
    plot(tspan_plot, total_sim, 'k-',  'LineWidth', 2);
    plot(tspan_plot, S_sim,     'b--', 'LineWidth', 1.5);
    plot(tspan_plot, R_sim,     'r:',  'LineWidth', 1.5);
    
    xlabel('Time (hours)', 'FontSize', 10);
    ylabel('CFU/mL',       'FontSize', 10);
    title(sprintf('%s (P_0=%.3f)\nR^2=%.3f', ...
        treatment_names{i}, P0_values(i), r2_label), 'FontSize', 11);
    
    if i == 1
        legend('Data', 'S+R', 'S', 'R', 'Location', 'best', 'FontSize', 8);
    end
    
    grid on;
    set(gca, 'FontSize', 9);
    xlim([0, max(time)]);
end

% Panel 5: All treatments overlaid
nexttile;
plot(time, data_low_mean,     'b-', 'LineWidth', 1.5, 'DisplayName', 'Low');
hold on;
plot(time, data_med_mean,     'g-', 'LineWidth', 1.5, 'DisplayName', 'Medium');
plot(time, data_medhigh_mean, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Med-High');
plot(time, data_high_mean,    'm-', 'LineWidth', 1.5, 'DisplayName', 'High');
xlabel('Time (hours)', 'FontSize', 10);
ylabel('CFU/mL',       'FontSize', 10);
title('All Treatments Comparison', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 9);

% Panel 6: Summary text
nexttile;
axis off;

summary_text = {
    '\bf{Results (All 4 Treatments)}', ...
    '', ...
    '\bf{Fixed:}', ...
    sprintf('  r=%.4f, a_{ss}=%.4f', r_fixed, a_ss_fixed), ...
    sprintf('  S_0=%.2e', S0_fixed), ...
    '', ...
    '\bf{Shared:}', ...
    sprintf('  a_{rr}=%.2f, \\beta=%.1f, m=%.3f', a_rr_fit, beta_fit, m_fit), ...
    '', ...
    '\bf{R^2 Values:}', ...
    sprintf('  Low: %.3f',      R_squared_all(1)), ...
    sprintf('  Medium: %.3f',   R_squared_all(2)), ...
    sprintf('  Med-High: %.3f', R_squared_all(3)), ...
    sprintf('  High: %.3f',     R_squared_all(4)), ...
    sprintf('  \\bf{Mean: %.3f}', mean(R_squared_all)), ...
};

text(0.02, 0.98, summary_text, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 9, 'Interpreter', 'tex');
box on;
title('Summary', 'FontSize', 11);

exportgraphics(fig1, 'Y3650_Fern_4treatment_fit_3shared.pdf', ...
    'ContentType', 'vector', 'Resolution', 300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting: Treatment-Specific Parameter Trends
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig2 = figure('Units', 'inches', 'Position', [1, 1, 14, 8]);

t_layout2 = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_layout2, 'Y3650 Fern: Treatment-Specific Parameters vs MOI (3 Shared)', ...
    'FontSize', 14, 'FontWeight', 'bold');

params_to_plot = {a_rs_fit, a_sr_fit, eps_fit, gamma_fit, delta_fit, phi_s_fit, R0_fit};
param_labels   = {'a_{rs}', 'a_{sr}', '\epsilon', '\gamma', '\delta', '\phi_s', 'R_0'};
param_titles   = {'Interspecific (R→S)', 'Interspecific (S→R)', ...
                  'Resistance Efficiency', 'Fitness Cost', ...
                  'Resistance Acquisition', 'Adsorption Rate', 'Initial Resistant'};

for i = 1:7
    nexttile;
    loglog(P0_values, params_to_plot{i}, 'ko-', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'b', 'LineWidth', 2);
    xlabel('P_0 (PFU/mL)',   'FontSize', 10);
    ylabel(param_labels{i},  'FontSize', 10);
    title(param_titles{i},   'FontSize', 11);
    grid on;
    set(gca, 'FontSize', 9);
end

% Panel 8: R-squared
nexttile;
semilogx(P0_values, R_squared_all, 'ko-', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'g', 'LineWidth', 2);
xlabel('P_0 (PFU/mL)', 'FontSize', 10);
ylabel('R^2',          'FontSize', 10);
title('Fit Quality',   'FontSize', 11);
ylim([0, 1]);
yline(0.9, 'r--', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 9);

exportgraphics(fig2, 'Y3650_Fern_4treatment_params_3shared.pdf', ...
    'ContentType', 'vector', 'Resolution', 300);

fprintf('\nFigures saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results.fixed.r    = r_fixed;
results.fixed.a_ss = a_ss_fixed;
results.fixed.S0   = S0_fixed;
results.fixed.K    = K;

results.shared.a_rr = a_rr_fit;
results.shared.beta = beta_fit;
results.shared.m    = m_fit;

results.treatment_specific.a_rs   = a_rs_fit;
results.treatment_specific.a_sr   = a_sr_fit;
results.treatment_specific.epsilon = eps_fit;
results.treatment_specific.gamma  = gamma_fit;
results.treatment_specific.delta  = delta_fit;
results.treatment_specific.phi_s  = phi_s_fit;
results.treatment_specific.R0     = R0_fit;

results.treatment_names = treatment_names;
results.P0_values       = P0_values;
results.R_squared       = R_squared_all;
results.RMSE            = RMSE_all;
results.AIC             = AIC;
results.BIC             = BIC;
results.note            = 'All 4 treatments; 3 shared (a_rr, beta, m); 7 treatment-specific';

save('fitting_results_4treatments_Y3650F_3shared.mat', 'results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manuscript Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s\n', repmat('=', 1, 100));
fprintf('TABLE FOR MANUSCRIPT\n');
fprintf('%s\n', repmat('=', 1, 100));
fprintf('\nShared parameters: a_rr=%.4f, beta=%.2f, m=%.4f\n\n', ...
    a_rr_fit, beta_fit, m_fit);
fprintf('Table X: Treatment-specific parameters (Y-3650 + Fern) - All 4 Treatments\n\n');
fprintf('%-10s | %8s | %8s | %8s | %8s | %8s | %10s | %12s | %12s | %6s\n', ...
    'Treatment', 'P0', 'a_rs', 'a_sr', 'epsilon', 'gamma', 'delta', 'phi_s', 'R0', 'R^2');
fprintf('%s\n', repmat('-', 1, 110));
for i = 1:n_treatments
    fprintf('%-10s | %8.2f | %8.4f | %8.4f | %8.4f | %8.4f | %10.2e | %12.2e | %12.2e | %6.3f\n', ...
        treatment_names{i}, P0_values(i), a_rs_fit(i), a_sr_fit(i), ...
        eps_fit(i), gamma_fit(i), delta_fit(i), phi_s_fit(i), R0_fit(i), R_squared_all(i));
end
fprintf('%s\n', repmat('-', 1, 110));
fprintf('Mean R^2: %.3f\n', mean(R_squared_all));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function residuals = compute_residuals(params, time, all_data_mean, ...
    P0_values, r_fixed, a_ss_fixed, S0_fixed, K)

    n_treatments = length(all_data_mean);
    n_timepoints = length(time);
    residuals    = zeros(n_treatments * n_timepoints, 1);

    for i = 1:n_treatments
        [~, S_sim, R_sim, ~] = simulate_model(params, time, P0_values(i), ...
            i, r_fixed, a_ss_fixed, S0_fixed, K);
        total_sim = S_sim + R_sim;

        idx_start = (i-1) * n_timepoints + 1;
        idx_end   = i     * n_timepoints;
        residuals(idx_start:idx_end) = total_sim - all_data_mean{i};
    end
end

function [t_out, S_out, R_out, P_out] = simulate_model(params, time, P0, ...
    treatment_idx, r_fixed, a_ss_fixed, S0_fixed, K)

    % params(1)     = a_rr   [SHARED]
    % params(2)     = beta   [SHARED]
    % params(3)     = m      [SHARED]
    % params(4:7)   = a_rs   [treatment-specific]
    % params(8:11)  = a_sr   [treatment-specific]
    % params(12:15) = epsilon [treatment-specific]
    % params(16:19) = gamma  [treatment-specific]
    % params(20:23) = delta  [treatment-specific]
    % params(24:27) = phi_s  [treatment-specific]
    % params(28:31) = R0     [treatment-specific]

    a_rr    = params(1);
    beta    = params(2);
    m       = params(3);

    a_rs    = params(3  + treatment_idx);   % indices 4-7
    a_sr    = params(7  + treatment_idx);   % indices 8-11
    epsilon = params(11 + treatment_idx);   % indices 12-15
    gamma   = params(15 + treatment_idx);   % indices 16-19
    delta   = params(19 + treatment_idx);   % indices 20-23
    phi_s   = params(23 + treatment_idx);   % indices 24-27
    R0      = params(27 + treatment_idx);   % indices 28-31

    phi_r = epsilon * phi_s;
    y0    = [S0_fixed; R0; P0];

    ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NonNegative', [1, 2, 3]);

    odefun = @(t, y) [
        r_fixed * y(1) * (1 - (a_ss_fixed * y(1) + a_rs * y(2)) / K) ...
            - delta * y(1) - phi_s * y(3) * y(1);
        gamma * r_fixed * y(2) * (1 - (a_sr * y(1) + a_rr * y(2)) / K) ...
            + delta * y(1) - phi_r * y(3) * y(2);
        beta * y(3) * (phi_s * y(1) + phi_r * y(2)) - m * y(3)
    ];

    try
        [t_out, Y] = ode45(odefun, time, y0, ode_opts);
        S_out = Y(:, 1);
        R_out = Y(:, 2);
        P_out = Y(:, 3);
    catch
        t_out = time;
        S_out = 1e15 * ones(size(time));
        R_out = 1e15 * ones(size(time));
        P_out = 1e15 * ones(size(time));
    end
end