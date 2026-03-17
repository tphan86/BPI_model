%% Fit Control Data for PA103 Rd3
% This script fits control data (no phage) to estimate:
%   - r: intrinsic growth rate of susceptible bacteria
%   - a_ss: intraspecific competition coefficient
%
% Model for control (P=0, R0=0):
%   dS/dt = r*S*(1 - a_ss*S/K) - delta*S
%
% Since delta << r, we approximate:
%   dS/dt = r*S*(1 - a_ss*S/K)   [Logistic growth]
%
% Conversion: 1 OD600 ≈ 1E9 CFU/mL
% Fixed: K = 1E9 CFU/mL

clear; clc; close all;

% Note: For actual fitting, load the full dataset from Excel:
data_table = readtable('dataset_PA103_Rd3.xlsx', 'Sheet', 'PA103_Rd3_Grouped');
time_hr = data_table.Time; % time in minutes
od_mean = data_table.Control_Group; % OD600
% od_std = data_table.Ctrl_Group_Std;   % OD600

% Convert OD600 to CFU/mL (1 OD600 ≈ 1E9 CFU/mL)
OD_to_CFU = 1e9;            % Conversion factor
S_data = od_mean * OD_to_CFU;   % CFU/mL
% S_std = od_std * OD_to_CFU;     % Std in CFU/mL

%% Fixed Parameters
K = 1e9;                    % Carrying capacity (CFU/mL)
S0 = 6e5;   % can be estimated using S0 = S(t)/exp(r*t)                 

%% Model Definition
% Logistic growth model: dS/dt = r*S*(1 - a_ss*S/K)
% Parameters to fit: [r, a_ss]

ode_model = @(t, S, params) logistic_growth(t, S, params, K);

%% Initial Parameter Guesses
% From literature: r ∈ [0.462, 1.65] h^-1
% From data: bacteria reach ~1E9 (OD~1) at ~15-20 hours
% Estimate r from exponential phase (first few hours)
% ln(S/S0) = r*t during exponential phase

% Use early exponential phase data (t < 5 hours) for initial r estimate
idx_exp = time_hr < 8;
if sum(idx_exp) > 2
    p = polyfit(time_hr(idx_exp), log(S_data(idx_exp)/S0), 1);
    r_init = p(1);
else
    r_init = 0.5;  % Default guess
end

% Estimate a_ss from stationary phase
% At stationary phase: S_max ≈ K/a_ss
S_max = max(S_data);
a_ss_init = K / S_max;

fprintf('Initial parameter estimates:\n');
fprintf('  r_init = %.4f h^-1\n', r_init);
fprintf('  a_ss_init = %.4f (dimensionless)\n\n', a_ss_init);

params_init = [r_init, a_ss_init];

%% Parameter Bounds
% r: [0.1, 3.0] h^-1 (reasonable growth rates)
% a_ss: [0.1, 10] (dimensionless, allows flexibility around 1)
lb = [0.462, 0.1];
ub = [1.65, 10.0];

%% Fitting Options
options = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-10);

%% Objective Function
% Minimize sum of squared residuals between model and data
objective = @(params) compute_residuals(params, time_hr, S_data, S0, K);

%% Run Optimization
fprintf('Starting parameter optimization...\n');
fprintf('%s\n', repmat('=', 1, 50));

[params_opt, resnorm, residuals, exitflag, output, lambda, jacobian] = ...
    lsqnonlin(objective, params_init, lb, ub, options);

%% Extract Fitted Parameters
r_fit = params_opt(1);
a_ss_fit = params_opt(2);

fprintf('\n');
fprintf('%s\n', repmat('=', 1, 50));
fprintf('FITTING RESULTS\n');
fprintf('%s\n', repmat('=', 1, 50));
fprintf('Fitted parameters:\n');
fprintf('  r = %.4f h^-1 (intrinsic growth rate)\n', r_fit);
fprintf('  a_ss = %.4f (intraspecific competition coefficient)\n', a_ss_fit);
fprintf('\nDerived quantities:\n');
fprintf('  Effective carrying capacity K/a_ss = %.2e CFU/mL\n', K/a_ss_fit);
fprintf('  Doubling time = %.2f hours\n', log(2)/r_fit);
fprintf('\nFit quality:\n');
fprintf('  Sum of squared residuals: %.4e\n', resnorm);
fprintf('  Exit flag: %d\n', exitflag);

%% Compute Confidence Intervals (approximate)
% Using Jacobian from lsqnonlin
try
    % Estimate variance of residuals
    n = length(S_data);
    np = length(params_opt);
    mse = resnorm / (n - np);

    % Covariance matrix
    J = jacobian;
    cov_matrix = mse * inv(J' * J);

    % Standard errors
    se = sqrt(diag(cov_matrix));

    % 95% confidence intervals
    t_crit = tinv(0.975, n - np);
    ci_r = [r_fit - t_crit*se(1), r_fit + t_crit*se(1)];
    ci_ass = [a_ss_fit - t_crit*se(2), a_ss_fit + t_crit*se(2)];

    fprintf('\n95%% Confidence Intervals:\n');
    fprintf('  r: [%.4f, %.4f] h^-1\n', ci_r(1), ci_r(2));
    fprintf('  a_ss: [%.4f, %.4f]\n', ci_ass(1), ci_ass(2));
catch
    fprintf('\nNote: Could not compute confidence intervals.\n');
    ci_r = [NaN, NaN];
    ci_ass = [NaN, NaN];
end

%% Compute R-squared
SS_res = resnorm;
SS_tot = sum((S_data - mean(S_data)).^2);
R_squared = 1 - SS_res / SS_tot;
fprintf('\nR-squared: %.4f\n', R_squared);

%% Simulate Model with Fitted Parameters
tspan = [2, max(time_hr)];
[t_sim, S_sim] = ode45(@(t, S) logistic_growth(t, S, params_opt, K), tspan, S0);

%% Plot Results using tiledlayout
fig = figure('Units', 'inches', 'Position', [1, 1, 10, 8], ...
    'PaperUnits', 'inches', 'PaperSize', [10, 8], 'PaperPosition', [0, 0, 10, 8]);

% Create tiled layout with better spacing
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Title for entire figure
title(t, 'PA103 Rd3 Control Data Fitting (No Phage)', ...
    'FontSize', 18, 'FontWeight', 'bold');

% Panel A: CFU/mL scale
nexttile;
ax = gca;
plot(time_hr, S_data, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'Control data', 'LineWidth', 1);
hold on;
plot(t_sim, S_sim, 'b-', 'LineWidth', 3, 'DisplayName', 'Fitted model');
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Bacterial density (CFU/mL)', 'FontSize', 14);
title('(A) CFU/mL Scale', 'FontSize', 16);
ax.TitleHorizontalAlignment = 'left';
legend('Location', 'southeast', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
xlim([0, max(time_hr)]);

% Panel B: OD600 scale
nexttile;
ax = gca;
plot(time_hr, od_mean, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'Control data', 'LineWidth', 1);
hold on;
plot(t_sim, S_sim / OD_to_CFU, 'b-', 'LineWidth', 3, 'DisplayName', 'Fitted model');
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Cell density (OD_{600})', 'FontSize', 14);
title('(B) OD_{600} Scale', 'FontSize', 18);
ax.TitleHorizontalAlignment = 'left';
legend('Location', 'southeast', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
xlim([0, max(time_hr)]);

% Panel C: Residuals plot
nexttile;
ax = gca;
plot(time_hr, residuals / OD_to_CFU, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
hold on;
yline(0, 'r--', 'LineWidth', 2);
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Residuals (OD_{600})', 'FontSize', 14);
title('(C) Residuals', 'FontSize', 18);
ax.TitleHorizontalAlignment = 'left';
grid on;
set(gca, 'FontSize', 12);
xlim([0, max(time_hr)]);

% Panel D: Results summary (text)
nexttile;
axis off;
ax = gca;

% Create text summary
summary_text = {
    '\bf{Fitted Parameters:}', ...
    sprintf('  r = %.4f h^{-1}', r_fit), ...
    sprintf('  a_{ss} = %.4f', a_ss_fit), ...
    '', ...
    '\bf{Fixed Parameters:}', ...
    sprintf('  K = %.2e CFU/mL', K), ...
    sprintf('  S_0 = %.2e CFU/mL', S0), ...
    '', ...
    '\bf{Fit Quality:}', ...
    sprintf('  R^2 = %.4f', R_squared), ...
    sprintf('  SSR = %.4e', resnorm), ...
    '', ...
    '\bf{95% Confidence Intervals:}', ...
    sprintf('  r: [%.4f, %.4f] h^{-1}', ci_r(1), ci_r(2)), ...
    sprintf('  a_{ss}: [%.4f, %.4f]', ci_ass(1), ci_ass(2)), ...
};

text(0.01, 0.99, summary_text, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 12, 'FontName', 'Arial', ...
    'FontName', 'FixedWidth', 'Interpreter', 'tex');

title('(D) Summary', 'FontSize', 16);
ax.TitleHorizontalAlignment = 'left';
box on;

% Save figure as PDF
exportgraphics(fig, 'control_fitting_PA103_Rd3.pdf', 'ContentType', 'vector', 'Resolution', 300);
fprintf('\nFigure saved: control_fitting_PA103_Rd3.pdf\n');

% Also save as PNG for quick viewing
% exportgraphics(fig, 'control_fitting_Y3650.png', 'Resolution', 300);
% fprintf('Figure saved: control_fitting_Y3650.png\n');

%% Save Results
results.r = r_fit;
results.a_ss = a_ss_fit;
results.K = K;
results.S0 = S0;
results.R_squared = R_squared;
results.resnorm = resnorm;
results.params_init = params_init;
results.ci_r = ci_r;
results.ci_ass = ci_ass;
results.time_data = time_hr;
results.S_data = S_data;
results.residuals = residuals;
results.t_sim = t_sim;
results.S_sim = S_sim;

save('control_fitting_results_PA103_Rd3.mat', 'results');
fprintf('Results saved: control_fitting_results_PA103_Rd3.mat\n');

%% Summary for use in global fitting
fprintf('\n');
fprintf('%s\n', repmat('=', 1, 50));
fprintf('SUMMARY: Initial estimates for global fitting\n');
fprintf('%s\n', repmat('=', 1, 50));
fprintf('From control data fitting:\n');
fprintf('  r = %.4f h^-1\n', r_fit);
fprintf('  a_ss = %.4f\n', a_ss_fit);
fprintf('  K = %.2e CFU/mL (fixed)\n', K);
fprintf('  S0 = %.2e CFU/mL\n', S0);
fprintf('\nThese values can be used as initial guesses for\n');
fprintf('fitting the full model to all MOI treatment data.\n');

%% ========== Helper Functions ==========

function dSdt = logistic_growth(~, S, params, K)
    % Logistic growth model
    % dS/dt = r*S*(1 - a_ss*S/K)
    %
    % Inputs:
    %   t: time (not used, autonomous system)
    %   S: bacterial density (CFU/mL)
    %   params: [r, a_ss]
    %   K: carrying capacity (CFU/mL)
    %
    % Output:
    %   dSdt: rate of change of S

    r = params(1);
    a_ss = params(2);

    dSdt = r * S * (1 - a_ss * S / K);
end

function residuals = compute_residuals(params, time_data, S_data, S0, K)
    % Compute residuals between model prediction and data
    %
    % Inputs:
    %   params: [r, a_ss]
    %   time_data: time points (hours)
    %   S_data: observed bacterial densities (CFU/mL)
    %   S0: initial condition
    %   K: carrying capacity
    %
    % Output:
    %   residuals: vector of residuals (model - data)

    % Solve ODE
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NonNegative', 1);

    try
        [~, S_model] = ode45(@(t, S) logistic_growth(t, S, params, K), ...
            time_data, S0, options);

        % Compute residuals
        residuals = S_model - S_data;

    catch
        % If ODE solver fails, return large residuals
        residuals = 1e10 * ones(size(S_data));
    end
end