%% Fit Control Data for P. larvae 25747
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
%
% REVISION NOTES:
%   - Fixed strain name inconsistency (Y-25747 vs Y-3650)
%   - Added weighted least squares option to account for heteroscedasticity
%   - Added option to fit in OD600 scale directly (avoids numerical issues)
%   - Improved initial parameter estimation
%   - Added AIC/BIC for model comparison
%   - Better handling of output file names

clear; clc; close all;

%% ==================== USER SETTINGS ====================
% Choose strain name (for file naming and plot titles)
STRAIN_NAME = '25747';  % Options: 'Y-25747' or 'Y-3650'

% Time window for fitting (hours)
% Set to Inf to use all data, or specify max time (e.g., 30)
MAX_FIT_TIME = 30;  % Only use data up to this time for fitting
                    % Recommended: 30 hours (before decline phase)

% Choose fitting scale
FIT_IN_OD_SCALE = true;  % true: fit in OD600 scale (recommended)
                          % false: fit in CFU/mL scale

% Choose weighting scheme
USE_WEIGHTED_FIT = false;  % true: weight by 1/variance (if std available)
                           % false: unweighted least squares

%% ==================== LOAD DATA ====================
fprintf('Loading data for P. larvae %s...\n', STRAIN_NAME);

data_table = readtable('dataset_25747F.xlsx', 'Sheet', 'Y25747_Fern_Grouped');
time_min = data_table.Time;           % time in minutes
od_mean = data_table.Ctrl_Group_Mean; % OD600
od_std = data_table.Ctrl_Group_Std;   % OD600

% Convert time to hours
time_hr = time_min / 60;

% Conversion factor
OD_to_CFU = 1e9;  % 1 OD600 ≈ 1E9 CFU/mL

% Convert to CFU/mL (for reference and plotting)
S_data_CFU = od_mean * OD_to_CFU;   % CFU/mL
S_std_CFU = od_std * OD_to_CFU;     % Std in CFU/mL

fprintf('Data loaded: %d time points over %.1f hours\n', length(time_hr), max(time_hr));
fprintf('OD600 range: [%.4f, %.4f]\n', min(od_mean), max(od_mean));

%% ==================== APPLY TIME WINDOW ====================
if isfinite(MAX_FIT_TIME)
    idx_fit = time_hr <= MAX_FIT_TIME;
    time_hr_fit = time_hr(idx_fit);
    od_mean_fit = od_mean(idx_fit);
    od_std_fit = od_std(idx_fit);
    
    fprintf('\nTime window applied: using data from 0 to %.1f hours\n', MAX_FIT_TIME);
    fprintf('  Points for fitting: %d (out of %d total)\n', sum(idx_fit), length(time_hr));
    fprintf('  OD600 range in window: [%.4f, %.4f]\n', min(od_mean_fit), max(od_mean_fit));
else
    time_hr_fit = time_hr;
    od_mean_fit = od_mean;
    od_std_fit = od_std;
    fprintf('\nUsing all data (no time window)\n');
end

% Keep full data for plotting
time_hr_full = time_hr;
od_mean_full = od_mean;
od_std_full = od_std;
S_data_CFU_full = S_data_CFU;
S_std_CFU_full = S_std_CFU;

%% ==================== FIXED PARAMETERS ====================
K_CFU = 1e9;                % Carrying capacity in CFU/mL
K_OD = K_CFU / OD_to_CFU;   % Carrying capacity in OD600 (= 1.0)

% Initial condition
S0_OD = od_mean(1);
S0_CFU = S0_OD * OD_to_CFU;

fprintf('\nInitial conditions:\n');
fprintf('  S0 = %.4f OD600 = %.2e CFU/mL\n', S0_OD, S0_CFU);
fprintf('  (Manuscript reports ~2.4E6 CFU/mL for 25747)\n');

%% ==================== SET UP FITTING ====================
% Use filtered data for fitting
if FIT_IN_OD_SCALE
    fprintf('\nFitting in OD600 scale (recommended for numerical stability)\n');
    y_data = od_mean_fit;
    y_std = od_std_fit;
    y0 = S0_OD;
    K = K_OD;
    scale_label = 'OD_{600}';
else
    fprintf('\nFitting in CFU/mL scale\n');
    y_data = od_mean_fit * OD_to_CFU;
    y_std = od_std_fit * OD_to_CFU;
    y0 = S0_CFU;
    K = K_CFU;
    scale_label = 'CFU/mL';
end

% Time for fitting
time_fit = time_hr_fit;

%% ==================== INITIAL PARAMETER ESTIMATES ====================
% Estimate r from exponential phase (first few hours)
% During exponential phase: ln(S/S0) ≈ r*t

% Find exponential phase (before significant slowdown)
idx_exp = time_fit < 10 & time_fit > 0;  % Use first 10 hours
if sum(idx_exp) > 3
    % Linear regression on log-transformed data
    log_ratio = log(y_data(idx_exp) ./ y0);
    p = polyfit(time_fit(idx_exp), log_ratio, 1);
    r_init = max(0.1, p(1));  % Ensure positive
else
    r_init = 0.5;  % Default guess
end

% Estimate a_ss from stationary phase
% At stationary phase: S_max ≈ K/a_ss, so a_ss ≈ K/S_max
y_max = max(y_data);
a_ss_init = K / y_max;

% Ensure initial guess is within bounds
r_init = min(max(r_init, 0.1), 2.0);
a_ss_init = min(max(a_ss_init, 0.5), 5.0);

fprintf('\nInitial parameter estimates:\n');
fprintf('  r_init = %.4f h^-1\n', r_init);
fprintf('  a_ss_init = %.4f\n', a_ss_init);

params_init = [r_init, a_ss_init];

%% ==================== PARAMETER BOUNDS ====================
% r: [0.1, 3.0] h^-1 (literature range for P. larvae: 0.462-1.65)
% a_ss: [0.1, 10] (dimensionless)
lb = [0.1, 0.1];
ub = [3.0, 10.0];

%% ==================== WEIGHTING ====================
if USE_WEIGHTED_FIT && any(y_std > 0)
    % Weights inversely proportional to variance
    % Avoid division by zero
    weights = 1 ./ max(y_std.^2, 1e-10);
    weights = weights / mean(weights);  % Normalize
    fprintf('\nUsing weighted least squares (weights based on 1/variance)\n');
else
    weights = ones(size(y_data));
    fprintf('\nUsing unweighted least squares\n');
end

%% ==================== FITTING OPTIONS ====================
options = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12);

%% ==================== OBJECTIVE FUNCTION ====================
objective = @(params) compute_weighted_residuals(params, time_fit, y_data, y0, K, weights);

%% ==================== RUN OPTIMIZATION ====================
fprintf('\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('Starting parameter optimization...\n');
fprintf('%s\n', repmat('=', 1, 60));

[params_opt, resnorm, residuals, exitflag, output, lambda, jacobian] = ...
    lsqnonlin(objective, params_init, lb, ub, options);

%% ==================== EXTRACT FITTED PARAMETERS ====================
r_fit = params_opt(1);
a_ss_fit = params_opt(2);

fprintf('\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('FITTING RESULTS FOR P. larvae %s\n', STRAIN_NAME);
fprintf('%s\n', repmat('=', 1, 60));
fprintf('\nFitted parameters:\n');
fprintf('  r = %.4f h^-1 (intrinsic growth rate)\n', r_fit);
fprintf('  a_ss = %.4f (intraspecific competition coefficient)\n', a_ss_fit);
fprintf('\nDerived quantities:\n');
fprintf('  Effective carrying capacity K/a_ss = %.4f %s\n', K/a_ss_fit, scale_label);
if FIT_IN_OD_SCALE
    fprintf('  Effective carrying capacity K/a_ss = %.2e CFU/mL\n', (K/a_ss_fit)*OD_to_CFU);
end
fprintf('  Doubling time = %.2f hours\n', log(2)/r_fit);
fprintf('\nFit quality:\n');
fprintf('  Sum of squared residuals: %.4e\n', resnorm);
fprintf('  Exit flag: %d (%s)\n', exitflag, get_exitflag_message(exitflag));

%% ==================== CONFIDENCE INTERVALS ====================
% Method 1: From Jacobian (may fail if parameters are correlated)
ci_method = 'none';
use_bootstrap = false;

try
    n = length(y_data);
    np = length(params_opt);
    dof = n - np;  % degrees of freedom
    
    % Mean squared error
    mse = resnorm / dof;
    
    % Covariance matrix from Jacobian
    J = full(jacobian);  % Convert sparse to full matrix if needed
    
    % Detailed Jacobian diagnostics
    fprintf('\nJacobian diagnostics:\n');
    fprintf('  Size of Jacobian: %d x %d\n', size(J, 1), size(J, 2));
    fprintf('  Number of data points: %d\n', n);
    fprintf('  Number of parameters: %d\n', np);
    fprintf('  Contains NaN: %s\n', mat2str(any(isnan(J(:)))));
    fprintf('  Contains Inf: %s\n', mat2str(any(isinf(J(:)))));
    fprintf('  Min value: %.4e\n', min(J(:)));
    fprintf('  Max value: %.4e\n', max(J(:)));
    
    % Check if Jacobian has correct size
    if size(J, 1) ~= n || size(J, 2) ~= np
        fprintf('  WARNING: Jacobian size mismatch!\n');
        fprintf('  Expected: %d x %d, Got: %d x %d\n', n, np, size(J,1), size(J,2));
        use_bootstrap = true;
    end
    
    % Check for NaN or Inf
    if any(isnan(J(:))) || any(isinf(J(:)))
        fprintf('  WARNING: Jacobian contains NaN or Inf values\n');
        use_bootstrap = true;
    end
    
    if ~use_bootstrap
        % Check Jacobian condition
        JtJ = J' * J;
        cond_num = cond(JtJ);
        fprintf('  Condition number of J''*J: %.2e\n', cond_num);
        
        % Show JtJ matrix for 2x2 case
        fprintf('  J''*J matrix:\n');
        fprintf('    [%.4e, %.4e]\n', JtJ(1,1), JtJ(1,2));
        fprintf('    [%.4e, %.4e]\n', JtJ(2,1), JtJ(2,2));
        
        if cond_num > 1e10
            fprintf('  WARNING: Jacobian is poorly conditioned (>1e10)\n');
            fprintf('  Parameters r and a_ss are likely highly correlated.\n');
            use_bootstrap = true;
        else
            % Compute covariance matrix
            cov_matrix = mse * pinv(JtJ);
            
            fprintf('  MSE: %.4e\n', mse);
            fprintf('  Covariance matrix:\n');
            fprintf('    [%.4e, %.4e]\n', cov_matrix(1,1), cov_matrix(1,2));
            fprintf('    [%.4e, %.4e]\n', cov_matrix(2,1), cov_matrix(2,2));
            
            % Check for negative variances (indicates numerical issues)
            var_diag = diag(cov_matrix);
            fprintf('  Variances: [%.4e, %.4e]\n', var_diag(1), var_diag(2));
            
            if any(var_diag < 0)
                fprintf('  WARNING: Negative variance detected, switching to bootstrap\n');
                use_bootstrap = true;
            elseif any(isnan(var_diag)) || any(isinf(var_diag))
                fprintf('  WARNING: NaN/Inf in variances, switching to bootstrap\n');
                use_bootstrap = true;
            else
                % Standard errors
                se = sqrt(var_diag);
                fprintf('  Standard errors: [%.4e, %.4e]\n', se(1), se(2));
                
                % 95% confidence intervals (t-distribution)
                t_crit = tinv(0.975, dof);
                fprintf('  t-critical (df=%d): %.4f\n', dof, t_crit);
                
                ci_r = [r_fit - t_crit*se(1), r_fit + t_crit*se(1)];
                ci_ass = [a_ss_fit - t_crit*se(2), a_ss_fit + t_crit*se(2)];
                
                % Check if CIs are valid
                if any(isnan(ci_r)) || any(isnan(ci_ass))
                    fprintf('  WARNING: NaN in confidence intervals, switching to bootstrap\n');
                    use_bootstrap = true;
                else
                    ci_method = 'jacobian';
                    fprintf('\n95%% Confidence Intervals (from Jacobian):\n');
                    fprintf('  r: [%.4f, %.4f] h^-1\n', ci_r(1), ci_r(2));
                    fprintf('  a_ss: [%.4f, %.4f]\n', ci_ass(1), ci_ass(2));
                    
                    % Parameter correlation
                    corr_val = cov_matrix(1,2) / (se(1) * se(2));
                    fprintf('\nParameter correlation (r, a_ss): %.4f\n', corr_val);
                    if abs(corr_val) > 0.95
                        fprintf('  Note: Parameters are highly correlated.\n');
                    end
                end
            end
        end
    end
    
catch ME
    fprintf('\nJacobian method failed with error:\n');
    fprintf('  %s\n', ME.message);
    fprintf('  at line %d in %s\n', ME.stack(1).line, ME.stack(1).name);
    use_bootstrap = true;
end

% Method 2: Bootstrap (more robust for correlated parameters)
if use_bootstrap
    fprintf('\n--- Computing Bootstrap Confidence Intervals ---\n');
    n_boot = 1000;  % Number of bootstrap samples
    fprintf('  Running %d bootstrap iterations...\n', n_boot);
    
    boot_params = zeros(n_boot, 2);
    boot_success = 0;
    
    % Suppress optimization output for bootstrap
    boot_options = optimoptions('lsqnonlin', ...
        'Display', 'off', ...
        'MaxIterations', 500, ...
        'MaxFunctionEvaluations', 2000, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8);
    
    rng(42);  % For reproducibility
    
    for i = 1:n_boot
        % Resample with replacement
        boot_idx = randsample(n, n, true);
        y_boot = y_data(boot_idx);
        t_boot = time_fit(boot_idx);
        w_boot = weights(boot_idx);
        
        % Sort by time for ODE solver
        [t_boot_sorted, sort_idx] = sort(t_boot);
        y_boot_sorted = y_boot(sort_idx);
        w_boot_sorted = w_boot(sort_idx);
        
        % Fit to bootstrap sample
        try
            boot_obj = @(p) compute_weighted_residuals(p, t_boot_sorted, y_boot_sorted, y0, K, w_boot_sorted);
            [p_boot, ~, ~, exitflag_boot] = lsqnonlin(boot_obj, params_opt, lb, ub, boot_options);
            
            if exitflag_boot > 0
                boot_success = boot_success + 1;
                boot_params(boot_success, :) = p_boot;
            end
        catch
            % Skip failed bootstrap iterations
        end
        
        % Progress indicator
        if mod(i, 200) == 0
            fprintf('    Progress: %d/%d (%.0f%% successful)\n', i, n_boot, 100*boot_success/i);
        end
    end
    
    % Compute bootstrap confidence intervals
    if boot_success >= 100
        boot_params = boot_params(1:boot_success, :);
        
        % Percentile method (2.5th and 97.5th percentiles)
        ci_r = prctile(boot_params(:,1), [2.5, 97.5]);
        ci_ass = prctile(boot_params(:,2), [2.5, 97.5]);
        
        % Standard errors from bootstrap
        se = std(boot_params);
        
        ci_method = 'bootstrap';
        fprintf('\n95%% Confidence Intervals (Bootstrap, n=%d):\n', boot_success);
        fprintf('  r: [%.4f, %.4f] h^-1\n', ci_r(1), ci_r(2));
        fprintf('  a_ss: [%.4f, %.4f]\n', ci_ass(1), ci_ass(2));
        
        fprintf('\nBootstrap Standard Errors:\n');
        fprintf('  SE(r) = %.4f\n', se(1));
        fprintf('  SE(a_ss) = %.4f\n', se(2));
        
        % Bootstrap correlation
        boot_corr = corr(boot_params(:,1), boot_params(:,2));
        fprintf('\nParameter correlation (bootstrap): %.4f\n', boot_corr);
    else
        fprintf('\n  WARNING: Too few successful bootstrap iterations (%d)\n', boot_success);
        ci_r = [NaN, NaN];
        ci_ass = [NaN, NaN];
        se = [NaN, NaN];
    end
end

%% ==================== GOODNESS OF FIT METRICS ====================
% R-squared
SS_res = resnorm;
SS_tot = sum((y_data - mean(y_data)).^2);
R_squared = 1 - SS_res / SS_tot;

% Adjusted R-squared
R_squared_adj = 1 - (1 - R_squared) * (n - 1) / (n - np - 1);

% AIC and BIC (for model comparison)
% Assuming Gaussian errors
log_likelihood = -n/2 * log(2*pi) - n/2 * log(resnorm/n) - n/2;
AIC = -2 * log_likelihood + 2 * np;
BIC = -2 * log_likelihood + np * log(n);

fprintf('\nGoodness of Fit:\n');
fprintf('  R-squared: %.4f\n', R_squared);
fprintf('  Adjusted R-squared: %.4f\n', R_squared_adj);
fprintf('  AIC: %.2f\n', AIC);
fprintf('  BIC: %.2f\n', BIC);

%% ==================== SIMULATE FITTED MODEL ====================
% Simulate over full time range (for comparison with all data)
tspan = linspace(0, max(time_hr_full), 500);
[t_sim, y_sim] = ode45(@(t, y) logistic_growth(t, y, params_opt, K), tspan, y0);

%% ==================== PLOTTING ====================
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 9], ...
    'PaperUnits', 'inches', 'PaperSize', [12, 9], 'PaperPosition', [0, 0, 12, 9]);

t_layout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Main title
if isfinite(MAX_FIT_TIME)
    title_str = sprintf('P. larvae %s Control Data Fitting (0-%.0f hours)', STRAIN_NAME, MAX_FIT_TIME);
else
    title_str = sprintf('P. larvae %s Control Data Fitting (No Phage)', STRAIN_NAME);
end
title(t_layout, title_str, 'FontSize', 18, 'FontWeight', 'bold');

% Panel A: Data and fit in OD600 scale
nexttile;
ax = gca;
% Plot all data (light color for data outside fitting window)
if isfinite(MAX_FIT_TIME)
    idx_outside = time_hr_full > MAX_FIT_TIME;
    errorbar(time_hr_full(idx_outside), od_mean_full(idx_outside), od_std_full(idx_outside), ...
        'o', 'MarkerSize', 4, 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.8 0.8 0.8], ...
        'DisplayName', 'Data (not used for fit)', 'LineWidth', 1, 'CapSize', 2);
    hold on;
end
errorbar(time_hr_fit, od_mean_fit, od_std_fit, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'Data used for fitting', 'LineWidth', 1, 'CapSize', 3);
hold on;
if FIT_IN_OD_SCALE
    plot(t_sim, y_sim, 'b-', 'LineWidth', 3, 'DisplayName', 'Fitted model');
else
    plot(t_sim, y_sim / OD_to_CFU, 'b-', 'LineWidth', 3, 'DisplayName', 'Fitted model');
end
% Add vertical line at fitting boundary
if isfinite(MAX_FIT_TIME)
    xline(MAX_FIT_TIME, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Fit boundary');
end
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Cell density (OD_{600})', 'FontSize', 14);
title('(A) OD_{600} Scale', 'FontSize', 16);
ax.TitleHorizontalAlignment = 'left';
legend('Location', 'southeast', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
xlim([0, max(time_hr_full)]);
ylim([0, max(od_mean_full)*1.1]);

% Panel B: Data and fit in CFU/mL scale
nexttile;
ax = gca;
if isfinite(MAX_FIT_TIME)
    errorbar(time_hr_full(idx_outside), S_data_CFU_full(idx_outside), S_std_CFU_full(idx_outside), ...
        'o', 'MarkerSize', 4, 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.8 0.8 0.8], ...
        'DisplayName', 'Data (not used for fit)', 'LineWidth', 1, 'CapSize', 2);
    hold on;
end
errorbar(time_hr_fit, od_mean_fit * OD_to_CFU, od_std_fit * OD_to_CFU, 'ro', ...
    'MarkerSize', 4, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'Data used for fitting', 'LineWidth', 1, 'CapSize', 3);
hold on;
if FIT_IN_OD_SCALE
    plot(t_sim, y_sim * OD_to_CFU, 'b-', 'LineWidth', 3, 'DisplayName', 'Fitted model');
else
    plot(t_sim, y_sim, 'b-', 'LineWidth', 3, 'DisplayName', 'Fitted model');
end
if isfinite(MAX_FIT_TIME)
    xline(MAX_FIT_TIME, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Fit boundary');
end
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Bacterial density (CFU/mL)', 'FontSize', 14);
title('(B) CFU/mL Scale', 'FontSize', 16);
ax.TitleHorizontalAlignment = 'left';
legend('Location', 'southeast', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
xlim([0, max(time_hr_full)]);

% Panel C: Residuals (only for fitting window)
nexttile;
ax = gca;
% Unweight residuals for plotting
unweighted_residuals = residuals ./ sqrt(weights);
if FIT_IN_OD_SCALE
    plot(time_fit, unweighted_residuals, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    ylabel_str = 'Residuals (OD_{600})';
else
    plot(time_fit, unweighted_residuals / OD_to_CFU, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    ylabel_str = 'Residuals (OD_{600})';
end
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
xlabel('Time (hours)', 'FontSize', 14);
ylabel(ylabel_str, 'FontSize', 14);
if isfinite(MAX_FIT_TIME)
    title(sprintf('(C) Residuals (0-%.0f hours)', MAX_FIT_TIME), 'FontSize', 16, 'HorizontalAlignment', 'left');
else
    title('(C) Residuals', 'FontSize', 16);
end
ax.TitleHorizontalAlignment = 'left';
grid on;
set(gca, 'FontSize', 12);
xlim([0, max(time_fit)*1.05]);

% Panel D: Summary text
nexttile;
ax = gca;
axis off;
box on;

if isfinite(MAX_FIT_TIME)
    fit_window_str = sprintf('Fit window: 0-%.0f hours', MAX_FIT_TIME);
else
    fit_window_str = 'Fit window: all data';
end

% CI method string
if exist('ci_method', 'var') && strcmp(ci_method, 'bootstrap')
    ci_label = '95% CI (bootstrap):';
else
    ci_label = '95% CI (Jacobian):';
end

summary_text = {
    '\bf{Fitted Parameters:}', ...
    sprintf('  r = %.4f h^{-1}', r_fit), ...
    sprintf('  a_{ss} = %.4f', a_ss_fit), ...
    '', ...
    '\bf{Fixed Parameters:}', ...
    sprintf('  K = %.2e CFU/mL', K_CFU), ...
    sprintf('  S_0 = %.2e CFU/mL', S0_CFU), ...
    '', ...
    '\bf{Fit Quality:}', ...
    sprintf('  R^2 = %.4f', R_squared), ...
    sprintf('  R^2_{adj} = %.4f', R_squared_adj), ...
    '', ...
    ['\bf{' ci_label '}'], ...
    sprintf('  r: [%.4f, %.4f]', ci_r(1), ci_r(2)), ...
    sprintf('  a_{ss}: [%.4f, %.4f]', ci_ass(1), ci_ass(2)), ...
};

text(0.05, 0.95, summary_text, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 12, 'FontName', 'FixedWidth', ...
    'Interpreter', 'tex');
title('(D) Summary', 'FontSize', 16);
ax.TitleHorizontalAlignment = 'left';

%% ==================== SAVE OUTPUTS ====================
% Generate file names based on strain
strain_tag = strrep(STRAIN_NAME, '-', '');  % Remove hyphen for filenames

% Save figure
if isfinite(MAX_FIT_TIME)
    fig_filename = sprintf('control_fitting_%s_%dhr.pdf', strain_tag, MAX_FIT_TIME);
else
    fig_filename = sprintf('control_fitting_%s_full.pdf', strain_tag);
end
exportgraphics(fig, fig_filename, 'ContentType', 'vector', 'Resolution', 300);
fprintf('\nFigure saved: %s\n', fig_filename);

% Save results structure
results = struct();
results.strain = STRAIN_NAME;
results.max_fit_time = MAX_FIT_TIME;
results.r = r_fit;
results.a_ss = a_ss_fit;
results.K_CFU = K_CFU;
results.S0_CFU = S0_CFU;
results.S0_OD = S0_OD;
results.R_squared = R_squared;
results.R_squared_adj = R_squared_adj;
results.AIC = AIC;
results.BIC = BIC;
results.resnorm = resnorm;
results.exitflag = exitflag;
results.params_init = params_init;
results.ci_r = ci_r;
results.ci_ass = ci_ass;
results.se = se;
results.time_data_hr = time_fit;
results.od_data = od_mean_fit;
results.od_std = od_std_fit;
results.time_full_hr = time_hr_full;
results.od_full = od_mean_full;
results.residuals = residuals;
results.t_sim = t_sim;
results.y_sim = y_sim;
results.fit_scale = ifelse(FIT_IN_OD_SCALE, 'OD600', 'CFU_mL');
results.weighted = USE_WEIGHTED_FIT;

if isfinite(MAX_FIT_TIME)
    mat_filename = sprintf('control_fitting_results_%s_%dhr.mat', strain_tag, MAX_FIT_TIME);
else
    mat_filename = sprintf('control_fitting_results_%s_full.mat', strain_tag);
end
save(mat_filename, 'results');
fprintf('Results saved: %s\n', mat_filename);

%% ==================== SUMMARY FOR GLOBAL FITTING ====================
fprintf('\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('SUMMARY: Parameters for global fitting of %s\n', STRAIN_NAME);
fprintf('%s\n', repmat('=', 1, 60));
fprintf('From control data fitting:\n');
fprintf('  r = %.4f h^-1  (95%% CI: [%.4f, %.4f])\n', r_fit, ci_r(1), ci_r(2));
fprintf('  a_ss = %.4f    (95%% CI: [%.4f, %.4f])\n', a_ss_fit, ci_ass(1), ci_ass(2));
fprintf('  K = %.2e CFU/mL (fixed)\n', K_CFU);
fprintf('  S0 = %.2e CFU/mL (from OD600 = %.4f)\n', S0_CFU, S0_OD);
fprintf('\nThese values can be used as initial guesses for\n');
fprintf('fitting the full model to all MOI treatment data.\n');
fprintf('\nNote: The effective carrying capacity is K/a_ss = %.2e CFU/mL\n', K_CFU/a_ss_fit);
fprintf('      which corresponds to OD600 = %.4f\n', K_CFU/a_ss_fit/OD_to_CFU);

%% ==================== HELPER FUNCTIONS ====================

function dSdt = logistic_growth(~, S, params, K)
    % Logistic growth model: dS/dt = r*S*(1 - a_ss*S/K)
    %
    % Inputs:
    %   t: time (not used, autonomous system)
    %   S: bacterial density (in fitting units)
    %   params: [r, a_ss]
    %   K: carrying capacity (in fitting units)
    %
    % Output:
    %   dSdt: rate of change of S
    
    r = params(1);
    a_ss = params(2);
    
    % Ensure non-negative
    S = max(S, 0);
    
    dSdt = r * S * (1 - a_ss * S / K);
end

function residuals = compute_weighted_residuals(params, time_data, y_data, y0, K, weights)
    % Compute weighted residuals between model prediction and data
    %
    % Inputs:
    %   params: [r, a_ss]
    %   time_data: time points (hours)
    %   y_data: observed data (in fitting units)
    %   y0: initial condition
    %   K: carrying capacity
    %   weights: weight vector (typically 1/variance)
    %
    % Output:
    %   residuals: vector of weighted residuals
    
    % ODE solver options
    ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NonNegative', 1);
    
    try
        [~, y_model] = ode45(@(t, y) logistic_growth(t, y, params, K), ...
            time_data, y0, ode_opts);
        
        % Weighted residuals
        residuals = sqrt(weights) .* (y_model - y_data);
        
    catch
        % If ODE solver fails, return large residuals
        residuals = 1e10 * ones(size(y_data));
    end
end

function msg = get_exitflag_message(flag)
    % Return human-readable message for lsqnonlin exit flag
    switch flag
        case 1
            msg = 'Function converged to a solution';
        case 2
            msg = 'Change in x smaller than tolerance';
        case 3
            msg = 'Change in residual smaller than tolerance';
        case 4
            msg = 'Magnitude of search direction smaller than tolerance';
        case 0
            msg = 'Maximum iterations or function evaluations reached';
        case -1
            msg = 'Output function terminated algorithm';
        case -2
            msg = 'Bounds are inconsistent';
        otherwise
            msg = 'Unknown exit condition';
    end
end

function result = ifelse(condition, true_val, false_val)
    % Simple if-else helper for inline conditionals
    if condition
        result = true_val;
    else
        result = false_val;
    end
end