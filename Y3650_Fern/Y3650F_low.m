%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Loading and Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the excel file name
excelFileName = 'dataset_Y3650F.xlsx';
sheetName = 'Y3650_Fern_Grouped';

% Read the excel file using readcell to preserve original column names
raw_data = readcell(excelFileName, 'Sheet', sheetName);

% Extract headers (first row) and data (remaining rows)
headers = raw_data(1, :); % Extract headers from the first row
data = raw_data(2:end, :); % Extract data from the remaining rows
y3650_grouped = cell2table(data, 'VariableNames', headers); % Convert to table for easier manipulation

% Data for model fitting
time = downsample(y3650_grouped.Time, 8) / 60;  % increment of 2 hours
data_low_mean = downsample(y3650_grouped.Low_Group_Mean * 1E8, 8);
data_low_std = downsample(y3650_grouped.Low_Group_Std * 1E8, 8);

fprintf('Data loaded successfully:\n'); %[output:3359f4b0]
fprintf('  Time points: %d\n', length(time)); %[output:2996ab42]
fprintf('  Time range: %.2f to %.2f hours\n', min(time), max(time)); %[output:44add7d4]
fprintf('  Data range: %.2e to %.2e CFU/mL', min(data_low_mean), max(data_low_mean)); %[output:80f01d8f]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Parameter Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter names for reference
param_names = {'r', 'a', 'a_{ss}', 'a_{rr}', '\delta', '\gamma', '\phi', '\beta', 'm', 'S_0', 'R_0', 'P_0'};

% Initial values
r = 0.7457;      % growth rate of susceptible bacteria S
a = 4.7349;      % interspecific competition rate between S and R
a_ss = 49.0987;  % intraspecific competition rate of S
a_rr = 9.5935;   % intraspecific competition rate of R
delta = 9.8149e-04;   % rate of resistance acquisition
gamma = 0.8575;  % fitness cost of phage resistance
phi = 1.01e-8;  % adsorption rate of phage P
beta = 22.44;       % burst size of phage P
m = 1e-1;      % decay rate of phage P
S0 = 1.09472e5;     % initial value of S
R0 = 1.00032e5;        % initial value of R
P0 = 0.2996;      % initial value of P

ini_para_fit = [r; a; a_ss; a_rr; delta; gamma; phi; beta; m; S0; R0; P0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lower_bound = [0.462;    % growth rate of susceptible bacteria S
               0;        % interspecific competition rate between S and R
               0;        % intraspecific competition rate of S
               0;        % intraspecific competition rate of R
               10^(-6);  % rate of resistance acquisition
               0;        % fitness cost of phage resistance
               10^(-8);  % adsorption rate of phage P
               10;       % burst size of phage P
               1.5E-3;   % decay rate of phage P
               0;        % initial value of S
               0;        % initial value of R
               1e-3      % initial value of P
               ];

upper_bound = [1.65;     % growth rate of susceptible bacteria S
               Inf;      % interspecific competition rate between S and R
               Inf;      % intraspecific competition rate of S
               Inf;      % intraspecific competition rate of R
               10^(-3);  % rate of resistance acquisition
               1;        % fitness cost of phage resistance
               10^(-6);  % adsorption rate of phage P
               100;      % burst size of phage P
               0.1;      % decay rate of phage P
               Inf;      % initial value of S
               Inf;      % initial value of R
               0.347       % initial value of P
               ];

% Check if bounds are reasonable
if any(ini_para_fit < lower_bound) || any(ini_para_fit > upper_bound)
    warning('Initial parameters are outside bounds! Adjusting...');
    ini_para_fit = max(ini_para_fit, lower_bound);
    ini_para_fit = min(ini_para_fit, upper_bound);
end

% Enhanced optimization options
opts = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'Cache', 'on', ...                    % Enable caching for repeated evaluations
    'UseCompletePoll', true, ...          % Use complete polling for better exploration
    'UseParallel', false, ...              % Enable parallel evaluation
    'CompletePoll', 'on', ...             % Complete polling at each iteration
    'InitialMeshSize', 0.1, ...           % Initial mesh size
    'MaxIterations', 500, ...             % Maximum iterations
    'MaxFunctionEvaluations', 1e5, ...    % Maximum function evaluations
    'ScaleMesh', 'off', ...               % Mesh scaling
    'MeshTolerance', 1e-10, ...           % Mesh tolerance
    'FunctionTolerance', 1e-12, ...       % Function tolerance for accuracy
    'StepTolerance', 1e-12, ...           % Step tolerance
    'PlotFcn', @psplotbestf, ...          % Plot function
    'PollMethod', 'MADSPositiveBasis2N', ... % More efficient polling method
    'SearchFcn', @searchlhs, ...          % Latin Hypercube search
    'MeshExpansionFactor', 2.0, ...       % Mesh expansion factor
    'MeshContractionFactor', 0.5, ...     % Mesh contraction factor
    'AccelerateMesh', true);           % Enable mesh acceleration
    % 'OutputFcn', @(x,opt,state) optimization_output(x,opt,state,param_names)); % Custom output function

% Enhanced objective function with robust error handling
func = @(para_fit) robust_objective_function(para_fit, data_low_mean, time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enhanced Optimization Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmse_target = 1.5E7;
rmse = 1E8;
iteration_count = 0;
max_restarts = 5;
best_rmse = Inf;
best_params = ini_para_fit;

fprintf('\n=== Starting Enhanced Pattern Search Optimization ===\n');
fprintf('Target RMSE: %.2e\n', rmse_target);
fprintf('Initial RMSE: %.2e\n', func(ini_para_fit));
fprintf('Number of parameters: %d\n', length(ini_para_fit));

tic; % start timing

while rmse > rmse_target && iteration_count < max_restarts
    iteration_count = iteration_count + 1;
    fprintf('\n--- Optimization Run %d/%d ---\n', iteration_count, max_restarts);

    % Run pattern search
    [fitted_para, fval, exitflag, output] = patternsearch(func, ini_para_fit,...
        [], [], [], [], lower_bound, upper_bound, [], opts);

    % Update results
    rmse = fval;

    % Keep track of best solution
    if rmse < best_rmse
        best_rmse = rmse;
        best_params = fitted_para;
    end

    fprintf('Run %d Results:\n', iteration_count);
    fprintf('  RMSE: %.4e (Target: %.4e)\n', rmse, rmse_target);
    fprintf('  Exit flag: %d (%s)\n', exitflag, get_exit_flag_description(exitflag));
    fprintf('  Function evaluations: %d\n', output.funccount);
    fprintf('  Iterations: %d\n', output.iterations);

    % Prepare for next iteration if needed
    if iteration_count < max_restarts
        if exitflag < 0
            fprintf(' Poor convergence detected, perturbing parameters...\n');
            % Add smart perturbation based on parameter sensitivity
            perturbation_factor = 0.05 * (1 + iteration_count * 0.02); % Increase perturbation with iterations
            perturbation = perturbation_factor * randn(size(fitted_para)) .* (upper_bound - lower_bound);
            ini_para_fit = fitted_para + perturbation;
        else
            % Use current best as starting point
            ini_para_fit = fitted_para;
        end

        % Ensure bounds are respected
        ini_para_fit = max(ini_para_fit, lower_bound + 1e-10);
        ini_para_fit = min(ini_para_fit, upper_bound - 1e-10);
    end

end

optimization_time = toc;

% Use best solution found
fitted_para = best_params;
rmse = best_rmse;

fprintf('\n=== Optimization Summary ===\n');
fprintf('Total optimization time: %.2f seconds\n', optimization_time);
fprintf('Number of restarts: %d\n', iteration_count);

if rmse <= rmse_target
    fprintf('SUCCESS: Target RMSE achieved!\n');
else
    fprintf('Optimization completed without reaching target.\n');
end
fprintf('Final RMSE: %.4e\n', rmse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results Analysis and Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final RMSE calculation
final_rmse = func(fitted_para);
fprintf('\nFinal RMSE = %.4f\n\n', final_rmse);

% Print fitted parameters with better formatting
fprintf('=== Fitted Parameters ===');
fprintf('%-8s %12s %12s %12s %12s\n', 'Param', 'Initial', 'Fitted', 'Lower', 'Upper');
fprintf('%-8s %12s %12s %12s %12s\n', '-----', '-------', '------', '-----', '-----');
for i = 1:length(fitted_para)
    fprintf('%-8s %12.5e %12.5e %12.5e %12.5e\n', ...
        param_names{i}, ini_para_fit(i), fitted_para(i), lower_bound(i), upper_bound(i));
end

% Parameter change analysis
param_changes = abs((fitted_para - ini_para_fit) ./ ini_para_fit) * 100;
fprintf('\nParameter Changes (%%):');
for i = 1:length(fitted_para)
    if param_changes(i) > 10
        fprintf('\n %s: %.1f%% (LARGE CHANGE)', param_names{i}, param_changes(i));
    end
end
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model validation and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute fitted curve with higher resolution
tspan = linspace(0, time(end), 200);
fitted_total_host_cells = total_host_cells(fitted_para, tspan);

% Calculate R-squared and other metrics
model_at_data_points = total_host_cells(fitted_para, time);
ss_res = sum((data_low_mean - model_at_data_points).^2);
ss_tot = sum((data_low_mean - mean(data_low_mean)).^2);
r_squared = 1 - ss_res/ss_tot;

fprintf('\n=== Model Quality Metrics ===');
fprintf('R-squared: %.4f\n', r_squared);
fprintf('RMSE: %.2e\n', sqrt(mean((data_low_mean - model_at_data_points).^2)));
fprintf('Mean Absolute Error: %.2e\n', mean(abs(data_low_mean - model_at_data_points)));

% Enhanced plotting
clf;
h = axes('Position', [0 0 1 1], 'Visible', 'off');
axes('Position', [.1 .1 .62 .8]);

ax = gca;
% Plot data with error bars if available
if exist('data_low_std', 'var') && ~isempty(data_low_std)
    errorbar(time, data_low_mean, data_low_std, '.', 'MarkerSize', 25, ...
        'MarkerEdgeColor', 'red', 'Color', 'red', 'LineWidth', 1.5);
else
    plot(time, data_low_mean, '.', 'MarkerSize', 25, ...
        'MarkerEdgeColor', 'red', 'Color', 'red');    
end
hold on;

% Plot fitted curve
plot(tspan, fitted_total_host_cells, 'k-', 'LineWidth', 3);

% Formatting
ylim([0 1.2*10^8]);
ax.FontSize = 18;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
xlabel('Hours', 'FontSize', 18);
ylabel('CFU/mL', 'FontSize', 18)
legend('Dataset', 'Fitted curve', 'FontSize', 16, 'Location', 'best');
title('Y3650 Fern - Low Group', 'FontSize', 18);
ax.TitleHorizontalAlignment = 'left';

% Enhanced parameter display with better formatting
str = cell(16, 1);
str{1} = 'Fitted parameter set:';
str{2} = sprintf('$R^2$ = %.4f', r_squared);
str{3} = sprintf('RMSE = %.2e', final_rmse);
for i = 1:12
    str{i+3} = sprintf('$%s$ = %.4e', param_names{i}, fitted_para(i));
end

set(gcf, 'CurrentAxes', h);
text(.75, .6, str, 'Interpreter', 'latex', 'FontSize', 14);

% Save results
set(0, 'Units', 'normalized');
set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8]);
print('Y3650_Fern_low', '-dpdf')

% Optional: Save fitted parameters
fitted_Y3650F_low = fitted_para;
save('fitted_param_Y3650F_low.mat', 'fitted_Y3650F_low', 'r_squared', 'final_rmse', 'optimization_time');

fprintf('\nOptimization complete! Results saved to Y3650_Fern_low.pdf\n');


