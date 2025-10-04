%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the groups and corresponding files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('fitted_param_Y3650F_low.mat','fitted_Y3650F_low')
load('fitted_param_Y3650F_medium.mat','fitted_Y3650F_medium')
load('fitted_param_Y3650F_medium_high.mat','fitted_Y3650F_medium_high')
load('fitted_param_Y3650F_high.mat','fitted_Y3650F_high')

%%%%%%%%%%%%%%%%%%%%%%
% Load the dataset
%%%%%%%%%%%%%%%%%%%%%%

% Specify the excel file name
excelFileName = 'dataset_Y3650F.xlsx';
sheetName = 'Y3650_Fern_Grouped';

% Read the excel file using readcell to preserve original column names
raw_data = readcell(excelFileName, 'Sheet', sheetName);

% Extract headers (first row) and data (remaining rows)
headers = raw_data(1, :); % Extract headers from the first row
data = raw_data(2:end, :); % Extract data from the remaining rows
y3650_grouped = cell2table(data, 'VariableNames', headers); % Convert to table for easier manipulation

% Extract time data
time = downsample(y3650_grouped.Time, 8) / 60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the figure for resolution matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('Position', [100, 100, 1200, 900]);

% A: Low Group
h1 = subplot(2, 2, 1);
J = calculateJacobian(@TSRP, fitted_Y3650F_low, time);
R = pinv(J' * J) * (J' * J);
imagesc(R);
colorbar;
title('Y3650F Low Group');
xlabel('Parameters');
ylabel('Parameters');
xticks(1:12);
yticks(1:12);
axis square;
set(gca, 'FontSize', 18);
text(-1.5, 0.5, 'A.', 'FontSize', 32, 'FontWeight', 'bold');

% B: Medium Group
h2 = subplot(2, 2, 2);
J = calculateJacobian(@TSRP, fitted_Y3650F_medium, time);
R = pinv(J' * J) * (J' * J);
imagesc(R);
colorbar;
title('Y3650F Medium Group');
xlabel('Parameters');
ylabel('Parameters');
xticks(1:12);
yticks(1:12);
axis square;
set(gca, 'FontSize', 18);
text(-1.5, 0.5, 'B.', 'FontSize', 32, 'FontWeight', 'bold');

% C: Medium High Group
h3 = subplot(2, 2, 3);
J = calculateJacobian(@TSRP, fitted_Y3650F_medium_high, time);
R = pinv(J' * J) * (J' * J);
imagesc(R);
colorbar;
title('Y3650F Medium High Group');
xlabel('Parameters');
ylabel('Parameters');
xticks(1:12);
yticks(1:12);
axis square;
set(gca, 'FontSize', 18);
text(-1.5, 0.5, 'C.', 'FontSize', 32, 'FontWeight', 'bold');

% D: High Group
h4 = subplot(2, 2, 4);
J = calculateJacobian(@TSRP, fitted_Y3650F_high, time);
R = pinv(J' * J) * (J' * J);
imagesc(R);
colorbar;
title('Y3650F High Group');
xlabel('Parameters');
ylabel('Parameters');
xticks(1:12);
yticks(1:12);
axis square;
set(gca, 'FontSize', 18);
text(-1.5, 0.5, 'D.', 'FontSize', 32, 'FontWeight', 'bold');

% Adjust subplot positions
set(0, 'Units', 'normalized');
set(h1, 'Position', [0.06, 0.57, 0.38, 0.37]);
set(h2, 'Position', [0.56, 0.57, 0.38, 0.37]);
set(h3, 'Position', [0.06, 0.08, 0.38, 0.37]);
set(h4, 'Position', [0.56, 0.08, 0.38, 0.37]);

% Set paper size and save the figure
set(gcf, 'PaperSize', [14, 14], 'PaperPosition', [0, 0, 14, 14]);
print('Y3650F_res_matrices', '-djpeg', '-r300');

