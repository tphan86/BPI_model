%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load fitted parameters Y3650F from mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
load('fitted_param_Y3650F_low.mat','fitted_Y3650F_low')
load('fitted_param_Y3650F_medium.mat','fitted_Y3650F_medium')
load('fitted_param_Y3650F_medium_high.mat','fitted_Y3650F_medium_high')
load('fitted_param_Y3650F_high.mat','fitted_Y3650F_high')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the dataset from excel file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the excel file name
excelFileName = 'dataset_Y3650F.xlsx';
sheetName = 'Y3650_Fern_Grouped';

% Read the excel file using readcell to preserve original column names
raw_data = readcell(excelFileName, 'Sheet', sheetName);

% Extract headers (first row) and data (remaining rows)
headers = raw_data(1, :); % Extract headers from the first row
data = raw_data(2:end, :); % Extract data from the remaining rows
y3650_grouped = cell2table(data, 'VariableNames', headers); % Convert to table for easier manipulation

% Extract growth data for 4 groups
time = downsample(y3650_grouped.Time,8)/60; % increment of 2 hours

data_low_mean = downsample(y3650_grouped.Low_Group_Mean*1E8,8);
data_low_std = downsample(y3650_grouped.Low_Group_Std*1E8,8);

data_medium_mean = downsample(y3650_grouped.Medium_Group_Mean*1E8,8);
data_medium_std = downsample(y3650_grouped.Medium_Group_Std*1E8,8);

data_medium_high_mean = downsample(y3650_grouped.Medium_High_Group_Mean*1E8,8);
data_medium_high_std = downsample(y3650_grouped.Medium_High_Group_Std*1E8,8);

data_high_mean = downsample(y3650_grouped.High_Group_Mean*1E8,8);
data_high_std = downsample(y3650_grouped.High_Group_Std*1E8,8);

% Calculated fitted curves for low group
tspan = linspace(time(1), time(end));
[T1, S1, R1, P1] = TSRP(fitted_Y3650F_low, tspan);
[T2, S2, R2, P2] = TSRP(fitted_Y3650F_medium, tspan);
[T3, S3, R3, P3] = TSRP(fitted_Y3650F_medium_high, tspan);
[T4, S4, R4, P4] = TSRP(fitted_Y3650F_high, tspan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = subplot(2,2,1);

set(gca,'FontSize', 20, 'FontName', 'Arial')
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
yyaxis left

plot(time, data_low_mean,'.', 'MarkerSize',40,...
     'MarkerEdgeColor','red','Color','red'); hold on
plot(tspan, T1, 'k-', 'LineWidth', 4.5) % the fitted curve S + R
plot(tspan, S1, 'b--', 'LineWidth', 4) % fitted susceptibles S
plot(tspan, R1, 'g:', 'LineWidth', 4) % fitted resistants R
% ylim([0 1.2*10^8]);
xlabel('Hours','FontSize',22)
ylabel('CFU/mL','FontSize',22)
% set(gca, 'YScale', 'log')

yyaxis right
plot(tspan, P1, 'm-', 'LineWidth', 4) % fitted phage P
% ylim([1E0 1E12]);
% set(gca, 'YScale', 'log')
ylabel('PFU/mL','FontSize',22)

title('Y3650F Low group', 'FontSize', 24)
% Add 'A' label to top-left outside the plot
text(-0.1, 1.07, 'A.', 'Units', 'normalized', 'FontSize', 32,...
                                'FontWeight', 'bold', 'FontName', 'Arial')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h2 = subplot(2,2,2);

set(gca,'FontSize', 20, 'FontName', 'Arial')
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
yyaxis left

plot(time, data_medium_mean,'.', 'MarkerSize',40,...
     'MarkerEdgeColor','red','Color','red'); hold on
plot(tspan, T2, 'k-', 'LineWidth', 4.5) % the fitted curve S + R
plot(tspan, S2, 'b--', 'LineWidth', 4) % fitted susceptibles S
plot(tspan, R2, 'g:', 'LineWidth', 4) % fitted resistants R
% ylim([0 1.2*10^8]);
xlabel('Hours','FontSize',22)
ylabel('CFU/mL','FontSize',22)
% set(gca, 'YScale', 'log')

yyaxis right
plot(tspan, P2, 'm-', 'LineWidth', 4) % fitted phage P
% ylim([1E0 1E12]);
% set(gca, 'YScale', 'log')
ylabel('PFU/mL','FontSize', 22)

title('Y3650F Medium group', 'FontSize', 24)
% Add 'B' label to top-left outside the plot
text(-0.1, 1.07, 'B.', 'Units', 'normalized', 'FontSize', 32,...
                               'FontWeight', 'bold', 'FontName', 'Arial')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h3 = subplot(2,2,3);

set(gca,'FontSize', 20, 'FontName', 'Arial')
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
yyaxis left

plot(time, data_medium_high_mean,'.', 'MarkerSize',40,...
     'MarkerEdgeColor','red','Color','red'); hold on
plot(tspan, T3, 'k-', 'LineWidth', 4.5) % the fitted curve S + R
plot(tspan, S3, 'b--', 'LineWidth', 4) % fitted susceptibles S
plot(tspan, R3, 'g:', 'LineWidth', 4) % fitted resistants R
ylim([0 1.2*10^8]);
xlabel('Hours','FontSize',22)
ylabel('CFU/mL','FontSize',22)
% set(gca, 'YScale', 'log')

yyaxis right
plot(tspan, P3, 'm-', 'LineWidth', 4) % fitted phage P
% ylim([1E0 1E12]);
% set(gca, 'YScale', 'log')
ylabel('PFU/mL','FontSize', 22)

title('Y3650F Medium High group', 'FontSize', 24)
% Add 'C' label to top-left outside the plot
text(-0.1, 1.07, 'C.', 'Units', 'normalized', 'FontSize', 32,...
                               'FontWeight', 'bold', 'FontName', 'Arial')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h4 = subplot(2,2,4);

set(gca,'FontSize', 20, 'FontName', 'Arial')
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
yyaxis left

plot(time, data_high_mean,'.', 'MarkerSize',40,...
     'MarkerEdgeColor','red','Color','red'); hold on
plot(tspan, T4, 'k-', 'LineWidth', 4.5) % the fitted curve S + R
plot(tspan, S4, 'b--', 'LineWidth', 4) % fitted susceptibles S
plot(tspan, R4, 'g:', 'LineWidth', 4) % fitted resistants R
ylim([0 1.2*10^8]);
xlabel('Hours','FontSize',22)
ylabel('CFU/mL','FontSize',22)
% set(gca, 'YScale', 'log')

yyaxis right
plot(tspan, P4, 'm-', 'LineWidth', 4) % fitted phage P
% ylim([1E0 1E12]);
% set(gca, 'YScale', 'log')
ylabel('PFU/mL','FontSize', 22)

title('Y3650F High group', 'FontSize', 24)
% Add 'D' label to top-left outside the plot
text(-0.1, 1.07, 'D.', 'Units', 'normalized', 'FontSize', 32,...
                               'FontWeight', 'bold', 'FontName', 'Arial')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add a common legend at the bottom of the figure
lgd = legend('growth dataset','fitted total cells', 'fitted susceptibles',...
       'fitted resistants','fitted phage');
lgd.FontSize = 22;
lgd.Orientation = 'horizontal';
lgd.Position = [0.5, 0.02, 0, 0];
lgd.Units = 'normalized';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','normalized')
set(h4,'position',[.56 .1 .38 .37])
set(h3,'position',[.05 .1 .38 .37])
set(h2,'position',[.56 .57 .38 .37])
set(h1,'position',[.05 .57 .38 .37])
set(gcf, 'PaperSize', [20 16], 'PaperPosition', [0 0 20 16])
print('Y3650F_fitted_curves', '-djpeg', '-r300')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%