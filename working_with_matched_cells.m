%%
% Initialize an empty array to store the concatenated values
concatenated_means_session_1 = [];
concatenated_means_session_2 = [];
% concatenated_means_session_3 = [];
% Initialize an empty array to store the concatenated values
concatenated_sems_session_1 = [];
concatenated_sems_session_2 = [];
% concatenated_sems_session_3 = [];


% Loop through each row of neuron_mean
for ii = 1:size(neuron_mean_array, 1)
    % Get the cell in column 1 for the current row
    cell_mean_data_1 = neuron_mean_array{ii, 1};
    cell_mean_data_2 = neuron_mean_array{ii, 2};
%     cell_mean_data_3 = neuron_mean_array{ii, 3};
    cell_sem_data_1 = neuron_sem_array{ii, 1};
    cell_sem_data_2 = neuron_sem_array{ii, 2};
%     cell_sem_data_3 = neuron_sem_array{ii, 3};
    % Concatenate the values from the cell to the overall array
    concatenated_means_session_1 = [concatenated_means_session_1; cell_mean_data_1];
    concatenated_means_session_2 = [concatenated_means_session_2; cell_mean_data_2];
%     concatenated_means_session_3 = [concatenated_means_session_3; cell_mean_data_3];
    concatenated_sems_session_1 = [concatenated_sems_session_1; cell_sem_data_1];
    concatenated_sems_session_2 = [concatenated_sems_session_2; cell_sem_data_2];
%     concatenated_sems_session_3 = [concatenated_sems_session_3; cell_sem_data_3];
end

figure; plot(ts1, mean(concatenated_means_session_1(respClass_all_array{1, 1} == 1,:)));
ylim([-0.6 0.6])
hold on; plot(ts1, mean(concatenated_means_session_2(respClass_all_array{1, 2} == 1,:)));
% hold on; plot(ts1, mean(concatenated_means_session_3(respClass_all(3,:) == 1,:)));

figure; plot(ts1, mean(concatenated_means_session_1(respClass_all_array{1, 1} == 2,:)));
ylim([-0.6 0.6])
hold on; plot(ts1, mean(concatenated_means_session_2(respClass_all_array{1, 2} == 2,:)));

figure; plot(ts1, mean(concatenated_means_session_1(respClass_all_array{1, 1} == 3,:)));
ylim([-0.6 0.6])
hold on; plot(ts1, mean(concatenated_means_session_2(respClass_all_array{1, 2} == 3,:)));

%% plot activated neurons

load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
colormap(batlowW); % c1 = colorbar; 



figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(concatenated_means_session_1(respClass_all_array{1, 1} == 1,:)), mean(concatenated_sems_session_1(respClass_all_array{1, 1} == 1,:)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(concatenated_means_session_2(respClass_all_array{1, 2} == 1,:)), mean(concatenated_sems_session_2(respClass_all_array{1, 2} == 1,:)), 'lineProps', {'color', batlowW(100,:)});
% shadedErrorBar(ts1, mean(concatenated_means_session_3(respClass_all(3,:) == 1,:)), mean(concatenated_sems_session_3(respClass_all(3,:) == 1,:)), 'lineProps', {'color', batlowW(200,:)});
hold off





figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(concatenated_means_session_1(respClass_all_array{1, 1} == 2,:)), mean(concatenated_sems_session_1(respClass_all_array{1, 1} == 2,:)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(concatenated_means_session_2(respClass_all_array{1, 2} == 2,:)), mean(concatenated_sems_session_2(respClass_all_array{1, 2} == 2,:)), 'lineProps', {'color', batlowW(100,:)});
hold off

figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(concatenated_means_session_1(respClass_all_array{1, 1} == 3,:)), mean(concatenated_sems_session_1(respClass_all_array{1, 1} == 3,:)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(concatenated_means_session_2(respClass_all_array{1, 2} == 3,:)), mean(concatenated_sems_session_2(respClass_all_array{1, 2} == 3,:)), 'lineProps', {'color', batlowW(100,:)});
hold off



%% get correlation coefficients across all cells
% Assuming you have two arrays: concatenated_values_session_1 and concatenated_values_session_2

% Initialize an array to store the correlations
correlations = zeros(size(concatenated_means_session_1, 1), 1);

% Loop through each row and calculate the correlation
for i = 1:size(concatenated_means_session_1, 1)
    row1 = concatenated_means_session_1(i, :);
    row2 = concatenated_means_session_2(i, :);
    
    % Calculate the correlation coefficient
    correlations(i) = corr(row1', row2');
end

% 'correlations' now contains the correlation coefficients for each pair of rows

mean_corr = mean(correlations(:,1));


%% Create a scatterplot to show relationship between session 1 and session 2 activity
start_time = 0; % sub-window start time
end_time = 2; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Extract the corresponding columns from neuron_mean


sub_window_activity_session_1 = concatenated_means_session_1(:, sub_window_idx);
sub_window_activity_session_2 = concatenated_means_session_2(:, sub_window_idx);

mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);
mean_sub_window_activity_session_2 = mean(sub_window_activity_session_2, 2);

x = mean_sub_window_activity_session_1;
y = mean_sub_window_activity_session_2;


% Create a scatter plot
figure;
set(gcf,'Position',[100 100 200 500])
% Group 1: respClass_all(1,:) == 1 (Orange)
idx_group_1 = (respClass_all_array{1, 1} == 1);
scatter(x(idx_group_1), y(idx_group_1), 'o', 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); % Orange

hold on;

% % Group 2: respClass_all(1,:) == 2 (Light Blue)
% idx_group_2 = (respClass_all_array{1, 1} == 2);
% scatter(x(idx_group_2), y(idx_group_2), 'o', 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k'); % Light Blue

% % Group 3: respClass_all(1,:) == 3 (Light Grey)
% idx_group_3 = (respClass_all_array{1, 1} == 3);
% scatter(x(idx_group_3), y(idx_group_3), 'o', 'filled', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'k'); % Light Grey

% Add a regression line (You can keep this part unchanged)
coefficients = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);
plot(x_fit, y_fit, 'r');

% Calculate R-squared value (You can keep this part unchanged)
y_pred = polyval(coefficients, x);
ssr = sum((y_pred - mean(y)).^2);
sst = sum((y - mean(y)).^2);
r_squared = ssr / sst;

% Add R-squared value to the plot (You can keep this part unchanged)
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared)], 'FontSize', 12);

% Add labels and a legend (You can keep this part unchanged)
% xlabel('X-axis Label');
% ylabel('Y-axis Label');
% title('Scatter Plot with Regression Line and R^2 Value');
% legend('Group 1', 'Group 2', 'Group 3', 'Regression Line');
ylim([-1 1.5])
xlim([-1 1.5])
hold off;

%% Sort activity for heatmaps


% Extract the corresponding columns from neuron_mean

[~, sorted_idx_session_1] = sort(mean_sub_window_activity_session_1, 'descend');
[~, sorted_idx_session_2] = sort(mean_sub_window_activity_session_2, 'descend');

% Sort neuron_mean based on the sorted indices
ranked_neurons_session_1 = concatenated_means_session_1(sorted_idx_session_1, :);
ranked_neurons_session_2 = concatenated_means_session_2(sorted_idx_session_2, :);

figure;
% Generate the heatmap
imagesc(ts1, 1, ranked_neurons_session_1);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');

figure;
% Generate the heatmap
imagesc(ts1, 1, ranked_neurons_session_2);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');

%%

% Define the categories
category_1 = 1;
category_2 = 2;
category_3 = 3;

% Define the time window of interest

start_time = -4; % sub-window start time
end_time = 0; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Calculate the mean activity of each neuron within the time window
mean_activity = mean(concatenated_means_session_1(:,sub_window_idx), 2);


% Find the indices of the neurons in each category
indices_1 = find(respClass_all_array{1, 1} == category_1);
indices_2 = find(respClass_all_array{1, 1} == category_2);
indices_3 = find(respClass_all_array{1, 1} == category_3);

% Select the three most representative neurons from each category
[~, sorted_indices_1] = sort(mean_activity(indices_1), 'descend');
[~, sorted_indices_2] = sort(mean_activity(indices_2), 'descend');
[~, sorted_indices_3] = sort(mean_activity(indices_3), 'descend');


midway_for_neutral = size(sorted_indices_3, 1)/2;

representative_neurons_1 = indices_1(sorted_indices_1(1:3));
representative_neurons_2 = indices_2(sorted_indices_2(end-3:end));
representative_neurons_3 = indices_3(sorted_indices_3(midway_for_neutral-3:midway_for_neutral));

representative_neurons_all = [concatenated_means_session_1(representative_neurons_1,:); concatenated_means_session_1(representative_neurons_2,:); concatenated_means_session_1(representative_neurons_3,:)];

% Create a figure
figure;

% Plot the activity of each representative neuron as a separate line
hold on;
plot(sgolayfilt((representative_neurons_all(:,:))', 9, 21));
plot(sgolayfilt(concatenated_means_session_1(representative_neurons_1,:)', 9, 21), 'LineWidth', 2)';
plot(sgolayfilt(concatenated_means_session_1(representative_neurons_2,:)', 9, 21), 'LineWidth', 2)';
plot(sgolayfilt(concatenated_means_session_1(representative_neurons_3,:)', 9, 21), 'LineWidth', 2)';

% Add labels and a legend
xlabel('Time (s)');
ylabel('Activity');
legend('Category 1', 'Category 2', 'Category 3');

%%
% Create a figure
figure;

% Define an offset for each category
offset_1 = max(concatenated_means_session_1(representative_neurons_1,:),[],2);
offset_2 = max(concatenated_means_session_1(representative_neurons_2,:),[],2) + 0.5;
offset_3 = max(concatenated_means_session_1(representative_neurons_3,:),[],2) + 1;

offset_all = [offset_1; offset_2; offset_3];

% Plot the activity of each representative neuron as a separate line with offset
hold on;
for i = 1:size(representative_neurons_all, 1)
    plot(ts1, representative_neurons_all(i, :) - offset_all(i)+ sum(offset_all(1:i)), 'LineWidth', 2);
end



%%
figure;
hold on;
for i = 1:length(representative_neurons_1)
    plot(ts1, concatenated_means_session_1(representative_neurons_1(i),:), 'LineWidth', 2);
end
ylim([-1.5 2])
xlim([-1.5 2])
hold off;
figure;
hold on;
for i = 1:length(representative_neurons_2)
    plot(ts1, concatenated_means_session_1(representative_neurons_2(i),:), 'LineWidth', 2);
end
ylim([-1.5 2])
xlim([-1.5 2])
hold off;
figure;
hold on;
for i = 1:length(representative_neurons_3)
    plot(ts1, concatenated_means_session_1(representative_neurons_3(i),:), 'LineWidth', 2);
end
ylim([-1.5 2])
xlim([-1.5 2])
% Add labels and a legend

xlabel('Time (s)');
ylabel('Activity');
legend('Category 1', 'Category 2', 'Category 3');
hold off;





%% Create a scatterplot to show relationship between session 1 and session 2 activity


%CREATE SCATTER PLOT BASED ON SPECIFIC EVENTS - ASSUMING THEY ARE IN PAIRS.
%CHECK AND UPDATE START & END TIME DEPENDING ON EVENT OF INTEREST
paired_neurons = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 1;
start_time = 0; % sub-window start time
end_time = 2; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Extract the corresponding columns from neuron_mean


sub_window_activity_session_1 = concatenated_means_session_1(paired_neurons, sub_window_idx);
sub_window_activity_session_2 = concatenated_means_session_2(paired_neurons, sub_window_idx);

mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);
mean_sub_window_activity_session_2 = mean(sub_window_activity_session_2, 2);

x = mean_sub_window_activity_session_1;
y = mean_sub_window_activity_session_2;


% Create a scatter plot
figure;
set(gcf,'Position',[100 100 200 500])
% Group 1: respClass_all(1,:) == 1 (Orange)
% idx_group_1 = (respClass_all_array{1, 1} == 1);
scatter(x, y, 'o', 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); % Orange

hold on;

% % Group 2: respClass_all(1,:) == 2 (Light Blue)
% idx_group_2 = (respClass_all_array{1, 1} == 2);
% scatter(x(idx_group_2), y(idx_group_2), 'o', 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k'); % Light Blue

% % Group 3: respClass_all(1,:) == 3 (Light Grey)
% idx_group_3 = (respClass_all_array{1, 1} == 3);
% scatter(x(idx_group_3), y(idx_group_3), 'o', 'filled', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'k'); % Light Grey

% Add a regression line (You can keep this part unchanged)
coefficients = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);
plot(x_fit, y_fit, 'r');

% Calculate R-squared value (You can keep this part unchanged)
y_pred = polyval(coefficients, x);
ssr = sum((y_pred - mean(y)).^2);
sst = sum((y - mean(y)).^2);
r_squared = ssr / sst;

% Add R-squared value to the plot (You can keep this part unchanged)
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared)], 'FontSize', 12);

% Add labels and a legend (You can keep this part unchanged)
% xlabel('X-axis Label');
% ylabel('Y-axis Label');
% title('Scatter Plot with Regression Line and R^2 Value');
% legend('Group 1', 'Group 2', 'Group 3', 'Regression Line');
% ylim([0 1.1])
% xlim([0 1.1])
hold off;

%% Create a scatterplot to show relationship between session 1 and session 2 activity


%CREATE SCATTER PLOT BASED ON SPECIFIC EVENTS - ASSUMING THEY ARE IN PAIRS.
%CHECK AND UPDATE START & END TIME DEPENDING ON EVENT OF INTEREST
paired_neurons = respClass_all_array{1, 1} == 1 | respClass_all_array{1, 2} == 1;
session_1_neurons = find(respClass_all_array{1, 1} == 1);
session_2_neurons = find(respClass_all_array{1, 2} == 1);
both_sessions = intersect(session_1_neurons, session_2_neurons);
start_time = -4; % sub-window start time
end_time = 0; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Extract the corresponding columns from neuron_mean


sub_window_activity_session_1 = concatenated_means_session_1(:, sub_window_idx);
sub_window_activity_session_2 = concatenated_means_session_2(:, sub_window_idx);

mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);
mean_sub_window_activity_session_2 = mean(sub_window_activity_session_2, 2);

x = mean_sub_window_activity_session_1;
y = mean_sub_window_activity_session_2;


% Create a scatter plot
figure;
set(gcf,'Position',[100 100 200 500])
% Group 1: respClass_all(1,:) == 1 (Orange)
idx_group_1 = session_1_neurons;
scatter(x(idx_group_1), y(idx_group_1), 'o', 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); % Orange

hold on;

% Group 2: respClass_all(1,:) == 2 (Light Blue)
idx_group_2 = session_2_neurons;
scatter(x(idx_group_2), y(idx_group_2), 'o', 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k'); % Light Blue

% Group 3: respClass_all(1,:) == 3 (Light Grey)
idx_group_3 = both_sessions;
scatter(x(idx_group_3), y(idx_group_3), 'o', 'filled', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'k'); % Light Grey

% Add a regression line (You can keep this part unchanged)
coefficients = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);
plot(x_fit, y_fit, 'r');

% Calculate R-squared value (You can keep this part unchanged)
y_pred = polyval(coefficients, x);
ssr = sum((y_pred - mean(y)).^2);
sst = sum((y - mean(y)).^2);
r_squared = ssr / sst;

% Add R-squared value to the plot (You can keep this part unchanged)
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared)], 'FontSize', 12);

% Add labels and a legend (You can keep this part unchanged)
% xlabel('X-axis Label');
% ylabel('Y-axis Label');
% title('Scatter Plot with Regression Line and R^2 Value');
% legend('Group 1', 'Group 2', 'Group 3', 'Regression Line');
% ylim([0 1.1])
% xlim([0 1.1])
hold off;