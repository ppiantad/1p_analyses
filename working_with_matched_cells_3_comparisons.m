%%

%% Create a scatterplot to show relationship between session 1 and session 2 activity
start_time = -5; % sub-window start time
end_time = 0; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Extract the corresponding columns from neuron_mean


sub_window_activity_session_1 = neuron_mean_array{1, 1}(:, sub_window_idx);
sub_window_activity_session_2 = neuron_mean_array{1, 2}(:, sub_window_idx);
sub_window_activity_session_3 = neuron_mean_array{1, 3}(:, sub_window_idx);

mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);
mean_sub_window_activity_session_2 = mean(sub_window_activity_session_2, 2);
mean_sub_window_activity_session_3 = mean(sub_window_activity_session_3, 2);
%%
% Calculate AUC for each session's sub-window
auc_session_1 = trapz(sub_window_activity_session_1, 2);
auc_session_2 = trapz(sub_window_activity_session_2, 2);
auc_session_3 = trapz(sub_window_activity_session_3, 2);

% Perform a statistical test to compare the AUCs (e.g., ANOVA or Kruskal-Wallis)
alpha = 0.05; % Set your significance level
[p_value, ~, stats] = kruskalwallis([auc_session_1, auc_session_2, auc_session_3]);
if p_value < alpha
    fprintf('There is a significant difference in AUCs between sessions.\n');
    % You can also perform post-hoc tests if the Kruskal-Wallis test is significant.
    % For example, pairwise comparisons using the multcompare function.
    c = multcompare(stats, 'CType', 'bonferroni');
    disp(c); % Display pairwise comparisons
else
    fprintf('No significant difference in AUCs between sessions.\n');
end

%%
% Calculate the slope (coefficient) of a linear fit for each session
for ii = 1:size(sub_window_activity_session_1, 1)
    coeff_session_1(ii,:) = polyfit(ts1(sub_window_idx), sub_window_activity_session_1(ii,:), 1);
    coeff_session_2(ii,:) = polyfit(ts1(sub_window_idx), sub_window_activity_session_2(ii,:), 1);
    coeff_session_3(ii,:) = polyfit(ts1(sub_window_idx), sub_window_activity_session_3(ii,:), 1);
end
% coeff_session_1 = polyfit(ts1(sub_window_idx), sub_window_activity_session_1, 1);
% coeff_session_2 = polyfit(ts1(sub_window_idx), sub_window_activity_session_2, 1);
% coeff_session_3 = polyfit(ts1(sub_window_idx), sub_window_activity_session_3, 1);

% Extract the slope (coefficient) values
slope_session_1 = coeff_session_1(:,1);
slope_session_2 = coeff_session_2(:,1);
slope_session_3 = coeff_session_3(:,1);

% Perform statistical tests to compare the slopes
alpha = 0.05; % Set your significance level

% For example, you can use a t-test to compare slopes between sessions.
[h_12,p_value_12,ci_12,stats_12] = ttest2(slope_session_1, slope_session_2, 'Alpha', alpha);
[h_13,p_value_13,ci_13,stats_13] = ttest2(slope_session_1, slope_session_3, 'Alpha', alpha);
[h_23,p_value_23,ci_23,stats_23] = ttest2(slope_session_2, slope_session_3, 'Alpha', alpha);

% Check if any of the comparisons are statistically significant
if p_value_12 < alpha
    fprintf('Slope of session 1 is significantly different from session 2.\n');
end
if p_value_13 < alpha
    fprintf('Slope of session 1 is significantly different from session 3.\n');
end
if p_value_23 < alpha
    fprintf('Slope of session 2 is significantly different from session 3.\n');
end


% Calculate the mean and SEM for each session's slopes
mean_slope_1 = mean(slope_session_1);
sem_slope_1 = std(slope_session_1) / sqrt(length(slope_session_1));

mean_slope_2 = mean(slope_session_2);
sem_slope_2 = std(slope_session_2) / sqrt(length(slope_session_2));

mean_slope_3 = mean(slope_session_3);
sem_slope_3 = std(slope_session_3) / sqrt(length(slope_session_3));

% Create a bar graph with error bars
figure;
bar([mean_slope_1, mean_slope_2, mean_slope_3]);
hold on;
errorbar([1, 2, 3], [mean_slope_1, mean_slope_2, mean_slope_3], [sem_slope_1, sem_slope_2, sem_slope_3], 'r.', 'LineWidth', 1.5);
hold off;

% Customize the plot
xticklabels({'Session 1', 'Session 2', 'Session 3'});
ylabel('Mean Slope');
title('Mean Slope for Each Session with SEM');
grid on;

% Optionally, you can save the figure to a file
% saveas(gcf, 'slope_comparison.png');

%%

load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
% colormap(batlowW); % c1 = colorbar; 
% Create label vector y (corresponding to trial blocks)
numNeuronsPerCondition = neuron_num;
% change depending on the number of behaviors to decode!



numConditions = size(neuron_mean_concat, 2)/numMeasurements;


y = repelem(1:numConditions, numMeasurements);



figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(neuron_mean_concat(:,y == 1)), mean(neuron_sem_concat(:,y == 1)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(neuron_mean_concat(:,y == 2)), mean(neuron_sem_concat(:,y == 2)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(neuron_mean_concat(:,y == 3)), mean(neuron_sem_concat(:,y == 3)), 'lineProps', {'color', batlowW(1,:)});
hold off

%% NEED TO RUN eventRelated first! 

figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(concatenated_means_session_1(respClass_all_array{1,1} == 1,:)), mean(concatenated_sems_session_1(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(concatenated_means_session_2(respClass_all_array{1,2} == 1,:)), mean(concatenated_sems_session_2(respClass_all_array{1,2} == 1,:)), 'lineProps', {'color', batlowW(100,:)});
shadedErrorBar(ts1, mean(concatenated_means_session_3(respClass_all_array{1,3} == 1,:)), mean(concatenated_sems_session_3(respClass_all_array{1,3} == 1,:)), 'lineProps', {'color', batlowW(200,:)});
hold off


figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(concatenated_means_session_1(respClass_all(1,:) == 2,:)), mean(concatenated_sems_session_1(respClass_all(1,:) == 2,:)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(concatenated_means_session_2(respClass_all(2,:) == 2,:)), mean(concatenated_sems_session_2(respClass_all(2,:) == 2,:)), 'lineProps', {'color', batlowW(100,:)});
shadedErrorBar(ts1, mean(concatenated_means_session_3(respClass_all(3,:) == 2,:)), mean(concatenated_sems_session_3(respClass_all(3,:) == 2,:)), 'lineProps', {'color', batlowW(200,:)});
hold off


figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(concatenated_means_session_1(respClass_all(1,:) == 3,:)), mean(concatenated_sems_session_1(respClass_all(1,:) == 3,:)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(concatenated_means_session_2(respClass_all(2,:) == 3,:)), mean(concatenated_sems_session_2(respClass_all(2,:) == 3,:)), 'lineProps', {'color', batlowW(100,:)});
shadedErrorBar(ts1, mean(concatenated_means_session_3(respClass_all(3,:) == 3,:)), mean(concatenated_sems_session_3(respClass_all(3,:) == 3,:)), 'lineProps', {'color', batlowW(200,:)});
hold off



%%
figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(neuron_mean_array{1, 1}(respClass_all_array{1,1} == 1,:)), mean(neuron_sem_array{1, 1}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(respClass_all_array{1,2} == 1,:)), mean(neuron_sem_array{1, 2}(respClass_all_array{1,2} == 1,:)), 'lineProps', {'color', batlowW(100,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 3}(respClass_all_array{1,3} == 1,:)), mean(neuron_sem_array{1, 3}(respClass_all_array{1,3} == 1,:)), 'lineProps', {'color', batlowW(200,:)});
hold off


% plot a comparison where you filter on 3 different events, but you want to
% look at how the activity pattern in the first comparison looks based on
% data from the 2nd and 3rd comparison
figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(neuron_mean_array{1, 1}(respClass_all_array{1,1} == 1,:)), mean(neuron_sem_array{1, 1}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(respClass_all_array{1,1} == 1,:)), mean(neuron_sem_array{1, 2}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(100,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), mean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(200,:)});
hold off


% use if you want to plot data aligned to diff events based on 1 timeseries
% (e.g., trial start neurons, pre-choice neurons, consumption neurons. if
% you want to have the figure aligned to CHOICE, use the 2nd
% neuron_mean_array, if TrialFilter was used in the above order!)
figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(respClass_all_array{1,1} == 1,:)), mean(neuron_sem_array{1, 2}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(respClass_all_array{1,2} == 1,:)), mean(neuron_sem_array{1, 2}(respClass_all_array{1,2} == 1,:)), 'lineProps', {'color', batlowW(100,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(respClass_all_array{1,3} == 1,:)), mean(neuron_sem_array{1, 2}(respClass_all_array{1,3} == 1,:)), 'lineProps', {'color', batlowW(200,:)});
hold off

start_active = neuron_mean_array{1, 2}(respClass_all_array{1,1} == 1, :);
choice_active = neuron_mean_array{1, 2}(respClass_all_array{1,2} == 1,:);
consumption_active = neuron_mean_array{1, 2}(respClass_all_array{1,3} == 1,:);


start_time = -3; % sub-window start time
end_time = 0; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Extract the corresponding columns from neuron_mean
sub_window_activity = start_active(:, sub_window_idx);
mean_activity = mean(sub_window_activity, 2);

[~, sorted_idx] = sort(mean_activity, 'ascend');

% Sort neuron_mean based on the sorted indices
start_ranked_neurons = start_active(sorted_idx, :);
figure;
% Generate the heatmap
imagesc(ts1, 1, start_ranked_neurons);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
% set(gca, 'YDir', 'reverse');


start_time = -3; % sub-window start time
end_time = 0; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Extract the corresponding columns from neuron_mean
sub_window_activity = start_active(:, sub_window_idx);
mean_activity = mean(sub_window_activity, 2);

[~, sorted_idx] = sort(mean_activity, 'ascend');

% Sort neuron_mean based on the sorted indices
start_ranked_neurons = start_active(sorted_idx, :);
figure;
% Generate the heatmap
imagesc(ts1, 1, start_ranked_neurons);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');




array_for_heatmap = [neuron_mean_array{1, 2}(respClass_all_array{1,1} == 1,:); neuron_mean_array{1, 2}(respClass_all_array{1,2} == 1,:); neuron_mean_array{1, 2}(respClass_all_array{1,3} == 1,:)]
