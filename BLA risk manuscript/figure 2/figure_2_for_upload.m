

%% For Figure 2C

% Define the custom colormap from white to orange
% custom_colormap = [
%     1, 1, 1; % white
%     1, 0.9, 0.8;
%     1, 0.8, 0.6;
%     1, 0.7, 0.4;
%     1, 0.6, 0.2;
%     1, 0.5, 0; % orange
% ];


% custom_colormap = [
%     1, 1, 1;       % white
%     0.9, 0.95, 0.9;
%     0.8, 0.9, 0.8;
%     0.6, 0.8, 0.6;
%     0.4, 0.7, 0.4;
%     0.2, 0.6, 0.2;
%     0.13, 0.55, 0.13; % forest green
% ];

% custom_colormap = [
%     1, 1, 1;         % white
%     0.9, 0.95, 0.95;
%     0.8, 0.9, 0.9;
%     0.6, 0.85, 0.85;
%     0.4, 0.8, 0.8;
%     0.2, 0.8, 0.8;
%     0.0, 0.8, 0.8;   % robin's egg blue
% ];

custom_colormap = [
    1, 1, 1;         % white
    0.9, 0.9, 0.95;
    0.8, 0.8, 0.9;
    0.6, 0.6, 0.8;
    0.4, 0.4, 0.7;
    0.2, 0.2, 0.6;
    0.0, 0.0, 0.55;   % dark blue
];

% custom_colormap = [
%     1, 1, 1; % white
%     1, 0.9, 0.9;
%     1, 0.8, 0.8;
%     1, 0.7, 0.7;
%     1, 0.6, 0.6;
%     1, 0.5, 0.5;
%     1, 0.4, 0.4;
%     1, 0.3, 0.3;
%     1, 0.2, 0.2;
%     1, 0.1, 0.1;
%     1, 0, 0;   % red
% ];


% Generate more intermediate colors for a smoother transition
n = 256; % Number of colors
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));



% this is done for each ensemble

% for RDT D1 BLA_Insc_40:
%prechoice neuron num 12
%postchoice rew num 70
%consumption num 10
%shock num 11

pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [2 5];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent. if you want to show that consum neurons peak after collect, use something like 2-5

windows = {pre_choice_window, post_choice_window, consumption_window};
plot_num = [15, 70, 10];

array_to_plot = [1, 1, 1]; % depends on the structure of zall

select_mouse = 'BLA_Insc_40';

% for RDT D1 BLA_Insc_25:
%prechoice neuron num 46
%postchoice rew num 38
%consumption num 39
%shock num 11


% for Pre RDT RM BLA_Insc_40:
%prechoice neuron num 12
%postchoice rew num 70
%consumption num 10
%shock num 11


select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'Pre_RDT_RM';

second_session = 'RDT_D1';


% Get the number of neurons
num_neurons = length(plot_num);

% Create a new figure
figure;

for neuron_idx = 1:num_neurons
    % Get the current neuron and corresponding time window
    current_neuron = plot_num(neuron_idx);
    current_window = windows{neuron_idx};
    array_to_plot_current = array_to_plot(neuron_idx);
    % Get the data for the current neuron
    neuron_data = zall_mouse{select_mouse_index, array_to_plot_current}{1, current_neuron};
    neuron_data_s_array = caTraceTrials_spikes_mouse{select_mouse_index, array_to_plot_current}{1, current_neuron};
    % Restrict time indices based on the current window
    time_indices = ts1 >= current_window(1) & ts1 <= current_window(2);
    restricted_data = neuron_data(:, time_indices);
    
    % Compute the mean activity within the restricted window for each trial
    mean_activity = mean(restricted_data, 2);
    
    % Sort trials by mean activity in descending order
    [~, sorted_trial_indices] = sort(mean_activity, 'descend');
    
    % Select the top 5 trials with the highest mean activity
    selected_trials = sorted_trial_indices(1:min(5, size(neuron_data, 1)));
    
    % Determine Y-axis limits for this neuron
    y_limits = [min(neuron_data(:)), max(neuron_data(:))];
    
    for subplot_idx = 1:length(selected_trials)
        % Calculate subplot position
        subplot_idx_global = (neuron_idx - 1) * 5 + subplot_idx;
        subplot(num_neurons, 5, subplot_idx_global);
        hold on;
        
        % Get the selected trial
        trial = selected_trials(subplot_idx);
        
        % Plot the selected trial
        plot(ts1, neuron_data(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
        hold on; plot(ts1, neuron_data_s_array(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
        % % Plot the mean as a thick black line
        % meanData = mean(neuron_data);
        % plot(ts1, meanData, 'LineWidth', 2, 'Color', 'k');
        
        % Set axis limits and labels
        ylim(y_limits);
        xlim([-8 8]);
        set(gca, 'XTick', [-8, 0, 8]);
        
        % Add reference lines
        xline(0, 'k--');
        yline(0, 'k--');
        
        % % Add title for the subplot
        % if neuron_idx == 1
        %     title(['Trial ' num2str(trial)]);
        % end
        
        % Adjust font size
        set(gca, 'FontSize', 12);
        
        hold off;
    end
    
    % % Add a label to the Y-axis for the first column in each row
    % ylabel(['Neuron ' num2str(current_neuron)], 'FontSize', 14);
end


%% Fig. 2E


pre_choice_neurons = neuron_mean_array{1, 1}(prechoice_block_1, :);
post_choice_reward_neurons = neuron_mean_array{1, 1}(postchoice_reward_block_1, :);
consumption_neurons = neuron_mean_array{1, 1}(collect_block_1, :);

only_active_array_stacked = [pre_choice_neurons; post_choice_reward_neurons; consumption_neurons];

% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(pre_choice_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
pre_choice_neurons_sorted = pre_choice_neurons(sort_indices, :);



% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(post_choice_reward_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
post_choice_reward_neurons_sorted = post_choice_reward_neurons(sort_indices, :);


% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(consumption_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
consumption_neurons_sorted = consumption_neurons(sort_indices, :);

sorted_only_active_array_stacked = [pre_choice_neurons_sorted; post_choice_reward_neurons_sorted; consumption_neurons_sorted];

% Now, activated_neuron_mean_sorted contains the rows of neuron_mean filtered by respClass_all == 1
% and sorted by the time of peak activity.

figure;
% Generate the heatmap
imagesc(ts1, 1, sorted_only_active_array_stacked);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');
clim([-1 1])
xline(0);


% If you want to customize the color map, you can use colormap function
% For example, using a blue-white-red colormap:
colormap(gray);

% If you want to limit the color scale to the range [0, 1]
caxis([-1 1]); % Assuming correlations range from -1 to 1
% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 


%% Fig. 2F - separate decoding scripts


%% For Figure 2G 
% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav

figure;
pre_choice_active_ind = find(respClass_all_array{1,1} == 1);
consum_active_ind = find(respClass_all_array{1,3} == 1);
post_choice_active_ind = find(respClass_all_array{1,2} == 1);
% consum_inhibited_ind = find(all_consum_inhibited == 1);
setListData = {pre_choice_active_ind, consum_active_ind, post_choice_active_ind};
setLabels = ["Pre-choice excited", "Consumption excited", "Post-choice excited"];

h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

h.ShowIntersectionCounts = true;
h.ShowIntersectionAreas = true;
% h.SetLabels = [];


%% Fig. 2H
%These data can be used to plot the median or mean choice
% time on a PCA graph, for example

behav_tbl_iter_single = behav_tbl_iter(1);

% Initialize the concatenated table
concatenatedTable = table();

% Iterate through the 3x1 cell array
for i = 1:numel(behav_tbl_iter_single)
    % Assuming each cell contains a 12x1 cell array of tables
    twelveByOneCellArray = behav_tbl_iter_single{i};
    
    % Initialize a temporary table to store the concatenated tables for this cell
    tempTable = table();
    
    % Iterate through the 12x1 cell array
    for j = 1:numel(twelveByOneCellArray)
        % Assuming each cell in the 12x1 cell array contains a table
        currentTable = twelveByOneCellArray{j};
        
        % Concatenate the current table to the temporary table vertically
        tempTable = vertcat(tempTable, currentTable);
    end
    
    % Concatenate the temporary table to the overall concatenated table vertically
    concatenatedTable = vertcat(concatenatedTable, tempTable);
end

median_choice_time_block_1 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_choice_time_block_2 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_choice_time_block_3 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));

median_collect_time_block_1 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_collect_time_block_2 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_collect_time_block_3 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));





median_start_time_from_choice = median(concatenatedTable.stTime - concatenatedTable.choiceTime);
median_collect_time_from_choice = median(concatenatedTable.collectionTime - concatenatedTable.choiceTime);

median_start_time_from_choice_large = median(concatenatedTable.stTime(concatenatedTable.bigSmall == 1.2) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 1.2));
median_start_time_from_choice_small = median(concatenatedTable.stTime(concatenatedTable.bigSmall == 0.3) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 0.3));

median_collect_time_from_choice_large = median(concatenatedTable.collectionTime(concatenatedTable.bigSmall == 1.2) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 1.2));
median_collect_time_from_choice_small = median(concatenatedTable.collectionTime(concatenatedTable.bigSmall == 0.3) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 0.3));



% For Figure 2H (top)
figure;
hold on
% Create a histogram for allCorrelations

width = 350; % Width of the figure
height = 300; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8], 'YTick', [-0.6, 0, 0.6]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(prechoice_block_1==1, :)), std(neuron_mean_array{1, 1}(prechoice_block_1==1, :)/sqrt(size(neuron_mean_array{1, 1}(prechoice_block_1==1, :), 1))), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(postchoice_reward_block_1==1, :)), std(neuron_mean_array{1, 1}(postchoice_reward_block_1==1, :)/sqrt(size(neuron_mean_array{1, 1}(postchoice_reward_block_1==1, :), 1))), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(collect_block_1==1, :)), std(neuron_mean_array{1, 1}(collect_block_1==1, :)/sqrt(size(neuron_mean_array{1, 1}(collect_block_1==1, :), 1))), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
ylim([-0.6 0.6]);

hold off


% For Figure 2H (bottom)
% Get AUCs for the relevant periods for the 3 defined events
% Define time windows
pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [1 3];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

% Initialize arrays to store AUCs
% auc_pre_choice = zeros(size(neuron_mean_array));
% auc_post_choice = zeros(size(neuron_mean_array));
% auc_consumption = zeros(size(neuron_mean_array));
pre_choice_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,1}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if prechoice_block_1(i) == 1
        pre_choice_neuron_count = pre_choice_neuron_count+1;
        selected_data = neuron_mean_array{1, 1}(i, :);
        % % Extract time variable (assuming it's named 'ts1')
        % ts1_data = ts1{i}(exclusive_activated_session_1{i} == 1);

        % Find indices corresponding to each time window
        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        % Compute AUC for each time window
        % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
        action_auc_pre_choice(pre_choice_neuron_count) = trapz(selected_data(pre_choice_indices));
        action_auc_post_choice(pre_choice_neuron_count) = trapz(selected_data(post_choice_indices));
        action_auc_consumption(pre_choice_neuron_count) = trapz(selected_data(consumption_indices));
    else

    end
end

post_choice_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,1}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if postchoice_reward_block_1(i) == 1
        post_choice_neuron_count = post_choice_neuron_count+1;
        selected_data = neuron_mean_array{1, 1}(i, :);
        % % Extract time variable (assuming it's named 'ts1')
        % ts1_data = ts1{i}(exclusive_activated_session_1{i} == 1);

        % Find indices corresponding to each time window
        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        % Compute AUC for each time window
        % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
        post_choice_reward_auc_pre_choice(post_choice_neuron_count) = trapz(selected_data(pre_choice_indices));
        post_choice_reward_auc_post_choice(post_choice_neuron_count) = trapz(selected_data(post_choice_indices));
        post_choice_reward_auc_consumption(post_choice_neuron_count) = trapz(selected_data(consumption_indices));
    else

    end
end


consumption_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,3}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if collect_block_1(i) == 1
        consumption_neuron_count = consumption_neuron_count+1;
        selected_data = neuron_mean_array{1, 3}(i, :);
        % % Extract time variable (assuming it's named 'ts1')
        % ts1_data = ts1{i}(exclusive_activated_session_1{i} == 1);

        % Find indices corresponding to each time window
        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        % Compute AUC for each time window
        % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
        consumption_auc_pre_choice(consumption_neuron_count) = trapz(selected_data(pre_choice_indices));
        consumption_auc_post_choice(consumption_neuron_count) = trapz(selected_data(post_choice_indices));
        consumption_auc_consumption(consumption_neuron_count) = trapz(selected_data(consumption_indices));
    else

    end
end

% Calculate mean and SEM for pre-choice period
mean_pre_choice = [mean(action_auc_pre_choice(:)), mean(post_choice_reward_auc_pre_choice(:)), mean(consumption_auc_pre_choice(:))];
sem_pre_choice = [std(action_auc_pre_choice(:))/sqrt(numel(action_auc_pre_choice)), std(post_choice_reward_auc_pre_choice(:))/sqrt(numel(post_choice_reward_auc_pre_choice)), std(consumption_auc_pre_choice(:))/sqrt(numel(consumption_auc_pre_choice))];

% Calculate mean and SEM for post-choice period
mean_post_choice = [mean(action_auc_post_choice(:)), mean(post_choice_reward_auc_post_choice(:)), mean(consumption_auc_post_choice(:))];
sem_post_choice = [std(action_auc_post_choice(:))/sqrt(numel(action_auc_post_choice)), std(post_choice_reward_auc_post_choice(:))/sqrt(numel(post_choice_reward_auc_post_choice)), std(consumption_auc_post_choice(:))/sqrt(numel(consumption_auc_post_choice))];

% Calculate mean and SEM for consumption period
mean_consumption = [mean(action_auc_consumption(:)), mean(post_choice_reward_auc_consumption(:)), mean(consumption_auc_consumption(:))];
sem_consumption = [std(action_auc_consumption(:))/sqrt(numel(action_auc_consumption)), std(post_choice_reward_auc_consumption(:))/sqrt(numel(post_choice_reward_auc_consumption)), std(consumption_auc_consumption(:))/sqrt(numel(consumption_auc_consumption))];

% Plot the bar graph
figure;
width = 350; % Width of the figure
height = 200; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
set(gca, 'YTick', [-10, 0, 10]);
bar_groups = 1:3; % Number of groups
bar_width = 0.3; % Width of each bar (adjust as needed)
hold on;
bar(bar_groups - bar_width, mean_pre_choice, bar_width, 'b');
bar(bar_groups, mean_post_choice, bar_width, 'g');
bar(bar_groups + bar_width, mean_consumption, bar_width, 'r');

% Add error bars
errorbar(bar_groups - bar_width, mean_pre_choice, sem_pre_choice, 'k', 'LineStyle', 'none');
errorbar(bar_groups, mean_post_choice, sem_post_choice, 'k', 'LineStyle', 'none');
errorbar(bar_groups + bar_width, mean_consumption, sem_consumption, 'k', 'LineStyle', 'none');

% Customize the plot
xlabel('AUC Periods');
ylabel('Mean AUC');
title('Mean AUC for Different Periods');
legend('Pre-choice', 'Post-choice', 'Consumption');
xticks(bar_groups);
xticklabels({'Set 1', 'Set 2', 'Set 3'}); % Replace with actual labels
grid off; % Remove grid lines
hold off;

% Combine the data into a single matrix
action_data = [action_auc_pre_choice(:); action_auc_post_choice(:); action_auc_consumption(:)];

% Create group labels
groups = [ones(size(action_auc_consumption(:)));  % Group 1: action_auc_consumption
          2*ones(size(action_auc_post_choice(:))); % Group 2: action_auc_post_choice
          3*ones(size(action_auc_pre_choice(:)))]; % Group 3: action_auc_pre_choice

figure;
% Perform one-way ANOVA
[p_anova, tbl, stats] = anova1(action_data, groups, 'off');

% Perform post-hoc Tukey's HSD test
c = multcompare(stats, 'CType', 'tukey-kramer');

% Extract p-values and group means
p_values = c(:, 6);
group_means = c(:, 4);

% Set significance threshold
alpha = 0.05;

% Determine significant pairs
significant_pairs = c(p_values < alpha, 1:2);

% Display results
disp('Significant pairwise comparisons:');
for i = 1:size(significant_pairs, 1)
    fprintf('Group %d vs Group %d\n', significant_pairs(i, 1), significant_pairs(i, 2));
end

% Combine the data into a single matrix
postchoice_data = [post_choice_reward_auc_pre_choice(:); post_choice_reward_auc_post_choice(:); post_choice_reward_auc_consumption(:)];

% Create group labels
groups = [ones(size(post_choice_reward_auc_consumption(:)));  % Group 1: action_auc_consumption
          2*ones(size(post_choice_reward_auc_post_choice(:))); % Group 2: action_auc_post_choice
          3*ones(size(post_choice_reward_auc_pre_choice(:)))]; % Group 3: action_auc_pre_choice

figure;
% Perform one-way ANOVA
[p_anova, tbl, stats] = anova1(postchoice_data, groups, 'off');

% Perform post-hoc Tukey's HSD test
c = multcompare(stats, 'CType', 'tukey-kramer');

% Extract p-values and group means
p_values = c(:, 6);
group_means = c(:, 4);

% Set significance threshold
alpha = 0.05;

% Determine significant pairs
significant_pairs = c(p_values < alpha, 1:2);

% Display results
disp('Significant pairwise comparisons:');
for i = 1:size(significant_pairs, 1)
    fprintf('Group %d vs Group %d\n', significant_pairs(i, 1), significant_pairs(i, 2));
end






%% For Figure 2I (left)
% create correlation matrix heatmap with exclusively active cells
test = [];
test = [neuron_mean_array{1, 1}(prechoice_block_1 == 1, :)];
% test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 1}~=1 & respClass_all_array{1, 2}==1 & respClass_all_array{1, 3}~=1, :)];
test = [test; neuron_mean_array{1, 1}(postchoice_reward_block_1 == 1, :)];
test = [test; neuron_mean_array{1, 1}(collect_block_1 == 1,:)];
% test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 1}~=1 & respClass_all_array{1, 2}~=1 & respClass_all_array{1,3}~=1, :)];

pre_choice_index = [1:sum(prechoice_block_1)];
post_choice_index = [pre_choice_index(end)+1:pre_choice_index(end)+sum(postchoice_reward_block_1)];
consumption_index = [post_choice_index(end)+1:post_choice_index(end)+sum(collect_block_1)];
neutral_index = [consumption_index(end)+1:consumption_index(end)+sum(respClass_all_array{1, 2}~=1 & respClass_all_array{1, 1}~=1 & respClass_all_array{1,3}~=1)];


data = test;

alpha = 0.0001;


correlation_matrix = zeros(size(data, 1));
p_value_matrix = zeros(size(data, 1));


for i = 1:size(data, 1)
    for j = 1:size(data, 1)
        [corr_coeff, p_value] = corrcoef(data(i, :)', data(j, :)');
        correlation_matrix(i, j) = corr_coeff(1, 2); % Store correlation coefficient
        p_value_matrix(i, j) = p_value(1, 2); % Store p-value
    end
end


figure;
imagesc(correlation_matrix);

axis square; %

ylim([1  size(test, 1)])

xticks([1  size(test, 1)]);
yticks([1  size(test, 1)]);

colormap(gray);


caxis([-1 1]); % assuming correlations range from -1 to 1

c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

%% For Figure 2I (right)
action_p_value_matrix = p_value_matrix(pre_choice_index, pre_choice_index);
action_correl_matrix = correlation_matrix(pre_choice_index, pre_choice_index);

n = size(pre_choice_index, 2); 
k = 2;   

num_combinations = nchoosek(n, k);
disp(['Number of distinct pairwise combinations: ', num2str(num_combinations)]);



action_positive_count = 0;
action_negative_count = 0;
action_no_correlation_count = 0;


matrix_size = size(action_correl_matrix, 1);
uu = 1;

for i = 1:matrix_size
    for j = i+1:matrix_size 

        action_ensemble_corr_overall(uu) = action_correl_matrix(i, j);
        uu = uu+1;
        % 
        if action_p_value_matrix(i, j) < alpha
            if action_correl_matrix(i, j) > 0
                action_positive_count = action_positive_count + 1;
            elseif action_correl_matrix(i, j) < 0
                action_negative_count = action_negative_count + 1;
            end
        else
            action_no_correlation_count = action_no_correlation_count + 1;
        end
    end
end

action_comparisons_possible = [action_positive_count + action_negative_count + action_no_correlation_count];


disp(['Number of positive correlations: ', num2str(action_positive_count)]);
disp(['Number of negative correlations: ', num2str(action_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(action_no_correlation_count)]);


action_data = [(action_positive_count/action_comparisons_possible)*100, (action_negative_count/action_comparisons_possible)*100, (action_no_correlation_count/action_comparisons_possible)*100];



%
post_choice_p_value_matrix = p_value_matrix(post_choice_index, post_choice_index);
post_choice_correl_matrix = correlation_matrix(post_choice_index, post_choice_index);

n = size(post_choice_index, 2); 
k = 2;   

num_combinations = nchoosek(n, k);
disp(['Number of distinct pairwise combinations: ', num2str(num_combinations)]);



post_choice_positive_count = 0;
post_choice_negative_count = 0;
post_choice_no_correlation_count = 0;


matrix_size = size(post_choice_correl_matrix, 1);
uu = 1

for i = 1:matrix_size
    for j = i+1:matrix_size % Start from i+1 to exclude the diagonal
        post_choice_ensemble_corr_overall(uu) = post_choice_correl_matrix(i, j);
        uu = uu+1;
        % Check if p-value is less than 0.01
        if post_choice_p_value_matrix(i, j) < alpha
            if post_choice_correl_matrix(i, j) > 0
                post_choice_positive_count = post_choice_positive_count + 1;
            elseif post_choice_correl_matrix(i, j) < 0
                post_choice_negative_count = post_choice_negative_count + 1;
            end
        else
            post_choice_no_correlation_count = post_choice_no_correlation_count + 1;
        end
    end
end

post_choice_comparisons_possible = [post_choice_positive_count + post_choice_negative_count + post_choice_no_correlation_count];


disp(['Number of positive correlations: ', num2str(post_choice_positive_count)]);
disp(['Number of negative correlations: ', num2str(post_choice_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(post_choice_no_correlation_count)]);

post_choice_data = [(post_choice_positive_count/post_choice_comparisons_possible)*100, (post_choice_negative_count/post_choice_comparisons_possible)*100, (post_choice_no_correlation_count/post_choice_comparisons_possible)*100];

consumption_p_value_matrix = p_value_matrix(consumption_index, consumption_index);
consumption_correl_matrix = correlation_matrix(consumption_index, consumption_index);

n = size(consumption_index, 2); 
k = 2;  

num_combinations = nchoosek(n, k);
disp(['Number of distinct pairwise combinations: ', num2str(num_combinations)]);



consumption_positive_count = 0;
consumption_negative_count = 0;
consumption_no_correlation_count = 0;

% Get the size of the correlation matrix
matrix_size = size(consumption_correl_matrix, 1);
uu = 1

for i = 1:matrix_size
    for j = i+1:matrix_size 
        consumption_ensemble_corr_overall(uu) = consumption_correl_matrix(i, j);
        uu = uu+1;
        % 
        if consumption_p_value_matrix(i, j) < alpha
            if consumption_correl_matrix(i, j) > 0
                consumption_positive_count = consumption_positive_count + 1;
            elseif consumption_correl_matrix(i, j) < 0
                consumption_negative_count = consumption_negative_count + 1;
            end
        else
            consumption_no_correlation_count = consumption_no_correlation_count + 1;
        end
    end
end

consumption_comparisons_possible = [consumption_positive_count + consumption_negative_count + consumption_no_correlation_count];


disp(['Number of positive correlations: ', num2str(consumption_positive_count)]);
disp(['Number of negative correlations: ', num2str(consumption_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(consumption_no_correlation_count)]);


consumption_data = [(consumption_positive_count/consumption_comparisons_possible)*100, (consumption_negative_count/consumption_comparisons_possible)*100, (consumption_no_correlation_count/consumption_comparisons_possible)*100];

action_post_choice_p_value_matrix = p_value_matrix(pre_choice_index, post_choice_index);
action_post_choice_correl_matrix = correlation_matrix(pre_choice_index, post_choice_index);

n1 = size(action_post_choice_p_value_matrix, 1); 
n2 = size(action_post_choice_p_value_matrix, 2); 
k = 2;    

num_combinations = nchoosek(n1, k) * nchoosek(n2, k);
disp(['Number of unique pairwise combinations: ', num2str(num_combinations)]);


action_post_choice_positive_count = 0;
action_post_choice_negative_count = 0;
action_post_choice_no_correlation_count = 0;

[num_neurons_1, num_neurons_2] = size(action_post_choice_correl_matrix);


for i = 1:num_neurons_1
    for j = 1:num_neurons_2
        correlation = action_post_choice_correl_matrix(i, j);
        p_value = action_post_choice_p_value_matrix(i, j);
        if p_value < alpha
            if correlation > 0
                action_post_choice_positive_count = action_post_choice_positive_count + 1;
            elseif correlation < 0
                action_post_choice_negative_count = action_post_choice_negative_count + 1;
            end
        else
            action_post_choice_no_correlation_count = action_post_choice_no_correlation_count + 1;
        end
    end
end


disp(['Number of positive correlations: ', num2str(action_post_choice_positive_count)]);
disp(['Number of negative correlations: ', num2str(action_post_choice_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(action_post_choice_no_correlation_count)]);

action_post_choice_comparisons_possible = [action_post_choice_positive_count+action_post_choice_negative_count+action_post_choice_no_correlation_count];



action_post_choice_data = [(action_post_choice_positive_count/action_post_choice_comparisons_possible)*100, (action_post_choice_negative_count/action_post_choice_comparisons_possible)*100, (action_post_choice_no_correlation_count/action_post_choice_comparisons_possible)*100];

action_consumption_p_value_matrix = p_value_matrix(pre_choice_index, consumption_index);
action_consumption_correl_matrix = correlation_matrix(pre_choice_index, consumption_index);

n1 = size(action_consumption_p_value_matrix, 1); 
n2 = size(action_consumption_p_value_matrix, 2);
k = 2;   

num_combinations = nchoosek(n1, k) * nchoosek(n2, k);
disp(['Number of unique pairwise combinations: ', num2str(num_combinations)]);




action_consumption_positive_count = 0;
action_consumption_negative_count = 0;
action_consumption_no_correlation_count = 0;


[num_neurons_1, num_neurons_2] = size(action_consumption_correl_matrix);

for i = 1:num_neurons_1
    for j = 1:num_neurons_2
        correlation = action_consumption_correl_matrix(i, j);
        p_value = action_consumption_p_value_matrix(i, j);
        if p_value < alpha
            if correlation > 0
                action_consumption_positive_count = action_consumption_positive_count + 1;
            elseif correlation < 0
                action_consumption_negative_count = action_consumption_negative_count + 1;
            end
        else
            action_consumption_no_correlation_count = action_consumption_no_correlation_count + 1;
        end
    end
end


disp(['Number of positive correlations: ', num2str(action_consumption_positive_count)]);
disp(['Number of negative correlations: ', num2str(action_consumption_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(action_consumption_no_correlation_count)]);

action_consumption_comparisons_possible = [action_consumption_positive_count+action_consumption_negative_count+action_consumption_no_correlation_count];


action_consumption_data = [(action_consumption_positive_count/action_consumption_comparisons_possible)*100, (action_consumption_negative_count/action_consumption_comparisons_possible)*100, (action_consumption_no_correlation_count/action_consumption_comparisons_possible)*100];



post_choice_consumption_p_value_matrix = p_value_matrix(post_choice_index, consumption_index);
post_choice_consumption_correl_matrix = correlation_matrix(post_choice_index, consumption_index);

n1 = size(post_choice_consumption_p_value_matrix, 1);
n2 = size(post_choice_consumption_p_value_matrix, 2); 
k = 2;    

num_combinations = nchoosek(n1, k) * nchoosek(n2, k);
disp(['Number of unique pairwise combinations: ', num2str(num_combinations)]);


post_choice_consumption_positive_count = 0;
post_choice_consumption_negative_count = 0;
post_choice_consumption_no_correlation_count = 0;


[num_neurons_1, num_neurons_2] = size(post_choice_consumption_correl_matrix);

for i = 1:num_neurons_1
    for j = 1:num_neurons_2
        correlation = post_choice_consumption_correl_matrix(i, j);
        p_value = post_choice_consumption_p_value_matrix(i, j);
        if p_value < alpha
            if correlation > 0
                post_choice_consumption_positive_count = post_choice_consumption_positive_count + 1;
            elseif correlation < 0
                post_choice_consumption_negative_count = post_choice_consumption_negative_count + 1;
            end
        else
            post_choice_consumption_no_correlation_count = post_choice_consumption_no_correlation_count + 1;
        end
    end
end


disp(['Number of positive correlations: ', num2str(post_choice_consumption_positive_count)]);
disp(['Number of negative correlations: ', num2str(post_choice_consumption_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(post_choice_consumption_no_correlation_count)]);

post_choice_consumption_comparisons_possible = [post_choice_consumption_positive_count+post_choice_consumption_negative_count+post_choice_consumption_no_correlation_count];



post_choice_consumption_data = [(post_choice_consumption_positive_count/post_choice_consumption_comparisons_possible)*100, (post_choice_consumption_negative_count/post_choice_consumption_comparisons_possible)*100, (post_choice_consumption_no_correlation_count/post_choice_consumption_comparisons_possible)*100];
%

x = 1:6;
figure;
barh(x, [action_data; post_choice_data; consumption_data; action_post_choice_data; action_consumption_data; post_choice_consumption_data], 'stacked');
legend('Positive correlation', 'Negative correlation', 'No sig correlation');



%% Figure 2J, Figure 2K, and Figure 2L

array_for_means = 3; 


for q = 1:length (behav_tbl_iter{array_for_means, 1})
    nestedCellArray_1 = behav_tbl_iter{array_for_means, 1}{q};
    if ~isempty(nestedCellArray_1)
        % nestedCellArray_2 = behav_tbl_iter{2, 1}{q};
        % if size(nestedCellArray_1, 1) > size(nestedCellArray_2, 1)
        %     delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime(1:end-1,:);
        % else
        %     delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
        % end
        %

        for zz = 1:size(nestedCellArray_1, 1)
            valid_start_times = nestedCellArray_1.stTime(2:end);
            valid_choice_times = nestedCellArray_1.choiceTime(1:end-1);
            delay_to_initiation = valid_start_times - valid_choice_times;
            trial_types = nestedCellArray_1.bigSmall;
        end

        trial_choice_times = nestedCellArray_1.choiceTime - nestedCellArray_1.stTime;
        % delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
        delay_to_collect_post_shk = nestedCellArray_1.collectionTime - nestedCellArray_1.choiceTime;
        
        trial_choice_times_by_mouse{q} = trial_choice_times;
        delay_to_initiation_by_mouse{q} = delay_to_initiation;
        delay_to_collect_post_shk_by_mouse{q} = delay_to_collect_post_shk;
        trial_types_by_mouse{q} = trial_types;
        consum_times = nestedCellArray_1.collectionTime_end - nestedCellArray_1.collectionTime;
        consum_times_by_mouse{q} = consum_times;
        clear trial_choice_times delay_to_initiation delay_to_collect_post_shk trial_types consum_times
    end


end

trial_types_concat = cat(1, trial_types_by_mouse{:});
trial_choice_times_concat = cat(1, trial_choice_times_by_mouse{:});
rew_collect_times_concat = cat(1, delay_to_collect_post_shk_by_mouse{:});
consum_times_concat = cat(1, consum_times_by_mouse{:});

bar_separation_value = 3;

figure;

% Adjust figure dimensions for subplots
width = 400; % Adjusted width for 3 subplots
height = 500; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]

% Subplot 1: trial_choice_times_concat
subplot(1, 3, 1); % First subplot
hold on;
% Define colors based on trial_types_concat
colors = repmat([0.5, 0.5, 0.5], length(trial_choice_times_concat), 1); % Default to gray
colors(trial_types_concat == 1.2, :) = repmat([0, 0, 1], sum(trial_types_concat == 1.2), 1); % Blue for trial_types_concat == 1.2
colors(trial_types_concat == 0.3, :) = repmat([1, 0, 0], sum(trial_types_concat == 0.3), 1); % Red for trial_types_concat == 0.3

% Create the swarmchart with colors
s = swarmchart(ones(1, length(trial_choice_times_concat)), trial_choice_times_concat, 36, colors, 'filled');
s.MarkerEdgeColor = 'k'; % Black stroke for all circles

% Calculate and plot means for each group
mean_12 = mean(trial_choice_times_concat(trial_types_concat == 1.2));
mean_03 = mean(trial_choice_times_concat(trial_types_concat == 0.3));
plot([0.8, 1.2], [mean_12, mean_12], 'b-', 'LineWidth', 2); % Blue line for trial_types_concat == 1.2
plot([0.8, 1.2], [mean_03, mean_03], 'r-', 'LineWidth', 2); % Red line for trial_types_concat == 0.3


[h_choice_times p_choice_times, ~, stats_choice_times] = ttest2(trial_choice_times_concat(trial_types_concat == 1.2), trial_choice_times_concat(trial_types_concat == 0.3))

% tests for normality for reviewers
[h_choice_times_large_normality, p_choicetimes_large_normality] = kstest(trial_choice_times_concat(trial_types_concat == 1.2))
[h,p,adstat,cv] = adtest(trial_choice_times_concat(trial_types_concat == 1.2))
[h_choice_times_small_normality, p_choicetimes_small_normality] = kstest(trial_choice_times_concat(trial_types_concat == 0.3))
[h,p,adstat,cv] = adtest(trial_choice_times_concat(trial_types_concat == 0.3))

normalitytest(trial_choice_times_concat(trial_types_concat == 1.2)')
normalitytest(trial_choice_times_concat(trial_types_concat == 0.3)')

[p_rank_sum_choice_times h_rank_sum_choice_times stats_rank_sum_choice_times] = ranksum(trial_choice_times_concat(trial_types_concat == 1.2), trial_choice_times_concat(trial_types_concat == 0.3))

U_choice_times = (stats_rank_sum_choice_times.ranksum  ) - [size(trial_choice_times_concat(trial_types_concat == 1.2), 1)*(size(rew_collect_times_concat(trial_types_concat == 1.2), 1)+1)/2];


normalitytest(rew_collect_times_concat(trial_types_concat == 1.2)')
normalitytest(rew_collect_times_concat(trial_types_concat == 0.3)')

[p_rank_sum_collect_times h_rank_sum_collect_times stats_rank_sum_collect_times] = ranksum(rew_collect_times_concat(trial_types_concat == 1.2), rew_collect_times_concat(trial_types_concat == 0.3))

U_collect_times = (stats_rank_sum_collect_times.ranksum  ) - [size(rew_collect_times_concat(trial_types_concat == 1.2), 1)*(size(rew_collect_times_concat(trial_types_concat == 1.2), 1)+1)/2];



hold off;
title('Trial Choice Times');
xlabel('Trial');
ylabel('Time');
yline(0); % Add yline if needed
xtickformat('%.1f');
ytickformat('%.1f');

% Subplot 2: rew_collect_times_concat
subplot(1, 3, 2); % Second subplot
hold on;
% Define colors based on trial_types_concat
colors = repmat([0.5, 0.5, 0.5], length(rew_collect_times_concat), 1); % Default to gray
colors(trial_types_concat == 1.2, :) = repmat([0, 0, 1], sum(trial_types_concat == 1.2), 1); % Blue for trial_types_concat == 1.2
colors(trial_types_concat == 0.3, :) = repmat([1, 0, 0], sum(trial_types_concat == 0.3), 1); % Red for trial_types_concat == 0.3

% Create the swarmchart with colors
s = swarmchart(ones(1, length(rew_collect_times_concat)) * bar_separation_value, rew_collect_times_concat, 36, colors, 'filled');
s.MarkerEdgeColor = 'k'; % Black stroke for all circles

% Calculate and plot means for each group
mean_12 = mean(rew_collect_times_concat(trial_types_concat == 1.2));
mean_03 = mean(rew_collect_times_concat(trial_types_concat == 0.3));
plot([bar_separation_value - 0.2, bar_separation_value + 0.2], [mean_12, mean_12], 'b-', 'LineWidth', 2); % Blue line for trial_types_concat == 1.2
plot([bar_separation_value - 0.2, bar_separation_value + 0.2], [mean_03, mean_03], 'r-', 'LineWidth', 2); % Red line for trial_types_concat == 0.3

[h_rew_times p_rew_times, ~, stats_rew_times] = ttest2(rew_collect_times_concat(trial_types_concat == 1.2), rew_collect_times_concat(trial_types_concat == 0.3))

hold off;
title('Reward Collection Times');
xlabel('Reward Collection');
ylabel('Time');
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');


variable_to_correlate = consum_times_by_mouse;


%



meanZallMouse = cell(size(zall_mouse, 2), 1);

% Define the time range
% timeRange = (ts1 >= -4) & (ts1 <= 0);
% timeRange = (ts1 >= 0) & (ts1 <= 2);
timeRange = (ts1 >= 1) & (ts1 <= 3);



for i = 1:length(zall_mouse)

    nestedCellArray_1 = zall_mouse{i, array_for_means};
    nestedCellArray_2 = zall_mouse{i, 1};

    meanNestedCellArray = cell(size(nestedCellArray_1));
    
    for j = 1:length(nestedCellArray_1)
     
        currentArray = nestedCellArray_1{j};
        comparisonArray_for_size = nestedCellArray_2{j};
        
        if isequal(variable_to_correlate, delay_to_collect_post_shk_by_mouse)
            currentArray = nestedCellArray_1{j};
        else
            if size(currentArray, 1) > size(comparisonArray_for_size, 1)
                currentArray = currentArray(1:end-1,:);
            else

            end
        end
        % uncomment below if you want to mean center
        % currentArray_mean = mean(currentArray, 2);
        % currentArray = currentArray-currentArray_mean;
        % Compute the mean activity for each row in the time range 0 to 2 seconds
        meanValues = mean(currentArray(:, timeRange), 2);
        % meanValues = max(currentArray(:, timeRange), [], 2);
        
        meanNestedCellArray{j} = meanValues;
    end
    
    meanZallMouse{i} = meanNestedCellArray;
end





%


correlationResults = cell(size(meanZallMouse));
correlationResults_sig = cell(size(meanZallMouse));



for i = 1:length(meanZallMouse)

    meanNestedCellArray = meanZallMouse{i};
    

    correlationNestedArray = zeros(size(meanNestedCellArray));
    corr_sig_NestedArray = zeros(size(meanNestedCellArray));
    trialIndex = mod(i-1, length(variable_to_correlate)) + 1;
    

    trialChoiceTimes = variable_to_correlate{i}';
    % trialChoiceTimes = variable_to_correlate{i};

    % trialChoiceTimes =  trialChoiceTimes(trial_types_by_mouse{1, i} == 1.2);
        

    for j = 1:length(meanNestedCellArray)

        meanValues = meanNestedCellArray{j};
        
        % meanValues = meanValues(trial_types_by_mouse{1, i} == 1.2)
        



        if length(trialChoiceTimes) == length(meanValues)

            [correlationCoeff, corr_sig_vals] = corr(meanValues, trialChoiceTimes(:));
        elseif length(trialChoiceTimes) < length(meanValues)
            [correlationCoeff, corr_sig_vals] = corr(meanValues(1:end-1), trialChoiceTimes(:));
        else
          
            [correlationCoeff, corr_sig_vals] = NaN;
        end
        
        
        correlationNestedArray(j) = correlationCoeff;
        corr_sig_NestedArray(j) = corr_sig_vals;
    end
    clear meanValues

    correlationResults{i} = correlationNestedArray;
    correlationResults_sig{i} = corr_sig_NestedArray;
end


%
% Assuming correlationResults is defined and contains the correlation coefficients

% Initialize an empty array to collect all correlation coefficients
allCorrelations = [];

% Iterate through each level of correlationResults
for i = 1:length(correlationResults)
    % Get the current nested cell array of correlation coefficients
    correlationNestedArray = correlationResults{i};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(correlationNestedArray)
        % Get the current correlation coefficient
        correlationCoeff = correlationNestedArray(j);
        
        % Check if the coefficient is not NaN (if applicable)
        if ~isnan(correlationCoeff)
            % Append the coefficient to the allCorrelations array
            allCorrelations = [allCorrelations; correlationCoeff];
        end
    end
end

% Now, allCorrelations contains all the correlation coefficients
% Create a histogram of the correlation coefficients
figure;
histogram(allCorrelations);
xlabel('Correlation Coefficient');
ylabel('Frequency');
title('Histogram of Correlation Coefficients');

% Optionally, you can add a vertical line at 0 for reference
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;

% SHK responsive neurons assumed to be stored in respClass_all_array{1, 1} for this purpose - change as necessary
% only_shk_responsive_corrs = allCorrelations(kmeans_idx' == 3);
% only_shk_responsive_corrs = allCorrelations(prechoice_block_1 == 1);
% only_shk_responsive_corrs = allCorrelations(postchoice_reward_block_1 == 1);
only_shk_responsive_corrs = allCorrelations(collect_block_1 == 1);
% only_shk_responsive_corrs = allCorrelations(prechoice_blocks_2_and_3 == 1);
% not_shk_responsive_corrs = allCorrelations(prechoice_block_1 ~=1);
% not_shk_responsive_corrs = allCorrelations(kmeans_idx' ~= 3);c
not_shk_responsive_corrs = allCorrelations(true_neutral ==1);

% Now, allCorrelations contains all the correlation coefficients
% Create a histogram of the correlation coefficients
figure;
histogram(only_shk_responsive_corrs);
xlabel('Correlation Coefficient');
ylabel('Frequency');
title('Histogram of Correlation Coefficients');

% Optionally, you can add a vertical line at 0 for reference
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;

%
% Assuming the following variables are defined:
% allCorrelations: array containing correlation coefficients from correlationResults
% only_shk_responsive_corrs: array containing correlation coefficients from a different variable

% Calculate means
mean_only_shk = mean(only_shk_responsive_corrs);
mean_not_shk = mean(not_shk_responsive_corrs);

% Create a histogram for allCorrelations
figure;
width = 250; % Width of the figure
height = 250; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
histogram(not_shk_responsive_corrs , 'Normalization', 'probability', 'FaceColor', 'blue','BinWidth', 0.05,'LineStyle','none');
hold on;

% Create a histogram for only_shk_responsive_corrs on the same figure
histogram(only_shk_responsive_corrs, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05, 'LineStyle','none');
xline(mean_only_shk, 'r')
xline(mean_not_shk, 'g')
% Add labels and title
xlabel('Correlation Coefficient');
ylabel('Probability');
% title('Histograms of Correlation Coefficients');
% legend('All Correlations', 'Only SHK Responsive Correlations');

% Optionally, you can add a vertical line at 0 for reference
yLimits = [0 0.15];
plot([0 0], yLimits, 'k', 'LineWidth', 2);
xtickformat('%.2f');
ytickformat('%.2f');
hold off;

% Perform a Kolmogorov-Smirnov test to compare the two distributions
[h, p, k] = kstest2(not_shk_responsive_corrs , only_shk_responsive_corrs)

% Display the results of the statistical test
fprintf('Kolmogorov-Smirnov test result:\n');
fprintf('h = %d (0 means the null hypothesis cannot be rejected, 1 means it can be rejected)\n', h);
fprintf('p-value = %.4f\n', p);

[h,p,ci,stats] = ttest2(not_shk_responsive_corrs , only_shk_responsive_corrs)





%

bar_separation_value = 3;

figure;
width = 250; % Width of the figure
height = 250; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
swarmchart(ones(1, length(only_shk_responsive_corrs)), only_shk_responsive_corrs)
hold on
swarmchart(ones(1, length(not_shk_responsive_corrs))*bar_separation_value, not_shk_responsive_corrs)

% yline(mean(only_shk_responsive_corrs), ones(length(only_shk_responsive_corrs)))
plot([0.5; 1.5], [mean(only_shk_responsive_corrs); mean(only_shk_responsive_corrs)], 'LineWidth',3)
plot([bar_separation_value-.5; bar_separation_value+.5], [mean(not_shk_responsive_corrs); mean(not_shk_responsive_corrs)], 'LineWidth',3)
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off

%% Fig. 2K and Fig. 2L
% plot scatters for individual neurons & behav variables



% % prechoice representative: 
% find(correlationResults{5, 1} < -0.3)
% start_time = -4;% sub-window start time
% end_time = 0; % sub-window end time
% sub_window_idx = ts1 >= start_time & ts1 <= end_time;
% sub_window_activity_session_1 = zall_mouse{5, 1}{1, 44}(:, sub_window_idx);
% r_val_for_representative = correlationResults{5, 1}(1, 44)
% p_val_for_representative = correlationResults_sig{5, 1}(1, 44)
% choice_times_mouse = trial_choice_times_by_mouse{1, 5};
% trial_types = trial_types_by_mouse{1, 5};

% postchoice representative:
find(correlationResults{5, 1} < -0.3)
start_time = 0;% sub-window start time
end_time = 2; % sub-window end time
sub_window_idx = ts1 >= start_time & ts1 <= end_time;
sub_window_activity_session_1 = zall_mouse{5, 1}{1, 65}(:, sub_window_idx);
r_val_for_representative = correlationResults{5, 1}(1, 65)
p_val_for_representative = correlationResults_sig{5, 1}(1, 65)
choice_times_mouse = trial_choice_times_by_mouse{1, 5};
trial_types = trial_types_by_mouse{1, 5};

%consumption representative:
% find(correlationResults{5, 1} > 0.3)
% start_time = 1;% sub-window start time
% end_time = 3; % sub-window end time
% sub_window_idx = ts1 >= start_time & ts1 <= end_time;
% sub_window_activity_session_1 = zall_mouse{5, 3}{1, 1}(:, sub_window_idx);
% choice_times_mouse = consum_times_by_mouse{1, 5};
% trial_types = trial_types_by_mouse{1, 5};



%CREATE SCATTER PLOT BASED ON SPECIFIC EVENTS - ASSUMING THEY ARE IN PAIRS.
%CHECK AND UPDATE START & END TIME DEPENDING ON EVENT OF INTEREST
paired_neurons = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 1;
% start_time = 0;% sub-window start time
% end_time = 2; % sub-window end time

% Extract the corresponding columns from neuron_mean


% sub_window_activity_session_1 = zall_mouse{5, 1}{1, 104}(:, sub_window_idx);
% choice_times_mouse = trial_choice_times_by_mouse{1, 5};
% trial_types = trial_types_by_mouse{1, 5};

% % Assume A and B are your 143x21 arrays
% correlation_coefficients = arrayfun(@(i) corr(sub_window_activity_session_1 (i, :)', sub_window_activity_session_2 (i, :)'), 1:size(sub_window_activity_session_1 , 1));


mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);


x = mean_sub_window_activity_session_1;
y = choice_times_mouse;
size(y)

% x = mean_sub_window_activity_session_1(trial_types == 1.2);
% y = choice_times_mouse(trial_types == 1.2);
% trial_types = trial_types(trial_types == 1.2);


% Define colors based on trial types
colors = repmat([0.5, 0.5, 0.5], length(trial_types), 1); % Default to gray
colors(trial_types == 1.2, :) = repmat([0, 0, 1], sum(trial_types == 1.2), 1); % Blue for trial_types == 1.2
colors(trial_types == 0.3, :) = repmat([1, 0, 0], sum(trial_types == 0.3), 1); % Red for trial_types == 0.3

% Create scatter plot
figure;
set(gcf, 'Position', [100, 100, 200, 200]); % Adjust figure position and size
scatter(x, y, 36, colors, 'filled', 'MarkerEdgeColor', 'k'); % Use 'colors' for MarkerFaceColor

hold on;

% Add separate regression lines for each trial type
% Large reward trials (1.2) - Blue dashed line
x_large = x(trial_types == 1.2);
y_large = y(trial_types == 1.2);
coefficients_large = polyfit(x_large, y_large, 1);
x_fit_large = linspace(min(x_large), max(x_large), 100);
y_fit_large = polyval(coefficients_large, x_fit_large);
plot(x_fit_large, y_fit_large, 'b', 'LineWidth', 2);

% Calculate R-squared for large reward trials
y_pred_large = polyval(coefficients_large, x_large);
ssr_large = sum((y_pred_large - mean(y_large)).^2);
sst_large = sum((y_large - mean(y_large)).^2);
r_squared_large = ssr_large / sst_large;

% Small reward trials (0.3) - Red dashed line
x_small = x(trial_types == 0.3);
y_small = y(trial_types == 0.3);
coefficients_small = polyfit(x_small, y_small, 1);
x_fit_small = linspace(min(x_small), max(x_small), 100);
y_fit_small = polyval(coefficients_small, x_fit_small);
plot(x_fit_small, y_fit_small, 'r', 'LineWidth', 2);

% Calculate R-squared for small reward trials
y_pred_small = polyval(coefficients_small, x_small);
ssr_small = sum((y_pred_small - mean(y_small)).^2);
sst_small = sum((y_small - mean(y_small)).^2);
r_squared_small = ssr_small / sst_small;

% Add R-squared values to the plot
text(min(x) + 0.1, max(y) - 0.1, ['Large R^2 = ' num2str(r_squared_large, '%.3f')], 'FontSize', 10, 'Color', 'b');
text(min(x) + 0.1, max(y) - 0.3, ['Small R^2 = ' num2str(r_squared_small, '%.3f')], 'FontSize', 10, 'Color', 'r');

% Add labels and a legend if needed
xlabel('Mean Sub-window Activity Session 1');
ylabel('Choice Times Mouse');
title('Scatter Plot with Regression Lines by Trial Type');
hold off;


