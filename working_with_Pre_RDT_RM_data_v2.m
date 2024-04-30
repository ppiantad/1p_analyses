


excited_to_excited_all = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} == 1;
excited_to_excited_sum = sum(excited_to_excited_all);

excited_to_excited_1_to_2 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} ~= 1;
excited_to_excited_1_to_2_sum = sum(excited_to_excited_1_to_2);


excited_to_excited_1_to_3 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} == 1;
excited_to_excited_1_to_3_sum = sum(excited_to_excited_1_to_3);

excited_to_excited_2_to_3 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} == 1;
excited_to_excited_2_to_3_sum = sum(excited_to_excited_2_to_3);


co_excited_three_events = excited_to_excited_1_to_2_sum+excited_to_excited_1_to_3_sum+excited_to_excited_2_to_3_sum;


exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} ~= 1;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1);
exclusive_activated_session_2 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} ~= 1;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2);
exclusive_activated_session_3 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} == 1;
exclusive_activated_session_3_sum = sum(exclusive_activated_session_3);

% more conservative approach below
% exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 3;
% exclusive_activated_session_1_sum = sum(exclusive_activated_session_1);
% exclusive_activated_session_2 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} == 3;
% exclusive_activated_session_2_sum = sum(exclusive_activated_session_2);
% exclusive_activated_session_3 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 1;
% exclusive_activated_session_3_sum = sum(exclusive_activated_session_3);


not_active = neuron_num - (exclusive_activated_session_1_sum + exclusive_activated_session_2_sum + exclusive_activated_session_3_sum);

neutral = respClass_all_array{1,1} ~=exclusive_activated_session_1...
    & respClass_all_array{1,1} ~=exclusive_activated_session_2


test_array = zeros(1, size(exclusive_activated_session_1, 2))

test_array = respClass_all_array{1,1} ~= exclusive_activated_session_1 & respClass_all_array{1,2} ~= exclusive_activated_session_2;


not_active = neuron_num - (exclusive_activated_session_1_sum + exclusive_activated_session_2_sum + exclusive_activated_session_3_sum + excited_to_excited_sum+ excited_to_excited_1_to_2_sum+excited_to_excited_1_to_3_sum+excited_to_excited_2_to_3_sum);


% get cells on a mouse x mouse basis (good for decoding).
% should build this out
% Initialize respClass_all_array_mouse_true_neutral
respClass_all_array_mouse_true_neutral = cell(size(respClass_all_array_mouse, 1), 1);
respClass_all_array_mouse_pre_choice_active = cell(size(respClass_all_array_mouse, 1), 1);
respClass_all_array_mouse_post_choice_reward = cell(size(respClass_all_array_mouse, 1), 1);
respClass_all_array_mouse_consumption = cell(size(respClass_all_array_mouse, 1), 1);

% Loop through each row
for row = 1:size(respClass_all_array_mouse, 1)
% Initialize the row result
    row_cells = respClass_all_array_mouse(row, :);
    D = vertcat(row_cells{:});
    % Compare each column of D to check if all values are equal to 3
    col_comparison_true_neutral = all(D == 3);
    col_comparison_pre_choice_active = D(1,:) == 1 & D(2,:) ~= 1 & D(3,:) ~= 1;
    col_comparison_post_choice_reward = D(1,:) ~= 1 & D(2,:) == 1 & D(3,:) ~= 1;
    col_comparison_consumption = D(1,:) ~= 1 & D(2,:) ~= 1 & D(3,:) == 1;
    % Assign 1 or 0 based on the comparison result
    respClass_all_array_mouse_true_neutral{row} = col_comparison_true_neutral;
    respClass_all_array_mouse_pre_choice_active{row} = col_comparison_pre_choice_active;
    respClass_all_array_mouse_post_choice_reward{row} = col_comparison_post_choice_reward;
    respClass_all_array_mouse_consumption{row} = col_comparison_consumption;
end



%%
% Example 2: Nested pie chart with custom colors for each wedge

% Initialize data points
inner_pie = [exclusive_activated_session_1_sum/neuron_num,...
            
            exclusive_activated_session_2_sum/neuron_num,...
            
            exclusive_activated_session_3_sum/neuron_num,...
           
            excited_to_excited_sum/neuron_num,...

            excited_to_excited_1_to_2_sum/neuron_num,...

            excited_to_excited_1_to_3_sum/neuron_num,...

            excited_to_excited_2_to_3_sum/neuron_num,...
            
            
            not_active/neuron_num];

figure; donutchart(inner_pie, 'InnerRadius', 0.5)
figure; pie(inner_pie)


%%
% Example 2: Nested pie chart with custom colors for each wedge






% Initialize data points
inner_pie = [exclusive_activated_session_1_sum/neuron_num,...
            
            exclusive_activated_session_2_sum/neuron_num,...
            
            exclusive_activated_session_3_sum/neuron_num,...
           

            co_excited_three_events/neuron_num,...
            
            
            not_active/neuron_num];

figure; donutchart(inner_pie, 'InnerRadius', 0.5)
figure; pie(inner_pie)



%%
% Example 2: Nested pie chart with custom colors for each wedge

clear inner_pie
not_active_alternate = neuron_num - (exclusive_activated_session_1_sum + exclusive_activated_session_2_sum + exclusive_activated_session_3_sum + excited_to_excited_sum);

% Initialize data points
inner_pie = [exclusive_activated_session_1_sum/neuron_num,...
            
            exclusive_activated_session_2_sum/neuron_num,...
            
            exclusive_activated_session_3_sum/neuron_num,...
          
            
            
            not_active_alternate/neuron_num];

figure; donutchart(inner_pie, 'InnerRadius', 0.5)
figure; pie(inner_pie)



%%
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


%%


figure;
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(exclusive_activated_session_1==1, :)), nanmean(neuron_sem_array{1, 1}(exclusive_activated_session_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(exclusive_activated_session_2==1, :)), nanmean(neuron_sem_array{1, 1}(exclusive_activated_session_2==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(exclusive_activated_session_3==1, :)), nanmean(neuron_sem_array{1, 1}(exclusive_activated_session_3==1, :)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')


%%
% Get AUCs for the relevant periods for the 3 defined events
% Define time windows
pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [1 3];     % Consumption period: 1 to 3 s

% Initialize arrays to store AUCs
% auc_pre_choice = zeros(size(neuron_mean_array));
% auc_post_choice = zeros(size(neuron_mean_array));
% auc_consumption = zeros(size(neuron_mean_array));
pre_choice_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,1}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if exclusive_activated_session_1(i) == 1
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
    if exclusive_activated_session_2(i) == 1
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
    if exclusive_activated_session_3(i) == 1
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

%%
%run data_loop with choiceTime large rew [-10 to 5]
large_pre_choice_ensemble_zall = zall_mean_all_array{1, 1}(exclusive_activated_session_1==1, :);
large_pre_choice_ensemble_sem = sem_all_array{1, 1}(exclusive_activated_session_1==1, :);

large_post_choice_ensemble_zall = zall_mean_all_array{1, 1}(exclusive_activated_session_2==1, :);
large_post_choice_ensemble_sem = sem_all_array{1, 1}(exclusive_activated_session_2==1, :);

large_consumption_ensemble_zall = zall_mean_all_array{1, 1}(exclusive_activated_session_3==1, :);
large_consumption_ensemble_sem = sem_all_array{1, 1}(exclusive_activated_session_3==1, :);




figure;
shadedErrorBar(ts1, nanmean(large_pre_choice_ensemble_zall), nanmean(large_pre_choice_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(large_post_choice_ensemble_zall), nanmean(large_post_choice_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(large_consumption_ensemble_zall), nanmean(large_consumption_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'consumption active', 'neutral'}, 'Location','northwest')






%%
%run data_loop with choiceTime large rew [-10 to 5]
small_pre_choice_ensemble_zall = zall_mean_all_array{1, 2}(exclusive_activated_session_1==1, :);
small_pre_choice_ensemble_sem = sem_all_array{1, 2}(exclusive_activated_session_1==1, :);

small_post_choice_ensemble_zall = zall_mean_all_array{1, 2}(exclusive_activated_session_2==1, :);
small_post_choice_ensemble_sem = sem_all_array{1, 2}(exclusive_activated_session_2==1, :);

small_consumption_ensemble_zall = zall_mean_all_array{1, 2}(exclusive_activated_session_3==1, :);
small_consumption_ensemble_sem = sem_all_array{1, 2}(exclusive_activated_session_3==1, :);




figure;
shadedErrorBar(ts1, nanmean(small_pre_choice_ensemble_zall), nanmean(small_pre_choice_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(small_post_choice_ensemble_zall), nanmean(small_post_choice_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(small_consumption_ensemble_zall), nanmean(small_consumption_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'consumption active', 'neutral'}, 'Location','northwest')





%% small vs large pre-choice ensemble
figure;
shadedErrorBar(ts1, nanmean(large_pre_choice_ensemble_zall), nanmean(large_pre_choice_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(small_pre_choice_ensemble_zall), nanmean(small_pre_choice_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});

xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'consumption active', 'neutral'}, 'Location','northwest')


%% small vs large post-choice reward ensemble
figure;
shadedErrorBar(ts1, nanmean(large_post_choice_ensemble_zall), nanmean(large_post_choice_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(small_post_choice_ensemble_zall), nanmean(small_post_choice_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});

xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'consumption active', 'neutral'}, 'Location','northwest')


%% small vs large post-choice reward ensemble
figure;
shadedErrorBar(ts1, nanmean(large_consumption_ensemble_zall), nanmean(large_consumption_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(small_consumption_ensemble_zall), nanmean(small_consumption_ensemble_sem), 'lineProps', {'color', batlowW(iter,:)});

xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'consumption active', 'neutral'}, 'Location','northwest')
