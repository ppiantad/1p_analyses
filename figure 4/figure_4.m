%%

figure;
hold on
% Create a histogram for allCorrelations

width = 200; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(prechoice_block_1==1, :)), nanmean(neuron_sem_array{1, 8}(prechoice_block_1==1, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(prechoice_block_1==1, :)), nanmean(neuron_sem_array{1, 8}(prechoice_block_1==1==1, :)), 'lineProps', {'color', 'k'});

xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.6 1.0]);
ytickformat('%.1f');
hold off


%%

figure;
hold on
% Create a histogram for allCorrelations

width = 200; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 2}(postchoice_reward_block_1==1, :)), nanmean(neuron_sem_array{1, 2}(postchoice_reward_block_1==1, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 9}(postchoice_reward_block_1==1, :)), nanmean(neuron_sem_array{1, 9}(postchoice_reward_block_1==1, :)), 'lineProps', {'color', 'k'});

xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.6 1.0]);
ytickformat('%.1f');
hold off
%%

figure;
hold on
% Create a histogram for allCorrelations

width = 200; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 3}(collect_block_1==1, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 10}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 10}(collect_block_1==1, :)), 'lineProps', {'color', 'k'});

xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.6 1.0]);
ytickformat('%.1f');
hold off




%%
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
        selected_data_block_1 = neuron_mean_array{1, 1}(i, :);
        selected_data_block_2_3 = neuron_mean_array{1, 8}(i, :);
        % % Extract time variable (assuming it's named 'ts1')
        % ts1_data = ts1{i}(exclusive_activated_session_1{i} == 1);

        % Find indices corresponding to each time window
        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        % Compute AUC for each time window
        % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
        action_auc_pre_choice_block_1(pre_choice_neuron_count) = trapz(selected_data_block_1(pre_choice_indices));
        action_auc_post_choice_block_1(pre_choice_neuron_count) = trapz(selected_data_block_1(post_choice_indices));
        action_auc_consumption_block_1(pre_choice_neuron_count) = trapz(selected_data_block_1(consumption_indices));

        action_auc_pre_choice_block_2_3(pre_choice_neuron_count) = trapz(selected_data_block_2_3(pre_choice_indices));
        action_auc_post_choice_block_2_3(pre_choice_neuron_count) = trapz(selected_data_block_2_3(post_choice_indices));
        action_auc_consumption_block_2_3(pre_choice_neuron_count) = trapz(selected_data_block_2_3(consumption_indices));
    else

    end
end

post_choice_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,1}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if postchoice_reward_block_1(i) == 1
        post_choice_neuron_count = post_choice_neuron_count+1;
        selected_data_block_1 = neuron_mean_array{1, 2}(i, :);
        selected_data_block_2_3 = neuron_mean_array{1, 9}(i, :);
        % % Extract time variable (assuming it's named 'ts1')
        % ts1_data = ts1{i}(exclusive_activated_session_1{i} == 1);

        % Find indices corresponding to each time window
        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        % Compute AUC for each time window
        % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice

        post_choice_reward_auc_pre_choice_block_1(post_choice_neuron_count) = trapz(selected_data_block_1(pre_choice_indices));
        post_choice_reward_auc_post_choice_block_1(post_choice_neuron_count) = trapz(selected_data_block_1(post_choice_indices));
        post_choice_reward_auc_consumption_block_1(post_choice_neuron_count) = trapz(selected_data_block_1(consumption_indices));

        post_choice_reward_auc_pre_choice_block_2_3(post_choice_neuron_count) = trapz(selected_data_block_2_3(pre_choice_indices));
        post_choice_reward_auc_post_choice_block_2_3(post_choice_neuron_count) = trapz(selected_data_block_2_3(post_choice_indices));
        post_choice_reward_auc_consumption_block_2_3(post_choice_neuron_count) = trapz(selected_data_block_2_3(consumption_indices));
    else

    end
end


consumption_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,3}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if collect_block_1(i) == 1
        consumption_neuron_count = consumption_neuron_count+1;

        selected_data_block_1 = neuron_mean_array{1, 3}(i, :);
        selected_data_block_2_3 = neuron_mean_array{1, 10}(i, :);
        % % Extract time variable (assuming it's named 'ts1')
        % ts1_data = ts1{i}(exclusive_activated_session_1{i} == 1);

        % Find indices corresponding to each time window
        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        % Compute AUC for each time window
        % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
        % consumption_auc_pre_choice(consumption_neuron_count) = trapz(selected_data(pre_choice_indices));
        % consumption_auc_post_choice(consumption_neuron_count) = trapz(selected_data(post_choice_indices));
        % consumption_auc_consumption(consumption_neuron_count) = trapz(selected_data(consumption_indices));

        consumption_auc_pre_choice_block_1(consumption_neuron_count) = trapz(selected_data_block_1(pre_choice_indices));
        consumption_auc_post_choice_block_1(consumption_neuron_count) = trapz(selected_data_block_1(post_choice_indices));
        consumption_auc_consumption_block_1(consumption_neuron_count) = trapz(selected_data_block_1(consumption_indices));

        consumption_auc_pre_choice_block_2_3(consumption_neuron_count) = trapz(selected_data_block_2_3(pre_choice_indices));
        consumption_auc_post_choice_block_2_3(consumption_neuron_count) = trapz(selected_data_block_2_3(post_choice_indices));
        consumption_auc_consumption_block_2_3(consumption_neuron_count) = trapz(selected_data_block_2_3(consumption_indices));
    else

    end
end

% calculate change in AUC and plot


% Calculate mean and SEM for pre-choice period
mean_pre_choice_block_1 = [mean(action_auc_pre_choice_block_1(:)), mean(post_choice_reward_auc_pre_choice_block_1(:)), mean(consumption_auc_pre_choice_block_1(:))];
sem_pre_choice_block_1 = [std(action_auc_pre_choice_block_1(:))/sqrt(numel(action_auc_pre_choice_block_1)), std(post_choice_reward_auc_pre_choice_block_1(:))/sqrt(numel(post_choice_reward_auc_pre_choice_block_1)), std(consumption_auc_pre_choice_block_1(:))/sqrt(numel(consumption_auc_pre_choice_block_1))];

mean_pre_choice_block_2_3 = [mean(action_auc_pre_choice_block_2_3(:)), mean(post_choice_reward_auc_pre_choice_block_2_3(:)), mean(consumption_auc_pre_choice_block_2_3(:))];
sem_pre_choice_block_2_3 = [std(action_auc_pre_choice_block_2_3(:))/sqrt(numel(action_auc_pre_choice_block_2_3)), std(post_choice_reward_auc_pre_choice_block_2_3(:))/sqrt(numel(post_choice_reward_auc_pre_choice_block_2_3)), std(consumption_auc_pre_choice_block_2_3(:))/sqrt(numel(consumption_auc_pre_choice_block_2_3))];

% Calculate mean and SEM for post-choice period
mean_post_choice_block_1 = [mean(action_auc_post_choice_block_1(:)), mean(post_choice_reward_auc_post_choice_block_1(:)), mean(consumption_auc_post_choice_block_1(:))];
sem_post_choice_block_1 = [std(action_auc_post_choice_block_1(:))/sqrt(numel(action_auc_post_choice_block_1)), std(post_choice_reward_auc_post_choice_block_1(:))/sqrt(numel(post_choice_reward_auc_post_choice_block_1)), std(consumption_auc_post_choice_block_1(:))/sqrt(numel(consumption_auc_post_choice_block_1))];

mean_post_choice_block_2_3 = [mean(action_auc_post_choice_block_2_3(:)), mean(post_choice_reward_auc_post_choice_block_2_3(:)), mean(consumption_auc_post_choice_block_2_3(:))];
sem_post_choice_block_2_3 = [std(action_auc_post_choice_block_2_3(:))/sqrt(numel(action_auc_post_choice_block_2_3)), std(post_choice_reward_auc_post_choice_block_2_3(:))/sqrt(numel(post_choice_reward_auc_post_choice_block_2_3)), std(consumption_auc_post_choice_block_2_3(:))/sqrt(numel(consumption_auc_post_choice_block_2_3))];


% Calculate mean and SEM for consumption period
mean_consumption_block_1 = [mean(action_auc_consumption_block_1(:)), mean(post_choice_reward_auc_consumption_block_1(:)), mean(consumption_auc_consumption_block_1(:))];
sem_consumption_block_1 = [std(action_auc_consumption_block_1(:))/sqrt(numel(action_auc_consumption_block_1)), std(post_choice_reward_auc_consumption_block_1(:))/sqrt(numel(post_choice_reward_auc_consumption_block_1)), std(consumption_auc_consumption_block_1(:))/sqrt(numel(consumption_auc_consumption_block_1))];

mean_consumption_block_2_3 = [mean(action_auc_consumption_block_2_3(:)), mean(post_choice_reward_auc_consumption_block_2_3(:)), mean(consumption_auc_consumption_block_2_3(:))];
sem_consumption_block_2_3 = [std(action_auc_consumption_block_2_3(:))/sqrt(numel(action_auc_consumption_block_2_3)), std(post_choice_reward_auc_consumption_block_2_3(:))/sqrt(numel(post_choice_reward_auc_consumption_block_2_3)), std(consumption_auc_consumption_block_2_3(:))/sqrt(numel(consumption_auc_consumption_block_2_3))];


mean_change_from_block_1_pre_choice = [sem_pre_choice_block_2_3 - sem_pre_choice_block_1];
mean_change_from_block_1_post_choice = [sem_post_choice_block_2_3 - sem_post_choice_block_1];
mean_change_from_block_1_consumption = [sem_consumption_block_2_3 - sem_consumption_block_1];

sem_change_from_block_1_pre_choice = [sem_pre_choice_block_2_3 - sem_pre_choice_block_1];
sem_change_from_block_1_post_choice = [sem_post_choice_block_2_3 - sem_post_choice_block_1];
sem_change_from_block_1_consumption = [sem_consumption_block_2_3 - sem_consumption_block_1];

% Plot the bar graph
figure;
width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
bar_groups = 1:3; % Number of groups
bar_width = 0.3; % Width of each bar (adjust as needed)
hold on;
bar(bar_groups - bar_width, mean_change_from_block_1_pre_choice, bar_width, 'b');
bar(bar_groups, mean_change_from_block_1_post_choice, bar_width, 'g');
bar(bar_groups + bar_width, mean_change_from_block_1_consumption, bar_width, 'r');

% Add error bars
errorbar(bar_groups - bar_width, mean_change_from_block_1_pre_choice, sem_change_from_block_1_pre_choice, 'k', 'LineStyle', 'none');
errorbar(bar_groups, mean_change_from_block_1_post_choice, sem_change_from_block_1_post_choice, 'k', 'LineStyle', 'none');
errorbar(bar_groups + bar_width, mean_change_from_block_1_consumption, sem_change_from_block_1_consumption, 'k', 'LineStyle', 'none');

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



