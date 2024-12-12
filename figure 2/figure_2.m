load('BLA_C_raw_no_additional_filtering_RDT_D1_only_completed_sessions_zall_window_base_workspace_10_categories.mat')
%% for heatmap, change "plot_num" (what neuron to plot) and array_to_plot (# corresponds to which dataset)

plot_num = 103 

array_to_plot = 1; % depends on the structure of zall

select_mouse = 'BLA_Insc_40';

% for RDT D1 BLA_Insc_25:
%prechoice neuron num 46
%postchoice rew num 38
%consumption num 39
%shock num 11


% for RDT D1 BLA_Insc_40:
%prechoice neuron num 12
%postchoice rew num 70
%consumption num 10
%shock num 11


select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'Pre_RDT_RM';

second_session = 'RDT_D1';



% in order to trim off excess (because calcium recording starts before
% behavior), you'll need the start time. unfortunately I don't save that
% variable in the "final" struct, but you can get it from the adjustment to
% some of the columns of BehavData. e.g., see below - but make sure to
% update the session etc as necessary! 

% BehavData = final.(select_mouse).(first_session).choiceTime.uv.BehavData;
BehavData = final.(select_mouse).(first_session).uv.BehavData;
% because the first trial possible is ALWAYS 60 seconds after ABET is
% issued, you can determine what adjustment has been made to this column
% (adding time to account for calcium recording starting first) by
% subtracting off 60 from the first element
stTime = BehavData.TrialPossible(1)-60; 



time2Collect = BehavData.collectionTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
trialStartTime = BehavData.stTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
median_trialStartTime = median(trialStartTime)
median_time2Collect = median(time2Collect)
xline(median_trialStartTime)
xline(median_time2Collect)
[numTrials, ~] = size(time2Collect);
Tris = [1:numTrials]';

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

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 200, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
ax1 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, 1:size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1), zall_mouse{select_mouse_index, array_to_plot}{1, plot_num});


% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar(ax1, 'eastoutside');
set(c, 'YTick', clim); % 

ylim([0.5, size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1) + 0.5]);

xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', []);
set(gca, 'YTick', [1, size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1)]);
xline(0)
scatter(time2Collect, Tris               , 'Marker', 'p')
scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;


% Plot the raw data in grey with transparency
for trial = 1:size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1)
    plot(ts1, zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
end

% Plot the mean as a thick black line
meanData = mean(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num});
plot(ts1, meanData, 'r', 'LineWidth', 2, 'Color', 'k');

ylim([-3 4]);
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
xline(0)
yline(0)
fontsize(18, 'points')
hold off;

%%
if size(respClass_all_array, 2) == 10 | size(respClass_all_array, 2) == 11
    comparison_arrays = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays = [1 2 3; 4 5 6]
elseif size(respClass_all_array, 2) == 3
    comparison_arrays = [1 2 3; 1 2 3]
end



arrays_to_examine = [1 8];

inhib_or_excite = 1;

event_for_figures = 1; 



% prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite;
prechoice_blocks_2_and_3 = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;

% postchoice_reward_block_1 = respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
postchoice_reward_block_1 = respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite;
postchoice_reward_blocks_2_and_3 = respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;

% collect_block_1 = respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
collect_block_1 = respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite;
collect_blocks_2_and_3 = respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite;

% shk_neurons = respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} == inhib_or_excite;
shk_neurons = 0;

rest_of_neurons = neuron_num - [sum(prechoice_block_1)+sum(postchoice_reward_block_1)+sum(collect_block_1)+sum(shk_neurons)];
figure; pie([sum(prechoice_block_1), sum(postchoice_reward_block_1), sum(collect_block_1), sum(shk_neurons), rest_of_neurons])

block_1_pre_and_post = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite;
sum(block_1_pre_and_post)
block_1_post_and_consumption = respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite;
sum(block_1_post_and_consumption)
block_1_pre_and_consumption = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite;
sum(block_1_pre_and_consumption)


block_2_and_3_pre_and_post = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;
sum(block_2_and_3_pre_and_post)
block_2_and_3_post_and_consumption = respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite;
sum(block_2_and_3_post_and_consumption)
block_2_and_3_pre_and_consumption = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite;
sum(block_2_and_3_pre_and_consumption)


%% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav
% this outputs a ever so slightly wonky diagram. a few nodes that do not
% actually overlap minimally overlap (but intersections are 0), and 1 node
% that has 1 overlap does not overlap at all. 
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

median_start_time_from_choice_large = median(concatenatedTable.stTime(concatenatedTable.bigSmall == 1.2) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 1.2));
median_start_time_from_choice_small = median(concatenatedTable.stTime(concatenatedTable.bigSmall == 0.3) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 0.3));

median_collect_time_from_choice_large = median(concatenatedTable.collectionTime(concatenatedTable.bigSmall == 1.2) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 1.2));
median_collect_time_from_choice_small = median(concatenatedTable.collectionTime(concatenatedTable.bigSmall == 0.3) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 0.3));


%%


figure;
hold on
% Create a histogram for allCorrelations

width = 600; % Width of the figure
height = 300; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(prechoice_block_1==1, :)), nanmean(neuron_sem_array{1, 1}(prechoice_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(postchoice_reward_block_1==1, :)), nanmean(neuron_sem_array{1, 1}(postchoice_reward_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 1}(collect_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.8 0.8]);
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
width = 600; % Width of the figure
height = 300; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
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






%% create correlation matrix heatmap with exclusively active cells
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

% Initialize matrices to store correlation coefficients and p-values
correlation_matrix = zeros(size(data, 1));
p_value_matrix = zeros(size(data, 1));

% Calculate correlation coefficients and p-values between rows
for i = 1:size(data, 1)
    for j = 1:size(data, 1)
        [corr_coeff, p_value] = corrcoef(data(i, :)', data(j, :)');
        correlation_matrix(i, j) = corr_coeff(1, 2); % Store correlation coefficient
        p_value_matrix(i, j) = p_value(1, 2); % Store p-value
    end
end

% Plot the correlation matrix
figure;
imagesc(correlation_matrix);

axis square; % Make the plot square for better visualization
% title('Correlation Matrix');
% xlabel('Neuron Number');
% ylabel('Neuron Number');
ylim([1  size(test, 1)])
% Show row and column indices on the plot
xticks([1  size(test, 1)]);
yticks([1  size(test, 1)]);

% If you want to customize the color map, you can use colormap function
% For example, using a blue-white-red colormap:
colormap(gray);

% If you want to limit the color scale to the range [0, 1]
caxis([-1 1]); % Assuming correlations range from -1 to 1
% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 
% % Display p-values as text on the plot
% for i = 1:size(data, 1)
%     for j = 1:size(data, 1)
%         text(j, i, sprintf('p = %.4f', p_value_matrix(i, j)), ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     end
% end


%%
action_p_value_matrix = p_value_matrix(pre_choice_index, pre_choice_index);
action_correl_matrix = correlation_matrix(pre_choice_index, pre_choice_index);

n = size(pre_choice_index, 2); % Total number of neurons
k = 2;   % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n, k);
disp(['Number of distinct pairwise combinations: ', num2str(num_combinations)]);


% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
action_positive_count = 0;
action_negative_count = 0;
action_no_correlation_count = 0;

% Get the size of the correlation matrix
matrix_size = size(action_correl_matrix, 1);
uu = 1;
% Loop over the upper triangular part of the correlation matrix
for i = 1:matrix_size
    for j = i+1:matrix_size % Start from i+1 to exclude the diagonal
        % get total pairwise correlations between all possible combos of
        % neurons
        action_ensemble_corr_overall(uu) = action_correl_matrix(i, j);
        uu = uu+1;
        % Check if p-value is less than 0.01
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

% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
action_data = [(action_positive_count/action_comparisons_possible)*100, (action_negative_count/action_comparisons_possible)*100, (action_no_correlation_count/action_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, action_data, 'stacked');

% Add labels and title
xlabel('Counts');
ylabel('Correlation Type');
title('Correlation Counts');

% Add legend
legend(labels);

% Adjust x-axis limits
xlim([0.5, 1.5]); % since we only have one set of data, we set the limits to center the bars

% Adjust y-axis limits if needed
% ylim([0, max(data) + 10]); % adjust ylim if needed for better visualization

% Optionally, you can add data labels on each bar
% text(1:length(data), data, num2str(data'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Optionally, you can rotate x-axis labels if needed
% xticklabels(labels);

% Optionally, you can change bar colors
% colormap([0.8 0.2 0.2; 0.2 0.8 0.2; 0.2 0.2 0.8]); % customize colors as needed


% uu = 1
% for i = 1:matrix_size
%     for j = i+1:matrix_size % Start from i+1 to exclude the diagonal
%         % Check if p-value is less than 0.01
%         action_ensemble_corr_overall(uu) = action_correl_matrix(i, j);
%         uu = uu+1;
%     end
% end


%%
post_choice_p_value_matrix = p_value_matrix(post_choice_index, post_choice_index);
post_choice_correl_matrix = correlation_matrix(post_choice_index, post_choice_index);

n = size(post_choice_index, 2); % Total number of neurons
k = 2;   % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n, k);
disp(['Number of distinct pairwise combinations: ', num2str(num_combinations)]);


% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
post_choice_positive_count = 0;
post_choice_negative_count = 0;
post_choice_no_correlation_count = 0;

% Get the size of the correlation matrix
matrix_size = size(post_choice_correl_matrix, 1);
uu = 1
% Loop over the upper triangular part of the correlation matrix
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
% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
post_choice_data = [(post_choice_positive_count/post_choice_comparisons_possible)*100, (post_choice_negative_count/post_choice_comparisons_possible)*100, (post_choice_no_correlation_count/post_choice_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, post_choice_data, 'stacked');

% Add labels and title
xlabel('Counts');
ylabel('Correlation Type');
title('Correlation Counts');

% Add legend
legend(labels);

% Adjust x-axis limits
xlim([0.5, 1.5]); % since we only have one set of data, we set the limits to center the bars

% Adjust y-axis limits if needed
% ylim([0, max(data) + 10]); % adjust ylim if needed for better visualization

% Optionally, you can add data labels on each bar
% text(1:length(data), data, num2str(data'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Optionally, you can rotate x-axis labels if needed
% xticklabels(labels);

% Optionally, you can change bar colors
% colormap([0.8 0.2 0.2; 0.2 0.8 0.2; 0.2 0.2 0.8]); % customize colors as needed





%%
consumption_p_value_matrix = p_value_matrix(consumption_index, consumption_index);
consumption_correl_matrix = correlation_matrix(consumption_index, consumption_index);

n = size(consumption_index, 2); % Total number of neurons
k = 2;   % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n, k);
disp(['Number of distinct pairwise combinations: ', num2str(num_combinations)]);


% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
consumption_positive_count = 0;
consumption_negative_count = 0;
consumption_no_correlation_count = 0;

% Get the size of the correlation matrix
matrix_size = size(consumption_correl_matrix, 1);
uu = 1
% Loop over the upper triangular part of the correlation matrix
for i = 1:matrix_size
    for j = i+1:matrix_size % Start from i+1 to exclude the diagonal
        consumption_ensemble_corr_overall(uu) = consumption_correl_matrix(i, j);
        uu = uu+1;
        % Check if p-value is less than 0.01
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

% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
consumption_data = [(consumption_positive_count/consumption_comparisons_possible)*100, (consumption_negative_count/consumption_comparisons_possible)*100, (consumption_no_correlation_count/consumption_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, consumption_data, 'stacked');

% Add labels and title
xlabel('Counts');
ylabel('Correlation Type');
title('Correlation Counts');

% Add legend
legend(labels);

% Adjust x-axis limits
xlim([0.5, 1.5]); % since we only have one set of data, we set the limits to center the bars

% Adjust y-axis limits if needed
% ylim([0, max(data) + 10]); % adjust ylim if needed for better visualization

% Optionally, you can add data labels on each bar
% text(1:length(data), data, num2str(data'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Optionally, you can rotate x-axis labels if needed
% xticklabels(labels);

% Optionally, you can change bar colors
% colormap([0.8 0.2 0.2; 0.2 0.8 0.2; 0.2 0.2 0.8]); % customize colors as needed



%%
action_post_choice_p_value_matrix = p_value_matrix(pre_choice_index, post_choice_index);
action_post_choice_correl_matrix = correlation_matrix(pre_choice_index, post_choice_index);

n1 = size(action_post_choice_p_value_matrix, 1); % Number of neurons in the first set
n2 = size(action_post_choice_p_value_matrix, 2); % Number of neurons in the second set
k = 2;    % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n1, k) * nchoosek(n2, k);
disp(['Number of unique pairwise combinations: ', num2str(num_combinations)]);



% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
action_post_choice_positive_count = 0;
action_post_choice_negative_count = 0;
action_post_choice_no_correlation_count = 0;

% % Initialize counters
% positive_count = 0;
% negative_count = 0;
% no_correlation_count = 0;

% Get the size of the correlation matrix
[num_neurons_1, num_neurons_2] = size(action_post_choice_correl_matrix);

% Loop through the matrices to count correlations
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


% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
action_post_choice_data = [(action_post_choice_positive_count/action_post_choice_comparisons_possible)*100, (action_post_choice_negative_count/action_post_choice_comparisons_possible)*100, (action_post_choice_no_correlation_count/action_post_choice_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, action_post_choice_data, 'stacked');

% Add labels and title
xlabel('Counts');
ylabel('Correlation Type');
title('Correlation Counts');

% Add legend
legend(labels);

% Adjust x-axis limits
xlim([0.5, 1.5]); % since we only have one set of data, we set the limits to center the bars

% Adjust y-axis limits if needed
% ylim([0, max(data) + 10]); % adjust ylim if needed for better visualization

% Optionally, you can add data labels on each bar
% text(1:length(data), data, num2str(data'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Optionally, you can rotate x-axis labels if needed
% xticklabels(labels);

% Optionally, you can change bar colors
% colormap([0.8 0.2 0.2; 0.2 0.8 0.2; 0.2 0.2 0.8]); % customize colors as needed
%%
action_consumption_p_value_matrix = p_value_matrix(pre_choice_index, consumption_index);
action_consumption_correl_matrix = correlation_matrix(pre_choice_index, consumption_index);

n1 = size(action_consumption_p_value_matrix, 1); % Number of neurons in the first set
n2 = size(action_consumption_p_value_matrix, 2); % Number of neurons in the second set
k = 2;    % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n1, k) * nchoosek(n2, k);
disp(['Number of unique pairwise combinations: ', num2str(num_combinations)]);



% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
action_consumption_positive_count = 0;
action_consumption_negative_count = 0;
action_consumption_no_correlation_count = 0;

% % Initialize counters
% positive_count = 0;
% negative_count = 0;
% no_correlation_count = 0;

% Get the size of the correlation matrix
[num_neurons_1, num_neurons_2] = size(action_consumption_correl_matrix);

% Loop through the matrices to count correlations
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


% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
action_consumption_data = [(action_consumption_positive_count/action_consumption_comparisons_possible)*100, (action_consumption_negative_count/action_consumption_comparisons_possible)*100, (action_consumption_no_correlation_count/action_consumption_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, action_consumption_data, 'stacked');

% Add labels and title
xlabel('Counts');
ylabel('Correlation Type');
title('Correlation Counts');

% Add legend
legend(labels);

% Adjust x-axis limits
xlim([0.5, 1.5]); % since we only have one set of data, we set the limits to center the bars

% Adjust y-axis limits if needed
% ylim([0, max(data) + 10]); % adjust ylim if needed for better visualization

% Optionally, you can add data labels on each bar
% text(1:length(data), data, num2str(data'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Optionally, you can rotate x-axis labels if needed
% xticklabels(labels);

% Optionally, you can change bar colors
% colormap([0.8 0.2 0.2; 0.2 0.8 0.2; 0.2 0.2 0.8]); % customize colors as needed

%%
post_choice_consumption_p_value_matrix = p_value_matrix(post_choice_index, consumption_index);
post_choice_consumption_correl_matrix = correlation_matrix(post_choice_index, consumption_index);

n1 = size(post_choice_consumption_p_value_matrix, 1); % Number of neurons in the first set
n2 = size(post_choice_consumption_p_value_matrix, 2); % Number of neurons in the second set
k = 2;    % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n1, k) * nchoosek(n2, k);
disp(['Number of unique pairwise combinations: ', num2str(num_combinations)]);



% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
post_choice_consumption_positive_count = 0;
post_choice_consumption_negative_count = 0;
post_choice_consumption_no_correlation_count = 0;

% % Initialize counters
% positive_count = 0;
% negative_count = 0;
% no_correlation_count = 0;

% Get the size of the correlation matrix
[num_neurons_1, num_neurons_2] = size(post_choice_consumption_correl_matrix);

% Loop through the matrices to count correlations
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


% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
post_choice_consumption_data = [(post_choice_consumption_positive_count/post_choice_consumption_comparisons_possible)*100, (post_choice_consumption_negative_count/post_choice_consumption_comparisons_possible)*100, (post_choice_consumption_no_correlation_count/post_choice_consumption_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, post_choice_consumption_data, 'stacked');

% Add labels and title
xlabel('Counts');
ylabel('Correlation Type');
title('Correlation Counts');

% Add legend
legend(labels);

% Adjust x-axis limits
xlim([0.5, 1.5]); % since we only have one set of data, we set the limits to center the bars

% Adjust y-axis limits if needed
% ylim([0, max(data) + 10]); % adjust ylim if needed for better visualization

% Optionally, you can add data labels on each bar
% text(1:length(data), data, num2str(data'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Optionally, you can rotate x-axis labels if needed
% xticklabels(labels);

% Optionally, you can change bar colors
% colormap([0.8 0.2 0.2; 0.2 0.8 0.2; 0.2 0.2 0.8]); % customize colors as needed



%%

% Assuming you have action_data, consumption_data, and action_consumption_data matrices

% Define x-values for each set of bars
x = 1:6;
% x = 1:3;
% Create a new figure
figure;

% Plot the stacked bar plots for each dataset
% bar(x, [action_data; post_choice_data; consumption_data; action_post_choice_data; action_consumption_data; post_choice_consumption_data], 'stacked');

barh(x, [action_data; post_choice_data; consumption_data; action_post_choice_data; action_consumption_data; post_choice_consumption_data], 'stacked');
% bar(x, [action_data; consumption_data; action_consumption_data], 'stacked');

% % Set x-axis tick locations and labels
% xticks(x);
% xticklabels({'aa', 'cc', 'ac'});
% ylabel('% Neuron Pairs');
% title('Title');

% Add legend
legend('Positive correlation', 'Negative correlation', 'No sig correlation');


%% use these data for mouse x mouse, which is likely better
% load a RDT dataset with the following variables filtered (in order):
    % "choiceTime.Outcome_Minus_4to0.OMITALL_0_BLANK_TOUCH_0_BLOCK_1"
    % "choiceTime.Outcome_0to2.OMITALL_0_BLANK_TOUCH_0_BLOCK_1"
    % "collectionTime.Outcome_1to3.OMITALL_0_BLANK_TOUCH_0_BLOCK_1"
    % "choiceTime.Outcome_0to2.SHK_1"
    % "choiceTime.Outcome_0to2.LOSS_PLUS_ONE_1"

%%

array_for_means = 1; 


for q = 1:length (behav_tbl_iter{1, 1})
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
        clear trial_choice_times delay_to_initiation delay_to_collect_post_shk trial_types
    end


end

trial_types_concat = cat(1, trial_types_by_mouse{:});
trial_choice_times_concat = cat(1, trial_choice_times_by_mouse{:});
rew_collect_times_concat = cat(1, delay_to_collect_post_shk_by_mouse{:});

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

hold off;
title('Reward Collection Times');
xlabel('Reward Collection');
ylabel('Time');
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');

% Subplot 3: path_length_concat
subplot(1, 3, 3); % Third subplot
hold on;
% Define colors based on trial_types_concat
colors = repmat([0.5, 0.5, 0.5], length(path_length_concat), 1); % Default to gray
colors(trial_types_concat == 1.2, :) = repmat([0, 0, 1], sum(trial_types_concat == 1.2), 1); % Blue for trial_types_concat == 1.2
colors(trial_types_concat == 0.3, :) = repmat([1, 0, 0], sum(trial_types_concat == 0.3), 1); % Red for trial_types_concat == 0.3

% Create the swarmchart with colors
s = swarmchart(ones(1, length(path_length_concat)) * bar_separation_value * 1.5, path_length_concat, 36, colors, 'filled');
s.MarkerEdgeColor = 'k'; % Black stroke for all circles

% Calculate and plot means for each group
mean_12 = mean(path_length_concat(trial_types_concat == 1.2));
mean_03 = mean(path_length_concat(trial_types_concat == 0.3));
plot([bar_separation_value * 1.5 - 0.2, bar_separation_value * 1.5 + 0.2], [mean_12, mean_12], 'b-', 'LineWidth', 2); % Blue line for trial_types_concat == 1.2
plot([bar_separation_value * 1.5 - 0.2, bar_separation_value * 1.5 + 0.2], [mean_03, mean_03], 'r-', 'LineWidth', 2); % Red line for trial_types_concat == 0.3

hold off;
title('Path Length');
xlabel('Path Length');
ylabel('Length (a.u.)');
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');

variable_to_correlate = trial_choice_times_by_mouse;


%%


% Initialize the new cell array to store the mean values
meanZallMouse = cell(size(zall_mouse, 2), 1);

% Define the time range for 0 to 2 seconds
% timeRange = (ts1 >= -4) & (ts1 <= 0);
timeRange = (ts1 >= 0) & (ts1 <= 2);
% timeRange = (ts1 >= 1) & (ts1 <= 3);


% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_mouse)
    % Get the current nested cell array
    nestedCellArray_1 = zall_mouse{i, array_for_means};
    nestedCellArray_2 = zall_mouse{i, 1};
    % Initialize the nested cell array for storing mean values
    meanNestedCellArray = cell(size(nestedCellArray_1));
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(nestedCellArray_1)
        % Get the current double array
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
        
        % Store the mean values in the corresponding cell of the nested cell array
        meanNestedCellArray{j} = meanValues;
    end
    
    % Store the nested cell array of mean values in the corresponding cell of the main cell array
    meanZallMouse{i} = meanNestedCellArray;
end

% Now, meanZallMouse contains the mean activity for each row in the time period 0 to 2 seconds
% Each cell in meanZallMouse contains a nested cell array with the



%% work in progress to do some correlations w/ block 1 activity.

%%

% Assuming the following variables are defined:
% meanZallMouse: 14x1 cell array where each cell contains another cell array with mean values
% trial_choice_times_by_mouse: 1x11 cell array containing values to correlate with

% Initialize the new cell array to store the correlation results
correlationResults = cell(size(meanZallMouse));



% Iterate through each level of meanZallMouse
for i = 1:length(meanZallMouse)
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{i};
    
    % Initialize the nested cell array for storing correlation results
    correlationNestedArray = zeros(size(meanNestedCellArray));
    
    % Determine the corresponding index in trial_choice_times_by_mouse
    % Adjust this logic based on how the indices are mapped
    trialIndex = mod(i-1, length(variable_to_correlate)) + 1;
    
    % Get the corresponding trial choice times
    trialChoiceTimes = variable_to_correlate{i}';
    % trialChoiceTimes = variable_to_correlate{i};

    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        % Get the current mean values array
        meanValues = meanNestedCellArray{j};
        
        % Check if trialChoiceTimes has the same length as meanValues
        if length(trialChoiceTimes) == length(meanValues)
            % Compute the correlation
            correlationCoeff = corr(meanValues, trialChoiceTimes(:));
        elseif length(trialChoiceTimes) < length(meanValues)
            correlationCoeff = corr(meanValues(1:end-1), trialChoiceTimes(:));
        else
            % If lengths do not match, handle the mismatch (e.g., set correlation to NaN)
            correlationCoeff = NaN;
        end
        
        % Store the correlation coefficient in the nested cell array
        correlationNestedArray(j) = correlationCoeff;
    end
    clear meanValues
    % Store the nested cell array of correlation coefficients in the main cell array
    correlationResults{i} = correlationNestedArray;
end

% Now, correlationResults contains the correlation coefficients for each nested structure in meanZallMouse
%%
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

%% SHK responsive neurons assumed to be stored in respClass_all_array{1, 1} for this purpose - change as necessary
only_shk_responsive_corrs = allCorrelations(postchoice_reward_block_1 == 1);
not_shk_responsive_corrs = allCorrelations(postchoice_reward_block_1 ~=1);
% not_shk_responsive_corrs = allCorrelations(true_neutral ==1);

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

%%
% Assuming the following variables are defined:
% allCorrelations: array containing correlation coefficients from correlationResults
% only_shk_responsive_corrs: array containing correlation coefficients from a different variable

% Calculate means
mean_only_shk = mean(only_shk_responsive_corrs);
mean_not_shk = mean(not_shk_responsive_corrs);

% Create a histogram for allCorrelations
figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
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
yLimits = ylim;
plot([0 0], yLimits, 'k', 'LineWidth', 2);
xtickformat('%.2f');
ytickformat('%.2f');
hold off;

% Perform a Kolmogorov-Smirnov test to compare the two distributions
[h, p] = kstest2(not_shk_responsive_corrs , only_shk_responsive_corrs);

% Display the results of the statistical test
fprintf('Kolmogorov-Smirnov test result:\n');
fprintf('h = %d (0 means the null hypothesis cannot be rejected, 1 means it can be rejected)\n', h);
fprintf('p-value = %.4f\n', p);

[h,p,ci,stats] = ttest2(not_shk_responsive_corrs , only_shk_responsive_corrs)


%%
for hh = 1:length(meanZallMouse)
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{hh};
    only_shk_meanNestedCellArray = meanNestedCellArray(:, respClass_all_array_mouse{hh, 4}==1);

    only_shk_meanNestedCellArray_mat = cell2mat(only_shk_meanNestedCellArray);
    mean_only_shk_meanNestedCellArray_mat = mean(only_shk_meanNestedCellArray_mat, 2);
    mean_mean_only_shk_meanNestedCellArray_mat(hh) = mean(mean_only_shk_meanNestedCellArray_mat);

end

scatter(mean_mean_only_shk_meanNestedCellArray_mat, riskiness)
corr(mean_mean_only_shk_meanNestedCellArray_mat', riskiness)


%%





% % Plot means as bars
% figure;
% hold on;
% bar(1, mean_only_shk, 'FaceColor', 'r'); % Red bar for 'only_shk_responsive_corrs'
% bar(2, mean_not_shk, 'FaceColor', 'b'); % Blue bar for 'not_shk_responsive_corrs'
% 
% swarmchart(1, only_shk_responsive_corrs, 5)
% 
% % Scatter individual data points
% scatter(ones(size(only_shk_responsive_corrs)), only_shk_responsive_corrs, 'r', 'filled', 'jitter', 'on', 'jitterAmount', 0.15); % Red points with jitter
% scatter(2 * ones(size(not_shk_responsive_corrs)), not_shk_responsive_corrs, 'b', 'filled', 'jitter', 'on', 'jitterAmount', 0.15); % Blue points with jitter
% 
% % Customize plot
% xlim([0.5, 2.5]);
% xticks([1 2]);
% xticklabels({'Only Shk Responsive', 'Not Shk Responsive'});
% ylabel('Correlation Values');
% title('Correlation Values and Means');
% grid on;
% 
% % Add legend
% legend({'Mean Only Shk', 'Mean Not Shk', 'Individual Only Shk', 'Individual Not Shk'}, 'Location', 'best');
% 
% hold off;

bar_separation_value = 3;

figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
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




%% minor attempts to correlate with MOTION (velocity after shock)

meanNestedCellArray = meanZallMouse{6};

for j = 1:length(meanNestedCellArray)
    % Get the current mean values array
    meanValues = meanNestedCellArray{j};

    % Check if trialChoiceTimes has the same length as meanValues
    if length(meanValues_motion) == length(meanValues)
        % Compute the correlation
        correlationCoeff = corr(meanValues, meanValues_motion(:));
    else
        % If lengths do not match, handle the mismatch (e.g., set correlation to NaN)
        correlationCoeff = NaN;
    end

    % Store the correlation coefficient in the nested cell array
    correlationNestedArray(j) = correlationCoeff;
end




%% USE DATA FROM Pre_RDT_RM (load data with 10 variables), then run data_loop on REW, 1.2 and REW, 0.3


figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
shadedErrorBar(ts1, nanmean(zall_mean_all_array{1, 11}(collect_block_1, :)), nanmean(sem_all_array{1, 11} (collect_block_1==1, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(zall_mean_all_array{1, 12}(collect_block_1==1, :)), nanmean(sem_all_array{1, 12}(collect_block_1==1, :)), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 1}(collect_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice_large, 'g', {'Median', 'start', 'time'})
xline(median_start_time_from_choice_small, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_large, 'r', {'Median', 'collect', 'latency'})
xline(median_collect_time_from_choice_small, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.8 0.8]);
hold off

%%

mean_data_array = {zall_mean_all_array{1, 11}(prechoice_block_1==1, :), zall_mean_all_array{1, 12}(prechoice_block_1==1, :)};
sem_data_array = {sem_all_array{1, 11}(prechoice_block_1==1, :), sem_all_array{1, 12}(prechoice_block_1==1, :)};

[comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1)

%%

%CREATE SCATTER PLOT BASED ON SPECIFIC EVENTS - ASSUMING THEY ARE IN PAIRS.
%CHECK AND UPDATE START & END TIME DEPENDING ON EVENT OF INTEREST
paired_neurons = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 1;
start_time = 1;% sub-window start time
end_time = 3; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Extract the corresponding columns from neuron_mean


sub_window_activity_session_1 = zall_mean_all_array{1, 11}(collect_block_1, sub_window_idx);
sub_window_activity_session_2 = zall_mean_all_array{1, 12}(collect_block_1, sub_window_idx);

% Assume A and B are your 143x21 arrays
correlation_coefficients = arrayfun(@(i) corr(sub_window_activity_session_1 (i, :)', sub_window_activity_session_2 (i, :)'), 1:size(sub_window_activity_session_1 , 1));


mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);
mean_sub_window_activity_session_2 = mean(sub_window_activity_session_2, 2);

x = mean_sub_window_activity_session_1;
y = mean_sub_window_activity_session_2;


% Create a scatter plot
figure;
set(gcf,'Position',[100 100 200 200])
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


%% Plot https://stackoverflow.com/questions/54528239/boxplot-for-paired-observations

combined_data = [mean_sub_window_activity_session_1 , mean_sub_window_activity_session_2];
figure();
set(gcf,'Position',[100 100 200 200])
coordLineStyle = 'k.';
boxplot(combined_data, 'Symbol', coordLineStyle); hold on;
parallelcoords(combined_data, 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
  'Marker', '.', 'MarkerSize', 10);

TF = isoutlier(combined_data, 'grubbs');

[h, p] = ttest2(mean_sub_window_activity_session_1, mean_sub_window_activity_session_2)


%%
ca_data_corresponding_to_large = zall_mouse(:, 11);

behav_data_corresponding_to_large = behav_tbl_iter{11, 1};
neuron_num = 1; 
slow_trial_means = []; 
fast_trial_means = []; 
for hh = 1:size(behav_data_corresponding_to_large, 1)
    current_behav = behav_data_corresponding_to_large{hh};
    current_ca = ca_data_corresponding_to_large{hh}; 
    current_choice_times = current_behav.choiceTime - current_behav.stTime;
    fast_trials = current_choice_times <= 2; 
    slow_trials = current_choice_times >= 4; 
    for jj = 1:size(current_ca, 2)
        
        current_current_ca = current_ca{1, jj}; 
        slow_trial_means(neuron_num, :) = mean(current_current_ca(fast_trials, :));
        fast_trial_means(neuron_num, :) = mean(current_current_ca(slow_trials, :));
        neuron_num = neuron_num + 1; 

    end





end

figure; plot(ts1, mean(slow_trial_means(prechoice_block_1, :)))
hold on; plot(ts1, mean(fast_trial_means(prechoice_block_1, :)))


%% plot scatters for individual neurons & behav variables

%CREATE SCATTER PLOT BASED ON SPECIFIC EVENTS - ASSUMING THEY ARE IN PAIRS.
%CHECK AND UPDATE START & END TIME DEPENDING ON EVENT OF INTEREST
paired_neurons = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 1;
start_time = -4;% sub-window start time
end_time = 0; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Extract the corresponding columns from neuron_mean


sub_window_activity_session_1 = zall_mouse{10, 1}{1, 12}(:, sub_window_idx);
choice_times_mouse = trial_choice_times_by_mouse{1, 10};
trial_types = trial_types_by_mouse{1, 10};

% % Assume A and B are your 143x21 arrays
% correlation_coefficients = arrayfun(@(i) corr(sub_window_activity_session_1 (i, :)', sub_window_activity_session_2 (i, :)'), 1:size(sub_window_activity_session_1 , 1));


mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);


x = mean_sub_window_activity_session_1;
y = choice_times_mouse;





% Define colors based on trial types
colors = repmat([0.5, 0.5, 0.5], length(trial_types), 1); % Default to gray
colors(trial_types == 1.2, :) = repmat([0, 0, 1], sum(trial_types == 1.2), 1); % Blue for trial_types == 1.2
colors(trial_types == 0.3, :) = repmat([1, 0, 0], sum(trial_types == 0.3), 1); % Red for trial_types == 0.3

% Create scatter plot
figure;
set(gcf, 'Position', [100, 100, 200, 200]); % Adjust figure position and size
scatter(x, y, 36, colors, 'filled', 'MarkerEdgeColor', 'k'); % Use 'colors' for MarkerFaceColor

hold on;

% Add a regression line
coefficients = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);
plot(x_fit, y_fit, 'r');

% Calculate R-squared value
y_pred = polyval(coefficients, x);
ssr = sum((y_pred - mean(y)).^2);
sst = sum((y - mean(y)).^2);
r_squared = ssr / sst;

% Add R-squared value to the plot
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared)], 'FontSize', 12);

% Add labels and a legend if needed
xlabel('Mean Sub-window Activity Session 1');
ylabel('Choice Times Mouse');
title('Scatter Plot with Regression Line and R^2 Value');
hold off;


%% these arrays are just the combined_data = [mean_sub_window_activity_session_1 , mean_sub_window_activity_session_2] arrays from above, but saved for Early & Late
load('Late_RM_arrays.mat')
load('RM_D1_arrays.mat')
small_means_prechoice = [mean(RM_D1_prechoice_combined_data(:, 2))  mean(Late_RM_prechoice_combined_data(:, 2))]
large_means_prechoice = [mean(RM_D1_prechoice_combined_data(:, 1))  mean(Late_RM_prechoice_combined_data(:, 1))]


% combined_data = [mean_sub_window_activity_session_1 , mean_sub_window_activity_session_2];
figure();
set(gcf,'Position',[100 100 200 600])
coordLineStyle = 'k.';
boxplot(Late_RM_collect_combined_data, 'Symbol', coordLineStyle); hold on;
parallelcoords(Late_RM_collect_combined_data, 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
  'Marker', '.', 'MarkerSize', 10);
ylim([-0.5 1.5]);
yticks([-0.5 0.0 0.5 1.0 1.5])
ytickformat('%.2f');
TF = isoutlier(Late_RM_collect_combined_data, 'grubbs');

[h, p] = ttest2(Late_RM_collect_combined_data(:, 1), Late_RM_collect_combined_data(:, 2))



%% to compare large vs small early, load 'BLA_RM_D1_small_vs_large_no_yoking_6_categories.mat'

consum_large_and_small_excited = respClass_all_array{1, 3} == 1 & respClass_all_array{1, 6} == 1;
sum(consum_large_and_small_excited)

consum_large_only = respClass_all_array{1, 3} == 1 & respClass_all_array{1, 6} == 3;
consum_small_only = respClass_all_array{1, 3} == 3 & respClass_all_array{1, 6} == 1;
sum(consum_large_only)
sum(consum_small_only)

consum_large_excited_small_inhibited = respClass_all_array{1, 3} == 1 & respClass_all_array{1, 6} == 2;
consum_small_excited_large_inhibited = respClass_all_array{1, 3} == 2 & respClass_all_array{1, 6} == 1;
sum(consum_large_excited_small_inhibited)
sum(consum_small_excited_large_inhibited)


figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(consum_large_and_small_excited, :)), nanmean(neuron_sem_array{1, 3}(consum_large_and_small_excited, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 6}(consum_large_and_small_excited, :)), nanmean(neuron_sem_array{1, 6}(consum_large_and_small_excited, :)), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 1}(collect_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice_large, 'g', {'Median', 'start', 'time'})
xline(median_start_time_from_choice_small, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_large, 'r', {'Median', 'collect', 'latency'})
xline(median_collect_time_from_choice_small, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.8 0.8]);
hold off


figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(consum_large_only, :)), nanmean(neuron_sem_array{1, 3}(consum_large_only, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 6}(consum_large_only, :)), nanmean(neuron_sem_array{1, 6}(consum_large_only, :)), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 1}(collect_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice_large, 'g', {'Median', 'start', 'time'})
xline(median_start_time_from_choice_small, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_large, 'r', {'Median', 'collect', 'latency'})
xline(median_collect_time_from_choice_small, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.8 0.8]);
hold off

figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(consum_small_only, :)), nanmean(neuron_sem_array{1, 3}(consum_small_only, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 6}(consum_small_only, :)), nanmean(neuron_sem_array{1, 6}(consum_small_only, :)), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 1}(collect_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice_large, 'g', {'Median', 'start', 'time'})
xline(median_start_time_from_choice_small, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_large, 'r', {'Median', 'collect', 'latency'})
xline(median_collect_time_from_choice_small, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.8 0.8]);
hold off


%% for heatmap, change "plot_num" (what neuron to plot) and array_to_plot (# corresponds to which dataset)

plot_num = 71

array_to_plot = 6; % depends on the structure of zall

select_mouse = 'BLA_Insc_40'; %34 looks promising for early RM small/large

% for RDT D1 BLA_Insc_25:
%prechoice neuron num 46
%postchoice rew num 38
%consumption num 39
%shock num 11


% for RDT D1 BLA_Insc_40:
%prechoice neuron num 12
%postchoice rew num 70
%consumption num 10
%shock num 11


select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RM_D1';

second_session = 'RDT_D1';



% in order to trim off excess (because calcium recording starts before
% behavior), you'll need the start time. unfortunately I don't save that
% variable in the "final" struct, but you can get it from the adjustment to
% some of the columns of BehavData. e.g., see below - but make sure to
% update the session etc as necessary! 

% BehavData = final.(select_mouse).(first_session).choiceTime.uv.BehavData;
BehavData = final.(select_mouse).(first_session).uv.BehavData;
% because the first trial possible is ALWAYS 60 seconds after ABET is
% issued, you can determine what adjustment has been made to this column
% (adding time to account for calcium recording starting first) by
% subtracting off 60 from the first element
stTime = BehavData.TrialPossible(1)-60; 



time2Collect = BehavData.collectionTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
trialStartTime = BehavData.stTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
median_trialStartTime = median(trialStartTime)
median_time2Collect = median(time2Collect)
xline(median_trialStartTime)
xline(median_time2Collect)
[numTrials, ~] = size(time2Collect);
Tris = [1:numTrials]';

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

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 200, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
ax1 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, 1:size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1), zall_mouse{select_mouse_index, array_to_plot}{1, plot_num});


% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar(ax1, 'eastoutside');
set(c, 'YTick', clim); % 

ylim([0.5, size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1) + 0.5]);

xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', []);
set(gca, 'YTick', [1, size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1)]);
xline(0)
scatter(time2Collect, Tris               , 'Marker', 'p')
scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;


% Plot the raw data in grey with transparency
for trial = 1:size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1)
    plot(ts1, zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
end

% Plot the mean as a thick black line
meanData = mean(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num});
plot(ts1, meanData, 'r', 'LineWidth', 2, 'Color', 'k');

ylim([-3 4]);
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
xline(0)
yline(0)
fontsize(18, 'points')
hold off;


%% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav
% this outputs a ever so slightly wonky diagram. a few nodes that do not
% actually overlap minimally overlap (but intersections are 0), and 1 node
% that has 1 overlap does not overlap at all. 
shk_ind = find(respClass_all_array{1,4} == 1);
pre_choice_active_ind = find(remapped_prechoice == 1);
% pre_choice_active_ind = find(prechoice_block_1 == 1);
pre_choice_active_block_2_3_ind = find(conserved_prechoice== 1);

% consum_inhibited_ind = find(all_consum_inhibited == 1);
setListData = {shk_ind, pre_choice_active_ind, pre_choice_active_block_2_3_ind};
setLabels = ["shk ind", "remapped prechoice", "conserved prechoice"];
figure;
ve_diagram = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

ve_diagram.ShowIntersectionCounts = true;
ve_diagram.ShowIntersectionAreas = true;
% h.SetLabels = [];


%% ONLY USE CODE BELOW IF YOU WANT TO COMBINE DATA ACROSS 2 FILTERS, E.G., BLOCK 1 & BLOCKS 2/3

array_for_means = 1; 
array_for_means_second = 5; 


for q = 1:length (behav_tbl_iter{1, 1})

    behav_part_1 = behav_tbl_iter{array_for_means, 1}{q, 1}
    behav_part_2 = behav_tbl_iter{array_for_means_second, 1}{q, 1}
    nestedCellArray_1 = vertcat(behav_part_1, behav_part_2 )



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
        clear trial_choice_times delay_to_initiation delay_to_collect_post_shk trial_types
    end


end

trial_types_concat = cat(1, trial_types_by_mouse{:});
trial_choice_times_concat = cat(1, trial_choice_times_by_mouse{:});
rew_collect_times_concat = cat(1, delay_to_collect_post_shk_by_mouse{:});

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

hold off;
title('Reward Collection Times');
xlabel('Reward Collection');
ylabel('Time');
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');

% Subplot 3: path_length_concat
subplot(1, 3, 3); % Third subplot
hold on;
% Define colors based on trial_types_concat
colors = repmat([0.5, 0.5, 0.5], length(path_length_concat), 1); % Default to gray
colors(trial_types_concat == 1.2, :) = repmat([0, 0, 1], sum(trial_types_concat == 1.2), 1); % Blue for trial_types_concat == 1.2
colors(trial_types_concat == 0.3, :) = repmat([1, 0, 0], sum(trial_types_concat == 0.3), 1); % Red for trial_types_concat == 0.3

% Create the swarmchart with colors
s = swarmchart(ones(1, length(path_length_concat)) * bar_separation_value * 1.5, path_length_concat, 36, colors, 'filled');
s.MarkerEdgeColor = 'k'; % Black stroke for all circles

% Calculate and plot means for each group
mean_12 = mean(path_length_concat(trial_types_concat == 1.2));
mean_03 = mean(path_length_concat(trial_types_concat == 0.3));
plot([bar_separation_value * 1.5 - 0.2, bar_separation_value * 1.5 + 0.2], [mean_12, mean_12], 'b-', 'LineWidth', 2); % Blue line for trial_types_concat == 1.2
plot([bar_separation_value * 1.5 - 0.2, bar_separation_value * 1.5 + 0.2], [mean_03, mean_03], 'r-', 'LineWidth', 2); % Red line for trial_types_concat == 0.3

hold off;
title('Path Length');
xlabel('Path Length');
ylabel('Length (a.u.)');
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');

variable_to_correlate = trial_choice_times_by_mouse;


%%


% Initialize the new cell array to store the mean values
meanZallMouse = cell(size(zall_mouse, 2), 1);

% Define the time range for 0 to 2 seconds
timeRange = (ts1 >= -4) & (ts1 <= 0);
% timeRange = (ts1 >= 0) & (ts1 <= 2);
% timeRange = (ts1 >= 1) & (ts1 <= 3);


% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_mouse)
    nestedCellArray_1 = {};
    nestedCellArray_2 = {}
    concatenated_columns = {};
    % Get the current nested cell array
    nestedCellArray_1 = zall_mouse{i, array_for_means};
    nestedCellArray_2 = zall_mouse{i, array_for_means_second};
    for qq = 1:size(nestedCellArray_1, 2)
        concatenated_columns{qq} = vertcat(nestedCellArray_1{qq}, ...
            nestedCellArray_2{qq});
    end
    nestedCellArray_1 = concatenated_columns; 
    % nestedCellArray_2 = zall_mouse{i, 1};
    % Initialize the nested cell array for storing mean values
    meanNestedCellArray = cell(size(nestedCellArray_1));
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(nestedCellArray_1)
        % Get the current double array
        currentArray = nestedCellArray_1{j};


        % uncomment below if you want to mean center
        % currentArray_mean = mean(currentArray, 2);
        % currentArray = currentArray-currentArray_mean;
        % Compute the mean activity for each row in the time range 0 to 2 seconds
        meanValues = mean(currentArray(:, timeRange), 2);
        % meanValues = max(currentArray(:, timeRange), [], 2);
        
        % Store the mean values in the corresponding cell of the nested cell array
        meanNestedCellArray{j} = meanValues;
    end
    
    % Store the nested cell array of mean values in the corresponding cell of the main cell array
    meanZallMouse{i} = meanNestedCellArray;
    concatenated_columns_array{i} = concatenated_columns;
    clear meanNestedCellArray
end

% Now, meanZallMouse contains the mean activity for each row in the time period 0 to 2 seconds
% Each cell in meanZallMouse contains a nested cell array with the







% Initialize the new cell arrays to store the correlation results
correlationResults_Type1_2 = cell(size(meanZallMouse)); % For trial_types_by_mouse == 1.2
correlationResults_Type0_3 = cell(size(meanZallMouse)); % For trial_types_by_mouse == 0.3

% Iterate through each level of meanZallMouse
for i = 1:length(meanZallMouse)
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{i};
    % Initialize the nested cell arrays for storing correlation results
    correlationNestedArray_Type1_2 = zeros(size(meanNestedCellArray));
    correlationNestedArray_Type0_3 = zeros(size(meanNestedCellArray));

    
    % Determine the corresponding index in variable_to_correlate
    trialIndex = mod(i-1, length(variable_to_correlate)) + 1;
    
    % Get the corresponding trial choice times and trial types
    trialChoiceTimes = variable_to_correlate{trialIndex}';
    trialTypes = trial_types_by_mouse{trialIndex};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        % Get the current mean values array
        meanValues = meanNestedCellArray{j}';

        % Check for matching lengths and extract subsets based on trial types
        if length(trialChoiceTimes) == length(meanValues)
            % Subset for trial_types_by_mouse == 1.2
            indices_Type1_2 = trialTypes == 1.2;

            if sum(indices_Type1_2) > 1
                correlationNestedArray_Type1_2(j) = corr(meanValues(indices_Type1_2)', trialChoiceTimes(indices_Type1_2)');
            else
                correlationNestedArray_Type1_2(j) = NaN; % Handle insufficient data
            end
            
            % Subset for trial_types_by_mouse == 0.3
            indices_Type0_3 = trialTypes == 0.3;
            if sum(indices_Type0_3) > 1
                correlationNestedArray_Type0_3(j) = corr(meanValues(indices_Type0_3)', trialChoiceTimes(indices_Type0_3)');
            else
                correlationNestedArray_Type0_3(j) = NaN; % Handle insufficient data
            end
        elseif length(trialChoiceTimes) < length(meanValues)
            % Adjust length mismatch scenario for trial_types_by_mouse == 1.2
            indices_Type1_2 = trialTypes == 1.2;
            if sum(indices_Type1_2) > 1
                correlationNestedArray_Type1_2(j) = corr(meanValues(1:end-1 & indices_Type1_2), trialChoiceTimes(indices_Type1_2));
            else
                correlationNestedArray_Type1_2(j) = NaN;
            end
            
            % Adjust length mismatch scenario for trial_types_by_mouse == 0.3
            indices_Type0_3 = trialTypes == 0.3;
            if sum(indices_Type0_3) > 1
                correlationNestedArray_Type0_3(j) = corr(meanValues(1:end-1 & indices_Type0_3), trialChoiceTimes(indices_Type0_3));
            else
                correlationNestedArray_Type0_3(j) = NaN;
            end
        else
            % If lengths do not match, handle the mismatch by setting correlation to NaN
            correlationNestedArray_Type1_2(j) = NaN;
            correlationNestedArray_Type0_3(j) = NaN;
        end
    end

    % Store the nested cell arrays of correlation coefficients in the main cell arrays
    correlationResults_Type1_2{i} = correlationNestedArray_Type1_2;
    correlationResults_Type0_3{i} = correlationNestedArray_Type0_3;
end

% The results are stored in correlationResults_Type1_2 and correlationResults_Type0_3

% Now, correlationResults contains the correlation coefficients for each nested structure in meanZallMouse


% Assuming correlationResults_Type1_2 and correlationResults_Type0_3 are defined and contain the correlation coefficients

% Initialize empty arrays to collect all correlation coefficients
allCorrelations_Type1_2 = [];
allCorrelations_Type0_3 = [];

% Iterate through each level of correlationResults_Type1_2 and correlationResults_Type0_3
for i = 1:length(correlationResults_Type1_2)
    % Get the current nested arrays of correlation coefficients
    correlationNestedArray_Type1_2 = correlationResults_Type1_2{i};
    correlationNestedArray_Type0_3 = correlationResults_Type0_3{i};
    
    % Iterate through each cell in the nested arrays
    for j = 1:length(correlationNestedArray_Type1_2)
        % Get the current correlation coefficients
        correlationCoeff_Type1_2 = correlationNestedArray_Type1_2(j);
        correlationCoeff_Type0_3 = correlationNestedArray_Type0_3(j);
        
        % Check if the coefficients are not NaN and append to respective arrays
        if ~isnan(correlationCoeff_Type1_2)
            allCorrelations_Type1_2 = [allCorrelations_Type1_2; correlationCoeff_Type1_2];
        end
        if ~isnan(correlationCoeff_Type0_3)
            allCorrelations_Type0_3 = [allCorrelations_Type0_3; correlationCoeff_Type0_3];
        end
    end
end

% Now, allCorrelations_Type1_2 and allCorrelations_Type0_3 contain all the correlation coefficients for the respective types

% Create histograms of the correlation coefficients for both types
figure;
subplot(2, 1, 1);
histogram(allCorrelations_Type1_2);
xlabel('Correlation Coefficient (Type 1.2)');
ylabel('Frequency');
title('Histogram of Correlation Coefficients (Type 1.2)');

subplot(2, 1, 2);
histogram(allCorrelations_Type0_3);
xlabel('Correlation Coefficient (Type 0.3)');
ylabel('Frequency');
title('Histogram of Correlation Coefficients (Type 0.3)');

% Optionally, you can add a vertical line at 0 for reference in both histograms
subplot(2, 1, 1);
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;

subplot(2, 1, 2);
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;


only_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(prechoice_block_1 == 1);
not_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(prechoice_block_1 ~= 1);

only_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(prechoice_block_1 == 1);
not_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(prechoice_block_1 ~= 1);

% Create histograms for SHK responsive correlations for both types
figure;
subplot(2, 2, 1);
histogram(only_shk_responsive_corrs_Type1_2);
xlabel('Correlation Coefficient (Type 1.2)');
ylabel('Frequency');
title('Histogram of SHK Responsive Correlations (Type 1.2)');

subplot(2, 2, 2);
histogram(not_shk_responsive_corrs_Type1_2);
xlabel('Correlation Coefficient (Type 1.2)');
ylabel('Frequency');
title('Histogram of Non-SHK Responsive Correlations (Type 1.2)');

subplot(2, 2, 3);
histogram(only_shk_responsive_corrs_Type0_3);
xlabel('Correlation Coefficient (Type 0.3)');
ylabel('Frequency');
title('Histogram of SHK Responsive Correlations (Type 0.3)');

subplot(2, 2, 4);
histogram(not_shk_responsive_corrs_Type0_3);
xlabel('Correlation Coefficient (Type 0.3)');
ylabel('Frequency');
title('Histogram of Non-SHK Responsive Correlations (Type 0.3)');

% Optionally, you can add vertical lines at 0 for reference in all histograms
for i = 1:4
    subplot(2, 2, i);
    hold on;
    yLimits = ylim;
    plot([0 0], yLimits, 'r--', 'LineWidth', 2);
    hold off;
end



% Assuming the following variables are defined:
% allCorrelations_Type1_2, allCorrelations_Type0_3: arrays containing correlation coefficients for Types 1.2 and 0.3
% only_shk_responsive_corrs_Type1_2, only_shk_responsive_corrs_Type0_3: arrays containing SHK responsive correlation coefficients for Types 1.2 and 0.3
% not_shk_responsive_corrs_Type1_2, not_shk_responsive_corrs_Type0_3: arrays containing non-SHK responsive correlation coefficients for Types 1.2 and 0.3

% Calculate means for Type 1.2
mean_only_shk_Type1_2 = mean(only_shk_responsive_corrs_Type1_2);
mean_not_shk_Type1_2 = mean(not_shk_responsive_corrs_Type1_2);

% Calculate means for Type 0.3
mean_only_shk_Type0_3 = mean(only_shk_responsive_corrs_Type0_3);
mean_not_shk_Type0_3 = mean(not_shk_responsive_corrs_Type0_3);

% Create histograms for both types
figure;
width = 250; % Width of the figure
height = 1000; % Height of the figure (double the width)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]

% Subplot for Type 1.2
subplot(2, 1, 1);
histogram(not_shk_responsive_corrs_Type1_2, 'Normalization', 'probability', 'FaceColor', 'blue', 'BinWidth', 0.05, 'LineStyle', 'none');
hold on;
histogram(only_shk_responsive_corrs_Type1_2, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05, 'LineStyle', 'none');
xline(mean_only_shk_Type1_2, 'r');
xline(mean_not_shk_Type1_2, 'g');
xlabel('Correlation Coefficient (Type 1.2)');
ylabel('Probability');
title('Histograms of Correlation Coefficients (Type 1.2)');
hold off;

% Subplot for Type 0.3
subplot(2, 1, 2);
histogram(not_shk_responsive_corrs_Type0_3, 'Normalization', 'probability', 'FaceColor', 'blue', 'BinWidth', 0.05, 'LineStyle', 'none');
hold on;
histogram(only_shk_responsive_corrs_Type0_3, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05, 'LineStyle', 'none');
xline(mean_only_shk_Type0_3, 'r');
xline(mean_not_shk_Type0_3, 'g');
xlabel('Correlation Coefficient (Type 0.3)');
ylabel('Probability');
title('Histograms of Correlation Coefficients (Type 0.3)');
hold off;

% Optionally, add vertical lines at 0 for reference
for i = 1:2
    subplot(2, 1, i);
    hold on;
    yLimits = ylim;
    plot([0 0], yLimits, 'k', 'LineWidth', 2);
    xtickformat('%.2f');
    ytickformat('%.2f');
    hold off;
end

% Perform Kolmogorov-Smirnov tests for both types
[h_Type1_2, p_Type1_2] = kstest2(not_shk_responsive_corrs_Type1_2, only_shk_responsive_corrs_Type1_2);
[h_Type0_3, p_Type0_3] = kstest2(not_shk_responsive_corrs_Type0_3, only_shk_responsive_corrs_Type0_3);

% Display the results of the statistical tests for Type 1.2
fprintf('Kolmogorov-Smirnov test result for Type 1.2:\n');
fprintf('h = %d (0 means the null hypothesis cannot be rejected, 1 means it can be rejected)\n', h_Type1_2);
fprintf('p-value = %.4f\n', p_Type1_2);

% Display the results of the statistical tests for Type 0.3
fprintf('Kolmogorov-Smirnov test result for Type 0.3:\n');
fprintf('h = %d (0 means the null hypothesis cannot be rejected, 1 means it can be rejected)\n', h_Type0_3);
fprintf('p-value = %.4f\n', p_Type0_3);

% Perform t-tests for both types
[h_ttest_Type1_2, p_ttest_Type1_2, ci_Type1_2, stats_Type1_2] = ttest2(not_shk_responsive_corrs_Type1_2, only_shk_responsive_corrs_Type1_2);
[h_ttest_Type0_3, p_ttest_Type0_3, ci_Type0_3, stats_Type0_3] = ttest2(not_shk_responsive_corrs_Type0_3, only_shk_responsive_corrs_Type0_3);

% Display t-test results for Type 1.2
fprintf('T-test result for Type 1.2:\n');
fprintf('h = %d\n', h_ttest_Type1_2);
fprintf('p-value = %.4f\n', p_ttest_Type1_2);

% Display t-test results for Type 0.3
fprintf('T-test result for Type 0.3:\n');
fprintf('h = %d\n', h_ttest_Type0_3);
fprintf('p-value = %.4f\n', p_ttest_Type0_3);








bar_separation_value = 3;

figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]

% Swarm chart for Type 1.2 correlations
swarmchart(ones(1, length(only_shk_responsive_corrs_Type1_2)), only_shk_responsive_corrs_Type1_2, 'o', 'MarkerFaceColor', 'red');
hold on;
swarmchart(ones(1, length(not_shk_responsive_corrs_Type1_2)) * bar_separation_value, not_shk_responsive_corrs_Type1_2, 'o', 'MarkerFaceColor', 'blue');

% Plot mean lines for Type 1.2 correlations
plot([0.5; 1.5], [mean(only_shk_responsive_corrs_Type1_2); mean(only_shk_responsive_corrs_Type1_2)], 'LineWidth', 3, 'Color', 'red');
plot([bar_separation_value - 0.5; bar_separation_value + 0.5], [mean(not_shk_responsive_corrs_Type1_2); mean(not_shk_responsive_corrs_Type1_2)], 'LineWidth', 3, 'Color', 'blue');

% Swarm chart for Type 0.3 correlations
bar_separation_value_0_3 = bar_separation_value + 3;
swarmchart(ones(1, length(only_shk_responsive_corrs_Type0_3)) + 2 * bar_separation_value, only_shk_responsive_corrs_Type0_3, 'o', 'MarkerFaceColor', 'magenta');
swarmchart(ones(1, length(not_shk_responsive_corrs_Type0_3)) + 3 * bar_separation_value, not_shk_responsive_corrs_Type0_3, 'o', 'MarkerFaceColor', 'cyan');

% Plot mean lines for Type 0.3 correlations
plot([bar_separation_value_0_3 - 0.5; bar_separation_value_0_3 + 0.5], [mean(only_shk_responsive_corrs_Type0_3); mean(only_shk_responsive_corrs_Type0_3)], 'LineWidth', 3, 'Color', 'magenta');
plot([bar_separation_value_0_3 + 2.5; bar_separation_value_0_3 + 3.5], [mean(not_shk_responsive_corrs_Type0_3); mean(not_shk_responsive_corrs_Type0_3)], 'LineWidth', 3, 'Color', 'cyan');

% General plot settings
yline(0, 'k--');
xticks([1, bar_separation_value, bar_separation_value + 2, bar_separation_value_0_3, bar_separation_value_0_3 + 3]);
xticklabels({'Type 1.2 SHK Resp.', 'Type 1.2 Not SHK Resp.', 'Type 0.3 SHK Resp.', 'Type 0.3 Not SHK Resp.'});
xtickformat('%.1f');
ytickformat('%.1f');
hold off;



%% plot scatters for individual neurons & behav variables

%CHECK AND UPDATE START & END TIME DEPENDING ON EVENT OF INTEREST

start_time = -4;% sub-window start time
end_time = 0; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% update mouse number {10, 1} etc and which zall_mouse arrays to combime
% (e.g., 10, 1 for mouse 10 prechoice block 1, and 10, 5 for mouse 10
% prechoice blocks 2/3
for qq = 1:size(zall_mouse{10, 1}, 2)
    concatenated_zall{qq} = vertcat(zall_mouse{10, 1}{qq}, ...
        zall_mouse{10, 5}{qq});
end

%update  concatenated_zall{1, 37} number after "1," and
sub_window_activity_session_1 = concatenated_zall{1, 85}(:, sub_window_idx);
choice_times_mouse = trial_choice_times_by_mouse{1, 10};
trial_types = trial_types_by_mouse{1, 10};

% % Assume A and B are your 143x21 arrays
% correlation_coefficients = arrayfun(@(i) corr(sub_window_activity_session_1 (i, :)', sub_window_activity_session_2 (i, :)'), 1:size(sub_window_activity_session_1 , 1));


mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);


x = mean_sub_window_activity_session_1;
y = choice_times_mouse;





% Define colors based on trial types
colors = repmat([0.5, 0.5, 0.5], length(trial_types), 1); % Default to gray
colors(trial_types == 1.2, :) = repmat([0, 0, 1], sum(trial_types == 1.2), 1); % Blue for trial_types == 1.2
colors(trial_types == 0.3, :) = repmat([1, 0, 0], sum(trial_types == 0.3), 1); % Red for trial_types == 0.3

% Create scatter plot
figure;
set(gcf, 'Position', [100, 100, 200, 200]); % Adjust figure position and size
scatter(x, y, 36, colors, 'filled', 'MarkerEdgeColor', 'k'); % Use 'colors' for MarkerFaceColor

hold on;

% Add a regression line
coefficients = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);
plot(x_fit, y_fit, 'r');

% Calculate R-squared value
y_pred = polyval(coefficients, x);
ssr = sum((y_pred - mean(y)).^2);
sst = sum((y - mean(y)).^2);
r_squared = ssr / sst;

% Add R-squared value to the plot
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared)], 'FontSize', 12);

% Add labels and a legend if needed
xlabel('Mean Sub-window Activity Session 1');
ylabel('Choice Times Mouse');
title('Scatter Plot with Regression Line and R^2 Value');
hold off;


%% these arrays are just the combined_data = [mean_sub_window_activity_session_1 , mean_sub_window_activity_session_2] arrays from above, but saved for Early & Late
load('Late_RM_arrays.mat')
load('RM_D1_arrays.mat')

% Calculate means
diff_score_prechoice_late = Late_RM_prechoice_combined_data(:, 1) - Late_RM_prechoice_combined_data(:, 2);

diff_score_prechoice_early = RM_D1_prechoice_combined_data(:, 1) - RM_D1_prechoice_combined_data(:, 2); 

% Calculate means
diff_score_postchoice_late = Late_RM_postchoice_combined_data(:, 1) - Late_RM_postchoice_combined_data(:, 2);

diff_score_postchoice_early = RM_D1_postchoice_combined_data(:, 1) - RM_D1_postchoice_combined_data(:, 2); 

diff_score_collect_late = Late_RM_collect_combined_data(:, 1) - Late_RM_collect_combined_data(:, 2);

diff_score_collect_early = RM_D1_collect_combined_data(:, 1) - RM_D1_collect_combined_data(:, 2); 

figure;

% Define figure size
width = 300; % Width of the figure
height = 500; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]

% Prechoice histograms
data_prechoice = {diff_score_prechoice_early, diff_score_prechoice_late};
data_postchoice = {diff_score_postchoice_early, diff_score_postchoice_late};
data_collect = {diff_score_collect_early, diff_score_collect_late};

data_all = {data_prechoice, data_postchoice, data_collect};
titles = {'Prechoice', 'Postchoice', 'Collect'};
colors = {'blue', 'red'};

num_rows = 3;
num_cols = 2;

for i = 1:num_rows
    for j = 1:num_cols
        idx = (i - 1) * num_cols + j;
        subplot(num_rows, num_cols, idx);
        histogram(data_all{i}{j}, 'Normalization', 'probability', 'FaceColor', colors{j}, 'BinWidth', 0.05, 'LineStyle', 'none');
        xline(0, 'k-', 'LineWidth', 1); % Add vertical line at x = 0
        if j == 1
            ylabel('Probability');
        end
        if i == num_rows
            xlabel('Difference Score');
        end
        if j == 1
            title([titles{i}, ' Early']);
        else
            title([titles{i}, ' Late']);
        end
    end
end

% Set shared limits for X and Y axes
all_axes = findall(gcf, 'Type', 'axes');
y_limits = [0 max(cellfun(@(ax) max(ylim(ax)), num2cell(all_axes)))];
x_limits = [min(cellfun(@(ax) min(xlim(ax)), num2cell(all_axes))), max(cellfun(@(ax) max(xlim(ax)), num2cell(all_axes)))];
for ax = all_axes'
    ylim(ax, y_limits);
    xlim(ax, x_limits);
end

% Perform a Kolmogorov-Smirnov test to compare the first set of prechoice distributions
[h, p] = kstest2(diff_score_prechoice_late, diff_score_prechoice_early);

% Display the results of the statistical test
fprintf('Kolmogorov-Smirnov test result:\n');
fprintf('h = %d (0 means the null hypothesis cannot be rejected, 1 means it can be rejected)\n', h);
fprintf('p-value = %.4f\n', p);

% Perform a t-test to compare the means of the two distributions
[h, p, ci, stats] = ttest2(diff_score_prechoice_late, diff_score_prechoice_early);

%%
array_for_means = 1; 
array_for_means_second = 5; 
% Initialize the new cell array to store the mean values
meanZallMouse = cell(size(zall_mouse, 2), 1);

% Define the time range for 0 to 2 seconds
timeRange = (ts1 >= -4) & (ts1 <= 0);
% timeRange = (ts1 >= 0) & (ts1 <= 2);
% timeRange = (ts1 >= 1) & (ts1 <= 3);

concatenated_zall = {};
% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_array)
    nestedCellArray_1 = [];
    nestedCellArray_2 = []

    % Get the current nested cell array
    nestedCellArray_1 = zall_array{array_for_means, i};
    nestedCellArray_2 = zall_array{array_for_means_second, i};

    concatenated_zall{i} = vertcat(nestedCellArray_1, ...
        nestedCellArray_2);


end

%%

prechoice_concatenated_zall = concatenated_zall(prechoice_block_1 == 1);
postchoice_concatenated_zall = concatenated_zall(postchoice_reward_block_1 == 1);
collect_concatenated_zall = concatenated_zall(collect_block_1 == 1);

% Get the number of cells and the number of columns in the arrays
num_cells = length(prechoice_concatenated_zall);
num_columns = size(prechoice_concatenated_zall{1}, 2); % Assuming all arrays have the same dimensions

% Initialize variables to store the means for each part
part1_means = zeros(num_cells, num_columns);
part2_means = zeros(num_cells, num_columns);
part3_means = zeros(num_cells, num_columns);
part4_means = zeros(num_cells, num_columns);
part5_means = zeros(num_cells, num_columns);

% Loop through each cell in the cell array
for i = 1:num_cells
    % Extract the current 90x160 double array
    current_array = prechoice_concatenated_zall{i};
    
    % Compute the mean for each part (18 rows each)
    part1_means(i, :) = mean(current_array(1:18, :), 1);
    part2_means(i, :) = mean(current_array(19:36, :), 1);
    part3_means(i, :) = mean(current_array(37:54, :), 1);
    part4_means(i, :) = mean(current_array(55:72, :), 1);
    part5_means(i, :) = mean(current_array(73:90, :), 1);
end


%% big heatmap for categories


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
