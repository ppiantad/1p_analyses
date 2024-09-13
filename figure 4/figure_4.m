%% for heatmap, change "plot_num" (what neuron to plot) and array_to_plot (# corresponds to which dataset)

plot_num = 178; % 81 / 31 or 58 or 70 / 2

array_to_plot = [1 8]; % depends on the structure of zall

select_mouse = 'BLA_Insc_35';

% for RDT D1 BLA_Insc_25:
%prechoice neuron num 46
%postchoice rew num 38
%consumption num 39
%shock num 11

select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';

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



time2Collect = BehavData.collectionTime(trials_per_mouse{select_mouse_index, array_to_plot(1)}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot(1)});
trialStartTime = BehavData.stTime(trials_per_mouse{select_mouse_index, array_to_plot(1)}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot(1)});
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

custom_colormap = [
    1, 1, 1;         % white
    0.9, 0.95, 0.95;
    0.8, 0.9, 0.9;
    0.6, 0.85, 0.85;
    0.4, 0.8, 0.8;
    0.2, 0.8, 0.8;
    0.0, 0.8, 0.8;   % robin's egg blue
];

% custom_colormap = [
%     1, 1, 1;         % white
%     0.9, 0.9, 0.95;
%     0.8, 0.8, 0.9;
%     0.6, 0.6, 0.8;
%     0.4, 0.4, 0.7;
%     0.2, 0.2, 0.6;
%     0.0, 0.0, 0.55;   % dark blue
% ];

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
figure('Position', [100, 100, 250, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile: Line plot
ax1 = nexttile(t);
hold on;
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', []);
shadedErrorBar(ts1, nanmean(zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num}), nanmean(neuron_sem_array{1, 1}(prechoice_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num}), nanmean(neuron_sem_array{1, 1}(postchoice_reward_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});

xline(0);
xline(median_trialStartTime, 'g')
xline(median_time2Collect, 'r')
% xlabel('Time from Large Rew Choice (s)');
% legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.8 0.8]);
hold off


% First tile (heatmap)
ax2 = nexttile(t);
hold on;

% Plot the heatmap
imagesc(ts1, 1:size(zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num}, 1), zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num});


% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);



ylim([0.5, size(zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num}, 1) + 0.5]);

xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', []);
set(gca, 'YTick', [1, size(zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num}, 1)]);
xline(0)
scatter(time2Collect, Tris               , 'Marker', 'p')
scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

% Second tile (mean and raw data)
ax3 = nexttile(t);
hold on;

time2Collect = BehavData.collectionTime(trials_per_mouse{select_mouse_index, array_to_plot(2)}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot(2)});
trialStartTime = BehavData.stTime(trials_per_mouse{select_mouse_index, array_to_plot(2)}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot(2)});
median_trialStartTime = median(trialStartTime)
median_time2Collect = median(time2Collect)
xline(median_trialStartTime)
xline(median_time2Collect)
[numTrials, ~] = size(time2Collect);
Tris = [1:numTrials]';



hold on;

% Plot the heatmap
imagesc(ts1, 1:size(zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num}, 1), zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num});


% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar(ax2, 'eastoutside');
set(c, 'YTick', clim); % 

ylim([0.5, size(zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num}, 1) + 0.5]);

xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', []);
set(gca, 'YTick', [1, size(zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num}, 1)]);
xline(0)
scatter(time2Collect, Tris               , 'Marker', 'p')
scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar(ax3, 'eastoutside');
set(c, 'YTick', clim); % 



%%

%These data can be used to plot the median or mean choice
% time on a PCA graph, for example

behav_tbl_iter_single = behav_tbl_iter(8);

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
hold on
% Create a histogram for allCorrelations

width = 200; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);

shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(prechoice_block_1==1, :)), nanmean(neuron_sem_array{1, 8}(prechoice_block_1==1, :)), 'lineProps', {'color', 'r'});
hold on; hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 9}(postchoice_reward_block_1==1, :)), nanmean(neuron_sem_array{1, 9}(postchoice_reward_block_1==1, :)), 'lineProps', {'color', 'k'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 8}(collect_block_1==1, :)), 'lineProps', {'color', 'b'});

xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.8 0.8]);
ytickformat('%.1f');
hold off




%%
% Get AUCs for the relevant periods for the 3 defined events
% Define time windows
pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [1 3];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

% Initialize arrays to store AUCs
action_auc_pre_choice_block_1 = [];
action_auc_post_choice_block_1 = [];
action_auc_consumption_block_1= [];

action_auc_pre_choice_block_2_3 = [];
action_auc_post_choice_block_2_3 = [];
action_auc_consumption_block_2_3 = [];
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

post_choice_auc_pre_choice_block_1 = []
post_choice_auc_post_choice_block_1 = [];
post_choice_auc_consumption_block_1= [];

post_choice_auc_pre_choice_block_2_3 = [];
post_choice_auc_post_choice_block_2_3 = [];
post_choice_auc_consumption_block_2_3 = [];
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

consumption_auc_pre_choice_block_1 = [];
consumption_auc_post_choice_block_1 = [];
consumption_auc_consumption_block_1= [];

consumption_auc_pre_choice_block_2_3 = [];
consumption_auc_post_choice_block_2_3 = [];
consumption_auc_consumption_block_2_3 = [];
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


mean_change_from_block_1_pre_choice = [mean_pre_choice_block_2_3 - mean_pre_choice_block_1];
mean_change_from_block_1_post_choice = [mean_post_choice_block_2_3 - mean_post_choice_block_1];
mean_change_from_block_1_consumption = [mean_consumption_block_2_3 - mean_consumption_block_1];

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

shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_prechoice==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_prechoice==1, :)), 'lineProps', {'color', 'r'});
hold on; hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 9}(conserved_postchoice==1, :)), nanmean(neuron_sem_array{1, 9}(conserved_postchoice==1, :)), 'lineProps', {'color', 'k'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});

xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.8 0.8]);
ytickformat('%.1f');
hold off

%%
% Get AUCs for the relevant periods for the 3 defined events
% Define time windows
pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [1 3];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

% Initialize arrays to store AUCs
action_auc_pre_choice_block_1 = [];
action_auc_post_choice_block_1 = [];
action_auc_consumption_block_1= [];

action_auc_pre_choice_block_2_3 = [];
action_auc_post_choice_block_2_3 = [];
action_auc_consumption_block_2_3 = [];
pre_choice_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,1}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if conserved_prechoice(i) == 1
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
end

post_choice_auc_pre_choice_block_1 = []
post_choice_auc_post_choice_block_1 = [];
post_choice_auc_consumption_block_1= [];

post_choice_auc_pre_choice_block_2_3 = [];
post_choice_auc_post_choice_block_2_3 = [];
post_choice_auc_consumption_block_2_3 = [];
post_choice_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,1}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if conserved_postchoice(i) == 1
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
end

consumption_auc_pre_choice_block_1 = [];
consumption_auc_post_choice_block_1 = [];
consumption_auc_consumption_block_1= [];

consumption_auc_pre_choice_block_2_3 = [];
consumption_auc_post_choice_block_2_3 = [];
consumption_auc_consumption_block_2_3 = [];
consumption_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,3}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if conserved_consumption(i) == 1
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


mean_change_from_block_1_pre_choice = [mean_pre_choice_block_2_3 - mean_pre_choice_block_1];
mean_change_from_block_1_post_choice = [mean_post_choice_block_2_3 - mean_post_choice_block_1];
mean_change_from_block_1_consumption = [mean_consumption_block_2_3 - mean_consumption_block_1];

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
ylim([-10 6]);
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

shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(remapped_prechoice==1, :)), nanmean(neuron_sem_array{1, 8}(remapped_prechoice==1, :)), 'lineProps', {'color', 'r'});
hold on; hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 9}(remapped_postchoice==1, :)), nanmean(neuron_sem_array{1, 9}(remapped_postchoice==1, :)), 'lineProps', {'color', 'k'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(remapped_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(remapped_consumption==1, :)), 'lineProps', {'color', 'b'});

xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.8 0.8]);
ytickformat('%.1f');
hold off

%%
% Get AUCs for the relevant periods for the 3 defined events
% Define time windows
pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [1 3];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

% Initialize arrays to store AUCs
action_auc_pre_choice_block_1 = [];
action_auc_post_choice_block_1 = [];
action_auc_consumption_block_1= [];

action_auc_pre_choice_block_2_3 = [];
action_auc_post_choice_block_2_3 = [];
action_auc_consumption_block_2_3 = [];
pre_choice_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,1}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if remapped_prechoice(i) == 1

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

post_choice_auc_pre_choice_block_1 = []
post_choice_auc_post_choice_block_1 = [];
post_choice_auc_consumption_block_1= [];

post_choice_auc_pre_choice_block_2_3 = [];
post_choice_auc_post_choice_block_2_3 = [];
post_choice_auc_consumption_block_2_3 = [];
post_choice_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,1}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if remapped_postchoice(i) == 1

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

consumption_auc_pre_choice_block_1 = [];
consumption_auc_post_choice_block_1 = [];
consumption_auc_consumption_block_1= [];

consumption_auc_pre_choice_block_2_3 = [];
consumption_auc_post_choice_block_2_3 = [];
consumption_auc_consumption_block_2_3 = [];
consumption_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,3}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if remapped_consumption(i) == 1

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
ylim([-10 6]);
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


%%
 mean_data_array = {neuron_mean_array{1, 1}(prechoice_block_1==1, :), neuron_mean_array{1, 8}(prechoice_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 1}(prechoice_block_1==1, :), neuron_sem_array{1, 8}(prechoice_block_1==1, :)};

 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1);



 %%
 mean_data_array = {neuron_mean_array{1, 2}(postchoice_reward_block_1==1, :), neuron_mean_array{1, 9}(postchoice_reward_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 2}(postchoice_reward_block_1==1, :), neuron_sem_array{1, 9}(postchoice_reward_block_1==1, :)};

 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1)

 %%
 if size(respClass_all_array, 2) == 10 | size(respClass_all_array, 2) == 11
    comparison_arrays = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays = [1 2 3; 4 5 6]
elseif size(respClass_all_array, 2) == 7
    comparison_arrays = [1 2 3; 5 6 7]
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


%%
% arrays_to_examine = [comparison_arrays(1, 3) comparison_arrays(2, 3)];
conserved_prechoice = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
conserved_prechoice_sum = sum(conserved_prechoice)

lost_prechoice = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 ~= event_for_figures;
lost_prechoice_sum = sum(lost_prechoice)

remapped_prechoice = prechoice_block_1 ~= event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
remapped_prechoice_sum = sum(remapped_prechoice)

% vars_to_use = {'collect_block_1', 'collect_blocks_2_and_3'};


mean_data_array = {neuron_mean_array{1, 1}(conserved_prechoice  ==1, :), neuron_mean_array{1, 8}(conserved_prechoice  ==1, :), neuron_mean_array{1, 8}(lost_prechoice  ==1, :)}
sem_data_array = {neuron_sem_array{1, 1}(conserved_prechoice  ==1, :), neuron_sem_array{1, 8}(conserved_prechoice  ==1, :), neuron_sem_array{1, 8}(lost_prechoice  ==1, :)}
[comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1)


%%
% arrays_to_examine = [comparison_arrays(1, 3) comparison_arrays(2, 3)];
conserved_postchoice_reward = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
conserved_postchoice_reward_sum = sum(conserved_postchoice_reward)

lost_postchoice_reward = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 ~= event_for_figures;
lost_postchoice_reward_sum = sum(lost_postchoice_reward)

remapped_postchoice_reward = postchoice_reward_block_1 ~= event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
remapped_postchoice_reward_sum = sum(remapped_postchoice_reward)

% vars_to_use = {'collect_block_1', 'collect_blocks_2_and_3'};


mean_data_array = {neuron_mean_array{1, 2}(conserved_postchoice_reward  ==1, :), neuron_mean_array{1, 6}(conserved_postchoice_reward  ==1, :), neuron_mean_array{1, 6}(lost_postchoice_reward  ==1, :)}
sem_data_array = {neuron_sem_array{1, 2}(conserved_postchoice_reward  ==1, :), neuron_sem_array{1, 6}(conserved_postchoice_reward  ==1, :), neuron_sem_array{1, 6}(lost_postchoice_reward  ==1, :)}
[comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1)

%%
% arrays_to_examine = [comparison_arrays(1, 3) comparison_arrays(2, 3)];
conserved_collection = collect_block_1 == event_for_figures & collect_blocks_2_and_3 == event_for_figures;
conserved_collection_sum = sum(conserved_collection)

lost_collection = collect_block_1 == event_for_figures & collect_blocks_2_and_3 ~= event_for_figures;
lost_collection_sum = sum(lost_collection)

remapped_collection = collect_block_1 ~= event_for_figures & collect_blocks_2_and_3 == event_for_figures;
remapped_collection_sum = sum(remapped_collection)

% vars_to_use = {'collect_block_1', 'collect_blocks_2_and_3'};


mean_data_array = {neuron_mean_array{1, 3}(conserved_collection  ==1, :), neuron_mean_array{1, 10}(conserved_collection  ==1, :), neuron_mean_array{1, 10}(lost_collection  ==1, :), neuron_mean_array{1, 10}(remapped_collection  ==1, :)}
sem_data_array = {neuron_sem_array{1, 3}(conserved_collection  ==1, :), neuron_sem_array{1, 10}(conserved_collection  ==1, :), neuron_sem_array{1, 10}(lost_collection  ==1, :), neuron_sem_array{1, 10}(remapped_collection  ==1, :)}
[comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1)


%%
lost_all = [lost_prechoice_sum, lost_postchoice_reward_sum, lost_collection_sum]
conserved_all = [conserved_prechoice_sum, conserved_postchoice_reward_sum, conserved_collection_sum]
remapped_all = [remapped_prechoice_sum, remapped_postchoice_reward_sum, remapped_collection_sum]

%% for plotting changes on a donut (specific) and pie (broad) charts
all_conserved_sum = sum(conserved_all)
all_lost_sum = sum(lost_all)
all_remapped_sum = sum(remapped_all)
remaining_neurons = neuron_num - (all_conserved_sum + all_lost_sum +all_remapped_sum);

figure;
piechart([all_conserved_sum/neuron_num, all_lost_sum/neuron_num, all_remapped_sum/neuron_num, remaining_neurons/neuron_num])

figure;
donutchart([conserved_all/neuron_num, lost_all/neuron_num, remapped_all/neuron_num, remaining_neurons/neuron_num])

% all_neurons = [conserved_sum; lost_sum; remapped_sum]
%all_neurons_2 = [conserved_sum; lost_sum; remapped_sum]


