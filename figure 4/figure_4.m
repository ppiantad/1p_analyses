%% for heatmap, change "plot_num" (what neuron to plot) and array_to_plot (# corresponds to which dataset)

plot_num = 2 % 81 / 31 or 58 or 70 / 2

array_to_plot = [1 8]; % depends on the structure of zall

select_mouse = 'BLA_Insc_34';

% for RDT D1 BLA_Insc_24:
%prechoice neuron num 1
%postchoice rew num 38
%consumption num 39
%shock num 11

% for RDT D1 BLA_Insc_27:
%prechoice neuron num 72
%postchoice rew num 4 or 47
%consumption num 34
%shock num 11

% for RDT D1 BLA_Insc_34:
%prechoice neuron num 18
%postchoice rew num 104
%consumption num 2


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

% xlabel('Time from Large Rew Choice (s)');
% legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-1 1]);
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
xline(median_trialStartTime, 'g')
xline(median_time2Collect, 'r')
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
xline(median_trialStartTime, 'g')
xline(median_time2Collect, 'r')
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

 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 1]);



 %%
 mean_data_array = {neuron_mean_array{1, 2}(postchoice_reward_block_1==1, :), neuron_mean_array{1, 9}(postchoice_reward_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 2}(postchoice_reward_block_1==1, :), neuron_sem_array{1, 9}(postchoice_reward_block_1==1, :)};

 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 1])

 %%
 mean_data_array = {neuron_mean_array{1, 3}(collect_block_1==1, :), neuron_mean_array{1, 10}(collect_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 3}(collect_block_1==1, :), neuron_sem_array{1, 10}(collect_block_1==1, :)};

 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 1]);



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


mean_data_array = {neuron_mean_array{1, 8}(conserved_prechoice  ==1, :), neuron_mean_array{1, 8}(remapped_prechoice  ==1, :), neuron_mean_array{1, 8}(lost_prechoice  ==1, :)}
sem_data_array = {neuron_sem_array{1, 8}(conserved_prechoice  ==1, :), neuron_sem_array{1, 8}(remapped_prechoice  ==1, :), neuron_sem_array{1, 8}(lost_prechoice  ==1, :)}

% [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.7], 3)
[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.7], 3)


%%
% arrays_to_examine = [comparison_arrays(1, 3) comparison_arrays(2, 3)];
conserved_postchoice_reward = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
conserved_postchoice_reward_sum = sum(conserved_postchoice_reward)

lost_postchoice_reward = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 ~= event_for_figures;
lost_postchoice_reward_sum = sum(lost_postchoice_reward)

remapped_postchoice_reward = postchoice_reward_block_1 ~= event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
remapped_postchoice_reward_sum = sum(remapped_postchoice_reward)

% vars_to_use = {'collect_block_1', 'collect_blocks_2_and_3'};


mean_data_array = {neuron_mean_array{1, 9}(conserved_postchoice_reward  ==1, :), neuron_mean_array{1, 9}(remapped_postchoice_reward  ==1, :), neuron_mean_array{1, 9}(lost_postchoice_reward  ==1, :)}
sem_data_array = {neuron_sem_array{1, 9}(conserved_postchoice_reward  ==1, :), neuron_sem_array{1, 9}(remapped_postchoice_reward  ==1, :), neuron_sem_array{1, 9}(lost_postchoice_reward  ==1, :)}
% [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.8], 3)

[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.8], 3)

%%
% arrays_to_examine = [comparison_arrays(1, 3) comparison_arrays(2, 3)];
conserved_collection = collect_block_1 == event_for_figures & collect_blocks_2_and_3 == event_for_figures;
conserved_collection_sum = sum(conserved_collection)

lost_collection = collect_block_1 == event_for_figures & collect_blocks_2_and_3 ~= event_for_figures;
lost_collection_sum = sum(lost_collection)

remapped_collection = collect_block_1 ~= event_for_figures & collect_blocks_2_and_3 == event_for_figures;
remapped_collection_sum = sum(remapped_collection)

% vars_to_use = {'collect_block_1', 'collect_blocks_2_and_3'};


mean_data_array = {neuron_mean_array{1, 10}(conserved_collection  ==1, :), neuron_mean_array{1, 10}(remapped_collection  ==1, :), neuron_mean_array{1, 10}(lost_collection  ==1, :), neuron_mean_array{1, 10}(remapped_collection  ==1, :)}
sem_data_array = {neuron_sem_array{1, 10}(conserved_collection  ==1, :), neuron_sem_array{1, 10}(remapped_collection  ==1, :), neuron_sem_array{1, 10}(lost_collection  ==1, :), neuron_sem_array{1, 10}(remapped_collection  ==1, :)}
% [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.8], 3)
[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.8], 3)




%%
% arrays_to_examine = [comparison_arrays(1, 3) comparison_arrays(2, 3)];
conserved_prechoice = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
conserved_prechoice_sum = sum(conserved_prechoice)

lost_prechoice = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 ~= event_for_figures;
lost_prechoice_sum = sum(lost_prechoice)

remapped_prechoice = prechoice_block_1 ~= event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
remapped_prechoice_sum = sum(remapped_prechoice)

% vars_to_use = {'collect_block_1', 'collect_blocks_2_and_3'};


mean_data_array = {neuron_mean_array{1, 8}(conserved_prechoice  ==1, :), neuron_mean_array{1, 8}(lost_prechoice  ==1, :), neuron_mean_array{1, 8}(remapped_prechoice  ==1, :)}
sem_data_array = {neuron_sem_array{1, 8}(conserved_prechoice  ==1, :), neuron_sem_array{1, 8}(lost_prechoice  ==1, :), neuron_sem_array{1, 8}(remapped_prechoice  ==1, :)}
[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.7])


%%
% arrays_to_examine = [comparison_arrays(1, 3) comparison_arrays(2, 3)];
conserved_postchoice_reward = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
conserved_postchoice_reward_sum = sum(conserved_postchoice_reward)

lost_postchoice_reward = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 ~= event_for_figures;
lost_postchoice_reward_sum = sum(lost_postchoice_reward)

remapped_postchoice_reward = postchoice_reward_block_1 ~= event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
remapped_postchoice_reward_sum = sum(remapped_postchoice_reward)

% vars_to_use = {'collect_block_1', 'collect_blocks_2_and_3'};


mean_data_array = { neuron_mean_array{1, 6}(conserved_postchoice_reward  ==1, :), neuron_mean_array{1, 9}(lost_postchoice_reward  ==1, :), neuron_mean_array{1, 9}(remapped_postchoice_reward  ==1, :)}
sem_data_array = {neuron_sem_array{1, 6}(conserved_postchoice_reward  ==1, :), neuron_sem_array{1, 9}(lost_postchoice_reward  ==1, :), neuron_sem_array{1, 9}(remapped_postchoice_reward  ==1, :)}
[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.8])


%%
% arrays_to_examine = [comparison_arrays(1, 3) comparison_arrays(2, 3)];
conserved_collection = collect_block_1 == event_for_figures & collect_blocks_2_and_3 == event_for_figures;
conserved_collection_sum = sum(conserved_collection)

lost_collection = collect_block_1 == event_for_figures & collect_blocks_2_and_3 ~= event_for_figures;
lost_collection_sum = sum(lost_collection)

remapped_collection = collect_block_1 ~= event_for_figures & collect_blocks_2_and_3 == event_for_figures;
remapped_collection_sum = sum(remapped_collection)

% vars_to_use = {'collect_block_1', 'collect_blocks_2_and_3'};


mean_data_array = {neuron_mean_array{1, 10}(conserved_collection  ==1, :), neuron_mean_array{1, 10}(lost_collection  ==1, :), neuron_mean_array{1, 10}(remapped_collection  ==1, :)}
sem_data_array = {neuron_sem_array{1, 10}(conserved_collection  ==1, :), neuron_sem_array{1, 10}(lost_collection  ==1, :), neuron_sem_array{1, 10}(remapped_collection  ==1, :)}
[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.8])


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


%%
lost_arrays = lost_prechoice == 1 | lost_postchoice_reward == 1 | lost_collection == 1; 
shk_event = respClass_all_array{1,4} == 1;


shk_and_lost = shk_event == 1 & lost_arrays == 1;

total_modulated = [(sum(shk_event)/neuron_num)*100 (sum(lost_arrays)/neuron_num)*100];
A = total_modulated;
I = (sum(shk_and_lost)/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end

%%

diff_all_mean_prechoice = neuron_mean_array{1, 8}(prechoice_block_1 == 1, :) - neuron_mean_array{1, 1}(prechoice_block_1 == 1, :)
diff_all_sem_prechoice = sqrt(neuron_sem_array{1, 8}(prechoice_block_1 == 1, :).^2 + neuron_sem_array{1, 1}(prechoice_block_1 == 1, :).^2)



diff_all_mean_postchoice = neuron_mean_array{1, 9}(postchoice_reward_block_1 == 1, :) - neuron_mean_array{1, 2}(postchoice_reward_block_1 == 1, :)
diff_all_sem_postchoice = sqrt(neuron_sem_array{1, 9}(postchoice_reward_block_1 == 1, :).^2 + neuron_sem_array{1, 2}(postchoice_reward_block_1 == 1, :).^2)


diff_all_mean_collect = neuron_mean_array{1, 10}(collect_block_1 == 1, :) - neuron_mean_array{1, 3}(collect_block_1 == 1, :)
diff_all_sem_collect = sqrt(neuron_sem_array{1, 10}(collect_block_1 == 1, :).^2 + neuron_sem_array{1, 3}(collect_block_1 == 1, :).^2)


trimmed_diff_all_mean_prechoice = diff_all_mean_prechoice(:, ts1 > -4 & ts1 <= 0);
trimmed_diff_all_sem_prechoice = diff_all_sem_prechoice(:, ts1 > -4 & ts1 <= 0);

trimmed_diff_all_mean_postchoice = diff_all_mean_postchoice(:, ts1 > 0 & ts1 <= 4);
trimmed_diff_all_sem_postchoice = diff_all_sem_postchoice(:, ts1 > 0 & ts1 <= 4);


trimmed_diff_all_mean_collect = diff_all_mean_collect(:, ts1 > 0 & ts1 <= 4);
trimmed_diff_all_sem_collect = diff_all_sem_collect(:, ts1 > 0 & ts1 <= 4);

trimmed_ts1 = ts1(:, ts1 > 0 & ts1 <= 4);

figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(trimmed_ts1 , mean(trimmed_diff_all_mean_prechoice), mean(trimmed_diff_all_sem_prechoice), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(trimmed_ts1 , mean(trimmed_diff_all_mean_postchoice), mean(trimmed_diff_all_sem_postchoice), 'lineProps', {'color', 'r'});
h(3) = shadedErrorBar(trimmed_ts1 , mean(trimmed_diff_all_mean_collect), mean(trimmed_diff_all_sem_collect), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([0 4]);
ylim([-0.4 0.2])
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');






figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(diff_all_mean_prechoice), mean(diff_all_sem_prechoice), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, mean(diff_all_mean_postchoice), mean(diff_all_sem_postchoice), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, mean(diff_all_mean_collect), mean(diff_all_sem_collect), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-5 1]);
ylim([-0.55 0.2])
set(gca, 'YTick', [-.5 -.4 -.3 -.2 -.1 0 .1 .2]);
ytickformat('%.1f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');


figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
% h(1) = shadedErrorBar(ts1, mean(diff_all_mean_prechoice), mean(diff_all_sem_prechoice), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, mean(diff_all_mean_postchoice), mean(diff_all_sem_postchoice), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, mean(diff_all_mean_collect), mean(diff_all_sem_collect), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-1 3]);
ylim([-0.55 0.2])
set(gca, 'YTick', [-.5 -.4 -.3 -.2 -.1 0 .1 .2]);
ytickformat('%.1f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
% h(1) = shadedErrorBar(ts1, mean(diff_all_mean_prechoice), mean(diff_all_sem_prechoice), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, mean(diff_all_mean_postchoice), mean(diff_all_sem_postchoice), 'lineProps', {'color', 'r'});
h(3) = shadedErrorBar(ts1, mean(diff_all_mean_collect), mean(diff_all_sem_collect), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([0 4]);
ylim([-0.55 0.2])
set(gca, 'YTick', [-.5 -.4 -.3 -.2 -.1 0 .1 .2]);
ytickformat('%.1f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

% plot differences with error & bCI for statistics
mean_data_array = {diff_all_mean_prechoice}
sem_data_array = {diff_all_sem_prechoice}
[comparison] = bCI_and_tCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 1], [-0.6 0.2], 3)

mean_data_array = {diff_all_mean_postchoice}
sem_data_array = {diff_all_sem_postchoice}
[comparison] = bCI_and_tCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 3], [-0.8 0.2])

mean_data_array = {diff_all_mean_collect}
sem_data_array = {diff_all_sem_collect}
[comparison] = bCI_and_tCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [0 4], [-0.6 0.2])


%% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav
% this outputs a ever so slightly wonky diagram. a few nodes that do not
% actually overlap minimally overlap (but intersections are 0), and 1 node
% that has 1 overlap does not overlap at all. 
shk_ind = find(respClass_all_array{1,4} == 1);
% pre_choice_active_ind = find(respClass_all_array{1,1} == 1);
lost_ind = find(lost_all == 1);
conserved_ind = find(conserved_all == 1);
remapped_ind = find(remapped_all == 1);
% aa_active_ind = find(respClass_all_array{1,11} == 1);
% consum_inhibited_ind = find(all_consum_inhibited == 1);
setListData = {shk_ind, lost_ind, conserved_ind, remapped_ind};
setLabels = ["Shk excited", "Lost", "Conserved", "Remapped"];
figure;
ve_diagram = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

ve_diagram.ShowIntersectionCounts = true;
ve_diagram.ShowIntersectionAreas = true;
% h.SetLabels = [];

% shk_alone = respClass_all_array{1,4} == 1  & prechoice_block_1 ~=1 & postchoice_reward_block_1 ~= 1 & collect_block_1 ~= 1 & respClass_all_array{1,11} ~= 1;


%%
figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(lost_prechoice == 1, :)), mean(neuron_sem_array{1,4}(lost_prechoice == 1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(lost_postchoice == 1, :)), mean(neuron_sem_array{1,4}(lost_postchoice == 1, :)), 'lineProps', {'color', 'r'});
h(3) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(lost_consumption == 1, :)), mean(neuron_sem_array{1,4}(lost_consumption == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-2 4]);
ylim([-0.5 0.6])
set(gca, 'YTick', [-.50 -.25 0 .25 0.5]);
ytickformat('%.2f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');


%%
figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(conserved_prechoice == 1, :)), mean(neuron_sem_array{1,4}(conserved_prechoice == 1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(conserved_postchoice == 1, :)), mean(neuron_sem_array{1,4}(conserved_postchoice == 1, :)), 'lineProps', {'color', 'r'});
h(3) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(conserved_consumption == 1, :)), mean(neuron_sem_array{1,4}(conserved_consumption == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-2 4]);
ylim([-0.5 0.6])
set(gca, 'YTick', [-.50 -.25 0 .25 0.5]);
ytickformat('%.2f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

%%
figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(remapped_prechoice == 1, :)), mean(neuron_sem_array{1,4}(remapped_prechoice == 1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(remapped_postchoice == 1, :)), mean(neuron_sem_array{1,4}(remapped_postchoice == 1, :)), 'lineProps', {'color', 'r'});
h(3) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(remapped_consumption == 1, :)), mean(neuron_sem_array{1,4}(remapped_consumption == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-2 4]);
ylim([-0.5 0.6])
set(gca, 'YTick', [-.50 -.25 0 .25 0.5]);
ytickformat('%.2f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

%%
figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(conserved_prechoice == 1, :)), mean(neuron_sem_array{1,4}(conserved_prechoice == 1, :)), 'lineProps', {'color', 'g'});
h(2) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(remapped_prechoice == 1, :)), mean(neuron_sem_array{1,4}(remapped_prechoice == 1, :)), 'lineProps', {'color', 'b'});
h(3) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(lost_prechoice == 1, :)), mean(neuron_sem_array{1,4}(lost_prechoice == 1, :)), 'lineProps', {'color', 'r'});

% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-2 4]);
ylim([-0.5 0.6])
set(gca, 'YTick', [-.50 -.25 0 .25 0.5]);
ytickformat('%.2f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');


%%
figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(conserved_postchoice == 1, :)), mean(neuron_sem_array{1,4}(conserved_postchoice == 1, :)), 'lineProps', {'color', 'g'});
h(2) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(remapped_postchoice == 1, :)), mean(neuron_sem_array{1,4}(remapped_postchoice == 1, :)), 'lineProps', {'color', 'b'});
h(3) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(lost_postchoice == 1, :)), mean(neuron_sem_array{1,4}(lost_postchoice == 1, :)), 'lineProps', {'color', 'r'});

% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-2 4]);
ylim([-0.5 0.6])
set(gca, 'YTick', [-.50 -.25 0 .25 0.5]);
ytickformat('%.2f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

%%
figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(conserved_consumption== 1, :)), mean(neuron_sem_array{1,4}(conserved_consumption == 1, :)), 'lineProps', {'color', 'g'});
h(2) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(remapped_consumption == 1, :)), mean(neuron_sem_array{1,4}(remapped_consumption == 1, :)), 'lineProps', {'color', 'b'});
h(3) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(lost_consumption == 1, :)), mean(neuron_sem_array{1,4}(lost_consumption == 1, :)), 'lineProps', {'color', 'r'});

% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-2 4]);
ylim([-0.5 0.6])
set(gca, 'YTick', [-.50 -.25 0 .25 0.5]);
ytickformat('%.2f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

%%
figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(lost_all == 1, :)), mean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(remapped_all == 1, :)), mean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'b'});
h(3) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(conserved_all == 1, :)), mean(neuron_sem_array{1,4}(conserved_all == 1, :)), 'lineProps', {'color', 'g'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-2 4]);
ylim([-0.5 0.6])
set(gca, 'YTick', [-.50 -.25 0 .25 0.5]);
ytickformat('%.2f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

%%
%% create correlation matrix heatmap with exclusively active cells
test = [];
test = [neuron_mean_array{1, 8}(prechoice_block_1 == 1, :)];
% test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 1}~=1 & respClass_all_array{1, 2}==1 & respClass_all_array{1, 3}~=1, :)];
test = [test; neuron_mean_array{1, 9}(postchoice_reward_block_1 == 1, :)];
test = [test; neuron_mean_array{1, 10}(collect_block_1 == 1,:)];
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
colormap(bluewhitered);

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
prechoice_conserved_and_pun = conserved_prechoice == 1 & respClass_all_array{1, 4} == 1; 
sum(prechoice_conserved_and_pun)

proportion_prechoice_conserved_and_pun = sum(prechoice_conserved_and_pun)/sum(conserved_prechoice)

prechoice_new_and_pun = remapped_prechoice == 1 & respClass_all_array{1, 4} == 1; 
sum(prechoice_new_and_pun)

proportion_prechoice_new_and_pun = sum(prechoice_new_and_pun)/sum(remapped_prechoice)


n1 = sum(prechoice_conserved_and_pun); N1 = sum(conserved_prechoice);

n2 = sum(prechoice_new_and_pun); N2 = sum(remapped_prechoice);

% n2 = postchoice_and_shk_sum; N2 = postchoice_and_shk_sum+postchoice_not_shk_sum;

% Pooled estimate of proportion

p0 = (n1+n2) / (N1+N2)

% Expected counts under H0 (null hypothesis)

n10 = N1 * p0;

n20 = N2 * p0;

% Chi-square test, by hand

observed = [n1 N1-n1 n2 N2-n2];

expected = [n10 N1-n10 n20 N2-n20];

chi2stat = sum((observed-expected).^2 ./ expected)

p = 1 - chi2cdf(chi2stat,1)

%% "new" fig 4 analyses below here





%%
 mean_data_array = {neuron_mean_array{1, 1}(prechoice_block_1==1, :), neuron_mean_array{1, 8}(prechoice_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 1}(prechoice_block_1==1, :), neuron_sem_array{1, 8}(prechoice_block_1==1, :)};

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 1], [-0.7 0.6]);



 %%
 mean_data_array = {neuron_mean_array{1, 2}(postchoice_reward_block_1==1, :), neuron_mean_array{1, 9}(postchoice_reward_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 2}(postchoice_reward_block_1==1, :), neuron_sem_array{1, 9}(postchoice_reward_block_1==1, :)};

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 3], [-0.7 0.6]);

 %%
 mean_data_array = {neuron_mean_array{1, 3}(collect_block_1==1, :), neuron_mean_array{1, 10}(collect_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 3}(collect_block_1==1, :), neuron_sem_array{1, 10}(collect_block_1==1, :)};

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [0 4], [-0.7 0.6]);



 %% below figures plot remapped neurons (in comparison to their block 1 activity

%%
 mean_data_array = {neuron_mean_array{1, 1}(remapped_prechoice==1, :), neuron_mean_array{1, 8}(remapped_prechoice==1, :)};
 sem_data_array = {neuron_sem_array{1, 1}(remapped_prechoice==1, :), neuron_sem_array{1, 8}(remapped_prechoice==1, :)};

  % [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 1], [-0.2 0.6], 3);

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 1], [-0.05 0.4], 3);



 %%
 mean_data_array = {neuron_mean_array{1, 2}(remapped_postchoice==1, :), neuron_mean_array{1, 9}(remapped_postchoice==1, :)};
 sem_data_array = {neuron_sem_array{1, 2}(remapped_postchoice==1, :), neuron_sem_array{1, 9}(remapped_postchoice==1, :)};

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 3], [-0.1 0.7]);

 %%
 mean_data_array = {neuron_mean_array{1, 3}(remapped_consumption==1, :), neuron_mean_array{1, 10}(remapped_consumption==1, :)};
 sem_data_array = {neuron_sem_array{1, 3}(remapped_consumption==1, :), neuron_sem_array{1, 10}(remapped_consumption==1, :)};

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [0 4], [-0.1 0.7]);







%% create correlation matrix heatmap with exclusively active cells
test = [];
test = [neuron_mean_array{1, 1}(prechoice_blocks_2_and_3 == 1, :)];
% test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 1}~=1 & respClass_all_array{1, 2}==1 & respClass_all_array{1, 3}~=1, :)];
test = [test; neuron_mean_array{1, 2}(postchoice_reward_blocks_2_and_3 == 1, :)];
test = [test; neuron_mean_array{1, 3}(collect_blocks_2_and_3 == 1,:)];
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

xline(pre_choice_index(end), '--g')
xline(post_choice_index(end), '--c')
xline(consumption_index(end), '--b')























%% create correlation matrix heatmap with exclusively active cells
test = [];
test = [neuron_mean_array{1, 8}(prechoice_block_1 == 1, :)];
% test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 1}~=1 & respClass_all_array{1, 2}==1 & respClass_all_array{1, 3}~=1, :)];
test = [test; neuron_mean_array{1, 9}(postchoice_reward_block_1 == 1, :)];
test = [test; neuron_mean_array{1, 10}(collect_block_1 == 1,:)];
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

xline(pre_choice_index(end), '--g')
xline(post_choice_index(end), '--c')
xline(consumption_index(end), '--b')
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
xline(pre_choice_index(end), '--g')
xline(post_choice_index(end), '--c')
xline(consumption_index(end), '--b')

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



%%
lost_prechoice_become_neutral = lost_prechoice == 1 & (collect_blocks_2_and_3 == 0 & postchoice_reward_blocks_2_and_3 == 0)
lost_prechoice_become_postchoice = lost_prechoice == 1 & postchoice_reward_blocks_2_and_3 == 1
lost_prechoice_become_collect = lost_prechoice == 1 & collect_blocks_2_and_3 == 1;

lost_postchoice_become_neutral = lost_postchoice == 1 & (prechoice_blocks_2_and_3 == 0 & collect_blocks_2_and_3 == 0);
lost_postchoice_become_prechoice = lost_postchoice == 1 & prechoice_blocks_2_and_3 == 1;
lost_postchoice_become_collect = lost_postchoice == 1 & collect_blocks_2_and_3 == 1;

lost_collect_become_neutral = lost_consumption == 1 & (prechoice_blocks_2_and_3 == 0 & postchoice_reward_blocks_2_and_3 == 0);
lost_collect_become_prechoice = lost_consumption == 1 & prechoice_blocks_2_and_3 == 1;
lost_collect_become_postchoice = lost_consumption == 1 & postchoice_reward_blocks_2_and_3 == 1;




% Calculate sums for each category
lost_prechoice_total = sum(lost_prechoice);
lost_postchoice_total = sum(lost_postchoice);
lost_consumption_total = sum(lost_consumption);



% Calculate stacks for each category
lost_prechoice_stacks = [sum(lost_prechoice_become_neutral),...
                        0,... %for prechoice
                        sum(lost_prechoice_become_postchoice), ...
                        sum(lost_prechoice_become_collect)];

lost_postchoice_stacks = [sum(lost_postchoice_become_neutral), ...
                          sum(lost_postchoice_become_prechoice), ...
                          0,... % for postchoice
                          sum(lost_postchoice_become_collect)];

lost_consumption_stacks = [sum(lost_collect_become_neutral), ...
                           sum(lost_collect_become_prechoice), ...
                           sum(lost_collect_become_postchoice),...
                           0]; %for collect



% Combine data into a matrix for plotting
stacked_data = [lost_prechoice_stacks;...
                lost_postchoice_stacks;...
                lost_consumption_stacks];

% Labels for categories
categories = {'Lost Prechoice', 'Lost Postchoice', 'Lost Consumption'};
subcategories = {'Neutral', 'Prechoice', 'Post-choice', 'Collect'};

% Create stacked bar plot
figure;
bar(stacked_data, 'stacked');
colormap(parula); % Optional: change the color map
legend(subcategories, 'Location', 'northeastoutside');
xticks(1:length(categories));
xticklabels(categories);
ylabel('Counts');
title('What are Lost neurons now');

% Add data labels (optional)
hold on;
for i = 1:size(stacked_data, 1)
    total = sum(stacked_data(i, :));
    text(i, total + max(total)*0.02, num2str(total), 'HorizontalAlignment', 'center', 'FontSize', 10);
end
hold off;


%%
remapped_prechoice_was_neutral = remapped_prechoice == 1 & (collect_block_1 == 0 & postchoice_reward_block_1 == 0);
remapped_prechoice_was_postchoice = remapped_prechoice == 1 & postchoice_reward_block_1 == 1;
remapped_prechoice_was_collect = remapped_prechoice == 1 & collect_block_1 == 1;

remapped_postchoice_was_neutral = remapped_postchoice == 1 & (collect_block_1 == 0 & prechoice_block_1 == 0);
remapped_postchoice_was_prechoice = remapped_postchoice == 1 & prechoice_block_1 == 1;
remapped_postchoice_was_collect = remapped_postchoice == 1 & collect_block_1 == 1;

remapped_consumption_was_neutral = remapped_consumption == 1 & (postchoice_reward_block_1 == 0 & prechoice_block_1 == 0);
remapped_consumption_was_prechoice = remapped_consumption == 1 & prechoice_block_1 == 1;
remapped_consumption_was_postchoice = remapped_consumption == 1 & postchoice_reward_block_1 == 1;

% new_prechoice_was_neutral = new_prechoice == 1 & (collect_block_1 == 0 & postchoice_reward_block_1 == 0);
% new_prechoice_was_postchoice = new_prechoice == 1 & postchoice_reward_block_1 == 1;
% new_prechoice_was_collect = new_prechoice == 1 & collect_block_1 == 1;
% 
% new_postchoice_was_neutral = new_postchoice == 1 & (collect_block_1 == 0 & prechoice_block_1 == 0);
% new_postchoice_was_prechoice = new_postchoice == 1 & prechoice_block_1 == 1;
% new_postchoice_was_collect = new_postchoice == 1 & collect_block_1 == 1;
% 
% new_consumption_was_neutral = new_consumption == 1 & (postchoice_reward_block_1 == 0 & prechoice_block_1 == 0);
% new_consumption_was_prechoice = new_consumption == 1 & prechoice_block_1 == 1;
% new_consumption_was_postchoice = new_consumption == 1 & postchoice_reward_block_1 == 1;



% % Retrospective stacks for "lost"
% new_prechoice_retrospective_stack = [sum(new_prechoice_was_neutral), ...
%                                       sum(new_prechoice_was_postchoice), ...
%                                       sum(new_prechoice_was_collect)];
% 
% new_postchoice_retrospective_stack = [sum(new_postchoice_was_neutral), ...
%                                        sum(new_postchoice_was_prechoice), ...
%                                        sum(new_postchoice_was_collect)];
% 
% new_consumption_retrospective_stack = [sum(new_consumption_was_neutral), ...
%                                         sum(new_consumption_was_prechoice), ...
%                                         sum(new_consumption_was_postchoice)];

% Retrospective stacks for "remapped"
remapped_prechoice_retrospective_stack = [sum(remapped_prechoice_was_neutral),...
                                          0,... %for prechoice
                                          sum(remapped_prechoice_was_postchoice), ...
                                          sum(remapped_prechoice_was_collect)];

remapped_postchoice_retrospective_stack = [sum(remapped_postchoice_was_neutral), ...
                                           sum(remapped_postchoice_was_prechoice), ...
                                           0,... %for postchoice
                                           sum(remapped_postchoice_was_collect)];

remapped_consumption_retrospective_stack = [sum(remapped_consumption_was_neutral), ...
                                            sum(remapped_consumption_was_prechoice), ...
                                            sum(remapped_consumption_was_postchoice),...
                                            0]; %for collect

% Combine all stacks into a single matrix for plotting
stacked_data = [remapped_prechoice_retrospective_stack; ...
                remapped_postchoice_retrospective_stack; ...
                remapped_consumption_retrospective_stack];

% Update labels for categories
categories = {'Remapped Prechoice Retrospective', 'Remapped Postchoice Retrospective', 'Remapped Consumption Retrospective'};
subcategories = {'Neutral', 'Prechoice', 'Post-choice', 'Collect'};

% Create the stacked bar plot
figure;
bar(stacked_data, 'stacked');
colormap(parula); % Optional: change the color map
legend(subcategories, 'Location', 'northeastoutside');
xticks(1:length(categories));
xticklabels(categories);
ylabel('Counts');
title('What cats do neurons remap to?');

% Add data labels for the total counts above each bar (optional)
hold on;
for i = 1:size(stacked_data, 1)
    total = sum(stacked_data(i, :));
    text(i, total + max(max(stacked_data)) * 0.02, num2str(total), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end
hold off;

%%
pre_active_choice_to_shk = prechoice_block_1 == 1 & shk_activated == 1;
sum(pre_active_choice_to_shk)
post_choice_reward_active_to_shk = postchoice_reward_block_1 == 1 & shk_activated == 1;
sum(post_choice_reward_active_to_shk)
consumption_active_to_shk = collect_block_1 == 1 & shk_activated == 1;
sum(consumption_active_to_shk)

pre_and_post_active_to_shk = block_1_pre_and_post == 1 & shk_activated == 1;
sum(pre_and_post_active_to_shk)
post_and_consum_active_to_shk =block_1_post_and_consumption == 1 & shk_activated == 1;
sum(post_and_consum_active_to_shk)
pre_and_post_and_consum_active_to_shk = block_1_pre_and_post_and_consumption == 1 & shk_activated == 1;
sum(pre_and_post_and_consum_active_to_shk)

pre_to_pre_and_shk = conserved_prechoice == 1 & shk_activated == 1;
sum(pre_to_pre_and_shk)
pre_to_post_and_shk = remapped_postchoice_was_prechoice == 1 & shk_activated == 1;
sum(pre_to_post_and_shk)
pre_to_consum_and_shk = remapped_consumption_was_prechoice == 1 & shk_activated == 1;
sum(pre_to_consum_and_shk)

post_to_pre_and_shk = remapped_prechoice_was_postchoice == 1 & shk_activated == 1;
sum(post_to_pre_and_shk)
post_to_post_and_shk = conserved_postchoice == 1 & shk_activated == 1;
sum(post_to_post_and_shk)
post_to_consum_and_shk = remapped_consumption_was_postchoice == 1 & shk_activated == 1;
sum(post_to_consum_and_shk)

consum_to_pre_and_shk = remapped_prechoice_was_collect == 1 & shk_activated == 1;
sum(consum_to_pre_and_shk)
consum_to_post_and_shk = remapped_postchoice_was_collect == 1 & shk_activated == 1;
sum(consum_to_post_and_shk)
consum_to_consum_and_shk = conserved_consumption == 1 & shk_activated == 1;
sum(consum_to_consum_and_shk)






neutral_to_shk = true_neutral == 1 & shk_activated == 1;
sum(neutral_to_shk)

pre_inhibited_choice_to_shk = respClass_all_array{1,1} == 2 & shk_activated == 1;
sum(pre_inhibited_choice_to_shk)

post_choice_reward_inhibited_to_shk = respClass_all_array{1,2} == 2 & shk_activated == 1;
sum(post_choice_reward_inhibited_to_shk)

consumption_inhibited_to_shk = respClass_all_array{1,3} == 2 & shk_activated == 1;
sum(consumption_inhibited_to_shk)


pie_slices = [sum(pre_active_choice_to_shk) sum(post_choice_reward_active_to_shk) sum(consumption_active_to_shk) sum(neutral_to_shk)]
figure;
pie(pie_slices)


%%

% Calculate sums
pre_active_choice_to_shk_sum = sum(prechoice_block_1 == 1 & shk_activated == 1);
post_choice_reward_active_to_shk_sum = sum(postchoice_reward_block_1 == 1 & shk_activated == 1);
consumption_active_to_shk_sum = sum(collect_block_1 == 1 & shk_activated == 1);

pre_and_post_active_to_shk_sum = sum(block_1_pre_and_post == 1 & shk_activated == 1);
post_and_consum_active_to_shk_sum = sum(block_1_post_and_consumption == 1 & shk_activated == 1);
pre_and_post_and_consum_active_to_shk_sum = sum(block_1_pre_and_post_and_consumption == 1 & shk_activated == 1);

pre_to_pre_and_shk_sum = sum(conserved_prechoice == 1 & shk_activated == 1);
pre_to_pre_and_shk_sum_redo = sum(prechoice_block_1 == 1 & prechoice_blocks_2_and_3 == 1 & shk_activated == 1);
pre_to_post_and_shk_sum = sum(remapped_postchoice_was_prechoice == 1 & shk_activated == 1);
pre_to_consum_and_shk_sum = sum(remapped_consumption_was_prechoice == 1 & shk_activated == 1);

post_to_pre_and_shk_sum = sum(remapped_prechoice_was_postchoice == 1 & shk_activated == 1);
post_to_post_and_shk_sum = sum(conserved_postchoice == 1 & shk_activated == 1);
post_to_consum_and_shk_sum = sum(remapped_consumption_was_postchoice == 1 & shk_activated == 1);

consum_to_pre_and_shk_sum = sum(remapped_prechoice_was_collect == 1 & shk_activated == 1);
consum_to_post_and_shk_sum = sum(remapped_postchoice_was_collect == 1 & shk_activated == 1);
consum_to_consum_and_shk_sum = sum(conserved_consumption == 1 & shk_activated == 1);

prechoice_to_consum_and_shk = prechoice_block_1 ==1 & collect_blocks_2_and_3 == 1 & shk_activated == 1;
postchoice_to_consum_and_shk = postchoice_reward_block_1 ==1 & collect_blocks_2_and_3 == 1 & shk_activated == 1;
consum_to_consum_and_shk = collect_block_1 ==1 & collect_blocks_2_and_3 == 1 & shk_activated == 1;
neutral_to_consum_and_shk = neutral_to_shk == 1 & collect_blocks_2_and_3 == 1 & shk_activated == 1;
pre_and_consumption_to_consum_and_shk = block_1_pre_and_consumption == 1 & collect_blocks_2_and_3 == 1 & shk_activated == 1;
pre_and_post_to_consum_and_shk = block_1_pre_and_post == 1 & collect_blocks_2_and_3 == 1 & shk_activated == 1;
post_and_consumption_to_consum_and_shk = block_1_post_and_consumption == 1 & collect_blocks_2_and_3 == 1 & shk_activated == 1;



pie_slices = [
    sum(prechoice_to_consum_and_shk) 
    sum(postchoice_to_consum_and_shk) 
    sum(consum_to_consum_and_shk) 
    sum(neutral_to_consum_and_shk)
    sum(pre_and_consumption_to_consum_and_shk)
    sum(pre_and_post_to_consum_and_shk)
    sum(post_and_consumption_to_consum_and_shk)];
figure;
pie(pie_slices)



neutral_to_shk_sum = sum(neutral_to_shk);

% Combine into array and calculate percentages
values = [
    neutral_to_shk_sum;
    pre_active_choice_to_shk_sum;
    post_choice_reward_active_to_shk_sum;
    consumption_active_to_shk_sum;
    pre_and_post_active_to_shk_sum;
    post_and_consum_active_to_shk_sum;
    pre_and_post_and_consum_active_to_shk_sum;
    pre_to_pre_and_shk_sum;
    pre_to_post_and_shk_sum;
    pre_to_consum_and_shk_sum;
    post_to_pre_and_shk_sum;
    post_to_post_and_shk_sum;
    post_to_consum_and_shk_sum;
    consum_to_pre_and_shk_sum;
    consum_to_post_and_shk_sum;
    consum_to_consum_and_shk_sum
];

% Convert to percentages
percentages = (values / neuron_num) * 100;

% Define bar positions
group_gap = 0.5; % Gap size between groups
bar_positions = zeros(size(percentages));
bar_positions(1) = 1; % Neutral bar starts at position 1

% Assign positions for the rest of the bars
for i = 2:length(bar_positions)
    if mod(i - 2, 3) == 0 % Start of a new group (after the Neutral bar)
        bar_positions(i) = bar_positions(i - 1) + group_gap + 1; % Add group gap
    else
        bar_positions(i) = bar_positions(i - 1) + 1; % Normal spacing within group
    end
end

% Bar plot
figure;
bar(bar_positions, percentages, 'FaceColor', 'b');

% Adjust x-ticks and labels
xticks(bar_positions); % Set x-ticks to match bar positions
xticklabels({'Neutral', 'PreChoice', 'PostchoiceRew', 'Consumption', ...
    'Pre & Post', 'Post & Consum', 'Pre, Post & Consum', ...
    'Pre to Pre', 'Pre to Post', 'Pre to Consum', ...
    'Post to Pre', 'Post to Post', 'Post to Consum', ...
    'Consum to Pre', 'Consum to Post', 'Consum to Consum'}); % Labels for bars
xtickangle(45); % Rotate x-axis labels for better readability

% Add labels and title
% xlabel('Conditions');
ylabel('Percentage (%)');
title('Percent of neurons activated by Pun for each category');
% grid on;

%% for below
load('BLA_RDT_6_variables_only_large_use_to_show_large_alone_decrease.mat')
% and run blockwise_changes

%%
 mean_data_array = {neuron_mean_array{1, 1}(prechoice_block_1==1, :), neuron_mean_array{1, 4}(prechoice_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 1}(prechoice_block_1==1, :), neuron_sem_array{1, 4}(prechoice_block_1==1, :)};

 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 1], [-0.5 0.5], 3);



 %%
 mean_data_array = {neuron_mean_array{1, 2}(postchoice_reward_block_1==1, :), neuron_mean_array{1, 5}(postchoice_reward_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 2}(postchoice_reward_block_1==1, :), neuron_sem_array{1, 5}(postchoice_reward_block_1==1, :)};

 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 1], [-0.5 .5])

 %%
 mean_data_array = {neuron_mean_array{1, 3}(collect_block_1==1, :), neuron_mean_array{1, 6}(collect_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 3}(collect_block_1==1, :), neuron_sem_array{1, 6}(collect_block_1==1, :)};

 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 1], [-0.5 1]);


 %% load('BLA_RDT_6_variables_only_large_use_to_show_large_alone_decrease.mat')
 % use below to show that decreases in activity are not driven solely by
 % more small reward trials!
 %%

diff_all_mean_prechoice = neuron_mean_array{1, 4}(prechoice_block_1 == 1, :) - neuron_mean_array{1, 1}(prechoice_block_1 == 1, :)
diff_all_sem_prechoice = sqrt(neuron_sem_array{1, 4}(prechoice_block_1 == 1, :).^2 + neuron_sem_array{1, 1}(prechoice_block_1 == 1, :).^2)



diff_all_mean_postchoice = neuron_mean_array{1, 5}(postchoice_reward_block_1 == 1, :) - neuron_mean_array{1, 2}(postchoice_reward_block_1 == 1, :)
diff_all_sem_postchoice = sqrt(neuron_sem_array{1, 5}(postchoice_reward_block_1 == 1, :).^2 + neuron_sem_array{1, 2}(postchoice_reward_block_1 == 1, :).^2)


diff_all_mean_collect = neuron_mean_array{1, 6}(collect_block_1 == 1, :) - neuron_mean_array{1, 3}(collect_block_1 == 1, :)
diff_all_sem_collect = sqrt(neuron_sem_array{1, 6}(collect_block_1 == 1, :).^2 + neuron_sem_array{1, 3}(collect_block_1 == 1, :).^2)


trimmed_diff_all_mean_prechoice = diff_all_mean_prechoice(:, ts1 > -4 & ts1 <= 0);
trimmed_diff_all_sem_prechoice = diff_all_sem_prechoice(:, ts1 > -4 & ts1 <= 0);

trimmed_diff_all_mean_postchoice = diff_all_mean_postchoice(:, ts1 > 0 & ts1 <= 4);
trimmed_diff_all_sem_postchoice = diff_all_sem_postchoice(:, ts1 > 0 & ts1 <= 4);


trimmed_diff_all_mean_collect = diff_all_mean_collect(:, ts1 > 0 & ts1 <= 4);
trimmed_diff_all_sem_collect = diff_all_sem_collect(:, ts1 > 0 & ts1 <= 4);

trimmed_ts1 = ts1(:, ts1 > 0 & ts1 <= 4);

figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(trimmed_ts1 , mean(trimmed_diff_all_mean_prechoice), mean(trimmed_diff_all_sem_prechoice), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(trimmed_ts1 , mean(trimmed_diff_all_mean_postchoice), mean(trimmed_diff_all_sem_postchoice), 'lineProps', {'color', 'r'});
h(3) = shadedErrorBar(trimmed_ts1 , mean(trimmed_diff_all_mean_collect), mean(trimmed_diff_all_sem_collect), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([0 4]);
ylim([-0.4 0.2])
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');






figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(diff_all_mean_prechoice), mean(diff_all_sem_prechoice), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, mean(diff_all_mean_postchoice), mean(diff_all_sem_postchoice), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, mean(diff_all_mean_collect), mean(diff_all_sem_collect), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-5 1]);
ylim([-0.55 0.2])
set(gca, 'YTick', [-.5 -.4 -.3 -.2 -.1 0 .1 .2]);
ytickformat('%.1f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');


figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
% h(1) = shadedErrorBar(ts1, mean(diff_all_mean_prechoice), mean(diff_all_sem_prechoice), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, mean(diff_all_mean_postchoice), mean(diff_all_sem_postchoice), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, mean(diff_all_mean_collect), mean(diff_all_sem_collect), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-1 3]);
ylim([-0.55 0.2])
set(gca, 'YTick', [-.5 -.4 -.3 -.2 -.1 0 .1 .2]);
ytickformat('%.1f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
% h(1) = shadedErrorBar(ts1, mean(diff_all_mean_prechoice), mean(diff_all_sem_prechoice), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, mean(diff_all_mean_postchoice), mean(diff_all_sem_postchoice), 'lineProps', {'color', 'r'});
h(3) = shadedErrorBar(ts1, mean(diff_all_mean_collect), mean(diff_all_sem_collect), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([0 4]);
ylim([-0.55 0.2])
set(gca, 'YTick', [-.5 -.4 -.3 -.2 -.1 0 .1 .2]);
ytickformat('%.1f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

% plot differences with error & bCI for statistics
mean_data_array = {diff_all_mean_prechoice}
sem_data_array = {diff_all_sem_prechoice}
[comparison] = bCI_and_tCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 1], [-0.5 0.2])

mean_data_array = {diff_all_mean_postchoice}
sem_data_array = {diff_all_sem_postchoice}
[comparison] = bCI_and_tCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 3], [-0.5 0.2])

mean_data_array = {diff_all_mean_collect}
sem_data_array = {diff_all_sem_collect}
[comparison] = bCI_and_tCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [0 4], [-0.5 0.2])

%% below for supplement 4 correlate pre-choice blocks 2/3 with choice
%% ONLY USE CODE BELOW IF YOU WANT TO COMBINE DATA ACROSS 2 FILTERS, E.G., BLOCK 1 & BLOCKS 2/3

array_for_means = 8; 
array_for_means_second = []; % 5


for q = 1:length (behav_tbl_iter{1, 1})

    behav_part_1 = behav_tbl_iter{array_for_means, 1}{q, 1};

    if isempty(array_for_means_second)
        nestedCellArray_1 = behav_part_1;
    else
        behav_part_2 = behav_tbl_iter{array_for_means_second, 1}{q, 1};
        nestedCellArray_1 = vertcat(behav_part_1, behav_part_2 );
    end




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
trial_choice_times_concat = trial_choice_times_concat(trial_types_concat == 1.2);
trial_types_concat = trial_types_concat(trial_types_concat == 1.2);
path_length_concat = path_length_concat(trial_types_concat == 1.2);
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

    if ~isempty(array_for_means_second)
        nestedCellArray_2 = zall_mouse{i, array_for_means_second};
        for qq = 1:size(nestedCellArray_1, 2)
            concatenated_columns{qq} = vertcat(nestedCellArray_1{qq}, ...
                nestedCellArray_2{qq});
        end
        nestedCellArray_1 = concatenated_columns;

    end
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
        % meanValues = currentArray(:, timeRange);

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





% only_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(prechoice_block_1 == 1 & postchoice_reward_block_1 ~=1);
% not_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(prechoice_block_1 ~= 1 & postchoice_reward_block_1 ~=1);
% 
% only_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(prechoice_block_1 == 1 & postchoice_reward_block_1 ~=1);
% not_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(prechoice_block_1 ~= 1 & postchoice_reward_block_1 ~=1);


% 
% only_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(prechoice_block_1 ~= 1 & postchoice_reward_block_1 ==1);
% not_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(prechoice_block_1 ~= 1 & postchoice_reward_block_1 ~=1);
% 
% only_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(prechoice_block_1 ~= 1 & postchoice_reward_block_1 ==1);
% not_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(prechoice_block_1 ~= 1 & postchoice_reward_block_1 ~=1);


only_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(collect_blocks_2_and_3 == 1);
not_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(true_neutral == 1);

only_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(collect_blocks_2_and_3 == 1);
not_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(true_neutral == 1);

% only_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(postchoice_reward_block_1 == 1);
% not_shk_responsive_corrs_Type1_2 = allCorrelations_Type1_2(postchoice_reward_block_1 ~= 1);
% 
% only_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(postchoice_reward_block_1 == 1);
% not_shk_responsive_corrs_Type0_3 = allCorrelations_Type0_3(postchoice_reward_block_1 ~= 1);



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


bar_separation_value = 3;

figure;
width = 200; % Width of the figure
height = 200; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
swarmchart(ones(1, length(only_shk_responsive_corrs_Type1_2)), only_shk_responsive_corrs_Type1_2)
hold on
swarmchart(ones(1, length(not_shk_responsive_corrs_Type1_2))*bar_separation_value, not_shk_responsive_corrs_Type1_2)

% yline(mean(only_shk_responsive_corrs), ones(length(only_shk_responsive_corrs)))
plot([0.5; 1.5], [mean(only_shk_responsive_corrs_Type1_2); mean(only_shk_responsive_corrs_Type1_2)], 'LineWidth',3)
plot([bar_separation_value-.5; bar_separation_value+.5], [mean(not_shk_responsive_corrs_Type1_2); mean(not_shk_responsive_corrs_Type1_2)], 'LineWidth',3)
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off


% Create a histogram for allCorrelations
figure;
width = 250; % Width of the figure
height = 250; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
histogram(not_shk_responsive_corrs_Type1_2 , 'Normalization', 'probability', 'FaceColor', 'blue','BinWidth', 0.05,'LineStyle','none');
hold on;

% Create a histogram for only_shk_responsive_corrs on the same figure
histogram(only_shk_responsive_corrs_Type1_2, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05, 'LineStyle','none');
% xline(mean_only_shk, 'r')
% xline(mean_not_shk, 'g')
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
