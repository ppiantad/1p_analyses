


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
