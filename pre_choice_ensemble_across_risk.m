
%use base_workspace_BLA_RDT_3_categories & then tack on the additional
%comparisons, which will add to neuron_mean_array and respClass_all in slots 4, 5, etc. 


for zz = 1:size(concatenatedTable_all, 2)
    temp_table = concatenatedTable_all{zz}; 
    % median_choice_time_all(zz) = median(temp_table.choiceTime - temp_table.stTime);
    median_collect_time_all(zz) = median(temp_table.collectionTime - temp_table.stTime);
    median_start_time_all(zz) = median(temp_table.stTime - temp_table.choiceTime);
    % median_choice_time_block_2_all(zz) = median(temp_table.choiceTime(temp_table.Block == 2) - temp_table.stTime(temp_table.Block == 2));
    % median_choice_time_block_3_all(zz) = median(temp_table.choiceTime(temp_table.Block == 3) - temp_table.stTime(temp_table.Block == 3));
    clear temp_table
end

figure;
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(exclusive_activated_session_1 == 1, :)), nanmean(neuron_sem_array{1, 1}(exclusive_activated_session_1 == 1, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(exclusive_activated_session_1 == 1, :)), nanmean(neuron_sem_array{1, 4}(exclusive_activated_session_1 == 1, :)), 'lineProps', {'color', 'g'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 5}(exclusive_activated_session_1 == 1, :)), nanmean(neuron_sem_array{1, 5}(exclusive_activated_session_1 == 1, :)), 'lineProps', {'color', 'b'});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_all(1, 3:5), 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active no risk', 'pre-choice active b2', 'pre-choice active b3'}, 'Location','northwest')

%% FOR COLLECTION RESPONSIVE NEURONS
for zz = 1:size(concatenatedTable_all, 2)
    temp_table = concatenatedTable_all{zz}; 
    % median_choice_time_all(zz) = median(temp_table.choiceTime - temp_table.stTime);
    median_collect_time_all(zz) = median(temp_table.collectionTime - temp_table.stTime);
    median_start_time_all(zz) = median(temp_table.stTime - temp_table.choiceTime);
    % median_choice_time_block_2_all(zz) = median(temp_table.choiceTime(temp_table.Block == 2) - temp_table.stTime(temp_table.Block == 2));
    % median_choice_time_block_3_all(zz) = median(temp_table.choiceTime(temp_table.Block == 3) - temp_table.stTime(temp_table.Block == 3));
    clear temp_table
end

figure;
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(exclusive_activated_session_3 == 1, :)), nanmean(neuron_sem_array{1, 3}(exclusive_activated_session_3 == 1, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(exclusive_activated_session_3 == 1, :)), nanmean(neuron_sem_array{1, 4}(exclusive_activated_session_3 == 1, :)), 'lineProps', {'color', 'g'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 5}(exclusive_activated_session_3 == 1, :)), nanmean(neuron_sem_array{1, 5}(exclusive_activated_session_3 == 1, :)), 'lineProps', {'color', 'b'});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_all(1, 3:5), 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active no risk', 'pre-choice active b2', 'pre-choice active b3'}, 'Location','northwest')

