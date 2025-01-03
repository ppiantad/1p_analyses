% load('BLA_RM_D1_3_categories.mat')
load('BLA_Pre_RDT_RM_3_categories.mat')
% run block_wise_changes_v1.m relevant sections first

%%



session_to_analyze = 'Pre_RDT_RM'; % set to whatever session to examine

fieldnames(final_behavior);
for gg = 1:size(fieldnames(final_behavior), 1)
    currentanimal = char(animalIDs(gg));
    session_string = session_to_analyze;
    behav_table = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
    all_consum_duration{gg} = behav_table.collectionTime_end(behav_table.bigSmall == 1.2 | behav_table.bigSmall == 0.3) - behav_table.collectionTime(behav_table.bigSmall == 1.2 | behav_table.bigSmall == 0.3); 
    large_consum_duration{gg} = behav_table.collectionTime_end(behav_table.bigSmall == 1.2) - behav_table.collectionTime(behav_table.bigSmall == 1.2);
    small_consum_duration{gg} = behav_table.collectionTime_end(behav_table.bigSmall == 0.3) - behav_table.collectionTime(behav_table.bigSmall == 0.3);
end


consumption_neuron_count = 0;
for qq = 1:size(zall_mouse, 1)
    zall_mouse_row = zall_mouse{qq, 3};
    all_consum_duration_current = all_consum_duration{qq};
    large_consum_duration_current = large_consum_duration{qq};
    small_consum_duration_current = small_consum_duration{qq};
    unnormalized_mouse_row = caTraceTrials_mouse{qq, 1};
    for ff = 1:size(zall_mouse_row, 2)
        zall_mouse_current = zall_mouse_row{1, ff};
        
        zall_mouse_current_large = zall_mouse_current(behav_tbl_iter{1, 1}{qq, 1}.bigSmall == 1.2, :);
        zall_mouse_current_small = zall_mouse_current(behav_tbl_iter{1, 1}{qq, 1}.bigSmall == 0.3, :);
        zall_mouse_current_all_subwindow =  mean(zall_mouse_current(:, ts1 >= 0 & ts1 <= 5), 2);
        zall_mouse_current_large_subwindow = mean(zall_mouse_current_large(:, ts1 >= 0 & ts1 <= 5), 2);
        zall_mouse_current_small_subwindow = mean(zall_mouse_current_small(:, ts1 >= 0 & ts1 <= 5), 2);
        % unnormalized_mouse_current = unnormalized_mouse_row{1, ff};
        % unnormalized_mouse_current = unnormalized_mouse_current(behav_tbl_iter{1, 1}{qq, 1}.bigSmall == 1.2, :);
        if collect_block_1_mouse{qq, 1}(ff) == 1
            consumption_neuron_count = consumption_neuron_count + 1;
            collect_array_all_correlation(consumption_neuron_count) = corr(zall_mouse_current_all_subwindow, all_consum_duration_current);
            collect_array_large_correlation(consumption_neuron_count) = corr(zall_mouse_current_large_subwindow, large_consum_duration_current);

            collect_array_small_correlation(consumption_neuron_count) = corr(zall_mouse_current_small_subwindow, small_consum_duration_current);

            % prechoice_unnormalized_array(prechoice_neuron_count, :) = unnormalized_mouse_current(1, :);
        end
    end
end
