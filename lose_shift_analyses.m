% lose_shift analyses - lose_shift across risky blocks vs. non-loseshift
% small choices
% choiceTime.Outcome_0to2.LOSESHIFT_1_BLOCK_2_BLOCK_3
% choiceTime.Outcome_0to2.REW_Small_LOSESHIFT_0_BLOCK_2_BLOCK_3

lose_shift_ensemble_data = neuron_mean_array{1, 1}(respClass_all_array{1,1}==1,:);
lose_shift_ensemble_sem = neuron_sem_array{1, 1}(respClass_all_array{1,1}==1,:);
non_lose_shift_small_choice_data = neuron_mean_array{1, 2}(respClass_all_array{1,2}==1,:);
non_lose_shift_small_choice_sem = neuron_sem_array{1, 2}(respClass_all_array{1,2}==1,:);

[concatenatedTable_all, concatenate_all_tables] = get_median_choice_and_collect_fn(behav_tbl_iter);

median_start_time_lose_shifts = median(concatenatedTable_all{1, 1}.stTime - concatenatedTable_all{1, 1}.choiceTime);
median_start_time_non_lose_shift_small = median(concatenatedTable_all{1, 2}.stTime - concatenatedTable_all{1, 2}.choiceTime);

median_start_times = [median_start_time_lose_shifts, median_start_time_non_lose_shift_small];

median_collect_time_lose_shifts = median(concatenatedTable_all{1, 1}.collectionTime - concatenatedTable_all{1, 1}.choiceTime);
median_collect_time_non_lose_shift_small = median(concatenatedTable_all{1, 2}.collectionTime - concatenatedTable_all{1, 2}.choiceTime);

median_collect_times = [median_collect_time_lose_shifts, median_collect_time_non_lose_shift_small];

figure;
shadedErrorBar(ts1, nanmean(lose_shift_ensemble_data), nanmean(lose_shift_ensemble_sem), 'lineProps', {'color', 'b'});
hold on;shadedErrorBar(ts1, nanmean(non_lose_shift_small_choice_data), nanmean(non_lose_shift_small_choice_sem), 'lineProps', {'color','g'});
% xtickformat('%.2f');
ytickformat('%.2f');
xline(0);
xline(median_start_times, 'g')
xline(median_collect_times, 'r')
% xlabel('Time from Large Rew Choice (s)');
% legend({'large_pre_choice_ensemble_block_1', 'consumption active', 'neutral'}, 'Location','northwest')