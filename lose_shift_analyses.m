% lose_shift analyses - lose_shift across risky blocks vs. non-loseshift
% small choices
% choiceTime.Outcome_0to2.LOSESHIFT_1_BLOCK_2_BLOCK_3
% choiceTime.Outcome_0to2.REW_Small_LOSESHIFT_0_BLOCK_2_BLOCK_3
% can also add lose_stay compoarison, which doesn't differ from loseshift
% with 'LOSESTAY', 1, 'BLOCK', 2, 'BLOCK', 3, 'SHK', 0

% CAN ALSO DO PRECHOICE
% choiceTime.prechoice_Minus_4to0.LOSESHIFT_1_REW_Small_BLOCK_2_BLOCK_3_TYPE_0
% choiceTime.prechoice_Minus_4to0.LOSESHIFT_0_REW_Small_BLOCK_2_BLOCK_3_TYPE_0

lose_shift_ensemble_data = neuron_mean_array{1, 1}(respClass_all_array{1,1}==1,:);
lose_shift_ensemble_sem = neuron_sem_array{1, 1}(respClass_all_array{1,1}==1,:);

non_lose_shift_small_choice_data = neuron_mean_array{1, 1}(respClass_all_array{1,2}==1,:);
non_lose_shift_small_choice_sem = neuron_sem_array{1, 1}(respClass_all_array{1,2}==1,:);

% lose_stay_ensemble_data = neuron_mean_array{1, 3}(respClass_all_array{1,3}==1,:);
% lose_stay_ensemble_sem = neuron_sem_array{1, 3}(respClass_all_array{1,3}==1,:);

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
% hold on;shadedErrorBar(ts1, nanmean(lose_stay_ensemble_data), nanmean(lose_stay_ensemble_sem), 'lineProps', {'color','g'});
% xtickformat('%.2f');
ytickformat('%.2f');
xline(0);
% xline(median_start_times, 'g')
% xline(median_collect_times, 'r')
% xlabel('Time from Large Rew Choice (s)');
% legend({'large_pre_choice_ensemble_block_1', 'consumption active', 'neutral'}, 'Location','northwest')

mean_data_array = {lose_shift_ensemble_data, non_lose_shift_small_choice_data};
sem_data_array = {lose_shift_ensemble_sem, non_lose_shift_small_choice_sem};
prechoice_x_limits = [-2 6];

[comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, prechoice_x_limits)


%%

inhib_or_excite = 1;


for kk = 1:size(animalIDs, 1)
    % prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
    loseshift_mouse{kk, :} = respClass_all_array_mouse{kk, 1} == inhib_or_excite;
    non_loseshift_mouse{kk, :} = respClass_all_array_mouse{kk, 2} == inhib_or_excite; 

  

    num_cells_mouse(kk) = size(respClass_all_array_mouse{kk, 1}, 2);
   
end


for kk = 1:size(animalIDs, 1)
    
    
    loseshift_sum(kk) = sum([sum(loseshift_mouse{kk, :})]);
    non_loseshift_sum(kk) = sum([sum(non_loseshift_mouse{kk, :})]);

end


%%
lose_shift_activated_num = 23+48;
other_cells = neuron_num - lose_shift_activated_num;

lose_shift_big_pie = [lose_shift_activated_num other_cells];

figure;
piechart(lose_shift_big_pie)
figure;
lose_shift_big_donutchart_data = [23 48 other_cells]
donutchart(lose_shift_big_donutchart_data)