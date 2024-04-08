pre_choice_cells = respClass_mouse.BLA_Insc_01.RDT_D1.choiceTime.Outcome_Minus_4to0.OMITALL_0_BLANK_TOUCH_0_BLOCK_1 == 1
pre_choice_data = zall_mouse{1, 1}(pre_choice_cells);

figure; plot(ts1, )

pre_choice_AUC_idx =  ts1 >= -4 & ts1 <= 0;

for qq = 1:size(pre_choice_data, 2)
    data = pre_choice_data{qq}

    pre_choice_mean(qq,:) = mean(data(:, pre_choice_AUC_idx))
end


large_pre_choice_ensemble_block_1_zall_auc = mean(large_pre_choice_ensemble_block_1_zall(:, pre_choice_AUC_idx))
large_pre_choice_ensemble_block_1_zall_auc_mean = mean(large_pre_choice_ensemble_block_1_zall_auc);
large_pre_choice_ensemble_block_1_zall_auc_sem = mean(nanstd(large_pre_choice_ensemble_block_1_zall(:, pre_choice_AUC_idx))/(sqrt(size(large_pre_choice_ensemble_block_1_zall_auc, 2))));

large_pre_choice_ensemble_block_2_zall_auc = mean(large_pre_choice_ensemble_block_2_zall(:, pre_choice_AUC_idx))
large_pre_choice_ensemble_block_2_zall_auc_mean = mean(large_pre_choice_ensemble_block_2_zall_auc);
large_pre_choice_ensemble_block_2_zall_auc_sem = mean(nanstd(large_pre_choice_ensemble_block_2_zall(:, pre_choice_AUC_idx))/(sqrt(size(large_pre_choice_ensemble_block_2_zall_auc, 2))));


large_pre_choice_ensemble_block_3_zall_auc = mean(large_pre_choice_ensemble_block_3_zall(:, pre_choice_AUC_idx))
large_pre_choice_ensemble_block_3_zall_auc_mean = mean(large_pre_choice_ensemble_block_3_zall_auc)
large_pre_choice_ensemble_block_3_zall_auc_sem = mean(nanstd(large_pre_choice_ensemble_block_3_zall(:, pre_choice_AUC_idx))/(sqrt(size(large_pre_choice_ensemble_block_3_zall_auc, 2))));

auc_data = [large_pre_choice_ensemble_block_1_zall_auc_mean large_pre_choice_ensemble_block_2_zall_auc_mean large_pre_choice_ensemble_block_3_zall_auc_mean]
auc_error_high = [large_pre_choice_ensemble_block_1_zall_auc_sem large_pre_choice_ensemble_block_2_zall_auc_sem large_pre_choice_ensemble_block_3_zall_auc_sem]
auc_error_low = -1* [large_pre_choice_ensemble_block_1_zall_auc_sem large_pre_choice_ensemble_block_2_zall_auc_sem large_pre_choice_ensemble_block_3_zall_auc_sem];

figure;
x = 1:3
bar(x,auc_data)                

hold on

er = errorbar(x,auc_data,auc_error_low,auc_error_high);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off