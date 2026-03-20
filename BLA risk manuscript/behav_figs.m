animalIDs = (fieldnames(final_behavior));

session_to_analyze = 'RDT_D1'

if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze) | strcmp('PR_D1', session_to_analyze)
    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39', 'BLA_Insc_13'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_behavior, fieldsToRemove{i})
            final_behavior = rmfield(final_behavior, fieldsToRemove{i});
        end
    end
elseif strcmp('RDT_D2', session_to_analyze)

    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_behavior, fieldsToRemove{i})
            final_behavior = rmfield(final_behavior, fieldsToRemove{i});
        end
    end
end



%%


animalIDs = (fieldnames(final_behavior));


risk_table = table;


if exist('hM4Di_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = hM4Di_IDs(valid_mice);

elseif exist('stGtACR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = stGtACR_IDs(valid_mice);
    valid_animalIDs = valid_animalIDs(~strcmp(valid_animalIDs, 'RRD391'));


elseif exist('PdCO_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = PdCO_IDs(valid_mice);

elseif exist('ChrimsonR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = ChrimsonR_IDs(valid_mice);

elseif exist('treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = mouse_IDs(valid_mice);

else
    valid_animalIDs = animalIDs;

end



for ii = 1:size(valid_animalIDs,1)
    currentanimal = char(valid_animalIDs(ii));
    
    if isfield(final_behavior.(currentanimal), session_to_analyze)



        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        for BehavDataRow = 1:size(BehavData,1)
            if BehavData.shock(BehavDataRow) == 1
                kk = 1;
                while true
                    if (BehavDataRow + kk) > size(BehavData, 1)
                        break;
                    end
                    if ~isnan(BehavData.bigSmall(BehavDataRow + kk)) & BehavData.ForceFree(BehavDataRow + kk) ~= 999
                        BehavData.trial_after_shk(BehavDataRow + kk) = 1;
                        break;
                    else
                        kk = kk + 1;
                    end
                end
            end

        end

        for BehavDataRow = 1:size(BehavData,1)
            if BehavData.bigSmall(BehavDataRow) ~= 999 & ~isnan(BehavData.bigSmall(BehavDataRow))
                kk = 1;
                while true
                    if (BehavDataRow + kk) > size(BehavData, 1)
                        break;
                    end
                    if ~isnan(BehavData.bigSmall(BehavDataRow + kk)) & BehavData.ForceFree(BehavDataRow + kk) ~= 999
                        BehavData.initiation_delay(BehavDataRow + kk) = BehavData.stTime(BehavDataRow + kk)- BehavData.collectionTime(BehavDataRow);
                        break;
                    else
                        BehavData.initiation_delay(BehavDataRow + kk) = nan;
                        kk = kk + 1;
                    end
                end
            else
                BehavData.initiation_delay(BehavDataRow) = nan;
            end


        end
        
        [BehavData,trials,varargin]=TrialFilter_test(BehavData,'ALL', 1); 
        large_small_trials_only = BehavData(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3, :);



        desired_rows = 90;

        large_trials_true = large_small_trials_only.bigSmall == 1.2;
        small_trials_true = large_small_trials_only.bigSmall == 0.3;

        if size(large_trials_true, 1) < desired_rows
            padding = zeros(desired_rows - size(large_trials_true, 1), size(large_trials_true, 2));
            large_trials_true = [large_trials_true; padding];
        end

        if size(small_trials_true, 1) < desired_rows
            padding = zeros(desired_rows - size(small_trials_true, 1), size(small_trials_true, 2));
            small_trials_true = [small_trials_true; padding];
        end
        
        large_sequences_mouse(ii, :) = large_trials_true;
        small_sequences_mouse(ii, :) = small_trials_true;

        block_1_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 1 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 1); 
        block_1_small_choice_percent = sum(BehavData.bigSmall == 0.3 & BehavData.Block == 1 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 1); 
        block_2_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 2 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 2);  
        block_2_small_choice_percent = sum(BehavData.bigSmall == 0.3 & BehavData.Block == 2 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 2);  
        block_3_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 3 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);  
        block_3_small_choice_percent = sum(BehavData.bigSmall == 0.3 & BehavData.Block == 3 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);  
        block_1_shock_total = sum(BehavData.shock == 1 & BehavData.Block == 1);
        block_2_shock_total = sum(BehavData.shock == 1 & BehavData.Block == 2);
        block_3_shock_total = sum(BehavData.shock == 1 & BehavData.Block == 3);

        block_1_omission_total = sum(BehavData.omissionALL == 1 & BehavData.Block == 1);
        block_2_omission_total = sum(BehavData.omissionALL == 1 & BehavData.Block == 2);
        block_3_omission_total = sum(BehavData.omissionALL == 1 & BehavData.Block == 3);
        
        lose_shift_percent = sum(BehavData.lose_shift == 1)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3 | BehavData.Block == 2);
        lose_omit_percent = sum(BehavData.lose_omit == 1)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3 | BehavData.Block == 2);
        lose_stay_percent = sum(BehavData.lose_stay == 1)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3 | BehavData.Block == 2);
        win_stay_percent = sum(BehavData.win_stay == 1)/sum(((BehavData.bigSmall == 1.2) & BehavData.ForceFree == 0) & BehavData.Block == 3 | BehavData.Block == 2);
        lose_shift_ratio = sum(BehavData.lose_shift == 1)/sum(BehavData.WL == 3 & BehavData.ForceFree == 0);
        lose_stay_ratio = sum(BehavData.lose_stay == 1)/sum(BehavData.WL == 3 & BehavData.ForceFree == 0);
        lose_omit_ratio = sum(BehavData.lose_omit == 1)/sum(BehavData.WL == 3 & BehavData.ForceFree == 0);
        win_stay_ratio = sum(BehavData.win_stay == 1)/sum(BehavData.WL == 1 & BehavData.ForceFree == 0 & BehavData.Block == 2 | BehavData.Block == 3);

        win_stay_ratio_block_1 = sum(BehavData.win_stay == 1 & BehavData.Block == 1)/sum(BehavData.WL == 1 & BehavData.ForceFree == 0 & BehavData.Block == 1);
        win_stay_ratio_block_2 = sum(BehavData.win_stay == 1 & BehavData.Block == 2)/sum(BehavData.WL == 1 & BehavData.ForceFree == 0 & BehavData.Block == 2);
        win_stay_ratio_block_3 = sum(BehavData.win_stay == 1 & BehavData.Block == 3)/sum(BehavData.WL == 1 & BehavData.ForceFree == 0 & BehavData.Block == 3);
        
        lose_shift_ratio_block_1 = sum(BehavData.lose_shift == 1 & BehavData.Block == 1)/sum(BehavData.WL == 3 & BehavData.ForceFree == 0 & BehavData.Block == 1);
        lose_shift_ratio_block_2 = sum(BehavData.lose_shift == 1 & BehavData.Block == 2)/sum(BehavData.WL == 3 & BehavData.ForceFree == 0 & BehavData.Block == 2);
        lose_shift_ratio_block_3 = sum(BehavData.lose_shift == 1 & BehavData.Block == 3)/sum(BehavData.WL == 3 & BehavData.ForceFree == 0 & BehavData.Block == 3);
        BehavData.choice_latency = BehavData.choiceTime - BehavData.stTime;
        BehavData.collect_latency = BehavData.collectionTime - BehavData.choiceTime;
        if ismember('collectionTime_end', BehavData.Properties.VariableNames)
            BehavData.consum_duration = BehavData.collectionTime_end - BehavData.collectionTime;
        end
        block_1_choice_latency_all = mean(BehavData.choice_latency(BehavData.Block == 1 & (BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3)));
        block_2_choice_latency_all = mean(BehavData.choice_latency(BehavData.Block == 2 & (BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3)));
        block_3_choice_latency_all = mean(BehavData.choice_latency(BehavData.Block == 3 & (BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3)));

        block_1_large_choice_latency = mean(BehavData.choice_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 1));
        block_2_large_choice_latency = mean(BehavData.choice_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 2));
        block_3_large_choice_latency = mean(BehavData.choice_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 3));

        block_1_small_choice_latency = mean(BehavData.choice_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 1));
        block_2_small_choice_latency = mean(BehavData.choice_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 2));
        block_3_small_choice_latency = mean(BehavData.choice_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 3));
        if ismember('trial_after_shk', BehavData.Properties.VariableNames)
            block_1_large_choice_latency_following_shk = mean(BehavData.choice_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 1 & BehavData.trial_after_shk == 1));
            block_2_large_choice_latency_following_shk = mean(BehavData.choice_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 2 & BehavData.trial_after_shk == 1));
            block_3_large_choice_latency_following_shk = mean(BehavData.choice_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 3 & BehavData.trial_after_shk == 1));

            block_1_small_choice_latency_following_shk = mean(BehavData.choice_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 1 & BehavData.trial_after_shk == 1));
            block_2_small_choice_latency_following_shk = mean(BehavData.choice_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 2 & BehavData.trial_after_shk == 1));
            block_3_small_choice_latency_following_shk = mean(BehavData.choice_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 3 & BehavData.trial_after_shk == 1));
        else
            block_1_large_choice_latency_following_shk = NaN;
            block_2_large_choice_latency_following_shk = NaN;
            block_3_large_choice_latency_following_shk = NaN;

            block_1_small_choice_latency_following_shk = NaN;
            block_2_small_choice_latency_following_shk = NaN;
            block_3_small_choice_latency_following_shk = NaN;    
        
        
        end
        block_1_collect_latency_all = mean(BehavData.collect_latency(BehavData.Block == 1 & [BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3]));
        block_2_collect_latency_all = mean(BehavData.collect_latency(BehavData.Block == 2 & [BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3]));
        block_3_collect_latency_all = mean(BehavData.collect_latency(BehavData.Block == 3 & [BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3]));

        block_1_large_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 1));
        block_2_large_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 2));
        block_3_large_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 3));

        block_1_large_collect_latency_no_shk_trials = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 1 & BehavData.shock == 0));
        block_2_large_collect_latency_no_shk_trials = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 2 & BehavData.shock == 0));
        block_3_large_collect_latency_no_shk_trials = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 3 & BehavData.shock == 0));

        block_1_small_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 1));
        block_2_small_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 2));
        block_3_small_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 3));
        
        [BehavData_filter_for_ITI_calc,trials,varargin]=TrialFilter_test(BehavData, 'OMITALL', 0, 'BLANK_TOUCH', 0); 
        ITI_times =  BehavData_filter_for_ITI_calc.stTime(2:end) - BehavData_filter_for_ITI_calc.collectionTime(1:end-1);
        blocks_for_ITI_times = BehavData_filter_for_ITI_calc.Block(2:end);
        
        block_1_mean_ITI_length = mean(ITI_times(blocks_for_ITI_times == 1));
        block_2_mean_ITI_length = mean(ITI_times(blocks_for_ITI_times == 2));
        block_3_mean_ITI_length = mean(ITI_times(blocks_for_ITI_times == 3));

        if ismember('type_binary', BehavData.Properties.VariableNames)
            large_aborts = sum(BehavData.type_binary == 1);
            small_aborts = sum(BehavData.type_binary == 2);
            large_aborts_block_1 = sum(BehavData.type_binary == 1 & BehavData.Block == 1);
            large_aborts_block_2 = sum(BehavData.type_binary == 1 & BehavData.Block == 2);
            large_aborts_block_3 = sum(BehavData.type_binary == 1 & BehavData.Block == 3);
            small_aborts_block_1 = sum(BehavData.type_binary == 2 & BehavData.Block == 1);
            small_aborts_block_2 = sum(BehavData.type_binary == 2 & BehavData.Block == 2);
            small_aborts_block_3 = sum(BehavData.type_binary == 2 & BehavData.Block == 3);
        else 
            large_aborts = 0;
            small_aborts = 0;
            large_aborts_block_1 = 0;
            large_aborts_block_2 = 0;
            large_aborts_block_3 = 0;
            small_aborts_block_1 = 0;
            small_aborts_block_2 = 0;
            small_aborts_block_3 = 0;
        end
        if ismember('collectionTime_end', BehavData.Properties.VariableNames)
            large_consum_duration_block_1 = nanmean(BehavData.consum_duration(BehavData.bigSmall == 1.2 & BehavData.Block == 1));
            large_consum_duration_block_2 = nanmean(BehavData.consum_duration(BehavData.bigSmall == 1.2 & BehavData.Block == 2));
            large_consum_duration_block_3 = nanmean(BehavData.consum_duration(BehavData.bigSmall == 1.2 & BehavData.Block == 3));
            small_consum_duration_block_1 = nanmean(BehavData.consum_duration(BehavData.bigSmall == 0.3 & BehavData.Block == 1));
            small_consum_duration_block_2 = nanmean(BehavData.consum_duration(BehavData.bigSmall == 0.3 & BehavData.Block == 2));
            small_consum_duration_block_3 = nanmean(BehavData.consum_duration(BehavData.bigSmall == 0.3 & BehavData.Block == 3));


        else
            large_consum_duration_block_1 = NaN;
            large_consum_duration_block_2 = NaN;
            large_consum_duration_block_3 = NaN;
            small_consum_duration_block_1 = NaN;
            small_consum_duration_block_2 = NaN;
            small_consum_duration_block_3 = NaN;

        end
        trials_completed = sum(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3);
        valid_trials = find(BehavData.Blank_Touch ~= 1 & BehavData.Blank_Touch ~= 2 & BehavData.omissionALL ~= 1);

        
        if contains(session_to_analyze, 'RDT')
            if ~isempty(valid_trials)
                if trials_completed >= 90
                    last_valid_idx = valid_trials(end);

                    session_length = BehavData.collectionTime(last_valid_idx);

                elseif trials_completed <= 90

                    session_length = 5460;
                end
            else

                session_length = NaN;
                warning('All trials are Blank Touch; no valid collection time found.');
            end
        else

            if ~isempty(valid_trials)

                last_valid_idx = valid_trials(end);


                session_length = BehavData.collectionTime(last_valid_idx);

            else

                session_length = NaN; 
                warning('All trials are Blank Touch; no valid collection time found.');
            end
        end


        risk_table(ii,:) = array2table([...
            block_1_large_choice_percent,...
            block_2_large_choice_percent,...
            block_3_large_choice_percent,... 
            block_1_small_choice_percent,...
            block_2_small_choice_percent,...
            block_3_small_choice_percent,...
            block_1_shock_total,...
            block_2_shock_total,...
            block_3_shock_total,...
            block_1_omission_total,...
            block_2_omission_total,...
            block_3_omission_total,...
            large_aborts,...
            small_aborts,...
            large_aborts_block_1,...
            large_aborts_block_2,...
            large_aborts_block_3,...
            small_aborts_block_1,...
            small_aborts_block_2,...
            small_aborts_block_3,...
            lose_shift_percent,...
            lose_omit_percent,...
            lose_stay_percent,...
            win_stay_percent,...
            lose_shift_ratio,...
            lose_stay_ratio,...
            lose_omit_ratio,...
            win_stay_ratio,...
            win_stay_ratio_block_1,...
            win_stay_ratio_block_2,...
            win_stay_ratio_block_3,...
            lose_shift_ratio_block_1,...
            lose_shift_ratio_block_2,...
            lose_shift_ratio_block_3,...
            trials_completed,...
            block_1_choice_latency_all,...
            block_2_choice_latency_all,...
            block_3_choice_latency_all,...
            block_1_large_choice_latency,...
            block_2_large_choice_latency,...
            block_3_large_choice_latency,...
            block_1_small_choice_latency,...
            block_2_small_choice_latency,...
            block_3_small_choice_latency,...
            block_1_collect_latency_all,...
            block_2_collect_latency_all,...
            block_3_collect_latency_all,...
            block_1_large_collect_latency,...
            block_2_large_collect_latency,...
            block_3_large_collect_latency,...
            block_1_large_collect_latency_no_shk_trials,...
            block_2_large_collect_latency_no_shk_trials,...
            block_3_large_collect_latency_no_shk_trials,...
            block_1_small_collect_latency,...
            block_2_small_collect_latency,...
            block_3_small_collect_latency,...
            large_consum_duration_block_1,...
            large_consum_duration_block_2,...
            large_consum_duration_block_3,...
            small_consum_duration_block_1,...
            small_consum_duration_block_2,...
            small_consum_duration_block_3,...
            block_1_large_choice_latency_following_shk,...
            block_2_large_choice_latency_following_shk,...
            block_3_large_choice_latency_following_shk,...
            block_1_small_choice_latency_following_shk,...
            block_2_small_choice_latency_following_shk,...
            block_3_small_choice_latency_following_shk,...
            session_length,...
            block_1_mean_ITI_length,...
            block_2_mean_ITI_length,...
            block_3_mean_ITI_length,...
            ]);

        variable_names = [...
            "block_1_large",...
            "block_2_large",...
            "block_3_large",...
            "block_1_small",...
            "block_2_small",...
            "block_3_small",...
            "block_1_shock_total",...
            "block_2_shock_total",...
            "block_3_shock_total",...
            "block_1_omission_total",...
            "block_2_omission_total",...
            "block_3_omission_total",...
            "large_abort",...
            "small_aborts",...
            "large_aborts_block_1",...
            "large_aborts_block_2",...
            "large_aborts_block_3",...
            "small_aborts_block_1",...
            "small_aborts_block_2",...
            "small_aborts_block_3",...
            "lose_shift",...
            "lose_omit",...
            "lose_stay",...
            "win_stay",...
            "lose_shift_ratio",...
            "lose_stay_ratio",...
            "lose_omit_ratio",...
            "win_stay_ratio",...
            "win_stay_ratio_block_1",...
            "win_stay_ratio_block_2",...
            "win_stay_ratio_block_3",...
            "lose_shift_ratio_block_1",...
            "lose_shift_ratio_block_2",...
            "lose_shift_ratio_block_3",...
            "trials_completed",...
            "block_1_choice_latency_all",...
            "block_2_choice_latency_all",...
            "block_3_choice_latency_all",...
            "block_1_large_choice_latency",...
            "block_2_large_choice_latency",...
            "block_3_large_choice_latency",...
            "block_1_small_choice_latency",...
            "block_2_small_choice_latency",...
            "block_3_small_choice_latency",...
            "block_1_collect_latency_all",...
            "block_2_collect_latency_all",...
            "block_3_collect_latency_all",...
            "block_1_large_collect_latency",...
            "block_2_large_collect_latency",...
            "block_3_large_collect_latency",...
            "block_1_large_collect_latency_no_shk_trials",...
            "block_2_large_collect_latency_no_shk_trials",...
            "block_3_large_collect_latency_no_shk_trials",...
            "block_1_small_collect_latency",...
            "block_2_small_collect_latency",...
            "block_3_small_collect_latency",...
            "large_consum_duration_block_1",...
            "large_consum_duration_block_2",...
            "large_consum_duration_block_3"...
            "small_consum_duration_block_1",...
            "small_consum_duration_block_2",...
            "small_consum_duration_block_3",...
            "block_1_large_choice_latency_following_shk",...
            "block_2_large_choice_latency_following_shk",...
            "block_3_large_choice_latency_following_shk",...
            "block_1_small_choice_latency_following_shk",...
            "block_2_small_choice_latency_following_shk",...
            "block_3_small_choice_latency_following_shk",...
            "session_length",...
            "block_1_mean_ITI_length",...
            "block_2_mean_ITI_length",...
            "block_3_mean_ITI_length",...
            ];

        if ismember('trial_after_shk', BehavData.Properties.VariableNames)
            mean_initiation_latency(ii,:) = [nanmean(BehavData.initiation_delay(BehavData.trial_after_shk == 1)); nanmean(BehavData.initiation_delay(BehavData.trial_after_shk == 0))];
        end
    elseif ~isfield(final_behavior.(currentanimal), session_to_analyze)
        risk_table{ii, :} = nan;
    end
    
end

risk_table.Properties.VariableNames = cellstr(variable_names);
animal_id_table =  array2table(valid_animalIDs);
risk_table = [animal_id_table risk_table];

row_means = nanmean(risk_table{:, 2:4}, 2);


risk_table.Mean_1_to_3 = row_means;
riskiness = risk_table.Mean_1_to_3;

mean_riskiness = mean([risk_table.block_1_large risk_table.block_2_large risk_table.block_3_large], 2);

median_riskiness = median(mean_riskiness);


risk_table.risky = mean_riskiness > median_riskiness;
risk_table.risky = double(risk_table.risky);

if exist('session_weight_raw_grams', 'var') == 1
    pretraining_sessions_sum = cellfun(@sum, pretraining_sessions);
    if strcmp('RDT_D1', session_to_analyze)
        session_weight_free_feed_raw = cell2mat(session_weight_raw_grams')
        session_weight_free_feed_filtered = session_weight_free_feed_raw(valid_mice, 1);
        session_weight_percent_free_feed = (session_weight_free_feed_filtered./freefeed_weight(valid_mice)')*100
        risk_table.freefeed_weight = freefeed_weight(valid_mice)';
        risk_table.session_weight_raw_grams = session_weight_free_feed_filtered;
        risk_table.session_weight_percent_free_feed = session_weight_percent_free_feed;
        risk_table.shock_sens = shock_sens(valid_mice)';
        risk_table.pre_training_sessions_sum = pretraining_sessions_sum(valid_mice)';
        risk_table.sessions_to_rm_criterion = sessions_to_rm_criterion(valid_mice)';
    elseif strcmp('RDT_D2', session_to_analyze)
        session_weight_free_feed_raw = cell2mat(session_weight_raw_grams')
        session_weight_free_feed_filtered = session_weight_free_feed_raw(valid_mice, 2);
        session_weight_percent_free_feed = (session_weight_free_feed_filtered./freefeed_weight(valid_mice)')*100
        risk_table.freefeed_weight = freefeed_weight(valid_mice)';
        risk_table.session_weight_raw_grams = session_weight_free_feed_filtered;
        risk_table.session_weight_percent_free_feed = session_weight_percent_free_feed;
        risk_table.shock_sens = shock_sens(valid_mice)';
        risk_table.pre_training_sessions_sum = pretraining_sessions_sum(valid_mice)';
        risk_table.sessions_to_rm_criterion = sessions_to_rm_criterion(valid_mice)';
    end
end

if exist('hM4Di_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = hM4Di_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);

    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, hM4Di_IDs);

    risk_table_trimmed.TreatmentCondition = hM4Di_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');
elseif exist('stGtACR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = stGtACR_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);

    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, stGtACR_IDs);

    risk_table_trimmed.TreatmentCondition = stGtACR_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');

elseif exist('PdCO_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = PdCO_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);

    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, PdCO_IDs);

    risk_table_trimmed.TreatmentCondition = PdCO_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');


elseif exist('ChrimsonR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = ChrimsonR_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);

    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, ChrimsonR_IDs);

    risk_table_trimmed.TreatmentCondition = ChrimsonR_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');


elseif exist('NAcSh_eNpHR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = mouse_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);

    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, mouse_IDs);

    risk_table_trimmed.TreatmentCondition = NAcSh_eNpHR_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');

elseif exist('treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = mouse_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);
    
    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, mouse_IDs);

    risk_table_trimmed.TreatmentCondition = treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');
end



% IF DOING PR ANALYSIS, DO NOT RUN MOST OF THE CODE ABOVE - ONLY RUN
% SECTION BELOW WITH FILTER ABOVE SET TO PR_D1
for ii = 1:size(valid_animalIDs,1)
    currentanimal = char(valid_animalIDs(ii));
    if isfield(final_behavior.(currentanimal), 'PR_D1')
        BehavData_PR = final_behavior.(currentanimal).PR_D1.uv.BehavData;
        total_PR_presses(ii) = size(BehavData_PR ,1);
        total_PR_rewards(ii) = sum(BehavData_PR.collectTrial == 1);
        PR_final_ratio_completed(ii) = max(BehavData_PR.Ratio - 1);
        PR_final_ratio_reached(ii) = max(BehavData_PR.Ratio);
    
    end



end


%% load Pre_RDT_RM 10 variable dataset, adjust the session_to_analyze = 'RM_D1', run top of script. save risk_table as RM_D1_risk_table.
% set session_to_analyze = 'Pre_RDT_RM', run top of script. then should be
% able to run code below this comment without issue

mean_large_RM_D1 = table2array(mean(RM_D1_risk_table(:, 2:4), 2));
mean_small_RM_D1 = table2array(mean(RM_D1_risk_table(:, 5:7), 2));
sem_large_RM_D1 = table2array(std(RM_D1_risk_table(:, 2:4), 0, 2) ./ sqrt(size(RM_D1_risk_table(:, 2:4), 2)));
sem_small_RM_D1 = table2array(std(RM_D1_risk_table(:, 5:7), 0, 2) ./ sqrt(size(RM_D1_risk_table(:, 5:7), 2)));

mean_large = table2array(mean(risk_table(:, 2:4), 2));
mean_small = table2array(mean(risk_table(:, 5:7), 2));
sem_large = table2array(std(risk_table(:, 2:4), 0, 2) ./ sqrt(size(risk_table(:, 2:4), 2)));
sem_small = table2array(std(risk_table(:, 5:7), 0, 2) ./ sqrt(size(risk_table(:, 5:7), 2)));

cross_sess_large_means = [mean(mean_large_RM_D1), mean(mean_large)]*100;
cross_sess_large_sems = [mean(sem_large_RM_D1), mean(sem_large)]*100;

cross_session_large_all = [mean_large_RM_D1, mean_large]*100;
cross_session_small_all = [mean_small_RM_D1, mean_small]*100;

cross_sess_small_means = [mean(mean_small_RM_D1), mean(mean_small)]*100;
cross_sess_small_sems = [mean(sem_small_RM_D1), mean(sem_small)]*100;


x_points = 1:size(cross_sess_large_means, 2);
figure;
hold on;



height = 450; 
set(gcf, 'Position', [50, 25, width, height]); 

for i = 1:size(cross_session_large_all, 1)
    plot(x_points, cross_session_large_all(i, :), '-', ...
        'Color', [0 0 1 0.6], ...
        'LineWidth', 1.2);
end


for i = 1:size(cross_session_small_all, 1)
    plot(x_points, cross_session_small_all(i, :), '-', ...
        'Color', [1 0 0 0.6], ...
        'LineWidth', 1.2);
end



errorbar(x_points, cross_sess_large_means, cross_sess_large_sems, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 18, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); 

errorbar(x_points, cross_sess_small_means, cross_sess_small_sems, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 18, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); 


xticks(x_points); 
xticklabels({'Early RM', 'Late RM'}); 
xlim([0.5, length(x_points) + 0.5]); 


ylim([0 1.1 * max([cross_sess_large_means + cross_sess_large_sems, ...
                   cross_sess_small_means + cross_sess_small_sems])]); 
set(gca, 'ytick', 0:25:100);


hold off;


num_subjects = size(cross_session_large_all, 1);

subject_id = (1:num_subjects)';
data_wide = table(subject_id, ...
    cross_session_large_all(:,1), cross_session_large_all(:,2), ...
    cross_session_small_all(:,1), cross_session_small_all(:,2), ...
    'VariableNames', {'Subject', 'Day1_Large', 'Day2_Large', 'Day1_Small', 'Day2_Small'});


withinDesign = table(categorical({'Day1'; 'Day2'; 'Day1'; 'Day2'}), ...
    categorical({'Large'; 'Large'; 'Small'; 'Small'}), ...
    'VariableNames', {'Day', 'RewardSize'});

rm = fitrm(data_wide, 'Day1_Large,Day2_Large,Day1_Small,Day2_Small ~ 1', 'WithinDesign', withinDesign);

ranova_table = ranova(rm, 'WithinModel', 'Day*RewardSize');

disp('Repeated Measures ANOVA Results:');
disp(ranova_table);

fprintf('Reward Size effect (RM ANOVA): F(%d,%d) = %.3f, p = %.3e\n', ...
    ranova_table{5,2}, ranova_table{6,2}, ranova_table{5,4}, ranova_table{5,5});
fprintf('Day effect (RM ANOVA): F(%d,%d) = %.3f, p = %.3e\n', ...
    ranova_table{3,2}, ranova_table{4,2}, ranova_table{3,4}, ranova_table{3,5});
fprintf('Interaction (RM ANOVA): F(%d,%d) = %.3f, p = %.3e\n', ...
    ranova_table{7,2}, ranova_table{8,2}, ranova_table{7,4}, ranova_table{7,5});

p_interaction = ranova_table{7,5};
if p_interaction < 0.05
    disp('Significant interaction found! Performing Tukey''s multiple comparisons...');
    fprintf('\n========== POST-HOC TESTS ==========\n');
    
    all_values = [];
    all_groups = {};
    
    all_values = [all_values; data_wide.Day1_Large];
    all_groups = [all_groups; repmat({'Day1_Large'}, num_subjects, 1)];
    
    all_values = [all_values; data_wide.Day2_Large];
    all_groups = [all_groups; repmat({'Day2_Large'}, num_subjects, 1)];
    
    all_values = [all_values; data_wide.Day1_Small];
    all_groups = [all_groups; repmat({'Day1_Small'}, num_subjects, 1)];
    
    all_values = [all_values; data_wide.Day2_Small];
    all_groups = [all_groups; repmat({'Day2_Small'}, num_subjects, 1)];
    
    [~, ~, stats_tukey] = anova1(all_values, all_groups, 'off');
    
    [c, m, h, gnames] = multcompare(stats_tukey, 'CType', 'tukey-kramer', 'Display', 'off');
    
    fprintf('\nTukey''s HSD Multiple Comparisons Results:\n');
    fprintf('%-12s %-12s %10s %10s %10s %10s\n', 'Group 1', 'Group 2', 'Diff', 'Lower CI', 'Upper CI', 'p-value');
    fprintf('%-12s %-12s %10s %10s %10s %10s\n', repmat('-', 1, 12), repmat('-', 1, 12), repmat('-', 1, 10), repmat('-', 1, 10), repmat('-', 1, 10), repmat('-', 1, 10));
    
    for i = 1:size(c, 1)
        group1_idx = c(i, 1);
        group2_idx = c(i, 2);
        group1_name = gnames{group1_idx};
        group2_name = gnames{group2_idx};
        diff = c(i, 4);
        lower_ci = c(i, 3);
        upper_ci = c(i, 5);
        p_val = c(i, 6);
        
        fprintf('%-12s %-12s %10.3f %10.3f %10.3f %10.3e', ...
            group1_name, group2_name, diff, lower_ci, upper_ci, p_val);
        
        if p_val < 0.05
            fprintf(' *');
        end
        fprintf('\n');
    end
    
    fprintf('\n* indicates significant difference (p < 0.05)\n');
    
    fprintf('\nGroup Means:\n');
    for i = 1:length(gnames)
        fprintf('%-12s: %.3f ± %.3f (SEM)\n', gnames{i}, m(i, 1), m(i, 2));
    end
    
    fprintf('\n--- Key Comparisons ---\n');
    fprintf('Reward Size comparisons within each day:\n');
    
    for i = 1:size(c, 1)
        group1_name = gnames{c(i, 1)};
        group2_name = gnames{c(i, 2)};
        
        if (strcmp(group1_name, 'Day1_Large') && strcmp(group2_name, 'Day1_Small')) || ...
           (strcmp(group1_name, 'Day1_Small') && strcmp(group2_name, 'Day1_Large'))
            diff = c(i, 4);
            lower_ci = c(i, 3);
            upper_ci = c(i, 5);
            p_val = c(i, 6);
            
            fprintf('Day 1 (Large vs Small): Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.3e', ...
                abs(diff), lower_ci, upper_ci, p_val);
            
            if p_val < 0.05
                fprintf(' *');
            end
            fprintf('\n');
            break;
        end
    end
    
    for i = 1:size(c, 1)
        group1_name = gnames{c(i, 1)};
        group2_name = gnames{c(i, 2)};
        
        if (strcmp(group1_name, 'Day2_Large') && strcmp(group2_name, 'Day2_Small')) || ...
           (strcmp(group1_name, 'Day2_Small') && strcmp(group2_name, 'Day2_Large'))
            diff = c(i, 4);
            lower_ci = c(i, 3);
            upper_ci = c(i, 5);
            p_val = c(i, 6);
            
            fprintf('Day 2 (Large vs Small): Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.3e', ...
                abs(diff), lower_ci, upper_ci, p_val);
            
            if p_val < 0.05
                fprintf(' *');
            end
            fprintf('\n');
            break;
        end
    end
    
    fprintf('\nDay comparisons within each reward size:\n');
    
    for i = 1:size(c, 1)
        group1_name = gnames{c(i, 1)};
        group2_name = gnames{c(i, 2)};
        
        if (strcmp(group1_name, 'Day1_Large') && strcmp(group2_name, 'Day2_Large')) || ...
           (strcmp(group1_name, 'Day2_Large') && strcmp(group2_name, 'Day1_Large'))
            diff = c(i, 4);
            lower_ci = c(i, 3);
            upper_ci = c(i, 5);
            p_val = c(i, 6);
            
            fprintf('Large Reward (Day 1 vs Day 2): Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.3e', ...
                abs(diff), lower_ci, upper_ci, p_val);
            
            if p_val < 0.05
                fprintf(' *');
            end
            fprintf('\n');
            break;
        end
    end
    
    for i = 1:size(c, 1)
        group1_name = gnames{c(i, 1)};
        group2_name = gnames{c(i, 2)};
        
        if (strcmp(group1_name, 'Day1_Small') && strcmp(group2_name, 'Day2_Small')) || ...
           (strcmp(group1_name, 'Day2_Small') && strcmp(group2_name, 'Day1_Small'))
            diff = c(i, 4);
            lower_ci = c(i, 3);
            upper_ci = c(i, 5);
            p_val = c(i, 6);
            
            fprintf('Small Reward (Day 1 vs Day 2): Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.3e', ...
                abs(diff), lower_ci, upper_ci, p_val);
            
            if p_val < 0.05
                fprintf(' *');
            end
            fprintf('\n');
            break;
        end
    end
    
else
    disp('No significant interaction effect found.');
    
    p_reward_size = ranova_table{5,5};
    p_day = ranova_table{3,5};
    
    if p_reward_size < 0.05
        fprintf('\n--- Main effect of Reward Size is significant (p = %.3e) ---\n', p_reward_size);
        large_mean = mean([data_wide.Day1_Large; data_wide.Day2_Large]);
        small_mean = mean([data_wide.Day1_Small; data_wide.Day2_Small]);
        fprintf('Large Reward mean: %.3f\n', large_mean);
        fprintf('Small Reward mean: %.3f\n', small_mean);
        fprintf('Difference: %.3f\n', large_mean - small_mean);
    end
    
    if p_day < 0.05
        fprintf('\n--- Main effect of Day is significant (p = %.3e) ---\n', p_day);
        day1_mean = mean([data_wide.Day1_Large; data_wide.Day1_Small]);
        day2_mean = mean([data_wide.Day2_Large; data_wide.Day2_Small]);
        fprintf('Day 1 mean: %.3f\n', day1_mean);
        fprintf('Day 2 mean: %.3f\n', day2_mean);
        fprintf('Difference: %.3f\n', day1_mean - day2_mean);
        
        dayComp = multcompare(rm, 'Day', 'ComparisonType', 'tukey-kramer');
        
        fprintf('Post-hoc comparisons for Day (across both reward sizes):\n');
        for i = 1:size(dayComp, 1)
            fprintf('Day %s vs Day %s: Mean Diff = %.3f, CI = [%.3f, %.3f], p = %.3e\n', ...
                dayComp.Day_1{i}, dayComp.Day_2{i}, ...
                dayComp.Difference(i), dayComp.Lower(i), ...
                dayComp.Upper(i), dayComp.pValue(i));
        end
    end
end

%%

large_choice = [risk_table.block_1_large, risk_table.block_2_large, risk_table.block_3_large]*100;
small_choice = [risk_table.block_1_small, risk_table.block_2_small, risk_table.block_3_small]*100;

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));

x_points = 1:size(large_choice, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ...
        'LineWidth', 1.2);
end

for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 1.1 * max([mean_large + sem_large, ...
                   mean_small + sem_small])]);
set(gca, 'ytick', 0:25:100);

hold off;

%%

large_choice = [risk_table.block_1_large_choice_latency, risk_table.block_2_large_choice_latency, risk_table.block_3_large_choice_latency];
small_choice = [risk_table.block_1_small_choice_latency, risk_table.block_2_small_choice_latency, risk_table.block_3_small_choice_latency];

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));

x_points = 1:size(large_choice, 2);

figure;
hold on;

width = 100;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ...
        'LineWidth', 1.2);
end

for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

set(gca, 'ytick', 0:5:20);
ylim([0 20])

hold off;

%%

large_choice = [risk_table.block_1_large_collect_latency, risk_table.block_2_large_collect_latency, risk_table.block_3_large_collect_latency];
small_choice = [risk_table.block_1_small_collect_latency, risk_table.block_2_small_collect_latency, risk_table.block_3_small_collect_latency];


mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));

x_points = 1:size(large_choice, 2);

figure;
hold on;

width = 100;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ...
        'LineWidth', 1.2);
end

for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

set(gca, 'ytick', 0:1:5);
ylim([0 4])

hold off;

%%

large_choice = [risk_table.large_consum_duration_block_1, risk_table.large_consum_duration_block_2, risk_table.large_consum_duration_block_3];
small_choice = [risk_table.small_consum_duration_block_1, risk_table.small_consum_duration_block_2, risk_table.small_consum_duration_block_3];

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));

x_points = 1:size(large_choice, 2);

figure;
hold on;

width = 100;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ...
        'LineWidth', 1.2);
end

for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

set(gca, 'ytick', 0:1:5);
ylim([0 5])

hold off;

%%

large_choice = [risk_table.large_aborts_block_1, risk_table.large_aborts_block_2, risk_table.large_aborts_block_3];
small_choice = [risk_table.small_aborts_block_1, risk_table.small_aborts_block_2, risk_table.small_aborts_block_3];

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));

x_points = 1:size(large_choice, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ...
        'LineWidth', 1.2);
end

for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

set(gca, 'ytick', 0:25:100);

hold off;




%% for plotting lose_shift ratio across blocks

large_choice = [risk_table.lose_shift_ratio_block_1, risk_table.lose_shift_ratio_block_2, risk_table.lose_shift_ratio_block_3];
small_choice = [risk_table.win_stay_ratio_block_1, risk_table.win_stay_ratio_block_2, risk_table.win_stay_ratio_block_3];

large_choice(isnan(large_choice))=0;
small_choice(isnan(small_choice))=0;

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));

x_points = 1:size(large_choice, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ...
        'LineWidth', 1.2);
end

for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

set(gca, 'ytick', 0:1);

hold off;




%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_large(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('hM4Di', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 100]);
set(gca, 'ytick', 0:25:100);

hold off;

%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_small(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('hM4Di', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 100]);
set(gca, 'ytick', 0:25:100);

hold off;

%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_aborts_block_1(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('hM4Di', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

set(gca, 'ytick', 0:25:100);

hold off;



%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_choice_latency(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('hM4Di', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 30]);
set(gca, 'ytick', 0:10:30);

hold off;

%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_collect_latency(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('hM4Di', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 15]);
set(gca, 'ytick', 0:5:15);

hold off;

%%
large_choice_mCherry = [risk_table.block_1_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_omission_total(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('hM4Di', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 60]);
set(gca, 'ytick', 0:20:60);

hold off;

%%
mCherry_session_length = [risk_table.session_length(strcmp('mCherry', risk_table.TreatmentCondition))]
hm4di_session_length = [risk_table.session_length(strcmp('hM4Di', risk_table.TreatmentCondition))]

bar_separation_value = 3;

figure;
width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

swarmchart(ones(1, length(mCherry_session_length)), mCherry_session_length, 90, 'o', 'MarkerFaceColor', mCherry_color);
hold on;
swarmchart(ones(1, length(hm4di_session_length)) * bar_separation_value, hm4di_session_length, 100, 'square', 'MarkerFaceColor', hM4Di_color);

plot([0.5; 1.5], [mean(mCherry_session_length); mean(mCherry_session_length)], 'LineWidth', 3, 'Color', 'red');
plot([bar_separation_value - 0.5; bar_separation_value + 0.5], [mean(hm4di_session_length); mean(hm4di_session_length)], 'LineWidth', 3, 'Color', 'blue');

yline(0, 'k--');
xticks([1, bar_separation_value, bar_separation_value + 2]);
xticklabels({'Type 1.2 SHK Resp.', 'Type 1.2 Not SHK Resp.'});
xtickformat('%.1f');
ytickformat('%.1f');
hold off;

[h p ci stats] = ttest2(mCherry_session_length, hm4di_session_length)


%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_mean_ITI_length(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_mean_ITI_length(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_mean_ITI_length(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_mean_ITI_length(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_mean_ITI_length(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_mean_ITI_length(strcmp('hM4Di', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 80]);
set(gca, 'ytick', 0:10:80);

yline(8)
hold off;

%%
large_choice_mCherry = [risk_table.large_consum_duration_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_consum_duration_block_1(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('hM4Di', risk_table.TreatmentCondition))];


mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_hM4Di = nanmean(large_choice_hM4Di, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_hM4Di = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_hM4Di, sem_hM4Di, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 6]);
set(gca, 'ytick', 0:2:6);

hold off;


%%
large_choice_mCherry = [risk_table.small_consum_duration_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.small_consum_duration_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.small_consum_duration_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.small_consum_duration_block_1(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.small_consum_duration_block_2(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.small_consum_duration_block_3(strcmp('hM4Di', risk_table.TreatmentCondition))];


mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_hM4Di = nanmean(large_choice_hM4Di, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_hM4Di = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_hM4Di, sem_hM4Di, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 6]);
set(gca, 'ytick', 0:2:6);

hold off;


%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.block_1_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_large(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('stGtACR', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 100]);
set(gca, 'ytick', 0:25:100);

hold off;

%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.block_1_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_small(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('stGtACR', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 100]);
set(gca, 'ytick', 0:25:100);

hold off;

%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.block_1_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_choice_latency(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('stGtACR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 17]);
set(gca, 'ytick', 0:5:15);

hold off;

%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.block_1_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_collect_latency(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('stGtACR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 17]);
set(gca, 'ytick', 0:5:15);

hold off;

%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_aborts_block_1(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('stGtACR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 250;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 125]);
set(gca, 'ytick', 0:25:125);

hold off;




%%
mCherry_session_length = [risk_table.session_length(strcmp('mCherry', risk_table.TreatmentCondition))]
hm4di_session_length = [risk_table.session_length(strcmp('stGtACR', risk_table.TreatmentCondition))]

bar_separation_value = 3;

figure;
width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

swarmchart(ones(1, length(mCherry_session_length)), mCherry_session_length, 90, 'o', 'MarkerFaceColor', mCherry_color);
hold on;
swarmchart(ones(1, length(hm4di_session_length)) * bar_separation_value, hm4di_session_length, 100, 'square', 'MarkerFaceColor', stGtACR_color);

plot([0.5; 1.5], [mean(mCherry_session_length); mean(mCherry_session_length)], 'LineWidth', 3, 'Color', 'red');
plot([bar_separation_value - 0.5; bar_separation_value + 0.5], [mean(hm4di_session_length); mean(hm4di_session_length)], 'LineWidth', 3, 'Color', 'blue');

yline(0, 'k--');
xticks([1, bar_separation_value, bar_separation_value + 2]);
xticklabels({'Type 1.2 SHK Resp.', 'Type 1.2 Not SHK Resp.'});
xtickformat('%.1f');
ytickformat('%.1f');
hold off;

[h p ci stats] = ttest2(mCherry_session_length, hm4di_session_length)


%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_mean_ITI_length(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_mean_ITI_length(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_mean_ITI_length(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_mean_ITI_length(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_2_mean_ITI_length(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_3_mean_ITI_length(strcmp('stGtACR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 150]);
set(gca, 'ytick', 0:25:150);

yline(8)
hold off;



%% for PdCO vs ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_PdCO = [risk_table.block_1_large(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('PdCO', risk_table.TreatmentCondition))]*100;
large_choice_ChrimsonR = [risk_table.block_1_large(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('ChrimsonR', risk_table.TreatmentCondition))]*100;

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 100]);
set(gca, 'ytick', 0:25:100);

hold off;

%% for PdCO vs ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_PdCO = [risk_table.block_1_small(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('PdCO', risk_table.TreatmentCondition))]*100;
large_choice_ChrimsonR = [risk_table.block_1_small(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('ChrimsonR', risk_table.TreatmentCondition))]*100;

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 100]);
set(gca, 'ytick', 0:25:100);

hold off;

%% for PdCO vs ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_choice_latency_all(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_choice_latency_all(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_choice_latency_all(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_PdCO = [risk_table.block_1_choice_latency_all(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_choice_latency_all(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_choice_latency_all(strcmp('PdCO', risk_table.TreatmentCondition))];
large_choice_ChrimsonR = [risk_table.block_1_choice_latency_all(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_choice_latency_all(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_choice_latency_all(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 17]);
set(gca, 'ytick', 0:5:15);

hold off;

%% for PdCO vs ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_PdCO = [risk_table.block_1_large_collect_latency(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('PdCO', risk_table.TreatmentCondition))];
large_choice_ChrimsonR = [risk_table.block_1_large_collect_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 17]);
set(gca, 'ytick', 0:5:15);

hold off;



%% for PdCO vs ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_PdCO = [risk_table.large_aborts_block_1(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('PdCO', risk_table.TreatmentCondition))];
large_choice_ChrimsonR = [risk_table.large_aborts_block_1(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 125]);
set(gca, 'ytick', 0:25:125);

hold off;





%%
mCherry_session_length = [risk_table.session_length(strcmp('mCherry', risk_table.TreatmentCondition))]
PdCO_session_length = [risk_table.session_length(strcmp('PdCO', risk_table.TreatmentCondition))]
ChrimsonR_session_length = [risk_table.session_length(strcmp('ChrimsonR', risk_table.TreatmentCondition))]

bar_separation_value = [3 5];

figure;
width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

swarmchart(ones(1, length(mCherry_session_length)), mCherry_session_length, 90, 'o', 'MarkerFaceColor', mCherry_color);
hold on;
swarmchart(ones(1, length(PdCO_session_length)) * bar_separation_value(1), PdCO_session_length, 100, '^', 'MarkerFaceColor', PdCO_color);
hold on;
swarmchart(ones(1, length(ChrimsonR_session_length)) * bar_separation_value(2), ChrimsonR_session_length, 100, 'square', 'MarkerFaceColor', ChrimsonR_color);

plot([0.5; 1.5], [mean(mCherry_session_length); mean(mCherry_session_length)], 'LineWidth', 3, 'Color', 'red');
plot([bar_separation_value(1) - 0.5; bar_separation_value(1) + 0.5], [mean(PdCO_session_length); mean(PdCO_session_length)], 'LineWidth', 3, 'Color', 'blue');
plot([bar_separation_value(2) - 0.5; bar_separation_value(2) + 0.5], [mean(ChrimsonR_session_length); mean(ChrimsonR_session_length)], 'LineWidth', 3, 'Color', 'black');

yline(0, 'k--');
xticks([1, bar_separation_value, bar_separation_value + 2]);
xticklabels({'Type 1.2 SHK Resp.', 'Type 1.2 Not SHK Resp.'});
xtickformat('%.1f');
ytickformat('%.1f');
hold off;

all_data = [mCherry_session_length; PdCO_session_length; ChrimsonR_session_length];

group = [ones(length(mCherry_session_length), 1); 
         2*ones(length(PdCO_session_length), 1); 
         3*ones(length(ChrimsonR_session_length), 1)];

[p, tbl, stats] = anova1(all_data, group);

if p < 0.05
    figure;
    [c, m, h, gnames] = multcompare(stats);
    
    comparison_table = array2table(c, ...
        'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDiff', 'UpperCI', 'pValue'});
    
    disp('Pairwise Comparisons:');
    disp(comparison_table);
    
    sig_comparisons = comparison_table(comparison_table.pValue < 0.05, :);
    disp('Significant pairwise comparisons (p < 0.05):');
    disp(sig_comparisons);
end


%% for PdCO vs ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_mean_ITI_length(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_mean_ITI_length(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_mean_ITI_length(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_PdCO = [risk_table.block_1_mean_ITI_length(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_mean_ITI_length(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_mean_ITI_length(strcmp('PdCO', risk_table.TreatmentCondition))];
large_choice_ChrimsonR = [risk_table.block_1_mean_ITI_length(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_mean_ITI_length(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_mean_ITI_length(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 120]);
set(gca, 'ytick', 0:20:120);

yline(8)
hold off;










%% for BLA-NAcSh 
large_choice_mCherry = [risk_table.block_1_large(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('Conting', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_large(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('Yoked', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 100]);
set(gca, 'ytick', 0:25:100);

hold off;

%% for BLA-NAcSh Contingent vs Yoked

large_choice_mCherry = [risk_table.block_1_small(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('Conting', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_small(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('Yoked', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 100]);
set(gca, 'ytick', 0:25:100);

hold off;

%% 

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('Conting', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_aborts_block_1(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('Yoked', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 115]);
set(gca, 'ytick', 0:25:100);

hold off;

%%

large_choice_mCherry = [risk_table.block_1_large_choice_latency(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('Conting', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_choice_latency(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('Yoked', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 20]);
set(gca, 'ytick', 0:5:20);

hold off;


%%

large_choice_mCherry = [risk_table.block_1_large_collect_latency(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('Conting', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('Conting', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_collect_latency(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('Yoked', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('Yoked', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));

x_points = 1:size(large_choice_mCherry, 2);

figure;
hold on;

width = 200;
height = 450;
set(gcf, 'Position', [50, 25, width, height]);

for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ...
        'LineWidth', 1.2);
end

for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ...
        'LineWidth', 1.2);
end

errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none');

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none');

xticks(x_points);
xticklabels({'0', '50', '75'});
xlim([0.5, length(x_points) + 0.5]);

ylim([0 6]);
set(gca, 'ytick', 0:1:5);

hold off;

%% For Figure 1G and Figure S3N
% animalIDs = (fieldnames(final_SLEAP));
session_to_analyze = 'RDT_D1';

if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
    fieldsToRemove = {'BLA_Insc_13', 'BLA_Insc_18', 'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_SLEAP, fieldsToRemove{i})
            final_SLEAP = rmfield(final_SLEAP, fieldsToRemove{i});
        end
    end
elseif strcmp('RDT_D2', session_to_analyze)

    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_SLEAP, fieldsToRemove{i})
            final_SLEAP = rmfield(final_SLEAP, fieldsToRemove{i});
        end
    end
end

animalIDs = (fieldnames(final_SLEAP));

array_for_means = 1; 
array_for_means_second = []; % 5

b1_large_path_length = [];
b2_large_path_length = [];
b3_large_path_length = [];

b1_small_path_length = [];
b2_small_path_length = [];
b3_small_path_length = [];

animals_with_sessions = {}; 

for dd = 1:size(animalIDs)

    select_mouse = animalIDs{dd};
    path_length_array = [];
    if isfield(final_SLEAP.(select_mouse), session_to_analyze)
        animals_with_sessions{dd} = select_mouse; 

        



        SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;

        X_data = SLEAP_data.corrected_x_pix;
        Y_data = SLEAP_data.corrected_y_pix;
        [X_data, Y_data] = correct_XY_outliers_v1(X_data, Y_data);

        % onset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.stTime';
        % choice_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.choiceTime';
        % offset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.collectionTime';
        fs_cam = 30; %set sampling rate according to camera, this is hard coded for now
        % time_ranges_trials = [onset_trials; choice_trials; offset_trials];


        % gcamp_samples = 1:1:size(Y_dF_all_session, 2);

        % gcamp_time = (0:length(F405_downsampled_data)-1)/fs_cam;

        SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;

        % velocity_data = final_SLEAP.(select_mouse).(session_to_analyze).zscored_SLEAP_data_velocity';

        velocity_data = zscore(SLEAP_data.vel_cm_s)';

        % SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data;

        
        BehavData = final_SLEAP.(select_mouse).(session_to_analyze).BehavData;

        % uncomment below if you want to run this on specific trials, eg
        % using the 10x variable data
        
        % BehavData = behav_tbl_iter{1, 1}{dd};
        % behav_part_1 = behav_tbl_iter{array_for_means, 1}{dd, 1};
        % if isempty(array_for_means_second)
        %     BehavData = behav_part_1;
        % 
        % else
        %     behav_part_2 = behav_tbl_iter{array_for_means_second, 1}{dd, 1};
        %     BehavData = vertcat(behav_part_1, behav_part_2);
        % end

        onset_trials = BehavData.stTime';
        choice_trials = BehavData.choiceTime';
        offset_trials = BehavData.collectionTime';
        time_ranges_trials = [onset_trials; choice_trials; offset_trials];
        adjusted_start_time = BehavData.TrialPossible(1)-60;
        SLEAP_data.idx_time = SLEAP_data.idx_time+adjusted_start_time;



        %%
        % FILTER ALL EXISTING DATA ON THESE TIME RANGES
        % filter streams
        if ~isempty(SLEAP_data)
            filtered_motion = [];
            max_ind = SLEAP_data.idx_frame(end);
            idx_frame_redo = 1:1:size(SLEAP_data, 1);
            good_index = 1;
            for j = 1:size(time_ranges_trials,2)
                onset = round(time_ranges_trials(1,j)*fs_cam)+1;
                choice = round(time_ranges_trials(2,j)*fs_cam)+1;
                offset = round(time_ranges_trials(3,j)*fs_cam)+1;

                % throw it away if onset or offset extends beyond recording window
                if isinf(offset)
                    if onset <= max_ind && onset > 0
                        filtered_motion{good_index} = SLEAP_data(onset:end);
                        break %return
                    end
                else
                    % if offset <= max_ind && offset > 0 && onset <= max_ind && onset > 0
                        % buffering this by adding +1 to the end time for now, for
                        % some reason the array seems too short without?
                        % after some extensive checking, it seems like the
                        % strangeness where the body @ start and @ end does not
                        % overlap often comes from the fact that the mouse's tail
                        % can trigger the IR beam in the food cup on a non-trivial
                        % # of trials
                        filtered_motion{j}= [X_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j))'; Y_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j))']; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
                        % choice_times{j} = [X_data(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'; Y_data(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'];
                        filtered_velocity{j}= velocity_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j));
                        % filtered_gcamp{j}= Y_dF_all_session(gcamp_samples(onset:offset));
                        % filtered{j} = Y_data_filtered(SLEAP_data.idx_frame(onset:offset))'; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
                        good_index = good_index + 1;
                    % end
                end
            end

        end

        for qq = 1:size(filtered_motion, 2)
            coordinates = filtered_motion{1, qq}';
            distances = pdist2(coordinates, coordinates);

            % Sum up the distances to get the total path length
            path_length = sum(diag(distances, 1));

            % disp(['Path Length: ', num2str(path_length)]);

            distances_matrix{qq} = distances;
            path_length_array(qq) = path_length;
            
            
            clear coordinates distances path_length
        end
        % filtered_motion_array{}
      
        b1_large_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 1));
        b2_large_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 2));
        b3_large_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 3));

        b1_small_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 1));
        b2_small_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 2));
        b3_small_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 3));

        large_path_length_all_blocks(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2));
        small_path_length_all_blocks(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3));

        b1_path_length_mouse{dd} = path_length_array(1, BehavData.Block == 1 &  BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0);
        b2_path_length_mouse{dd} = path_length_array(1, BehavData.Block == 2 &  BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0);
        b3_path_length_mouse{dd} = path_length_array(1, BehavData.Block == 3 &  BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0);
        path_length_array_mouse{dd} = path_length_array';

    elseif ~isfield(final_SLEAP.(select_mouse), session_to_analyze)
        b1_path_length_mouse{dd} = nan;
        b2_path_length_mouse{dd} = nan;
        b3_path_length_mouse{dd} = nan;
        path_length_array_mouse{dd} = nan';

    end


end

path_length_concat = cat(1, path_length_array_mouse{:});



large_path_length = ([b1_large_path_length; b2_large_path_length; b3_large_path_length])';
small_path_length = ([b1_small_path_length; b2_small_path_length; b3_small_path_length])';



mean_large_path = nanmean(large_path_length, 1);
mean_small_path = nanmean(small_path_length, 1);

sem_large = nanstd(large_path_length, 0, 1) ./ sqrt(size(large_path_length, 1));
sem_small = nanstd(small_path_length, 0, 1) ./ sqrt(size(small_path_length, 1));



x_points = 1:3;



figure;
hold on;


width = 100; 
height = 450; 
set(gcf, 'Position', [50, 25, width, height]); 

for i = 1:size(large_path_length, 1)
    plot(x_points, large_path_length(i, :), '-', ...
        'Color', [0 0 1 0.6], ...
        'LineWidth', 1.2);
end


% Plot individual lines for "Small" data
for i = 1:size(small_path_length, 1)
    plot(x_points, small_path_length(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % 
        'LineWidth', 1.2);
end


errorbar(x_points, mean_large_path, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large'); 

errorbar(x_points, mean_small_path, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); 


xticks(x_points); 
xticklabels({'0', '50%', '75%'}); 
xlim([0.5, length(x_points) + 0.5]); 


set(gca, 'ytick', 0:200:1000);
ylim([0 800])


hold off;


