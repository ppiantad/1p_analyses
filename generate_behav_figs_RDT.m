% iter = iter+1;
% neuron_num = 0;
% these are mice that did not complete the entire session - kinda have to
% toss them to do some comparisons during RDT

% final_behavior = final_SLEAP; % for hM4Di data;


animalIDs = (fieldnames(final_behavior));

session_to_analyze = 'RDT_OPTO_CHOICE'

if any(startsWith(string(animalIDs), 'RDT-F'))
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
end


%%


animalIDs = (fieldnames(final_behavior));



possible_effectors = {'stGtACR', 'PdCO', 'ChrimsonR', 'hM4Di'}
possible_controls = {'mCherry', 'EGFP'}


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

elseif exist('NAcSh_eNpHR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = mouse_IDs(valid_mice);

else
    valid_animalIDs = animalIDs;

end



for ii = 1:size(valid_animalIDs,1) % 1:size(fieldnames(final),1)
    currentanimal = char(valid_animalIDs(ii));
    
    if isfield(final_behavior.(currentanimal), session_to_analyze)
        % if contains(session_to_analyze, 'CNO')
        %     BehavData = final_behavior.(currentanimal).(session_to_analyze).BehavData;
        % else
        %     BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        % end




        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        for BehavDataRow = 1:size(BehavData,1)
            if BehavData.shock(BehavDataRow) == 1
                kk = 1;
                while true
                    if (BehavDataRow + kk) > size(BehavData, 1)  % Check if index exceeds the number of rows
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
            % if BehavDataRow > 1
            %     BehavData.initiation_delay(BehavDataRow+1) = BehavData.stTime(BehavDataRow+1)-BehavData.collectionTime(BehavDataRow); 
            % end

        end

        for BehavDataRow = 1:size(BehavData,1)
            if BehavData.bigSmall(BehavDataRow) ~= 999 & ~isnan(BehavData.bigSmall(BehavDataRow))
                kk = 1;
                while true
                    if (BehavDataRow + kk) > size(BehavData, 1)  % Check if index exceeds the number of rows
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
            % if BehavDataRow > 1
            %     BehavData.initiation_delay(BehavDataRow+1) = BehavData.stTime(BehavDataRow+1)-BehavData.collectionTime(BehavDataRow);
            % end

        end
        
        % only use rewarded trials for this, otherwise things get wonky
        [BehavData,trials,varargin]=TrialFilter_test(BehavData,'ALL', 1); 
        large_small_trials_only = BehavData(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3, :);
        % large_small_trials_only = BehavData((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0, :);


        % if a mouse finishes all trials, there should be 66 free choice
        % trials
        desired_rows = 90;

        large_trials_true = large_small_trials_only.bigSmall == 1.2;
        small_trials_true = large_small_trials_only.bigSmall == 0.3;

        % % if a mouse finishes all trials, there should be 66 free choice
        % % trials
        % desired_rows = 66;
        % 
        % 
        % large_trials_true = large_small_trials_only.bigSmall == 1.2 & large_small_trials_only.ForceFree == 0;
        % small_trials_true = large_small_trials_only.bigSmall == 0.3 & large_small_trials_only.ForceFree == 0;

        % Pad large_trials_true with zeros if it has fewer rows than desired
        if size(large_trials_true, 1) < desired_rows
            padding = zeros(desired_rows - size(large_trials_true, 1), size(large_trials_true, 2));
            large_trials_true = [large_trials_true; padding];
        end

        % Pad small_trials_true with zeros if it has fewer rows than desired
        if size(small_trials_true, 1) < desired_rows
            padding = zeros(desired_rows - size(small_trials_true, 1), size(small_trials_true, 2));
            small_trials_true = [small_trials_true; padding];
        end
        % Store padded arrays in the mouse sequences
        large_sequences_mouse(ii, :) = large_trials_true;
        small_sequences_mouse(ii, :) = small_trials_true;

        block_1_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 1 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 1); 
        block_1_small_choice_percent = sum(BehavData.bigSmall == 0.3 & BehavData.Block == 1 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 1); 
        % block_1_mouse(ii,:) = [block_1(1, 1) block_1(end, 2)];
        block_2_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 2 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 2);  
        block_2_small_choice_percent = sum(BehavData.bigSmall == 0.3 & BehavData.Block == 2 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 2);  
        % block_2_mouse(ii,:) = [block_2(1, 1) block_2(end, 2)];
        block_3_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 3 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);  
        block_3_small_choice_percent = sum(BehavData.bigSmall == 0.3 & BehavData.Block == 3 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);  
        block_1_omission_total = sum(BehavData.omissionALL == 1 & BehavData.Block == 1);
        block_2_omission_total = sum(BehavData.omissionALL == 1 & BehavData.Block == 2);
        block_3_omission_total = sum(BehavData.omissionALL == 1 & BehavData.Block == 3);
        % block_3_mouse(ii,:) = [block_3(1, 1) block_3(end, 2)];
        
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
        block_1_large_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 1));
        block_2_large_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 2));
        block_3_large_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 3));

        block_1_large_collect_latency_no_shk_trials = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 1 & BehavData.shock == 0));
        block_2_large_collect_latency_no_shk_trials = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 2 & BehavData.shock == 0));
        block_3_large_collect_latency_no_shk_trials = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 3 & BehavData.shock == 0));

        block_1_small_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 1));
        block_2_small_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 2));
        block_3_small_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 0.3 & BehavData.Block == 3));
        if ismember('type_binary', BehavData.Properties.VariableNames)
            large_aborts = sum(BehavData.type_binary == 1); %[] sum(BehavData.type_binary == 1)
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
        session_length = BehavData.collectionTime(end); 


        risk_table(ii,:) = array2table([...
            %{currentanimal},...
            block_1_large_choice_percent,...
            block_2_large_choice_percent,...
            block_3_large_choice_percent,... 
            block_1_small_choice_percent,...
            block_2_small_choice_percent,...
            block_3_small_choice_percent,...
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
            ]);

        variable_names = [...
            %"animalID",...
            "block_1_large",...
            "block_2_large",...
            "block_3_large",...
            "block_1_small",...
            "block_2_small",...
            "block_3_small",...
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
% some mice have NaNs if they didn't make it to this trial block. replace
% the NaNs with 0 because not making it to the trial block is basically
% being as risk averse as possible. 
% risk_table{:, :}(isnan(risk_table{:, :})) = 0;
row_means = nanmean(risk_table{:, 2:4}, 2);


risk_table.Mean_1_to_3 = row_means;
riskiness = risk_table.Mean_1_to_3;
% aborts = risk_table.Var4;
% lose_shift = risk_table.Var5;

mean_riskiness = mean([risk_table.block_1_large risk_table.block_2_large risk_table.block_3_large], 2);

median_riskiness = median(mean_riskiness);


risk_table.risky = mean_riskiness > median_riskiness;
risk_table.risky = double(risk_table.risky); % Convert logical to numeric (0 or 1)

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

    % Add a new column with the treatment condition
    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, hM4Di_IDs);

    % Create a new column with the corresponding treatment conditions
    risk_table_trimmed.TreatmentCondition = hM4Di_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');
elseif exist('stGtACR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = stGtACR_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);

    % Add a new column with the treatment condition
    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, stGtACR_IDs);

    % Create a new column with the corresponding treatment conditions
    risk_table_trimmed.TreatmentCondition = stGtACR_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');

elseif exist('PdCO_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = PdCO_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);

    % Add a new column with the treatment condition
    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, PdCO_IDs);

    % Create a new column with the corresponding treatment conditions
    risk_table_trimmed.TreatmentCondition = PdCO_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');


elseif exist('ChrimsonR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = ChrimsonR_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);

    % Add a new column with the treatment condition
    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, ChrimsonR_IDs);

    % Create a new column with the corresponding treatment conditions
    risk_table_trimmed.TreatmentCondition = ChrimsonR_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');


elseif exist('NAcSh_eNpHR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = mouse_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);

    % Add a new column with the treatment condition
    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, mouse_IDs);

    % Create a new column with the corresponding treatment conditions
    risk_table_trimmed.TreatmentCondition = NAcSh_eNpHR_treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');

elseif exist('treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = mouse_IDs(valid_mice);
    risk_table_trimmed = risk_table(ismember(risk_table.valid_animalIDs, valid_animalIDs), :);
    
        
    % Add a new column with the treatment condition
    [isMatch, matchIdx] = ismember(risk_table_trimmed.valid_animalIDs, mouse_IDs);

    % Create a new column with the corresponding treatment conditions
    risk_table_trimmed.TreatmentCondition = treatment_groups(matchIdx);
    risk_table = risk_table_trimmed;

    risk_table_sorted = sortrows(risk_table, 'TreatmentCondition');
end

% Plot the raw data in grey with transparency

figure;
% Define a colormap for unique colors
% colors = lines(size(large_sequences_mouse, 1)); % Generates a unique color for each trial
% 
% % Plot each trial with a unique color
% for trial = 1:size(large_sequences_mouse, 1)
%     plot(large_sequences_mouse(trial, :), 'Color', colors(trial, :));
%     hold on;
% end
% Plot the mean as a thick black line
meanData_large = mean(large_sequences_mouse);
plot(meanData_large, 'b', 'LineWidth', 2, 'Color', 'b');
hold on
meanData_small = mean(small_sequences_mouse);
plot(meanData_small, 'r', 'LineWidth', 2, 'Color', 'r');

ylim([-0.1 1.1]);
% xlim([-8 8]);
% Set X-axis ticks

xline(0)
yline(0)
fontsize(18, 'points')
hold off;


% Get the number of mice (rows)
numMice = size(large_sequences_mouse, 1);

% Create the figure and set its size
figure('Units', 'normalized', 'Position', [0.1 0.1 0.3 1]); % 1 column, height 3x width

% Loop through each mouse and create a subplot
for mouse = 1:numMice
    subplot(numMice, 1, mouse); % Create subplot
    
    % Plot large_sequences in red
    plot(large_sequences_mouse(mouse, :), 'b', 'LineWidth', 1.5); 
    hold on;
    
    % Plot small_sequences in blue
    plot(small_sequences_mouse(mouse, :), 'r', 'LineWidth', 1.5);
    
    % Add labels and a title
    title(['Mouse ', num2str(mouse)]);
    
    % Customize axis limits for consistency
    ylim([-0.1 1.1]);
    xlim([1 size(large_sequences_mouse, 2)]);
    xline(0, '--k'); % Optional: Add x=0 line
    yline(0, '--k'); % Optional: Add y=0 line
    
    % Remove x-axis labels for all but the bottom subplot
    if mouse < numMice
        set(gca, 'XTickLabel', []);
    end
    hold off;
end

% Add a global label for x-axis
xlabel('Trial');


% IF DOING PR ANALYSIS, DO NOT RUN MOST OF THE CODE ABOVE - ONLY RUN
% SECTION BELOW WITH FILTER ABOVE SET TO PR_D1
for ii = 1:size(valid_animalIDs,1) % 1:size(fieldnames(final),1)
    currentanimal = char(valid_animalIDs(ii));
    if isfield(final_behavior.(currentanimal), 'PR_D1')
        BehavData_PR = final_behavior.(currentanimal).PR_D1.uv.BehavData;
        total_PR_presses(ii) = size(BehavData_PR ,1);
        total_PR_rewards(ii) = sum(BehavData_PR.collectTrial == 1);
        PR_final_ratio_completed(ii) = max(BehavData_PR.Ratio - 1);
        PR_final_ratio_reached(ii) = max(BehavData_PR.Ratio);
    
    end



end

%%

% Calculate means and SEMs
mean_1_3 = table2array(mean(risk_table(:, 1:3), 1));
sem_1_3 = table2array(std(risk_table(:, 1:3), 0, 1) ./ sqrt(size(risk_table(:, 1:3), 1)));

mean_4_6 = table2array(mean(risk_table(:, 4:6), 1));
sem_4_6 = table2array(std(risk_table(:, 4:6), 0, 1) ./ sqrt(size(risk_table(:, 4:6), 1)));

% X-axis points
x_points = 1:size(mean_1_3, 2);

% Plotting
figure;
hold on;

% Plot lines for risk_table(:, 1:3) and risk_table(:, 4:6)
plot(mean_1_3, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'DisplayName', 'Risk 1:3');
plot(mean_4_6, 's-', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'DisplayName', 'Risk 4:6');

% Add error bars manually using "line"
for i = 1:length(x_points)
    % Error bars for risk_table(:, 1:3)
    line([x_points(i), x_points(i)], [mean_1_3(i) - sem_1_3(i), mean_1_3(i) + sem_1_3(i)], ...
        'Color', 'blue', 'LineWidth', 1.2);

    % Error bars for risk_table(:, 4:6)
    line([x_points(i), x_points(i)], [mean_4_6(i) - sem_4_6(i), mean_4_6(i) + sem_4_6(i)], ...
        'Color', 'red', 'LineWidth', 1.2);
end

% Formatting
ylim([0 1]);
xlabel('Time Points');
ylabel('Mean Value');
legend('Location', 'Best');
title('Mean Risk Values with SEM');
grid on;

%%
mean_large = table2array(mean(risk_table(:, 1:3), 2));
mean_small = table2array(mean(risk_table(:, 4:6), 2));
sem_large = table2array(std(risk_table(:, 1:3), 0, 2) ./ sqrt(size(risk_table(:, 1:3), 2)));
sem_small = table2array(std(risk_table(:, 4:6), 0, 2) ./ sqrt(size(risk_table(:, 4:6), 2)));

% Calculate means for the bar plot
group_means = [mean(mean_large), mean(mean_small)];
group_sems = [mean(sem_large), mean(sem_small)];

% Create a bar plot
figure;
hold on;
bar_handle = bar(group_means, 'BarWidth', 0.5); % Create bar plot with some transparency
errorbar(bar_handle.XEndPoints(1), group_means(1), group_sems(1)); % Create bar plot with some transparency
errorbar(bar_handle.XEndPoints(2), group_means(2), group_sems(2)); % Create bar plot with some transparency


colors = lines(2); % Generate distinct colors for each group

% Add scatter points for each group
x_locations = bar_handle.XEndPoints; % X locations of the bars
scatter_jitter = 0.1; % Jitter width for scatter points

% Scatter points for the first group (mean_large)
scatter(x_locations(1) + (rand(size(mean_large)) - 0.5) * scatter_jitter, mean_large, ...
    50, colors(1, :), 'filled');

% Scatter points for the second group (mean_small)
scatter(x_locations(2) + (rand(size(mean_small)) - 0.5) * scatter_jitter, mean_small, ...
    50, colors(2, :), 'filled');

% Customize plot
set(gca, 'XTick', 1:2, 'XTickLabel', {'Large', 'Small'});
ylabel('Values');
hold off;
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

% Calculate means for the bar plot
cross_sess_large_means = [mean(mean_large_RM_D1), mean(mean_large)]*100;
cross_sess_large_sems = [mean(sem_large_RM_D1), mean(sem_large)]*100;

cross_session_large_all = [mean_large_RM_D1, mean_large]*100;
cross_session_small_all = [mean_small_RM_D1, mean_small]*100;

% Calculate means for the bar plot
cross_sess_small_means = [mean(mean_small_RM_D1), mean(mean_small)]*100;
cross_sess_small_sems = [mean(sem_small_RM_D1), mean(sem_small)]*100;

% X-axis points
x_points = 1:size(cross_sess_large_means, 2);


% Plotting
figure;
hold on;


% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size


% Set figure size
% width = 550; % Width of the figure
% height = 450; % Height of the figure
% set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(cross_session_large_all, 1)
    plot(x_points, cross_session_large_all(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(cross_session_small_all, 1)
    plot(x_points, cross_session_small_all(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, cross_sess_large_means, cross_sess_large_sems, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 18, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, cross_sess_small_means, cross_sess_small_sems, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 18, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'Early RM', 'Late RM'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 1.1 * max([cross_sess_large_means + cross_sess_large_sems, ...
                   cross_sess_small_means + cross_sess_small_sems])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

% 2×2 Repeated Measures ANOVA for behavioral task data
% Factors: Day (Day 1, Day 2) × Reward Size (Large, Small)
% Assuming cross_session_large_all and cross_session_small_all are 10×2 matrices

% Step 1: Prepare data in the format required for repeated measures ANOVA
num_subjects = size(cross_session_large_all, 1);
% Create a table with subject IDs and the 4 condition measurements
subject_id = (1:num_subjects)';
data_wide = table(subject_id, ...
    cross_session_large_all(:,1), cross_session_large_all(:,2), ...
    cross_session_small_all(:,1), cross_session_small_all(:,2), ...
    'VariableNames', {'Subject', 'Day1_Large', 'Day2_Large', 'Day1_Small', 'Day2_Small'});

% Create the within-subjects design table
% This defines the structure of our repeated measures design
withinDesign = table(categorical({'Day1'; 'Day2'; 'Day1'; 'Day2'}), ...
    categorical({'Large'; 'Large'; 'Small'; 'Small'}), ...
    'VariableNames', {'Day', 'RewardSize'});

% Step 2: Conduct the 2×2 repeated measures ANOVA
% Define the repeated measures model - all measures with the same between-subjects model (~ 1)
rm = fitrm(data_wide, 'Day1_Large,Day2_Large,Day1_Small,Day2_Small ~ 1', 'WithinDesign', withinDesign);

% Run the repeated measures ANOVA
ranova_table = ranova(rm, 'WithinModel', 'Day*RewardSize');

% Display results
disp('Repeated Measures ANOVA Results:');
disp(ranova_table);

% Display ANOVA results
fprintf('Reward Size effect (RM ANOVA): F(%d,%d) = %.3f, p = %.3e\n', ...
    ranova_table{5,2}, ranova_table{6,2}, ranova_table{5,4}, ranova_table{5,5});
% Display ANOVA results
fprintf('Day effect (RM ANOVA): F(%d,%d) = %.3f, p = %.3e\n', ...
    ranova_table{3,2}, ranova_table{4,2}, ranova_table{3,4}, ranova_table{3,5});
fprintf('Interaction (RM ANOVA): F(%d,%d) = %.3f, p = %.3e\n', ...
    ranova_table{7,2}, ranova_table{8,2}, ranova_table{7,4}, ranova_table{7,5});

% Check if the interaction effect is significant
p_interaction = ranova_table{7,5}; % The interaction p-value
if p_interaction < 0.05
    disp('Significant interaction found! Performing Tukey''s multiple comparisons...');
    fprintf('\n========== POST-HOC TESTS ==========\n');
    
    % Reshape data for Tukey's multiple comparisons
    % Create long-format data for all Day-RewardSize combinations
    all_values = [];
    all_groups = {};
    
    % Day1 Large
    all_values = [all_values; data_wide.Day1_Large];
    all_groups = [all_groups; repmat({'Day1_Large'}, num_subjects, 1)];
    
    % Day2 Large
    all_values = [all_values; data_wide.Day2_Large];
    all_groups = [all_groups; repmat({'Day2_Large'}, num_subjects, 1)];
    
    % Day1 Small
    all_values = [all_values; data_wide.Day1_Small];
    all_groups = [all_groups; repmat({'Day1_Small'}, num_subjects, 1)];
    
    % Day2 Small
    all_values = [all_values; data_wide.Day2_Small];
    all_groups = [all_groups; repmat({'Day2_Small'}, num_subjects, 1)];
    
    % Use anova1 to get stats structure for multcompare
    [~, ~, stats_tukey] = anova1(all_values, all_groups, 'off');
    
    % Perform Tukey's multiple comparisons
    [c, m, h, gnames] = multcompare(stats_tukey, 'CType', 'tukey-kramer', 'Display', 'off');
    
    % Display Tukey's results in a readable format
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
    
    % Display group means for reference
    fprintf('\nGroup Means:\n');
    for i = 1:length(gnames)
        fprintf('%-12s: %.3f ± %.3f (SEM)\n', gnames{i}, m(i, 1), m(i, 2));
    end
    
    % Highlight specific comparisons of interest
    fprintf('\n--- Key Comparisons ---\n');
    fprintf('Reward Size comparisons within each day:\n');
    
    % Day 1: Large vs Small
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
    
    % Day 2: Large vs Small
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
    
    % Large: Day 1 vs Day 2
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
    
    % Small: Day 1 vs Day 2
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
    
    % Check main effects if interaction is not significant
    p_reward_size = ranova_table{5,5}; % p-value for RewardSize main effect
    p_day = ranova_table{3,5};  % p-value for Day main effect
    
    % Check main effect of Reward Size
    if p_reward_size < 0.05
        fprintf('\n--- Main effect of Reward Size is significant (p = %.3e) ---\n', p_reward_size);
        % Compare overall means (Large vs Small)
        large_mean = mean([data_wide.Day1_Large; data_wide.Day2_Large]);
        small_mean = mean([data_wide.Day1_Small; data_wide.Day2_Small]);
        fprintf('Large Reward mean: %.3f\n', large_mean);
        fprintf('Small Reward mean: %.3f\n', small_mean);
        fprintf('Difference: %.3f\n', large_mean - small_mean);
    end
    
    % Check main effect of Day
    if p_day < 0.05
        fprintf('\n--- Main effect of Day is significant (p = %.3e) ---\n', p_day);
        % Compare overall means (Day 1 vs Day 2)
        day1_mean = mean([data_wide.Day1_Large; data_wide.Day1_Small]);
        day2_mean = mean([data_wide.Day2_Large; data_wide.Day2_Small]);
        fprintf('Day 1 mean: %.3f\n', day1_mean);
        fprintf('Day 2 mean: %.3f\n', day2_mean);
        fprintf('Difference: %.3f\n', day1_mean - day2_mean);
        
        % Run pairwise comparisons for Day using multcompare
        dayComp = multcompare(rm, 'Day', 'ComparisonType', 'tukey-kramer');
        
        % Display results
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
mean_latency_large = mean([risk_table.block_1_large_choice_latency, risk_table.block_2_large_choice_latency, risk_table.block_3_large_choice_latency], 2)
mean_latency_small = mean([risk_table.block_1_small_choice_latency, risk_table.block_2_small_choice_latency, risk_table.block_3_small_choice_latency], 2)

latencies_for_boxplot = [mean_latency_large mean_latency_small]; 

mean_collect_latency_large = mean([risk_table.block_1_large_collect_latency, risk_table.block_2_large_collect_latency, risk_table.block_3_large_collect_latency], 2)
mean_collect_latency_small = mean([risk_table.block_1_small_collect_latency, risk_table.block_2_small_collect_latency, risk_table.block_3_small_collect_latency], 2)

collect_latencies_for_boxplot = [mean_collect_latency_large mean_collect_latency_small]; 

mean_consum_latency_large = mean([risk_table.large_consum_duration_block_1, risk_table.large_consum_duration_block_2, risk_table.large_consum_duration_block_3], 2)
mean_consum_latency_small = mean([risk_table.small_consum_duration_block_1, risk_table.small_consum_duration_block_2, risk_table.small_consum_duration_block_3], 2)

consum_latencies_for_boxplot = [mean_consum_latency_large mean_consum_latency_small]; 

hold on
figure;
tiledlayout(3, 1)

% First boxplot
nexttile;
boxplot(latencies_for_boxplot, 'Orientation', 'horizontal');
xlim([0 10]);
xticks([0 10]); % Set x-ticks to the first and last values
xticklabels({'0', '10'}); % Label the x-ticks

% Second boxplot
nexttile;
boxplot(collect_latencies_for_boxplot, 'Orientation', 'horizontal');
xlim([0 5]);
xticks([0 5]); % Set x-ticks to the first and last values
xticklabels({'0', '5'}); % Label the x-ticks

% Third boxplot
nexttile;
boxplot(consum_latencies_for_boxplot, 'Orientation', 'horizontal');
xlim([0 5]);
xticks([0 5]); % Set x-ticks to the first and last values
xticklabels({'0', '5'}); % Label the x-ticks

hold off

%%

large_choice = [risk_table.block_1_large, risk_table.block_2_large, risk_table.block_3_large]*100;
small_choice = [risk_table.block_1_small, risk_table.block_2_small, risk_table.block_3_small]*100;

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));







% X-axis points
x_points = 1:size(large_choice, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 1.1 * max([mean_large + sem_large, ...
                   mean_small + sem_small])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice = [risk_table.block_1_large_choice_latency, risk_table.block_2_large_choice_latency, risk_table.block_3_large_choice_latency];
small_choice = [risk_table.block_1_small_choice_latency, risk_table.block_2_small_choice_latency, risk_table.block_3_small_choice_latency];

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));







% X-axis points
x_points = 1:size(large_choice, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
% ylim([0 1.1 * max([mean_large + sem_large, ...
%                    mean_small + sem_small])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:20);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

% large_choice = [risk_table.block_1_large_collect_latency_no_shk_trials, risk_table.block_2_large_collect_latency_no_shk_trials, risk_table.block_3_large_collect_latency_no_shk_trials];
% small_choice = [risk_table.block_1_small_collect_latency, risk_table.block_2_small_collect_latency, risk_table.block_3_small_collect_latency];

large_choice = [risk_table.block_1_large_collect_latency, risk_table.block_2_large_collect_latency, risk_table.block_3_large_collect_latency];
small_choice = [risk_table.block_1_small_collect_latency, risk_table.block_2_small_collect_latency, risk_table.block_3_small_collect_latency];


mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));







% X-axis points
x_points = 1:size(large_choice, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
% ylim([0 1.1 * max([mean_large + sem_large, ...
%                    mean_small + sem_small])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:1:5);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice = [risk_table.large_consum_duration_block_1, risk_table.large_consum_duration_block_2, risk_table.large_consum_duration_block_3];
small_choice = [risk_table.small_consum_duration_block_1, risk_table.small_consum_duration_block_2, risk_table.small_consum_duration_block_3];

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));







% X-axis points
x_points = 1:size(large_choice, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
% ylim([0 1.1 * max([mean_large + sem_large, ...
%                    mean_small + sem_small])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:1:5);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice = [risk_table.large_aborts_block_1, risk_table.large_aborts_block_2, risk_table.large_aborts_block_3];
small_choice = [risk_table.small_aborts_block_1, risk_table.small_aborts_block_2, risk_table.small_aborts_block_3];

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));







% X-axis points
x_points = 1:size(large_choice, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
% ylim([0 1.1 * max([mean_large + sem_large, ...
%                    mean_small + sem_small])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;




%% for plotting lose_shift ratio across blocks

large_choice = [risk_table.lose_shift_ratio_block_1, risk_table.lose_shift_ratio_block_2, risk_table.lose_shift_ratio_block_3];
small_choice = [risk_table.win_stay_ratio_block_1, risk_table.win_stay_ratio_block_2, risk_table.win_stay_ratio_block_3];

% turn NaNs (trials where the denominator, number of losses, is 0) into 0s
large_choice(isnan(large_choice))=0;
small_choice(isnan(small_choice))=0;

mean_large = nanmean(large_choice, 1);
mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));







% X-axis points
x_points = 1:size(large_choice, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(small_choice, 1)
    plot(x_points, small_choice(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
% ylim([0 1.1 * max([mean_large + sem_large, ...
%                    mean_small + sem_small])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:1);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%%

large_choice = [risk_table.block_1_omission_total, risk_table.block_2_omission_total, risk_table.block_3_omission_total];
% small_choice = [risk_table.block_1_small, risk_table.block_2_small, risk_table.block_3_small]*100;

mean_large = nanmean(large_choice, 1);
% mean_small = nanmean(small_choice, 1);
sem_large = nanstd(large_choice, 0, 1) ./ sqrt(size(large_choice, 1));
% sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));







% X-axis points
x_points = 1:size(large_choice, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice, 1)
    plot(x_points, large_choice(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% % Plot individual lines for "Small" data
% for i = 1:size(small_choice, 1)
%     plot(x_points, small_choice(i, :), '-', ...
%         'Color', [1 0 0 0.6], ... % Red with 60% opacity
%         'LineWidth', 1.2);
% end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% errorbar(x_points, mean_small, sem_small, '^-', ...
%     'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
%     'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
% ylim([0 1.1 * max([mean_large + sem_large])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:25);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_large(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('hM4Di', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_small(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('hM4Di', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'
% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_aborts_block_1(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('hM4Di', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
% ylim([0 1.1 * max([mean_large + sem_large, ...
%                    mean_small + sem_small])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;



%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_choice_latency(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('hM4Di', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 15]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for hM4Di vs mCherry

large_choice_mCherry = [risk_table.block_1_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_collect_latency(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('hM4Di', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 15]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%
large_choice_mCherry = [risk_table.block_1_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_omission_total(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('hM4Di', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('hM4Di', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 250; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% % Plot individual lines for "Large" data
% for i = 1:size(large_choice_mCherry, 1)
%     plot(x_points, large_choice_mCherry(i, :), '-', ...
%         'Color', mCherry_color, ... % Blue with 60% opacity
%         'LineWidth', 1.2);
% end
% 
% % Plot individual lines for "Small" data
% for i = 1:size(large_choice_hM4Di, 1)
%     plot(x_points, large_choice_hM4Di(i, :), '-', ...
%         'Color', hM4Di_color, ... % Red with 60% opacity
%         'LineWidth', 1.2);
% end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 15]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;
%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.block_1_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_large(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('stGtACR', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.block_1_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_small(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('stGtACR', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.block_1_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_choice_latency(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('stGtACR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.block_1_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_collect_latency(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('stGtACR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for stGtACR vs mCherry

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_aborts_block_1(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('stGtACR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('stGtACR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 250; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 125]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:125);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for stGtACR vs mCherry

block_1_abort_sum = risk_table.large_aborts_block_1 + risk_table.small_aborts_block_1;
block_2_abort_sum = risk_table.large_aborts_block_2 + risk_table.small_aborts_block_2;
block_3_abort_sum = risk_table.large_aborts_block_3 + risk_table.small_aborts_block_3;

large_choice_mCherry = [block_1_abort_sum(strcmp('mCherry', risk_table.TreatmentCondition)), block_2_abort_sum(strcmp('mCherry', risk_table.TreatmentCondition)), block_3_abort_sum(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [block_1_abort_sum(strcmp('stGtACR', risk_table.TreatmentCondition)), block_2_abort_sum(strcmp('stGtACR', risk_table.TreatmentCondition)), block_3_abort_sum(strcmp('stGtACR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 125]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:125);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%
%% for PdCO vs mCherry

large_choice_mCherry = [risk_table.block_1_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_large(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('PdCO', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 100; % Width of the figure
height = 250; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size
fontsize(12, 'points')
% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, 'MarkerEdgeColor', 'none', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for PdCO vs mCherry

large_choice_mCherry = [risk_table.block_1_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_small(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('PdCO', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 100; % Width of the figure
height = 250; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size
fontsize(12, 'points')
% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, 'MarkerEdgeColor', 'none', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for PdCO vs mCherry

large_choice_mCherry = [risk_table.block_1_choice_latency_all(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_choice_latency_all(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_choice_latency_all(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_choice_latency_all(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_choice_latency_all(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_choice_latency_all(strcmp('PdCO', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 100; % Width of the figure
height = 250; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size
fontsize(12, 'points')
% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, 'MarkerEdgeColor', 'none', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for PdCO vs mCherry

large_choice_mCherry = [risk_table.block_1_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_omission_total(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('PdCO', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 250; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size
fontsize(12, 'points')
% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, 'MarkerEdgeColor', 'none', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;





%% for PdCO vs mCherry

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_aborts_block_1(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('PdCO', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 100; % Width of the figure
height = 250; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size
fontsize(12, 'points')
% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, 'MarkerEdgeColor', 'none', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 125]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:125);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;



%% for BLA-NAcSh ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_large(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('ChrimsonR', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for BLA-NAcSh ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('mCherry', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_small(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('ChrimsonR', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% 

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_aborts_block_1(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 125]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:125);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice_mCherry = [risk_table.block_1_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_choice_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice_mCherry = [risk_table.block_1_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_omission_total(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 60]); % Adjust ylim dynamically
set(gca, 'ytick', 0:10:60);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%%

large_choice_mCherry = [risk_table.block_1_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_collect_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for PdCO vs mCherry

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_aborts_block_1(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 125]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:125);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice = [risk_table.block_1_large, risk_table.block_2_large, risk_table.block_3_large]*100;
large_choice_risky = large_choice(risk_table.risky == 1, :);
large_choice_not_risky = large_choice(risk_table.risky == 0, :);
small_choice = [risk_table.block_1_small, risk_table.block_2_small, risk_table.block_3_small]*100;

mean_large_risky = nanmean(large_choice_risky, 1);
mean_large_not_risky = nanmean(large_choice_not_risky, 1);
mean_small = nanmean(small_choice, 1);

sem_large_risky = nanstd(large_choice_risky, 0, 1) ./ sqrt(size(large_choice_risky, 1));
sem_large_not_risky = nanstd(large_choice_not_risky, 0, 1) ./ sqrt(size(large_choice_not_risky, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));







% X-axis points
x_points = 1:size(large_choice_risky, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_risky, 1)
    plot(x_points, large_choice_risky(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_not_risky, 1)
    plot(x_points, large_choice_not_risky(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large_risky, sem_large_risky, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_large_not_risky, sem_large_not_risky, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 1.1 * max([mean_large_risky + sem_large_risky, ...
                   mean_large_not_risky + sem_large_not_risky])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice = [risk_table.large_aborts_block_1, risk_table.large_aborts_block_2, risk_table.large_aborts_block_3];
large_choice_risky = large_choice(risk_table.risky == 1, :);
large_choice_not_risky = large_choice(risk_table.risky == 0, :);
small_choice = [risk_table.block_1_small, risk_table.block_2_small, risk_table.block_3_small];

mean_large_risky = nanmean(large_choice_risky, 1);
mean_large_not_risky = nanmean(large_choice_not_risky, 1);
mean_small = nanmean(small_choice, 1);

sem_large_risky = nanstd(large_choice_risky, 0, 1) ./ sqrt(size(large_choice_risky, 1));
sem_large_not_risky = nanstd(large_choice_not_risky, 0, 1) ./ sqrt(size(large_choice_not_risky, 1));
sem_small = nanstd(small_choice, 0, 1) ./ sqrt(size(small_choice, 1));







% X-axis points
x_points = 1:size(large_choice_risky, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_risky, 1)
    plot(x_points, large_choice_risky(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_not_risky, 1)
    plot(x_points, large_choice_not_risky(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large_risky, sem_large_risky, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_large_not_risky, sem_large_not_risky, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 1.1 * max([mean_large_risky + sem_large_risky, ...
                   mean_large_not_risky + sem_large_not_risky])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%
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






% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

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






% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

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






% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

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






% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for PdCO vs ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_PdCO = [risk_table.block_1_omission_total(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('PdCO', risk_table.TreatmentCondition))];
large_choice_ChrimsonR = [risk_table.block_1_omission_total(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));






% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 60]); % Adjust ylim dynamically
set(gca, 'ytick', 0:10:60);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

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






% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 125]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:125);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for PdCO vs ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.small_aborts_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.small_aborts_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.small_aborts_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_PdCO = [risk_table.small_aborts_block_1(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.small_aborts_block_2(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.small_aborts_block_3(strcmp('PdCO', risk_table.TreatmentCondition))];
large_choice_ChrimsonR = [risk_table.small_aborts_block_1(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.small_aborts_block_2(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.small_aborts_block_3(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));






% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 125]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:125);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%
large_choice_mCherry = [risk_table.large_consum_duration_block_1(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_PdCO = [risk_table.large_consum_duration_block_1(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('PdCO', risk_table.TreatmentCondition))];
large_choice_ChrimsonR = [risk_table.large_consum_duration_block_1(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));






% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 10]); % Adjust ylim dynamically
set(gca, 'ytick', 0:10:10);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for PdCO vs ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_large_collect_latency(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('mCherry', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('mCherry', risk_table.TreatmentCondition))];
large_choice_PdCO = [risk_table.block_1_large_collect_latency(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('PdCO', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('PdCO', risk_table.TreatmentCondition))];
large_choice_ChrimsonR = [risk_table.block_1_large_collect_latency(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('ChrimsonR', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('ChrimsonR', risk_table.TreatmentCondition))];

mean_mCherry = nanmean(large_choice_mCherry, 1);
mean_PdCO = nanmean(large_choice_PdCO, 1);
mean_ChrimsonR = nanmean(large_choice_ChrimsonR, 1);

sem_mCherry = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_PdCO = nanstd(large_choice_PdCO, 0, 1) ./ sqrt(size(large_choice_PdCO, 1));
sem_ChrimsonR = nanstd(large_choice_ChrimsonR, 0, 1) ./ sqrt(size(large_choice_ChrimsonR, 1));






% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_PdCO, 1)
    plot(x_points, large_choice_PdCO(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_ChrimsonR, 1)
    plot(x_points, large_choice_ChrimsonR(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_mCherry, sem_mCherry, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_PdCO, sem_PdCO, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_ChrimsonR, sem_ChrimsonR, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for A2A vs D1 vs Control

large_choice_Control = [risk_table.block_1_large(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('Control', risk_table.TreatmentCondition))]*100;
large_choice_A2A = [risk_table.block_1_large(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('A2A', risk_table.TreatmentCondition))]*100;
large_choice_D1 = [risk_table.block_1_large(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('D1', risk_table.TreatmentCondition))]*100;

mean_Control = nanmean(large_choice_Control, 1);
mean_A2A = nanmean(large_choice_A2A, 1);
mean_D1 = nanmean(large_choice_D1, 1);

sem_Control = nanstd(large_choice_Control, 0, 1) ./ sqrt(size(large_choice_Control, 1));
sem_A2A = nanstd(large_choice_A2A, 0, 1) ./ sqrt(size(large_choice_A2A, 1));
sem_D1 = nanstd(large_choice_D1, 0, 1) ./ sqrt(size(large_choice_D1, 1));






% X-axis points
x_points = 1:size(large_choice_Control, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_Control, 1)
    plot(x_points, large_choice_Control(i, :), '-', ...
        'Color', Control_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_A2A, 1)
    plot(x_points, large_choice_A2A(i, :), '-', ...
        'Color', A2A_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_D1, 1)
    plot(x_points, large_choice_D1(i, :), '-', ...
        'Color', D1_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_Control, sem_Control, Control_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', Control_color, 'MarkerFaceColor', Control_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_A2A, sem_A2A, A2A_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', A2A_color, 'MarkerFaceColor', A2A_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_D1, sem_D1, D1_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', D1_color, 'MarkerFaceColor', D1_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for A2A vs D1 vs Control

large_choice_Control = [risk_table.block_1_small(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('Control', risk_table.TreatmentCondition))]*100;
large_choice_A2A = [risk_table.block_1_small(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('A2A', risk_table.TreatmentCondition))]*100;
large_choice_D1 = [risk_table.block_1_small(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('D1', risk_table.TreatmentCondition))]*100;

mean_Control = nanmean(large_choice_Control, 1);
mean_A2A = nanmean(large_choice_A2A, 1);
mean_D1 = nanmean(large_choice_D1, 1);

sem_Control = nanstd(large_choice_Control, 0, 1) ./ sqrt(size(large_choice_Control, 1));
sem_A2A = nanstd(large_choice_A2A, 0, 1) ./ sqrt(size(large_choice_A2A, 1));
sem_D1 = nanstd(large_choice_D1, 0, 1) ./ sqrt(size(large_choice_D1, 1));






% X-axis points
x_points = 1:size(large_choice_Control, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_Control, 1)
    plot(x_points, large_choice_Control(i, :), '-', ...
        'Color', Control_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_A2A, 1)
    plot(x_points, large_choice_A2A(i, :), '-', ...
        'Color', A2A_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_D1, 1)
    plot(x_points, large_choice_D1(i, :), '-', ...
        'Color', D1_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_Control, sem_Control, Control_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', Control_color, 'MarkerFaceColor', Control_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_A2A, sem_A2A, A2A_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', A2A_color, 'MarkerFaceColor', A2A_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_D1, sem_D1, D1_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', D1_color, 'MarkerFaceColor', D1_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for A2A vs D1 vs Control

large_choice_Control = [risk_table.block_1_choice_latency_all(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_2_choice_latency_all(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_3_choice_latency_all(strcmp('Control', risk_table.TreatmentCondition))];
large_choice_A2A = [risk_table.block_1_choice_latency_all(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_2_choice_latency_all(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_3_choice_latency_all(strcmp('A2A', risk_table.TreatmentCondition))];
large_choice_D1 = [risk_table.block_1_choice_latency_all(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_2_choice_latency_all(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_3_choice_latency_all(strcmp('D1', risk_table.TreatmentCondition))];

mean_Control = nanmean(large_choice_Control, 1);
mean_A2A = nanmean(large_choice_A2A, 1);
mean_D1 = nanmean(large_choice_D1, 1);

sem_Control = nanstd(large_choice_Control, 0, 1) ./ sqrt(size(large_choice_Control, 1));
sem_A2A = nanstd(large_choice_A2A, 0, 1) ./ sqrt(size(large_choice_A2A, 1));
sem_D1 = nanstd(large_choice_D1, 0, 1) ./ sqrt(size(large_choice_D1, 1));






% X-axis points
x_points = 1:size(large_choice_Control, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_Control, 1)
    plot(x_points, large_choice_Control(i, :), '-', ...
        'Color', Control_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_A2A, 1)
    plot(x_points, large_choice_A2A(i, :), '-', ...
        'Color', A2A_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_D1, 1)
    plot(x_points, large_choice_D1(i, :), '-', ...
        'Color', D1_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_Control, sem_Control, Control_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', Control_color, 'MarkerFaceColor', Control_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_A2A, sem_A2A, A2A_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', A2A_color, 'MarkerFaceColor', A2A_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_D1, sem_D1, D1_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', D1_color, 'MarkerFaceColor', D1_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for A2A vs D1 vs Control

large_choice_Control = [risk_table.block_1_large_collect_latency(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('Control', risk_table.TreatmentCondition))];
large_choice_A2A = [risk_table.block_1_large_collect_latency(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('A2A', risk_table.TreatmentCondition))];
large_choice_D1 = [risk_table.block_1_large_collect_latency(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('D1', risk_table.TreatmentCondition))];

mean_Control = nanmean(large_choice_Control, 1);
mean_A2A = nanmean(large_choice_A2A, 1);
mean_D1 = nanmean(large_choice_D1, 1);

sem_Control = nanstd(large_choice_Control, 0, 1) ./ sqrt(size(large_choice_Control, 1));
sem_A2A = nanstd(large_choice_A2A, 0, 1) ./ sqrt(size(large_choice_A2A, 1));
sem_D1 = nanstd(large_choice_D1, 0, 1) ./ sqrt(size(large_choice_D1, 1));






% X-axis points
x_points = 1:size(large_choice_Control, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_Control, 1)
    plot(x_points, large_choice_Control(i, :), '-', ...
        'Color', Control_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_A2A, 1)
    plot(x_points, large_choice_A2A(i, :), '-', ...
        'Color', A2A_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_D1, 1)
    plot(x_points, large_choice_D1(i, :), '-', ...
        'Color', D1_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_Control, sem_Control, Control_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', Control_color, 'MarkerFaceColor', Control_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_A2A, sem_A2A, A2A_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', A2A_color, 'MarkerFaceColor', A2A_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_D1, sem_D1, D1_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', D1_color, 'MarkerFaceColor', D1_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for A2A vs D1 vs Control

large_choice_Control = [risk_table.block_1_omission_total(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('Control', risk_table.TreatmentCondition))];
large_choice_A2A = [risk_table.block_1_omission_total(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('A2A', risk_table.TreatmentCondition))];
large_choice_D1 = [risk_table.block_1_omission_total(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('D1', risk_table.TreatmentCondition))];

mean_Control = nanmean(large_choice_Control, 1);
mean_A2A = nanmean(large_choice_A2A, 1);
mean_D1 = nanmean(large_choice_D1, 1);

sem_Control = nanstd(large_choice_Control, 0, 1) ./ sqrt(size(large_choice_Control, 1));
sem_A2A = nanstd(large_choice_A2A, 0, 1) ./ sqrt(size(large_choice_A2A, 1));
sem_D1 = nanstd(large_choice_D1, 0, 1) ./ sqrt(size(large_choice_D1, 1));






% X-axis points
x_points = 1:size(large_choice_Control, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_Control, 1)
    plot(x_points, large_choice_Control(i, :), '-', ...
        'Color', Control_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_A2A, 1)
    plot(x_points, large_choice_A2A(i, :), '-', ...
        'Color', A2A_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_D1, 1)
    plot(x_points, large_choice_D1(i, :), '-', ...
        'Color', D1_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_Control, sem_Control, Control_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', Control_color, 'MarkerFaceColor', Control_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_A2A, sem_A2A, A2A_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', A2A_color, 'MarkerFaceColor', A2A_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_D1, sem_D1, D1_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', D1_color, 'MarkerFaceColor', D1_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 45]); % Adjust ylim dynamically
set(gca, 'ytick', 0:10:40);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;
%% for A2A vs D1 vs Control

large_choice_Control = [risk_table.large_aborts_block_1(strcmp('Control', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('Control', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('Control', risk_table.TreatmentCondition))];
large_choice_A2A = [risk_table.large_aborts_block_1(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('A2A', risk_table.TreatmentCondition))];
large_choice_D1 = [risk_table.large_aborts_block_1(strcmp('D1', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('D1', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('D1', risk_table.TreatmentCondition))];

mean_Control = nanmean(large_choice_Control, 1);
mean_A2A = nanmean(large_choice_A2A, 1);
mean_D1 = nanmean(large_choice_D1, 1);

sem_Control = nanstd(large_choice_Control, 0, 1) ./ sqrt(size(large_choice_Control, 1));
sem_A2A = nanstd(large_choice_A2A, 0, 1) ./ sqrt(size(large_choice_A2A, 1));
sem_D1 = nanstd(large_choice_D1, 0, 1) ./ sqrt(size(large_choice_D1, 1));






% X-axis points
x_points = 1:size(large_choice_Control, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_Control, 1)
    plot(x_points, large_choice_Control(i, :), '-', ...
        'Color', Control_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_A2A, 1)
    plot(x_points, large_choice_A2A(i, :), '-', ...
        'Color', A2A_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_D1, 1)
    plot(x_points, large_choice_D1(i, :), '-', ...
        'Color', D1_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_Control, sem_Control, Control_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', Control_color, 'MarkerFaceColor', Control_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_A2A, sem_A2A, A2A_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', A2A_color, 'MarkerFaceColor', A2A_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_D1, sem_D1, D1_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', D1_color, 'MarkerFaceColor', D1_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 90]); % Adjust ylim dynamically
set(gca, 'ytick', 0:30:90);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for A2A vs D1 vs Control

large_choice_Control = [risk_table.small_aborts_block_1(strcmp('Control', risk_table.TreatmentCondition)), risk_table.small_aborts_block_2(strcmp('Control', risk_table.TreatmentCondition)), risk_table.small_aborts_block_3(strcmp('Control', risk_table.TreatmentCondition))];
large_choice_A2A = [risk_table.small_aborts_block_1(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.small_aborts_block_2(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.small_aborts_block_3(strcmp('A2A', risk_table.TreatmentCondition))];
large_choice_D1 = [risk_table.small_aborts_block_1(strcmp('D1', risk_table.TreatmentCondition)), risk_table.small_aborts_block_2(strcmp('D1', risk_table.TreatmentCondition)), risk_table.small_aborts_block_3(strcmp('D1', risk_table.TreatmentCondition))];

mean_Control = nanmean(large_choice_Control, 1);
mean_A2A = nanmean(large_choice_A2A, 1);
mean_D1 = nanmean(large_choice_D1, 1);

sem_Control = nanstd(large_choice_Control, 0, 1) ./ sqrt(size(large_choice_Control, 1));
sem_A2A = nanstd(large_choice_A2A, 0, 1) ./ sqrt(size(large_choice_A2A, 1));
sem_D1 = nanstd(large_choice_D1, 0, 1) ./ sqrt(size(large_choice_D1, 1));






% X-axis points
x_points = 1:size(large_choice_Control, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_Control, 1)
    plot(x_points, large_choice_Control(i, :), '-', ...
        'Color', Control_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_A2A, 1)
    plot(x_points, large_choice_A2A(i, :), '-', ...
        'Color', A2A_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_D1, 1)
    plot(x_points, large_choice_D1(i, :), '-', ...
        'Color', D1_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_Control, sem_Control, Control_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', Control_color, 'MarkerFaceColor', Control_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_A2A, sem_A2A, A2A_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', A2A_color, 'MarkerFaceColor', A2A_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_D1, sem_D1, D1_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', D1_color, 'MarkerFaceColor', D1_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 90]); % Adjust ylim dynamically
set(gca, 'ytick', 0:30:90);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%
large_choice_Control = [risk_table.large_consum_duration_block_1(strcmp('Control', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('Control', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('Control', risk_table.TreatmentCondition))];
large_choice_A2A = [risk_table.large_consum_duration_block_1(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('A2A', risk_table.TreatmentCondition))];
large_choice_D1 = [risk_table.large_consum_duration_block_1(strcmp('D1', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('D1', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('D1', risk_table.TreatmentCondition))];

mean_Control = nanmean(large_choice_Control, 1);
mean_A2A = nanmean(large_choice_A2A, 1);
mean_D1 = nanmean(large_choice_D1, 1);

sem_Control = nanstd(large_choice_Control, 0, 1) ./ sqrt(size(large_choice_Control, 1));
sem_A2A = nanstd(large_choice_A2A, 0, 1) ./ sqrt(size(large_choice_A2A, 1));
sem_D1 = nanstd(large_choice_D1, 0, 1) ./ sqrt(size(large_choice_D1, 1));






% X-axis points
x_points = 1:size(large_choice_Control, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_Control, 1)
    plot(x_points, large_choice_Control(i, :), '-', ...
        'Color', Control_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_A2A, 1)
    plot(x_points, large_choice_A2A(i, :), '-', ...
        'Color', A2A_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_D1, 1)
    plot(x_points, large_choice_D1(i, :), '-', ...
        'Color', D1_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_Control, sem_Control, Control_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', Control_color, 'MarkerFaceColor', Control_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_A2A, sem_A2A, A2A_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', A2A_color, 'MarkerFaceColor', A2A_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_D1, sem_D1, D1_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', D1_color, 'MarkerFaceColor', D1_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 10]); % Adjust ylim dynamically
set(gca, 'ytick', 0:10:10);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for A2A vs D1 vs Control

large_choice_Control = [risk_table.block_1_large_collect_latency(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('Control', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('Control', risk_table.TreatmentCondition))];
large_choice_A2A = [risk_table.block_1_large_collect_latency(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('A2A', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('A2A', risk_table.TreatmentCondition))];
large_choice_D1 = [risk_table.block_1_large_collect_latency(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('D1', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('D1', risk_table.TreatmentCondition))];

mean_Control = nanmean(large_choice_Control, 1);
mean_A2A = nanmean(large_choice_A2A, 1);
mean_D1 = nanmean(large_choice_D1, 1);

sem_Control = nanstd(large_choice_Control, 0, 1) ./ sqrt(size(large_choice_Control, 1));
sem_A2A = nanstd(large_choice_A2A, 0, 1) ./ sqrt(size(large_choice_A2A, 1));
sem_D1 = nanstd(large_choice_D1, 0, 1) ./ sqrt(size(large_choice_D1, 1));






% X-axis points
x_points = 1:size(large_choice_Control, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_Control, 1)
    plot(x_points, large_choice_Control(i, :), '-', ...
        'Color', Control_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_A2A, 1)
    plot(x_points, large_choice_A2A(i, :), '-', ...
        'Color', A2A_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_D1, 1)
    plot(x_points, large_choice_D1(i, :), '-', ...
        'Color', D1_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_Control, sem_Control, Control_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', Control_color, 'MarkerFaceColor', Control_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_A2A, sem_A2A, A2A_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', A2A_color, 'MarkerFaceColor', A2A_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_D1, sem_D1, D1_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', D1_color, 'MarkerFaceColor', D1_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 17]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:15);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;



%%
















































%%
%% for males vs females

large_choice_mCherry = [risk_table.block_1_large(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('Male', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_large(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_2_large(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_3_large(strcmp('Female', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% for BLA-NAcSh ChrimsonR vs mCherry

large_choice_mCherry = [risk_table.block_1_small(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('Male', risk_table.TreatmentCondition))]*100;
large_choice_hM4Di = [risk_table.block_1_small(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_2_small(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_3_small(strcmp('Female', risk_table.TreatmentCondition))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 100]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%% 

large_choice_mCherry = [risk_table.large_aborts_block_1(strcmp('Male', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('Male', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_aborts_block_1(strcmp('Female', risk_table.TreatmentCondition)), risk_table.large_aborts_block_2(strcmp('Female', risk_table.TreatmentCondition)), risk_table.large_aborts_block_3(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 125]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:125);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% 

large_choice_mCherry = [risk_table.small_aborts_block_1(strcmp('Male', risk_table.TreatmentCondition)), risk_table.small_aborts_block_2(strcmp('Male', risk_table.TreatmentCondition)), risk_table.small_aborts_block_3(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.small_aborts_block_1(strcmp('Female', risk_table.TreatmentCondition)), risk_table.small_aborts_block_2(strcmp('Female', risk_table.TreatmentCondition)), risk_table.small_aborts_block_3(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 125]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:125);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice_mCherry = [risk_table.block_1_large_choice_latency(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_choice_latency(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 20]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:20);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%%

large_choice_mCherry = [risk_table.block_1_small_choice_latency(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_2_small_choice_latency(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_3_small_choice_latency(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_small_choice_latency(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_2_small_choice_latency(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_3_small_choice_latency(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 20]); % Adjust ylim dynamically
set(gca, 'ytick', 0:5:20);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;
%%

large_choice_mCherry = [risk_table.block_1_omission_total(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_omission_total(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_2_omission_total(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_3_omission_total(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 30]); % Adjust ylim dynamically
set(gca, 'ytick', 0:10:30);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%%

large_choice_mCherry = [risk_table.block_1_large_collect_latency(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_collect_latency(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_2_large_collect_latency(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_3_large_collect_latency(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 8]); % Adjust ylim dynamically
set(gca, 'ytick', 0:2:8);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice_mCherry = [risk_table.block_1_small_collect_latency(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_2_small_collect_latency(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_3_small_collect_latency(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_small_collect_latency(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_2_small_collect_latency(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_3_small_collect_latency(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 8]); % Adjust ylim dynamically
set(gca, 'ytick', 0:2:8);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;


%%

large_choice_mCherry = [risk_table.large_consum_duration_block_1(strcmp('Male', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('Male', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.large_consum_duration_block_1(strcmp('Female', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_2(strcmp('Female', risk_table.TreatmentCondition)), risk_table.large_consum_duration_block_3(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 8]); % Adjust ylim dynamically
set(gca, 'ytick', 0:2:8);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice_mCherry = [risk_table.small_consum_duration_block_1(strcmp('Male', risk_table.TreatmentCondition)), risk_table.small_consum_duration_block_2(strcmp('Male', risk_table.TreatmentCondition)), risk_table.small_consum_duration_block_3(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.small_consum_duration_block_1(strcmp('Female', risk_table.TreatmentCondition)), risk_table.small_consum_duration_block_2(strcmp('Female', risk_table.TreatmentCondition)), risk_table.small_consum_duration_block_3(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 8]); % Adjust ylim dynamically
set(gca, 'ytick', 0:2:8);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice_mCherry = [risk_table.win_stay_ratio_block_1(strcmp('Male', risk_table.TreatmentCondition)), risk_table.win_stay_ratio_block_2(strcmp('Male', risk_table.TreatmentCondition)), risk_table.win_stay_ratio_block_3(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.win_stay_ratio_block_1(strcmp('Female', risk_table.TreatmentCondition)), risk_table.win_stay_ratio_block_2(strcmp('Female', risk_table.TreatmentCondition)), risk_table.win_stay_ratio_block_3(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 1]); % Adjust ylim dynamically
set(gca, 'ytick', 0:0.2:1);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice_mCherry = [risk_table.block_1_large_choice_latency_following_shk(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_choice_latency_following_shk(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 8]); % Adjust ylim dynamically
set(gca, 'ytick', 0:2:8);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%

large_choice_mCherry = [risk_table.block_1_large_choice_latency_following_shk(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('Male', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('Male', risk_table.TreatmentCondition))];
large_choice_hM4Di = [risk_table.block_1_large_choice_latency_following_shk(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_2_large_choice_latency_following_shk(strcmp('Female', risk_table.TreatmentCondition)), risk_table.block_3_large_choice_latency_following_shk(strcmp('Female', risk_table.TreatmentCondition))];

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small', 'MarkerEdgeColor', 'none'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 8]); % Adjust ylim dynamically
set(gca, 'ytick', 0:2:8);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;



%%

gray_color = [.8 .8 .8];

% x = risk_table.shock_sens;
x = risk_table.Mean_1_to_3;
% y = risk_table.session_weight_percent_free_feed;
% y = risk_table.freefeed_weight;
y = risk_table.shock_sens;
% treatment = risk_table.TreatmentCondition;

% Separate data by treatment group
female_idx = strcmp(treatment, 'Female');
male_idx = strcmp(treatment, 'Male');

x_female = x(female_idx);
y_female = y(female_idx);
x_male = x(male_idx);
y_male = y(male_idx);

% Calculate regression statistics for females
coefficients_female = polyfit(x_female, y_female, 1);
y_pred_female = polyval(coefficients_female, x_female);
ssr_female = sum((y_pred_female - mean(y_female)).^2);
sst_female = sum((y_female - mean(y_female)).^2);
r_squared_female = ssr_female / sst_female;
[r_female, pval_female] = corrcoef(x_female, y_female);
r_val_female = r_female(1, 2);
p_val_female = pval_female(1, 2);

% Calculate regression statistics for males
coefficients_male = polyfit(x_male, y_male, 1);
y_pred_male = polyval(coefficients_male, x_male);
ssr_male = sum((y_pred_male - mean(y_male)).^2);
sst_male = sum((y_male - mean(y_male)).^2);
r_squared_male = ssr_male / sst_male;
[r_male, pval_male] = corrcoef(x_male, y_male);
r_val_male = r_male(1, 2);
p_val_male = pval_male(1, 2);

% Calculate regression statistics for all
coefficients_all = polyfit(x, y, 1);
y_pred_all = polyval(coefficients_all, x);
ssr_all = sum((y_pred_all - mean(y)).^2);
sst_all = sum((y - mean(y)).^2);
r_squared_all = ssr_all / sst_all;
[r_all, pval_all] = corrcoef(x, y);
r_val_all = r_all(1, 2);
p_val_all = pval_all(1, 2);

% Create figure with specified dimensions
figure;
width = 350; % Width of the figure
height = 350; % Height of the figure
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;

% Plot scatter points for each group
scatter(x_female, y_female, 100, 'green', 'filled', "square");
scatter(x_male, y_male, 100, gray_color, 'filled');

% Set the axis labels to have 2 decimal places
xtickformat('%.2f');
ytickformat('%.2f');

% Create regression lines
x_fit_female = linspace(min(x_female), max(x_female), 100);
y_fit_female = polyval(coefficients_female, x_fit_female);
x_fit_male = linspace(min(x_male), max(x_male), 100);
y_fit_male = polyval(coefficients_male, x_fit_male);
x_fit_all = linspace(min(x), max(x), 100);
y_fit_all = polyval(coefficients_all, x_fit_all);

% Plot regression lines
plot(x_fit_female, y_fit_female, 'Color', 'green', 'LineWidth', 2); % Pink line
plot(x_fit_male, y_fit_male, 'Color', gray_color, 'LineWidth', 2); % Black line
plot(x_fit_all, y_fit_all, 'Color', 'black', 'LineWidth', 2); % Black line

% Get axes limits for text positioning
ax = gca;
xLimits = xlim(ax);
yLimits = ylim(ax);

% Position text for females (upper right area)
xPos_female = xLimits(2) - 0.05 * range(xLimits);
yPos_female = yLimits(2) - 0.05 * range(yLimits);

% Position text for males (lower right area)
xPos_male = xLimits(2) - 0.05 * range(xLimits);
yPos_male = yLimits(2) - 0.25 * range(yLimits);

% Add statistics text for females
text(xPos_female, yPos_female, ...
    {['Female:'], ...
     ['R^2 = ' num2str(r_squared_female, '%.2f')], ...
     ['R = ' num2str(r_val_female, '%.2f')], ...
     ['p = ' num2str(p_val_female, '%.2f')]}, ...
    'FontSize', 10, ...
    'Color', 'green', ... % Pink color
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold');

% Add statistics text for males
text(xPos_male, yPos_male, ...
    {['Male:'], ...
     ['R^2 = ' num2str(r_squared_male, '%.2f')], ...
     ['R = ' num2str(r_val_male, '%.2f')], ...
     ['p = ' num2str(p_val_male, '%.2f')]}, ...
    'FontSize', 10, ...
    'Color', gray_color, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold');

% Add statistics text for females
text(xPos_female, yPos_female, ...
    {['Female:'], ...
     ['R^2 = ' num2str(r_squared_all, '%.2f')], ...
     ['R = ' num2str(r_val_all, '%.2f')], ...
     ['p = ' num2str(p_val_all, '%.2f')]}, ...
    'FontSize', 10, ...
    'Color', 'black', ... % Pink color
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold');

% Add labels and legend
xlabel('Riskiness');
ylabel('Percent free-feeding weight');
legend({'Female', 'Male', 'Female fit', 'Male fit'}, 'Location', 'best');

hold off;

% Display summary statistics in command window
fprintf('Female group statistics:\n');
fprintf('  R-squared: %.3f\n', r_squared_female);
fprintf('  R: %.3f\n', r_val_female);
fprintf('  p-value: %.3f\n', p_val_female);
fprintf('  Regression equation: y = %.3fx + %.3f\n\n', coefficients_female(1), coefficients_female(2));

fprintf('Male group statistics:\n');
fprintf('  R-squared: %.3f\n', r_squared_male);
fprintf('  R: %.3f\n', r_val_male);
fprintf('  p-value: %.3f\n', p_val_male);
fprintf('  Regression equation: y = %.3fx + %.3f\n', coefficients_male(1), coefficients_male(2));

%%
% Extract data

% y = risk_table.freefeed_weight;
% y = risk_table.session_weight_raw_grams;
% y = risk_table.session_weight_percent_free_feed;
y = risk_table.shock_sens;
% y = risk_table.session_length;  
% y = risk_table.pre_training_sessions_sum;  
% y = risk_table.sessions_to_rm_criterion;  




treatment = risk_table.TreatmentCondition;

% Separate data by treatment group
female_idx = strcmp(treatment, 'Female');
male_idx = strcmp(treatment, 'Male');

y_female = y(female_idx);
y_male = y(male_idx);

% Calculate means for bar plot
mean_female = mean(y_female);
mean_male = mean(y_male);

% Calculate standard errors for error bars
se_female = std(y_female) / sqrt(length(y_female));
se_male = std(y_male) / sqrt(length(y_male));

% Create figure
figure;
width = 300; % Width of the figure
height = 400; % Height of the figure
set(gcf, 'Position', [100, 100, width, height]);

% Create bar plot
bar_data = [mean_female, mean_male];
bar_colors = [[0 1 0]; 0 0 0]; % Pink for female, black for male
b = bar(bar_data, 'FaceColor', 'flat');
b.CData = bar_colors;

hold on;

% Add error bars
errorbar([1, 2], [mean_female, mean_male], [se_female, se_male], ...
    'k', 'LineWidth', 1.5, 'CapSize', 10);

% Add individual data points as scatter overlay
% Female data points
x_female_scatter = ones(length(y_female), 1) + 0.1 * randn(length(y_female), 1); % Add jitter
scatter(x_female_scatter, y_female, 70, 'green', 'filled', "square");

% Male data points
x_male_scatter = 2 * ones(length(y_male), 1) + 0.1 * randn(length(y_male), 1); % Add jitter
scatter(x_male_scatter, y_male, 50, 'black', 'filled', "o");

% Customize the plot
set(gca, 'XTick', [1, 2]);
set(gca, 'XTickLabel', {'Female', 'Male'});
ylim([0 0.7]); % Adjust ylim dynamically
set(gca, 'ytick', 0:.2:.7);
xlabel('Treatment Condition');
ylabel('Weight');
% title('Footshock Sensitivity by Treatment Group');

% Add sample sizes to the plot
text(1, min(ylim) + 0.1 * range(ylim), ['n = ' num2str(length(y_female))], ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'green', 'FontWeight', 'bold');
text(2, min(ylim) + 0.1 * range(ylim), ['n = ' num2str(length(y_male))], ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'black', 'FontWeight', 'bold');

% Add legend
legend({'Mean ± SEM', '', 'Female', 'Male'}, 'Location', 'best');

hold off;

% Display summary statistics in command window
fprintf('Summary Statistics:\n');
fprintf('Female group:\n');
fprintf('  Mean: %.3f\n', mean_female);
fprintf('  Standard Error: %.3f\n', se_female);
fprintf('  Sample size: %d\n\n', length(y_female));

fprintf('Male group:\n');
fprintf('  Mean: %.3f\n', mean_male);
fprintf('  Standard Error: %.3f\n', se_male);
fprintf('  Sample size: %d\n', length(y_male));

% Perform t-test to compare groups
[h, p_ttest] = ttest2(y_female, y_male);
fprintf('\nStatistical comparison (two-sample t-test):\n');
fprintf('  p-value: %.4f\n', p_ttest);
if h == 1
    fprintf('  Result: Significantly different (p < 0.05)\n');
else
    fprintf('  Result: Not significantly different (p >= 0.05)\n');
end

%%

treatment = risk_table.TreatmentCondition;

% Separate data by treatment group
female_idx = strcmp(treatment, 'Female');
male_idx = strcmp(treatment, 'Male');

did_not_complete_sessions = risk_table.session_length > (91*60);


females_did_not_complete_percent = (sum(did_not_complete_sessions(female_idx))/size(did_not_complete_sessions(female_idx), 1))*100
males_did_not_complete_percent = (sum(did_not_complete_sessions(male_idx))/size(did_not_complete_sessions(male_idx), 1))*100

females_completed_percent = 100-females_did_not_complete_percent;
males_completed_percent = 100-males_did_not_complete_percent;

figure; donutchart([females_did_not_complete_percent females_completed_percent], 'InnerRadius', 0.7, 'StartAngle', 270)
figure; donutchart([males_did_not_complete_percent males_completed_percent], 'InnerRadius', 0.7, 'StartAngle', 270)
