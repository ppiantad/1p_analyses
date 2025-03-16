% iter = iter+1;
% neuron_num = 0;
% these are mice that did not complete the entire session - kinda have to
% toss them to do some comparisons during RDT

% final_behavior = final_SLEAP; % for hM4Di data;

session_to_analyze = 'RDT_D1'

if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
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
bar_handle = bar(group_means, 'FaceAlpha', 0.7, 'BarWidth', 0.5); % Create bar plot with some transparency
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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, cross_sess_small_means, cross_sess_small_sems, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 18, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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

large_choice = [risk_table.block_1_large_collect_latency_no_shk_trials, risk_table.block_2_large_collect_latency_no_shk_trials, risk_table.block_3_large_collect_latency_no_shk_trials];
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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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

%% need to run get_path_lengths

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
% small_choice = [risk_table.block_1_small_choice_latency, risk_table.block_2_small_choice_latency, risk_table.block_3_small_choice_latency];

% turn NaNs (trials where the denominator, number of losses, is 0) into 0s
large_choice(isnan(large_choice))=0;


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

% Plot individual lines for "Small" data
% for i = 1:size(small_choice, 1)
%     plot(x_points, small_choice(i, :), '-', ...
%         'Color', [1 0 0 0.6], ... % Red with 60% opacity
%         'LineWidth', 1.2);
% end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

% errorbar(x_points, mean_small, sem_small, '^-', ...
%     'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
%     'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'
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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 12, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, 'MarkerEdgeColor', 'none',...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_large_not_risky, sem_large_not_risky, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_large_not_risky, sem_large_not_risky, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

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

