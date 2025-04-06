% load dataset that contains data where events have been identified

clear similarityMatrix

data_to_load = {'BLA_C_raw_no_additional_filtering_RDT_D1_only_completed_sessions_zall_window_base_workspace_10_categories.mat', 'BLA_C_raw_no_additional_filtering_Pre_RDT_RM_only_completed_sessions_zall_window_base_workspace_10_categories.mat', 'BLA_RM_D1_6_categories_03212025_including_spike_prob_array_counts.mat'};

load(data_to_load{1})


if size(respClass_all_array, 2) == 10
    % comparison_arrays_full = [1 2 3; 8 9 10]
    comparison_arrays_full = [1 2 3; 5 6 7]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays_full = [1 2 3; 4 5 6]
end

% depending on analysis desired, may need to run upper part, then run
% block_wise_changes_v1.m, then run below

% some data need the following uncommented:
% ca_data_type = uv.ca_data_type;

%%
for aa = 1:size(comparison_arrays_full, 1)
    
    comparison_arrays = comparison_arrays_full(aa, :)
    first_session = session_to_analyze;
    

    for gg = 1:size(animalIDs, 1)
        select_mouse = animalIDs{gg};

        select_mouse_index = find(strcmp(animalIDs, select_mouse));

        BehavData = final.(select_mouse).(first_session).uv.BehavData;


        %%
        PV_prechoice_all_mouse = [];
        for ff = 1:size(caTraceTrials_mouse{select_mouse_index, comparison_arrays(1, 1)}, 2)
            for cc = 1:size(caTraceTrials_mouse{select_mouse_index,comparison_arrays(1, 1)}{1, ff}, 1)
                PV_prechoice_all_mouse(cc, ff) = mean(caTraceTrials_mouse{select_mouse_index, comparison_arrays(1, 1)}{1, ff}(cc, ts1 >= -4 & ts1 <= 0));

            end

        end

        % reorganize data for plotting the heatmap, with ensemble-specific
        % neurons on top, the rest below
        PV_prechoice_all_mouse_just_prechoice = PV_prechoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} == 1);
        PV_prechoice_all_mouse_just_prechoice = [PV_prechoice_all_mouse_just_prechoice PV_prechoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} ~= 1)];
        PV_prechoice_all_mouse_just_prechoice_array{aa} = PV_prechoice_all_mouse_just_prechoice;
        % data to use below, without reorg
        mean_PV_prechoice_all_mouse = mean(PV_prechoice_all_mouse)';
        
        % use the code below if you want to subsample the data by some
        % factor (change the demominator after
        % round(size(mean_PV_prechoice_all_mouse, 1)/4);)
        % set the desired number of rows
        % num_rows = round(size(mean_PV_prechoice_all_mouse, 1)/10);
        % selected_rows = datasample(1:size(mean_PV_prechoice_all_mouse, 1), num_rows, 'Replace', false);
        % mean_PV_prechoice_all_mouse  = mean_PV_prechoice_all_mouse(selected_rows, :);
        
        % mean_PV_prechoice_all_mouse = mean_PV_prechoice_all_mouse(prechoice_lost_mouse{gg, 1} == 1, :);
        % mean_PV_prechoice_all_mouse = mean_PV_prechoice_all_mouse(respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} == 1, :);
        % mean_PV_prechoice_all_mouse = mean_PV_prechoice_all_mouse(true_neutral_block_1_mouse{gg, 1} ~= 1, :);
        % mean_PV_prechoice_all_mouse = mean_PV_prechoice_all_mouse(prechoice_conserved_mouse{gg, 1} ~= 1 & prechoice_lost_mouse{gg, 1} ~= 1 & prechoice_remapped_mouse{gg, 1} ~= 1);
        
        % hold on;
        % figure; imagesc(1:size(PV_prechoice_all_mouse_just_prechoice, 1), [], PV_prechoice_all_mouse_just_prechoice')
        % colormap gray;
        % colorbar;
        % clim([0 1]);
        % hold off;
        % 
        % % Define the position for the first green line
        % y1_start = 1;  % Start at the first neuron
        % y1_end = sum(respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} == 1);  % End at the number of neurons in the first group
        % 
        % % Define the position for the second line
        % y2_start = y1_end + 1;  % Start right after the first group
        % y2_end = size(PV_prechoice_all_mouse_just_prechoice, 2);  % End at the last neuron
        % 
        % % Draw the first vertical green line
        % hold on;
        % line([1 1], [y1_start y1_end], 'Color', 'g', 'LineWidth', 2);
        % 
        % % Draw the second vertical line (can be any color, e.g., black here)
        % line([1 1], [y2_start y2_end], 'Color', 'r', 'LineWidth', 2);
        % title(['From animal ',select_mouse], 'Interpreter', 'none')
        % 
        % 
        % hold off;


        ca = final.(select_mouse).(first_session).CNMFe_data.(ca_data_type);
        
        % filter based on variable above
        % ca = ca(selected_rows, :);
        % ca = ca(prechoice_conserved_mouse{gg, 1} ~= 1 & prechoice_lost_mouse{gg, 1} ~= 1 & prechoice_remapped_mouse{gg, 1} ~= 1, :);
        % ca_zscored = zscore(ca, [], 2);
        ca_zscored = ca;

        % ca_zscored = ca(prechoice_lost_mouse{gg, 1} == 1, :);
        time_array = final.(select_mouse).(first_session).time;
        prechoice_similarityOverTime = [];
        for t = 1:size(ca_zscored, 2)
            % for i = 1:length(uniqueTypes)
            % typeIdx = (neuronTypes == uniqueTypes(i));
            % use specific subset of neurons
            activitySubset = ca_zscored(:, t);
            % uncomment if you want to use all neurons
            %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
            similarityMatrix = corrcoef(activitySubset, mean_PV_prechoice_all_mouse);
            prechoice_similarityOverTime(t) = similarityMatrix(2);
            % end
        end
        % 
        % figure; plot(time_array, prechoice_similarityOverTime(1, 1:size(time_array, 1)))
        % xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
        % xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
        % xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
        % xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
        % if strcmp('shock',BehavData.Properties.VariableNames)
        %     xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
        % end

        %%

        PV_postchoice_all_mouse = [];
        for ff = 1:size(caTraceTrials_mouse{select_mouse_index, comparison_arrays(1, 2)}, 2)
            for cc = 1:size(caTraceTrials_mouse{select_mouse_index, comparison_arrays(1, 2)}{1, ff}, 1)
                PV_postchoice_all_mouse(cc, ff) = mean(caTraceTrials_mouse{select_mouse_index, comparison_arrays(1, 2)}{1, ff}(cc, ts1 >= 0 & ts1 <= 2));

            end

        end


        PV_postchoice_all_mouse_just_postchoice = PV_postchoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 2)} == 1);
        PV_postchoice_all_mouse_just_postchoice = [PV_postchoice_all_mouse_just_postchoice PV_postchoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 2)} ~= 1)];
        mean_PV_postchoice_all_mouse = mean(PV_postchoice_all_mouse)';
        % hold on;
        % figure; imagesc(1:size(PV_postchoice_all_mouse_just_postchoice, 1), [], PV_postchoice_all_mouse_just_postchoice')
        % colormap gray;
        % colorbar;
        % clim([0 2]);
        % hold off;
        % 
        % hold off;
        % 
        % % Define the position for the first green line
        % y1_start = 1;  % Start at the first neuron
        % y1_end = sum(respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 2)} == 1);  % End at the number of neurons in the first group
        % 
        % % Define the position for the second line
        % y2_start = y1_end + 1;  % Start right after the first group
        % y2_end = size(PV_postchoice_all_mouse_just_postchoice, 2);  % End at the last neuron
        % 
        % % Draw the first vertical green line
        % hold on;
        % line([1 1], [y1_start y1_end], 'Color', 'g', 'LineWidth', 2);
        % 
        % % Draw the second vertical line (can be any color, e.g., black here)
        % line([1 1], [y2_start y2_end], 'Color', 'r', 'LineWidth', 2);
        % 
        % hold off;


        ca = final.(select_mouse).(session_to_analyze).CNMFe_data.(ca_data_type);
        % ca_zscored = zscore(ca, [], 2);
        ca_zscored = ca;
        time_array = final.(select_mouse).(first_session).time;
        postchoice_similarityOverTime = [];
        for t = 1:size(ca_zscored, 2)
            % for i = 1:length(uniqueTypes)
            % typeIdx = (neuronTypes == uniqueTypes(i));
            % use specific subset of neurons
            activitySubset = ca_zscored(:, t);
            % uncomment if you want to use all neurons
            %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
            similarityMatrix = corrcoef(activitySubset, mean_PV_postchoice_all_mouse);
            postchoice_similarityOverTime(t) = similarityMatrix(2);
            % end
        end

        % figure; plot(time_array, postchoice_similarityOverTime(1, 1:size(time_array, 1)))
        % xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
        % xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
        % xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
        % xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
        % if strcmp('shock',BehavData.Properties.VariableNames)
        %     xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
        % end


        %%

        PV_consumption_all_mouse = [];
        for ff = 1:size(caTraceTrials_mouse{select_mouse_index, comparison_arrays(1, 3)}, 2)
            for cc = 1:size(caTraceTrials_mouse{select_mouse_index, comparison_arrays(1, 3)}{1, ff}, 1)
                PV_consumption_all_mouse(cc, ff) = mean(caTraceTrials_mouse{select_mouse_index, comparison_arrays(1, 3)}{1, ff}(cc, ts1 >= 1 & ts1 <= 3));

            end

        end


        PV_consumption_all_mouse_just_consumption = PV_consumption_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 3)} == 1);
        PV_consumption_all_mouse_just_consumption = [PV_consumption_all_mouse_just_consumption PV_consumption_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 3)} ~= 1)];
        mean_PV_consumption_all_mouse = mean(PV_consumption_all_mouse)';
        % hold on;
        % figure; imagesc(1:size(PV_consumption_all_mouse_just_consumption, 1), [], PV_consumption_all_mouse_just_consumption')
        % colormap gray;
        % colorbar;
        % clim([0 2]);
        % hold off;
        % 
        % hold off;
        % 
        % % Define the position for the first green line
        % y1_start = 1;  % Start at the first neuron
        % y1_end = sum(respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 3)} == 1);  % End at the number of neurons in the first group
        % 
        % % Define the position for the second line
        % y2_start = y1_end + 1;  % Start right after the first group
        % y2_end = size(PV_consumption_all_mouse_just_consumption, 2);  % End at the last neuron
        % 
        % % Draw the first vertical green line
        % hold on;
        % line([1 1], [y1_start y1_end], 'Color', 'g', 'LineWidth', 2);
        % 
        % % Draw the second vertical line (can be any color, e.g., black here)
        % line([1 1], [y2_start y2_end], 'Color', 'r', 'LineWidth', 2);
        % 
        % hold off;


        ca = final.(select_mouse).(session_to_analyze).CNMFe_data.(ca_data_type);
        % ca_zscored = zscore(ca, [], 2);
        ca_zscored = ca;
        time_array = final.(select_mouse).(first_session).time;
        consumption_similarityOverTime = [];
        for t = 1:size(ca_zscored, 2)
            % for i = 1:length(uniqueTypes)
            % typeIdx = (neuronTypes == uniqueTypes(i));
            % use specific subset of neurons
            activitySubset = ca_zscored(:, t);
            % uncomment if you want to use all neurons
            %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
            similarityMatrix = corrcoef(activitySubset, mean_PV_consumption_all_mouse);
            consumption_similarityOverTime(t) = similarityMatrix(2);
            % end
        end

        % figure; plot(time_array, consumption_similarityOverTime(1, 1:size(time_array, 1)))
        % xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
        % xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
        % xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
        % xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
        % if strcmp('shock',BehavData.Properties.VariableNames)
        %     xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
        % end


        %%
        % figure; plot(time_array, prechoice_similarityOverTime(1, 1:size(time_array, 1)))
        % xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
        % xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
        % xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
        % xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
        % if strcmp('shock',BehavData.Properties.VariableNames)
        %     xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
        % end
        % hold on; plot(time_array, postchoice_similarityOverTime(1, 1:size(time_array, 1)))
        % hold on; plot(time_array, consumption_similarityOverTime(1, 1:size(time_array, 1)))

        %% attempting to get something similar to Courtin Fig 5e
        prechoice_similarityOverTimeTrials = [];
        postchoice_similarityOverTimeTrials = [];
        consumption_similarityOverTimeTrials = [];

        eTS = [BehavData.choiceTime BehavData.collectionTime];

        % Define time windows
        pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
        post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
        consumption_window = [1 3];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

        % Find indices corresponding to each time window
        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        for t = 1:size(eTS,1)

            timeWin_prechoice_post_choice = [eTS(t, 1)+uv.evtWin(1,1):uv.dt:eTS(t, 1)+uv.evtWin(1,2)];  %calculate time window around each event
            timeWin_consumption = [eTS(t, 2)+uv.evtWin(1,1):uv.dt:eTS(t, 2)+uv.evtWin(1,2)];  %calculate time window around each event
            BL_win = [eTS(t, 1)+uv.BLper(1,1):uv.dt:eTS(t, 1)+uv.BLper(1,2)];



            if min(timeWin_prechoice_post_choice) > min(time_array) && max(timeWin_prechoice_post_choice) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                % get unit event counts in trials
                % get unit ca traces in trials

                idx_prechoice_post_choice = time_array > min(timeWin_prechoice_post_choice) & time_array < max(timeWin_prechoice_post_choice);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                idx_consumption = time_array > min(timeWin_consumption) & time_array < max(timeWin_consumption);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                bl_idx = time_array > min(BL_win) & time_array < max(BL_win);
                %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                prechoice_similarityOverTimeTrials(t,1:sum(idx_prechoice_post_choice)) = prechoice_similarityOverTime(idx_prechoice_post_choice);
                postchoice_similarityOverTimeTrials(t, 1:sum(idx_prechoice_post_choice)) = postchoice_similarityOverTime(idx_prechoice_post_choice);
                consumption_similarityOverTimeTrials(t, 1:sum(idx_consumption)) = consumption_similarityOverTime(idx_consumption);
                % prechoice_similarityOverTimeTrials_zb(t,:) = nanmean(unitTrace(bl_idx)); %baseline mean
                % prechoice_similarityOverTimeTrials_zb_window(t,:) = nanmean(prechoice_similarityOverTimeTrials(t,:));
                % prechoice_similarityOverTimeTrials_zsd(t,:) = nanstd(unitTrace(bl_idx)); %baseline std
                % prechoice_similarityOverTimeTrials_zsd_window(t,:) = nanstd(prechoice_similarityOverTimeTrials(t,:));
                % tmp = 0;
                % for j = 1:size(prechoice_similarityOverTimeTrials,2)
                %     tmp = tmp+1;
                %     % prechoice_similarityOverTimeTrials_zall_baselined(t,tmp) = (prechoice_similarityOverTimeTrials(t,j) - prechoice_similarityOverTimeTrials_zb(t))/prechoice_similarityOverTimeTrials_zsd(t);
                %     prechoice_similarityOverTimeTrials_raw(t, tmp) = prechoice_similarityOverTimeTrials(t,j);
                %     % prechoice_similarityOverTimeTrials_zall_window(t,tmp) = (prechoice_similarityOverTimeTrials(t,j) - prechoice_similarityOverTimeTrials_zb_window(t))/prechoice_similarityOverTimeTrials_zsd_window(t);
                %     % prechoice_similarityOverTimeTrials_zall_session(t,tmp) = (prechoice_similarityOverTimeTrials(t,j) - zb_session(u))/zsd_session(u);
                % end
                % clear j;

            end


        end

        block_1_inds = BehavData.Block == 1 & BehavData.Blank_Touch == 0 & BehavData.omissionALL == 0;
        block_2_inds = BehavData.Block == 2 & BehavData.Blank_Touch == 0 & BehavData.omissionALL == 0;
        block_3_inds = BehavData.Block == 3 & BehavData.Blank_Touch == 0 & BehavData.omissionALL == 0;

        prechoice_period_b1_mean = mean(prechoice_similarityOverTimeTrials(block_1_inds, pre_choice_indices), 2);
        prechoice_period_b2_mean = mean(prechoice_similarityOverTimeTrials(block_2_inds, pre_choice_indices), 2);
        prechoice_period_b3_mean = mean(prechoice_similarityOverTimeTrials(block_3_inds, pre_choice_indices), 2);

        postchoice_period_b1_mean = mean(postchoice_similarityOverTimeTrials(block_1_inds, post_choice_indices), 2);
        postchoice_period_b2_mean = mean(postchoice_similarityOverTimeTrials(block_2_inds, post_choice_indices), 2);
        postchoice_period_b3_mean = mean(postchoice_similarityOverTimeTrials(block_3_inds, post_choice_indices), 2);

        consumption_period_b1_mean = mean(consumption_similarityOverTimeTrials(block_1_inds, consumption_indices), 2);
        consumption_period_b2_mean = mean(consumption_similarityOverTimeTrials(block_2_inds, consumption_indices), 2);
        consumption_period_b3_mean = mean(consumption_similarityOverTimeTrials(block_3_inds, consumption_indices), 2);


        prechoice_over_time_mean(:, gg) = [prechoice_period_b1_mean; prechoice_period_b2_mean; prechoice_period_b3_mean];
        postchoice_over_time_mean(:, gg) = [postchoice_period_b1_mean; postchoice_period_b2_mean; postchoice_period_b3_mean];
        consumption_over_time_mean(:, gg) = [consumption_period_b1_mean; consumption_period_b2_mean; consumption_period_b3_mean];


        prechoice_period_b1_sem = std(prechoice_similarityOverTimeTrials(block_1_inds, pre_choice_indices), 0, 2) / sqrt(length(find(block_1_inds == 1)));
        prechoice_period_b2_sem = std(prechoice_similarityOverTimeTrials(block_2_inds, pre_choice_indices), 0, 2) / sqrt(length(find(block_2_inds == 1)));
        prechoice_period_b3_sem = std(prechoice_similarityOverTimeTrials(block_3_inds, pre_choice_indices), 0, 2) / sqrt(length(find(block_3_inds == 1)));


        postchoice_period_b1_sem = std(postchoice_similarityOverTimeTrials(block_1_inds, post_choice_indices), 0, 2) / sqrt(length(find(block_1_inds == 1)));
        postchoice_period_b2_sem = std(postchoice_similarityOverTimeTrials(block_2_inds, post_choice_indices), 0, 2) / sqrt(length(find(block_2_inds == 1)));
        postchoice_period_b3_sem = std(postchoice_similarityOverTimeTrials(block_3_inds, post_choice_indices), 0, 2) / sqrt(length(find(block_3_inds == 1)));


        consumption_period_b1_sem = std(consumption_similarityOverTimeTrials(block_1_inds, consumption_indices), 0, 2) / sqrt(length(find(block_1_inds == 1)));
        consumption_period_b2_sem = std(consumption_similarityOverTimeTrials(block_2_inds, consumption_indices), 0, 2) / sqrt(length(find(block_2_inds == 1)));
        consumption_period_b3_sem = std(consumption_similarityOverTimeTrials(block_3_inds, consumption_indices), 0, 2) / sqrt(length(find(block_3_inds == 1)));



        prechoice_over_time_sem(:, gg) = [prechoice_period_b1_sem; prechoice_period_b2_sem; prechoice_period_b3_sem];
        postchoice_over_time_sem(:, gg) = [postchoice_period_b1_sem; postchoice_period_b2_sem; postchoice_period_b3_sem];
        consumption_over_time_sem(:, gg) = [consumption_period_b1_sem; consumption_period_b2_sem; consumption_period_b3_sem];

        prechoice_similarityOverTime_array{aa, gg} = prechoice_similarityOverTime;
        postchoice_similarityOverTime_array{aa, gg} = postchoice_similarityOverTime;
        consumption_similarityOverTime_array{aa, gg} = consumption_similarityOverTime;


        prechoice_similarityOverTimeTrials_array{aa, gg} = prechoice_similarityOverTimeTrials;
        postchoice_similarityOverTimeTrials_array{aa, gg} = postchoice_similarityOverTimeTrials;
        consumption_similarityOverTimeTrials_array{aa, gg} = consumption_similarityOverTimeTrials;
    end
    %%
    mean_prechoice_over_time(:, aa) = mean(prechoice_over_time_mean, 2);
    mean_postchoice_over_time(:, aa) = mean(postchoice_over_time_mean, 2);
    mean_consumption_over_time(:, aa) = mean(consumption_over_time_mean, 2);

    
    sem_prechoice_over_time(:, aa) = mean(prechoice_over_time_sem, 2); 
    sem_postchoice_over_time(:, aa) = mean(postchoice_over_time_sem, 2);
    sem_consumption_over_time(:, aa) = mean(consumption_over_time_sem, 2);

    prechoice_over_time_mean_iter{aa} = prechoice_over_time_mean;
    postchoice_over_time_mean_iter{aa} = postchoice_over_time_mean;
    consumption_over_time_mean_iter{aa} = consumption_over_time_mean;

    prechoice_over_time_sem_iter{aa} = prechoice_over_time_sem;
    postchoice_over_time_sem_iter{aa} = postchoice_over_time_sem;
    consumption_over_time_sem_iter{aa} = consumption_over_time_sem; 



    % % Calculate the standard deviation for each column
    % std_dev = std(prechoice_over_time_mean, 0, 2); % Use 2 for row-wise standard deviation
    % 
    % % Calculate the number of columns
    % n = size(prechoice_over_time_mean, 2);
    % 
    % % Calculate the SEM for each row
    % SEM = std_dev / sqrt(n);

    %
    % sem_prechoice_over_time = std(prechoice_over_time_mean, 0, 2)/sqrt(size(prechoice_over_time_mean, 2));
    % sem_postchoice_over_time = std(postchoice_over_time_mean, 0, 2)/sqrt(size(postchoice_over_time_mean, 2));
    % sem_consumption_over_time = std(consumption_over_time_mean, 0, 2)/sqrt(size(consumption_over_time_mean, 2));


    % figure; plot(mean(prechoice_over_time_mean, 2))
    % figure; plot(mean(postchoice_over_time_mean, 2))
    % figure; plot(mean(consumption_over_time_mean, 2))



end

SD_for_prechoice_over_time = prechoice_over_time_sem * sqrt(size(prechoice_over_time_sem, 2));
SD_for_postchoice_over_time = postchoice_over_time_sem * sqrt(size(postchoice_over_time_sem, 2));
SD_for_consumption_over_time = consumption_over_time_sem * sqrt(size(consumption_over_time_sem, 2));
SEM_for_prechoice_over_time = std(SD_for_prechoice_over_time, 0, 2) / sqrt(size(SD_for_prechoice_over_time, 2));

%%

figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 100; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
% xlim([-8 8]);
% % Set X-axis ticks
% set(gca, 'XTick', [-8, 0, 8]);

shadedErrorBar(1:90, mean_prechoice_over_time(:, 1), sem_prechoice_over_time(:, 1), 'lineProps', {'color', 'r'});
hold on; shadedErrorBar(1:90, mean_prechoice_over_time(:, 2), sem_prechoice_over_time(:, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});
xlim([1 90]);
% Set X-axis ticks
set(gca, 'XTick', [1, 30, 60, 90]);
xline([30 60])
ytickformat('%.2f');
ylim([0.1 0.7])
hold off

%%

figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 100; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
% xlim([-8 8]);
% % Set X-axis ticks
% set(gca, 'XTick', [-8, 0, 8]);

shadedErrorBar(1:90, mean_postchoice_over_time(:, 1), sem_postchoice_over_time(:, 1), 'lineProps', {'color', 'r'});
hold on; shadedErrorBar(1:90, mean_postchoice_over_time(:, 2), sem_postchoice_over_time(:, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});
xlim([1 90]);
% Set X-axis ticks
set(gca, 'XTick', [1, 30, 60, 90]);
xline([30 60])
ytickformat('%.2f');
ylim([0.1 0.7])
hold off
%%

figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 100; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
% xlim([-8 8]);
% % Set X-axis ticks
% set(gca, 'XTick', [-8, 0, 8]);

shadedErrorBar(1:90, mean_consumption_over_time(:, 1), sem_consumption_over_time(:, 1), 'lineProps', {'color', 'r'});
hold on; shadedErrorBar(1:90, mean_consumption_over_time(:, 2), sem_consumption_over_time(:, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});
xlim([1 90]);
% Set X-axis ticks
set(gca, 'XTick', [1, 30, 60, 90]);
xline([30 60])
ytickformat('%.2f');
ylim([0.1 0.7])
hold off

%%
% Assuming prechoice_over_time_mean_iter is a 1x2 cell array, each containing a 90x10 double
data1 = prechoice_over_time_mean_iter{1}; % First cell (90x10)
data2 = prechoice_over_time_mean_iter{2}; % Second cell (90x10)

% Define indices for splitting into 3 segments of 30 rows each
segments = [1, 30; 31, 60; 61, 90];

% Initialize mean matrices
mean_values = zeros(3,2); % 3 segments x 2 cell arrays
individual_means = cell(3,2); % To store individual column means

% Loop through each segment
for i = 1:3
    % Extract rows for this segment
    rows = segments(i,1):segments(i,2);
    
    % Compute means for each column and store
    mean_seg1 = mean(data1(rows, :), 1); % Mean over rows for first array
    mean_seg2 = mean(data2(rows, :), 1); % Mean over rows for second array
    
    % Store overall mean (for bar heights)
    mean_values(i,1) = mean(mean_seg1); 
    mean_values(i,2) = mean(mean_seg2);
    
    % Store individual column means (for scatter overlay)
    individual_means{i,1} = mean_seg1;
    individual_means{i,2} = mean_seg2;
end

% Plot bar graph
figure; hold on;

width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
bar_handle = bar(mean_values, 'grouped');

% Set bar colors
bar_handle(1).FaceColor = 'b'; % First cell (blue)
bar_handle(2).FaceColor = 'r'; % Second cell (red)

% Get X positions for scatter overlay
x_positions = zeros(3,2); % Store bar center positions
for i = 1:2
    x_positions(:,i) = bar_handle(i).XEndPoints; % Correct x positions
end

% Overlay individual data points correctly aligned
for i = 1:3
    scatter(x_positions(i,1) * ones(1,10), individual_means{i,1}, 'k', 'filled'); % First dataset
    scatter(x_positions(i,2) * ones(1,10), individual_means{i,2}, 'k', 'filled'); % Second dataset
end

% Labels and legend
xticks(mean(x_positions,2)); % Align x-ticks with bar groups
xticklabels({'0%', '50%', '75%'});
ylabel('Mean PV correlation');
legend({'Safe-identified', 'Risky-identified'}, 'Location', 'SouthEast');
ytickformat('%.2f');

hold off;

% can use ANOVA_test.m to analyze the results
large_choice = cell2mat(individual_means(1:3, 1))';
small_choice = cell2mat(individual_means(1:3, 2))';

%%
% Assuming prechoice_over_time_mean_iter is a 1x2 cell array, each containing a 90x10 double
data1 = postchoice_over_time_mean_iter{1}; % First cell (90x10)
data2 = postchoice_over_time_mean_iter{2}; % Second cell (90x10)

% Define indices for splitting into 3 segments of 30 rows each
segments = [1, 30; 31, 60; 61, 90];

% Initialize mean matrices
mean_values = zeros(3,2); % 3 segments x 2 cell arrays
individual_means = cell(3,2); % To store individual column means

% Loop through each segment
for i = 1:3
    % Extract rows for this segment
    rows = segments(i,1):segments(i,2);
    
    % Compute means for each column and store
    mean_seg1 = mean(data1(rows, :), 1); % Mean over rows for first array
    mean_seg2 = mean(data2(rows, :), 1); % Mean over rows for second array
    
    % Store overall mean (for bar heights)
    mean_values(i,1) = mean(mean_seg1); 
    mean_values(i,2) = mean(mean_seg2);
    
    % Store individual column means (for scatter overlay)
    individual_means{i,1} = mean_seg1;
    individual_means{i,2} = mean_seg2;
end

% Plot bar graph
figure; hold on;
bar_handle = bar(mean_values, 'grouped');

% Set bar colors
bar_handle(1).FaceColor = 'b'; % First cell (blue)
bar_handle(2).FaceColor = 'r'; % Second cell (red)

% Get X positions for scatter overlay
x_positions = zeros(3,2); % Store bar center positions
for i = 1:2
    x_positions(:,i) = bar_handle(i).XEndPoints; % Correct x positions
end

% Overlay individual data points correctly aligned
for i = 1:3
    scatter(x_positions(i,1) * ones(1,10), individual_means{i,1}, 'k', 'filled'); % First dataset
    scatter(x_positions(i,2) * ones(1,10), individual_means{i,2}, 'k', 'filled'); % Second dataset
end

% Labels and legend
xticks(mean(x_positions,2)); % Align x-ticks with bar groups
xticklabels({'0%', '50%', '75%'});
ylabel('Mean PV correlation');
legend({'Safe-identified', 'Risky-identified'}, 'Location', 'SouthEast');
ytickformat('%.2f');

hold off;

%%
% Assuming postchoice_over_time_mean_iter is a 1x2 cell array, each containing a 90x10 double
data1 = postchoice_over_time_mean_iter{1}; % First cell (90x10)
data2 = postchoice_over_time_mean_iter{2}; % Second cell (90x10)

% Define new grouping: first 30 rows, and remaining 60 rows
segments = {1:30, 31:90};

% Initialize mean matrices
mean_values = zeros(2,2); % 2 sets x 2 datasets
individual_means = cell(2,2); % To store individual column means

% Loop through each segment group
for i = 1:2
    rows = segments{i};
    
    % Compute means for each column and store
    mean_seg1 = mean(data1(rows, :), 1);
    mean_seg2 = mean(data2(rows, :), 1);
    
    % Store overall mean (for bar heights)
    mean_values(i,1) = mean(mean_seg1);
    mean_values(i,2) = mean(mean_seg2);
    
    % Store individual column means (for scatter overlay)
    individual_means{i,1} = mean_seg1;
    individual_means{i,2} = mean_seg2;
end

% Plot bar graph
figure; hold on;
bar_handle = bar(mean_values, 'grouped');

% Set bar colors
bar_handle(1).FaceColor = 'b'; % First cell (blue)
bar_handle(2).FaceColor = 'r'; % Second cell (red)

% Get X positions for scatter overlay
x_positions = zeros(2,2); % Store bar center positions
for i = 1:2
    x_positions(:,i) = bar_handle(i).XEndPoints;
end

% Overlay individual data points correctly aligned
for i = 1:2
    scatter(x_positions(i,1) * ones(1,10), individual_means{i,1}, 'k', 'filled'); % First dataset
    scatter(x_positions(i,2) * ones(1,10), individual_means{i,2}, 'k', 'filled'); % Second dataset
end

% Labels and legend
xticks(mean(x_positions,2));
xticklabels({'0%', '50–75%'}); % Adjusted labels
ylabel('Mean PV correlation');
legend({'Safe-identified', 'Risky-identified'}, 'Location', 'SouthEast');
ytickformat('%.2f');

hold off;

%%
% Assuming prechoice_over_time_mean_iter is a 1x2 cell array, each containing a 90x10 double
data1 = consumption_over_time_mean_iter{1}; % First cell (90x10)
data2 = consumption_over_time_mean_iter{2}; % Second cell (90x10)

% Define indices for splitting into 3 segments of 30 rows each
segments = [1, 30; 31, 60; 61, 90];

% Initialize mean matrices
mean_values = zeros(3,2); % 3 segments x 2 cell arrays
individual_means = cell(3,2); % To store individual column means

% Loop through each segment
for i = 1:3
    % Extract rows for this segment
    rows = segments(i,1):segments(i,2);
    
    % Compute means for each column and store
    mean_seg1 = mean(data1(rows, :), 1); % Mean over rows for first array
    mean_seg2 = mean(data2(rows, :), 1); % Mean over rows for second array
    
    % Store overall mean (for bar heights)
    mean_values(i,1) = mean(mean_seg1); 
    mean_values(i,2) = mean(mean_seg2);
    
    % Store individual column means (for scatter overlay)
    individual_means{i,1} = mean_seg1;
    individual_means{i,2} = mean_seg2;
end

% Plot bar graph
figure; hold on;
bar_handle = bar(mean_values, 'grouped');

% Set bar colors
bar_handle(1).FaceColor = 'b'; % First cell (blue)
bar_handle(2).FaceColor = 'r'; % Second cell (red)

% Get X positions for scatter overlay
x_positions = zeros(3,2); % Store bar center positions
for i = 1:2
    x_positions(:,i) = bar_handle(i).XEndPoints; % Correct x positions
end

% Overlay individual data points correctly aligned
for i = 1:3
    scatter(x_positions(i,1) * ones(1,10), individual_means{i,1}, 'k', 'filled'); % First dataset
    scatter(x_positions(i,2) * ones(1,10), individual_means{i,2}, 'k', 'filled'); % Second dataset
end

% Labels and legend
xticks(mean(x_positions,2)); % Align x-ticks with bar groups
xticklabels({'0%', '50%', '75%'});
ylabel('Mean PV correlation');
legend({'Safe-identified', 'Risky-identified'}, 'Location', 'SouthEast');
ytickformat('%.2f');

hold off;