% load dataset that contains data where events have been identified
% run raster_with_representative_neurons_v2.m for the mouse for whom you
% will analyze below
% also you should first run block_wise_changes_v1.m to get the blocktimes
clear similarityMatrix

data_to_load = {'BLA_C_raw_no_additional_filtering_RDT_D1_only_completed_sessions_zall_window_base_workspace_10_categories.mat', 'BLA_C_raw_no_additional_filtering_Pre_RDT_RM_only_completed_sessions_zall_window_base_workspace_10_categories.mat'};

load(data_to_load{1})

if size(respClass_all_array, 2) == 10
    comparison_arrays_full = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays_full = [1 2 3; 4 5 6]
end


%%
for aa = 1:size(comparison_arrays_full, 1)
    
    comparison_arrays = comparison_arrays_full(aa, :)
    first_session = session_to_analyze;
    

    for ii = 1:size(animalIDs,1)
        currentanimal = char(animalIDs(ii));
        if isfield(final.(currentanimal), session_to_analyze)
            BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
            % only use rewarded trials for this, otherwise things get wonky
            [BehavData,trials,varargin]=TrialFilter_test(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0);
            block_1 = [BehavData.stTime(BehavData.Block == 1) BehavData.collectionTime(BehavData.Block == 1)];
            block_1_mouse(ii,:) = [block_1(1, 1) block_1(end, 2)];
            block_2 = [BehavData.stTime(BehavData.Block == 2) BehavData.collectionTime(BehavData.Block == 2)];
            block_2_mouse(ii,:) = [block_2(1, 1) block_2(end, 2)];
            block_3 = [BehavData.stTime(BehavData.Block == 3) BehavData.collectionTime(BehavData.Block == 3)];
            block_3_mouse(ii,:) = [block_3(1, 1) block_3(end, 2)];

        end
    end
    for gg = 1:size(animalIDs, 1)
        select_mouse = animalIDs{gg};

        % for RDT D1 BLA_Insc_25:
        %prechoice neuron num 46
        %postchoice rew num 38
        %consumption num 39
        %shock num 11

        select_mouse_index = find(strcmp(animalIDs, select_mouse));

        % first_session = 'RDT_D1';


        neuronTypes = respClass_all_array_mouse{select_mouse_index, 3};
        neuralActivity = final.(select_mouse).(first_session).CNMFe_data.C_raw;
        BehavData = final.(select_mouse).(first_session).uv.BehavData;


        %%
        PV_prechoice_all_mouse = [];
        for ff = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 1)}, 2)
            for cc = 1:size(zall_mouse{select_mouse_index,comparison_arrays(1, 1)}{1, ff}, 1)
                PV_prechoice_all_mouse(cc, ff) = mean(zall_mouse{select_mouse_index, comparison_arrays(1, 1)}{1, ff}(cc, ts1 >= -4 & ts1 <= 0));

            end

        end


        PV_prechoice_all_mouse_just_prechoice = PV_prechoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} == 1);
        PV_prechoice_all_mouse_just_prechoice = [PV_prechoice_all_mouse_just_prechoice PV_prechoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} ~= 1)];
        mean_PV_prechoice_all_mouse = mean(PV_prechoice_all_mouse)';
        % mean_PV_prechoice_all_mouse = mean_PV_prechoice_all_mouse(respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} == 1, :);
        
        % hold on;
        % figure; imagesc(1:size(PV_prechoice_all_mouse_just_prechoice, 1), [], PV_prechoice_all_mouse_just_prechoice')
        % colormap gray;
        % colorbar;
        % clim([0 1.5]);
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
        % 
        % hold off;


        ca = final.(select_mouse).(first_session).CNMFe_data.(ca_data_type);
        % ca = ca(respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} == 1, :);
        ca_zscored = zscore(ca, [], 2);
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
        % 
        % b1_prechoice_corr = mean(prechoice_similarityOverTime(1, time_array > block_1_mouse(select_mouse_index, 1) & time_array <= block_1_mouse(select_mouse_index, 2)));
        % b2_prechoice_corr = mean(prechoice_similarityOverTime(1, time_array > block_2_mouse(select_mouse_index, 1) & time_array <= block_2_mouse(select_mouse_index, 2)));
        % b3_prechoice_corr = mean(prechoice_similarityOverTime(1, time_array > block_3_mouse(select_mouse_index, 1) & time_array <= block_3_mouse(select_mouse_index, 2)));
        % 

        %%

        PV_postchoice_all_mouse = [];
        for ff = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 2)}, 2)
            for cc = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 2)}{1, ff}, 1)
                PV_postchoice_all_mouse(cc, ff) = mean(zall_mouse{select_mouse_index, comparison_arrays(1, 2)}{1, ff}(cc, ts1 >= 0 & ts1 <= 2));

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
        ca_zscored = zscore(ca, [], 2);
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
        % 
        % 
        % 
        % b1_postchoice_corr = mean(postchoice_similarityOverTime(1, time_array > block_1_mouse(select_mouse_index, 1) & time_array <= block_1_mouse(select_mouse_index, 2)));
        % b2_postchoice_corr = mean(postchoice_similarityOverTime(1, time_array > block_2_mouse(select_mouse_index, 1) & time_array <= block_2_mouse(select_mouse_index, 2)));
        % b3_postchoice_corr = mean(postchoice_similarityOverTime(1, time_array > block_3_mouse(select_mouse_index, 1) & time_array <= block_3_mouse(select_mouse_index, 2)));
        % 
        % 


        %%

        PV_consumption_all_mouse = [];
        for ff = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 3)}, 2)
            for cc = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 3)}{1, ff}, 1)
                PV_consumption_all_mouse(cc, ff) = mean(zall_mouse{select_mouse_index, comparison_arrays(1, 3)}{1, ff}(cc, ts1 >= 1 & ts1 <= 3));

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
        ca_zscored = zscore(ca, [], 2);
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
        % 
        % 
        % b1_consumption_corr = mean(consumption_similarityOverTime(1, time_array > block_1_mouse(select_mouse_index, 1) & time_array <= block_1_mouse(select_mouse_index, 2)));
        % b2_consumption_corr = mean(consumption_similarityOverTime(1, time_array > block_2_mouse(select_mouse_index, 1) & time_array <= block_2_mouse(select_mouse_index, 2)));
        % b3_consumption_corr = mean(consumption_similarityOverTime(1, time_array > block_3_mouse(select_mouse_index, 1) & time_array <= block_3_mouse(select_mouse_index, 2)));





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

%%

figure;
hold on
% Create a histogram for allCorrelations

width = 200; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
% xlim([-8 8]);
% % Set X-axis ticks
% set(gca, 'XTick', [-8, 0, 8]);

shadedErrorBar(1:90, mean_prechoice_over_time(:, 1), sem_prechoice_over_time(:, 1), 'lineProps', {'color', 'r'});
hold on; shadedErrorBar(1:90, mean_prechoice_over_time(:, 2), sem_prechoice_over_time(:, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});
%
% xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
% xlabel('Time from Large Rew Choice (s)');
% legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
% ylim([-0.8 0.8]);
% ytickformat('%.1f');
hold off

%%

figure;
hold on
% Create a histogram for allCorrelations

width = 200; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
% xlim([-8 8]);
% % Set X-axis ticks
% set(gca, 'XTick', [-8, 0, 8]);

shadedErrorBar(1:90, mean_postchoice_over_time(:, 1), sem_postchoice_over_time(:, 1), 'lineProps', {'color', 'r'});
hold on; shadedErrorBar(1:90, mean_postchoice_over_time(:, 2), sem_postchoice_over_time(:, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});
%
% xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
% xlabel('Time from Large Rew Choice (s)');
% legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
% ylim([-0.8 0.8]);
% ytickformat('%.1f');
hold off
%%

figure;
hold on
% Create a histogram for allCorrelations

width = 200; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
% xlim([-8 8]);
% % Set X-axis ticks
% set(gca, 'XTick', [-8, 0, 8]);

shadedErrorBar(1:90, mean_consumption_over_time(:, 1), sem_consumption_over_time(:, 1), 'lineProps', {'color', 'r'});
hold on; shadedErrorBar(1:90, mean_consumption_over_time(:, 2), sem_consumption_over_time(:, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});
%
% xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
% xlabel('Time from Large Rew Choice (s)');
% legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
% ylim([-0.8 0.8]);
% ytickformat('%.1f');
hold off