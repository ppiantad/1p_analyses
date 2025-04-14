% load dataset that contains data where events have been identified

clear similarityMatrix

data_to_load = {'BLA_RDT_11_variables_incl_AA_use_for_AA_PV.mat', 'BLA_C_raw_no_additional_filtering_Pre_RDT_RM_only_completed_sessions_zall_window_base_workspace_10_categories.mat', 'BLA_RM_D1_6_categories_03212025_including_spike_prob_array_counts.mat'};

load(data_to_load{1})


if size(respClass_all_array, 2) == 10 | size(respClass_all_array, 2) == 11 | size(respClass_all_array, 2) == 12
    % comparison_arrays_full = [1 2 3; 8 9 10]
    comparison_arrays_full = [1 2 3; 5 6 7]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays_full = [1 2 3; 4 5 6]
end

% depending on analysis desired, may need to run upper part, then run
% block_wise_changes_v1.m, then run below

% some data need the following uncommented:
% ca_data_type = uv.ca_data_type;




for gg = 1:size(animalIDs, 1)
    select_mouse = animalIDs{gg};
    first_session = session_to_analyze;
    select_mouse_index = find(strcmp(animalIDs, select_mouse));

    BehavData = final.(select_mouse).(first_session).uv.BehavData;


    %%
    PV_AA_all_mouse = [];
    for ff = 1:size(caTraceTrials_mouse{select_mouse_index, 11}, 2)
        for cc = 1:size(caTraceTrials_mouse{select_mouse_index, 11}{1, ff}, 1)
            PV_AA_all_mouse(cc, ff) = mean(caTraceTrials_mouse{select_mouse_index, 11}{1, ff}(cc, ts1 >= 0 & ts1 <= 2));

        end

    end

    % reorganize data for plotting the heatmap, with ensemble-specific
    % neurons on top, the rest below
    PV_prechoice_all_mouse_just_shock = PV_AA_all_mouse(:, respClass_all_array_mouse{select_mouse_index, 11} == 1);
    PV_prechoice_all_mouse_just_shock = [PV_prechoice_all_mouse_just_shock PV_AA_all_mouse(:, respClass_all_array_mouse{select_mouse_index, 11} ~= 1)];
    PV_prechoice_all_mouse_just_prechoice_array{gg} = PV_prechoice_all_mouse_just_shock;
    % data to use below, without reorg
    mean_PV_AA_all_mouse = mean(PV_AA_all_mouse)';

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

    hold on;
    figure; imagesc(1:size(PV_prechoice_all_mouse_just_shock, 1), [], PV_prechoice_all_mouse_just_shock')
    colormap gray;
    colorbar;
    clim([0 1]);
    hold off;

    % Define the position for the first green line
    y1_start = 1;  % Start at the first neuron
    y1_end = sum(respClass_all_array_mouse{select_mouse_index, 11} == 1);  % End at the number of neurons in the first group

    % Define the position for the second line
    y2_start = y1_end + 1;  % Start right after the first group
    y2_end = size(PV_prechoice_all_mouse_just_shock, 2);  % End at the last neuron

    % Draw the first vertical green line
    hold on;
    line([1 1], [y1_start y1_end], 'Color', 'g', 'LineWidth', 2);

    % Draw the second vertical line (can be any color, e.g., black here)
    line([1 1], [y2_start y2_end], 'Color', 'r', 'LineWidth', 2);
    title(['From animal ',select_mouse], 'Interpreter', 'none')


    hold off;


    ca = final.(select_mouse).(first_session).CNMFe_data.(ca_data_type);

    % filter based on variable above
    % ca = ca(selected_rows, :);
    % ca = ca(prechoice_conserved_mouse{gg, 1} ~= 1 & prechoice_lost_mouse{gg, 1} ~= 1 & prechoice_remapped_mouse{gg, 1} ~= 1, :);
    % ca_zscored = zscore(ca, [], 2);
    ca_zscored = ca;

    % ca_zscored = ca(prechoice_lost_mouse{gg, 1} == 1, :);
    time_array = final.(select_mouse).(first_session).time;
    shock_similarityOverTime = [];
    for t = 1:size(ca_zscored, 2)
        % for i = 1:length(uniqueTypes)
        % typeIdx = (neuronTypes == uniqueTypes(i));
        % use specific subset of neurons
        activitySubset = ca_zscored(:, t);
        % uncomment if you want to use all neurons
        %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
        similarityMatrix = corrcoef(activitySubset, mean_PV_AA_all_mouse);
        shock_similarityOverTime(t) = similarityMatrix(2);
        % end
    end
    %
    figure; plot(time_array, shock_similarityOverTime(1, 1:size(time_array, 1)))
    xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
    xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
    xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
    xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
    if strcmp('shock',BehavData.Properties.VariableNames)
        xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
    end




    %% attempting to get something similar to Courtin Fig 5e
    shock_similarityOverTimeTrials = [];


    eTS = [BehavData.choiceTime BehavData.collectionTime];

    % Define time windows
    pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s


    % Find indices corresponding to each time window
    pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);


    for t = 1:size(eTS,1)

        timeWin_prechoice_post_choice = [eTS(t, 1)+uv.evtWin(1,1):uv.dt:eTS(t, 1)+uv.evtWin(1,2)];  %calculate time window around each event
        BL_win = [eTS(t, 1)+uv.BLper(1,1):uv.dt:eTS(t, 1)+uv.BLper(1,2)];



        if min(timeWin_prechoice_post_choice) > min(time_array) && max(timeWin_prechoice_post_choice) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
            % get unit event counts in trials
            % get unit ca traces in trials

            idx_prechoice_post_choice = time_array > min(timeWin_prechoice_post_choice) & time_array < max(timeWin_prechoice_post_choice);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
            bl_idx = time_array > min(BL_win) & time_array < max(BL_win);
            %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
            shock_similarityOverTimeTrials(t,1:sum(idx_prechoice_post_choice)) = shock_similarityOverTime(idx_prechoice_post_choice);
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

    block_2_3_inds = BehavData.Block ~= 1 & BehavData.Blank_Touch == 0 & BehavData.omissionALL == 0;

    prechoice_period_b1_mean = mean(shock_similarityOverTimeTrials(block_1_inds, pre_choice_indices), 2);
    prechoice_period_b2_mean = mean(shock_similarityOverTimeTrials(block_2_inds, pre_choice_indices), 2);
    prechoice_period_b3_mean = mean(shock_similarityOverTimeTrials(block_3_inds, pre_choice_indices), 2);

    prechoice_block_2_3_mean = mean(shock_similarityOverTimeTrials(block_2_3_inds, pre_choice_indices), 2);

    AA_inds = BehavData.type_binary == 1;

    AA_mean = mean(shock_similarityOverTimeTrials(AA_inds, pre_choice_indices), 2);



    prechoice_b2_3_over_time_mean{gg} = prechoice_block_2_3_mean;
    aa_over_time_mean{gg} = AA_mean;




    prechoice_period_b1_sem = std(shock_similarityOverTimeTrials(block_1_inds, pre_choice_indices), 0, 2) / sqrt(length(find(block_1_inds == 1)));
    prechoice_period_b2_sem = std(shock_similarityOverTimeTrials(block_2_inds, pre_choice_indices), 0, 2) / sqrt(length(find(block_2_inds == 1)));
    prechoice_period_b3_sem = std(shock_similarityOverTimeTrials(block_3_inds, pre_choice_indices), 0, 2) / sqrt(length(find(block_3_inds == 1)));


    prechoice_over_time_sem{gg} = [prechoice_period_b1_sem; prechoice_period_b2_sem; prechoice_period_b3_sem];


    prechoice_similarityOverTime_array{gg} = shock_similarityOverTime;



    prechoice_similarityOverTimeTrials_array{gg} = shock_similarityOverTimeTrials;


end
%%



mean_prechoice_on_AA = cellfun(@mean,aa_over_time_mean) 
mean_prechoice_on_choices = cellfun(@mean,prechoice_b2_3_over_time_mean) 



% Combine data into a matrix
data = [mean_prechoice_on_AA(:), mean_prechoice_on_choices(:)];
group_labels = {'AA', 'Block 2/3 choices'};

% Calculate mean and SEM
mean_vals = mean(data, 1);
sem_vals = std(data, 0, 1) ./ sqrt(size(data, 1));

% Create bar plot
figure;
hold on;
b = bar(mean_vals, 'FaceColor', 'flat');

% Add error bars
errorbar(1:2, mean_vals, sem_vals, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Overlay scatter plot of individual points with slight jitter
jitter = 0.1;
for i = 1:2
    x = i + (rand(size(data, 1), 1) - 0.5) * jitter;
    scatter(x, data(:, i), 30, 'filled', 'MarkerFaceAlpha', 0.6);
end

% Format plot
set(gca, 'XTick', 1:2, 'XTickLabel', group_labels);
ylabel('Mean PV corr.');





