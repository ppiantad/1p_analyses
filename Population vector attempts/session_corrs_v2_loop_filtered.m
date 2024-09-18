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
elseif size(respClass_all_array, 2) == 7
    comparison_arrays_full = [1 2 3; 5 6 7]
end

% ca_data_type = uv.ca_data_type

% Assume 'ca' is your calcium imaging data, sampled every 100 ms
% bin_size is the desired bin size in ms
bin_size = 100;  % Set this as needed
bin_factor = bin_size / (uv.dt*1000);  % Determine how many samples to bin together

ts1 = (uv.evtWin(1):bin_factor/10:uv.evtWin(2)-0.1);
%%
for aa = 1:size(comparison_arrays_full, 1) %size(comparison_arrays_full, 1)

    comparison_arrays = comparison_arrays_full(aa, :)
    first_session = session_to_analyze;

    for gg = 1:size(animalIDs, 1) %size(animalIDs, 1)
        % if gg ~= 2 && gg ~= 3 && gg ~= 6
        select_mouse = animalIDs{gg};

        select_mouse_index = find(strcmp(animalIDs, select_mouse));
        BehavData = final.(select_mouse).(first_session).uv.BehavData;


        ca = final.(select_mouse).(first_session).CNMFe_data.(ca_data_type);
        % ca = ca(prechoice_indices_for_PV{aa, select_mouse_index} == 1, :);
        % ca_zscored = zscore(ca, [], 2);
        ca_zscored = normalize(ca, 2);
        time_array = final.(select_mouse).(first_session).time;

        % Reshape the data and take the mean of each bin
        [n_neurons, n_samples] = size(ca_zscored);
        ca_binned = squeeze(mean(reshape(ca_zscored(:, 1:floor(n_samples/bin_factor)*bin_factor), n_neurons, bin_factor, []), 2));

        % Squeeze the data to remove the singleton dimension
        % ca_binned = squeeze(ca_binned);
        ca = ca_binned;
        % Assuming time_array is a 1D array of time points corresponding to each sample
        time_array_binned = mean(reshape(time_array(1:floor(length(time_array)/bin_factor)*bin_factor), bin_factor, []), 1);
        time_array = time_array_binned;



        
        
        %%
        PV_prechoice_all_mouse = [];
        for ff = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 1)}, 2)
            for cc = 1:size(zall_mouse{select_mouse_index,comparison_arrays(1, 1)}{1, ff}, 1)
                PV_prechoice_all_mouse(cc, ff) = mean(zall_mouse{select_mouse_index, comparison_arrays(1, 1)}{1, ff}(cc, ts1 >= -4 & ts1 <= 0));
                if isnan(PV_prechoice_all_mouse(cc, ff))
                    PV_prechoice_all_mouse(cc, ff) = 0;
                end

            end
        end


        % PV_prechoice_all_mouse_just_prechoice = PV_prechoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} == 1);
        % PV_prechoice_all_mouse_just_prechoice = [PV_prechoice_all_mouse_just_prechoice PV_prechoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 1)} ~= 1)];
        mean_PV_prechoice_all_mouse = mean(PV_prechoice_all_mouse)';
        % mean_PV_prechoice_all_mouse = mean_PV_prechoice_all_mouse(prechoice_indices_for_PV{aa, select_mouse_index} == 1, :);
        mean_PV_prechoice_all_mouse_array{aa, gg} = mean_PV_prechoice_all_mouse;

        



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


        %%

        PV_postchoice_all_mouse = [];
        for ff = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 2)}, 2)
            for cc = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 2)}{1, ff}, 1)
                PV_postchoice_all_mouse(cc, ff) = mean(zall_mouse{select_mouse_index, comparison_arrays(1, 2)}{1, ff}(cc, ts1 >= 0 & ts1 <= 2));
                if isnan(PV_postchoice_all_mouse(cc, ff))
                    PV_postchoice_all_mouse(cc, ff) = 0;
                end
            end

        end


        % PV_postchoice_all_mouse_just_postchoice = PV_postchoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 2)} == 1);
        % PV_postchoice_all_mouse_just_postchoice = [PV_postchoice_all_mouse_just_postchoice PV_postchoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 2)} ~= 1)];
        mean_PV_postchoice_all_mouse = mean(PV_postchoice_all_mouse)';
        % mean_PV_postchoice_all_mouse = mean_PV_postchoice_all_mouse(postchoice_indices_for_PV{aa, select_mouse_index} == 1, :);
        mean_PV_postchoice_all_mouse_array{aa, gg} = mean_PV_postchoice_all_mouse;



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

        %%

        PV_consumption_all_mouse = [];
        for ff = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 3)}, 2)
            for cc = 1:size(zall_mouse{select_mouse_index, comparison_arrays(1, 3)}{1, ff}, 1)
                PV_consumption_all_mouse(cc, ff) = mean(zall_mouse{select_mouse_index, comparison_arrays(1, 3)}{1, ff}(cc, ts1 >= 1 & ts1 <= 3));
                if isnan(PV_consumption_all_mouse(cc, ff))
                    PV_consumption_all_mouse(cc, ff) = 0;
                end
            end

        end


        % PV_consumption_all_mouse_just_consumption = PV_consumption_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 3)} == 1);
        % PV_consumption_all_mouse_just_consumption = [PV_consumption_all_mouse_just_consumption PV_consumption_all_mouse(:, respClass_all_array_mouse{select_mouse_index, comparison_arrays(1, 3)} ~= 1)];
        mean_PV_consumption_all_mouse = mean(PV_consumption_all_mouse)';
        % mean_PV_consumption_all_mouse = mean_PV_consumption_all_mouse(consumption_indices_for_PV{aa, select_mouse_index} == 1, :);
        mean_PV_consumption_all_mouse_array{aa, gg} = mean_PV_consumption_all_mouse;

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

            end


        end

        block_1_inds = [];
        block_2_inds = [];
        block_3_inds = [];
        prechoice_period_b1_mean = [];
        prechoice_period_b2_mean = [];
        prechoice_period_b3_mean = [];
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
        % else
        %     prechoice_over_time_sem(:, gg) = NaN;
        %     postchoice_over_time_sem(:, gg) = NaN;
        %     consumption_over_time_sem(:, gg) = NaN;
        %
        %     prechoice_over_time_mean(:, gg) = NaN;
        %     postchoice_over_time_mean(:, gg) = NaN;
        %     consumption_over_time_mean(:, gg) = NaN;
        %     continue
        % end
    end
    %%
    mean_prechoice_over_time(:, aa) = nanmean(prechoice_over_time_mean, 2);
    mean_postchoice_over_time(:, aa) = nanmean(postchoice_over_time_mean, 2);
    mean_consumption_over_time(:, aa) = nanmean(consumption_over_time_mean, 2);


    sem_prechoice_over_time(:, aa) = nanmean(prechoice_over_time_sem, 2);
    sem_postchoice_over_time(:, aa) = nanmean(postchoice_over_time_sem, 2);
    sem_consumption_over_time(:, aa) = nanmean(consumption_over_time_sem, 2);



end

%%

figure;
hold on
% Create a histogram for allCorrelations

width = 400; % Width of the figure
height = 250; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]

shadedErrorBar(1:90, mean_prechoice_over_time(:, 1), sem_prechoice_over_time(:, 1), 'lineProps', {'color', 'r'});
hold on; shadedErrorBar(1:90, mean_prechoice_over_time(:, 2), sem_prechoice_over_time(:, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});
xlim([1 90]);
% Set X-axis ticks
set(gca, 'XTick', [1, 30, 60, 90]);

xline(30);
xline(60)
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

width = 400; % Width of the figure
height = 250; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]


shadedErrorBar(1:90, mean_postchoice_over_time(:, 1), sem_postchoice_over_time(:, 1), 'lineProps', {'color', 'r'});
hold on; shadedErrorBar(1:90, mean_postchoice_over_time(:, 2), sem_postchoice_over_time(:, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});
%
xlim([1 90]);
% Set X-axis ticks
set(gca, 'XTick', [1, 30, 60, 90]);

xline(30);
xline(60)
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

width = 400; % Width of the figure
height = 250; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]


shadedErrorBar(1:90, mean_consumption_over_time(:, 1), sem_consumption_over_time(:, 1), 'lineProps', {'color', 'r'});
hold on; shadedErrorBar(1:90, mean_consumption_over_time(:, 2), sem_consumption_over_time(:, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(conserved_consumption==1, :)), nanmean(neuron_sem_array{1, 8}(conserved_consumption==1, :)), 'lineProps', {'color', 'b'});
xlim([1 90]);
% Set X-axis ticks
set(gca, 'XTick', [1, 30, 60, 90]);

xline(30);
xline(60)
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
% xlabel('Time from Large Rew Choice (s)');
% legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
% ylim([-0.8 0.8]);
% ytickformat('%.1f');
hold off