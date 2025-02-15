%% Edit these user variables with what you want to look at
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5]
% uv.BLper = [-10 -5];
uv.BLper = [uv.evtWin(1) uv.evtWin(1)+3];
uv.dt = 0.1; %what is your frame rate
% uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at

ca_data_type = "C_raw"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate

session_to_analyze = 'RDT_D1';
uv.yoke_data = 0; % set to 1 if you want to be prompted to yoke the number of trials analyzed, set to 0 otherwise

epoc_to_align = 'choiceTime'; % stTime choiceTime collectionTime
period_of_interest = 'postchoice';

if strcmp(epoc_to_align, 'stTime')
    period_of_interest = 'trial_start';
    uv.evtSigWin.outcome = [-1 1]; %for trial start
elseif strcmp(epoc_to_align, 'choiceTime')
    if strcmp(period_of_interest, 'prechoice')
        uv.evtSigWin.outcome = [-4 0]; %for pre-choice   [-4 0]    [-4 1]
    elseif strcmp(period_of_interest, 'postchoice')
        uv.evtSigWin.outcome = [0 2]; %for SHK or immediate post-choice [0 2]
    end
elseif strcmp(epoc_to_align, 'collectionTime')
    period_of_interest = 'reward_collection';
    uv.evtSigWin.outcome = [1 3]; %for REW collection [1 3]
end

ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
animalIDs = (fieldnames(final));
neuron_num = 0;
use_normalized_time = 0;
clear neuron_mean neuron_sem neuron_num zall_mean zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all 


% these are mice that did not complete the entire session - kinda have to
% toss them to do some comparisons during RDT


if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final, fieldsToRemove{i})
            final = rmfield(final, fieldsToRemove{i});
        end
    end
elseif strcmp('RDT_D2', session_to_analyze)

    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final, fieldsToRemove{i})
            final = rmfield(final, fieldsToRemove{i});
        end
    end
end


%% get trials that are aborts, plus the previous trial

neuron_num = 0;
animalIDs = (fieldnames(final));
for ii = 1:size(animalIDs, 1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        % Extract the table for easy reference
        data = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        % Initialize counters for trials where a shock is followed by an abort or no abort
        rows_for_new_table = 1;

        % Loop through each row
        i = 1;
        while i <= height(data) - 1
            % Check if the current trial is a shock (shock == 1)
            if data.type_binary(i) == 1
                % Initialize a flag to determine if we found a matching condition
                match_found = false;
                behav_data_extracted(rows_for_new_table, :) = data(i, :);
                rows_for_new_table = rows_for_new_table + 1;

                % Look for a matching condition in the subsequent rows
                j = i - 1; % Check previous row
                while j >= 1 && ~match_found
                    % Check if the previous trial is a real trial (not abort or blank touch)
                    if data.Blank_Touch(j) ~= 1 && data.type_binary(j) ~= 1 && data.type_binary(j) ~= 2
                        match_found = true;
                        behav_data_extracted(rows_for_new_table, :) = data(j, :);
                        rows_for_new_table = rows_for_new_table + 1;
                    else
                        j = j - 1; % Move further back
                    end
                end

                % Find the first subsequent row where bigSmall == 1.2 or 0.3
                found_bigSmall = false;
                k = i + 1;
                while k <= height(data) && ~found_bigSmall
                    if data.bigSmall(k) == 1.2 || data.bigSmall(k) == 0.3
                        found_bigSmall = true;
                        i = k; % Set i to this row and continue
                    else
                        k = k + 1;
                    end
                end

                % If no match was found, increment i normally
                if ~found_bigSmall
                    i = i + 1;
                end
            else
                % Move to the next row if no shock was found
                i = i + 1;
            end
        end

        % Store extracted data
        behav_data_extracted_array{ii} = behav_data_extracted;
        clear behav_data_extracted
    end
end



%% **use this code only if you want to select all aborts (not just the first abort following a punished or un-punished trials, which is what the above code does)**

for ii = 1:size(animalIDs, 1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        % Extract the table for easy reference
        data = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData  ;

        % Initialize counters for trials where a shock is followed by an abort or no abort
        total_shock_abort_trials = 0;
        total_shock_no_abort_trials = 0;
        rows_for_new_table = 1;
        % Loop through each row
        i = 1;
        while i <= height(data) - 1
            % Check if the current trial is a shock (shock == 1)
            if data.type_binary(i) == 1
                % Initialize a flag to determine if we found a matching condition
                match_found = false;
                behav_data_extracted(rows_for_new_table, :) = data(i, :);
                rows_for_new_table = rows_for_new_table + 1;
                % Look for a matching condition in the subsequent rows
                j = i + -1;
                while j <= height(data) && ~match_found
                    % check back in prior to trials, track trials where the
                    % previous trial was a trial proper (not an abort or
                    % blank touch
                    if data.Blank_Touch(j) ~= 1 && data.type_binary(j) ~= 1 && data.type_binary(j) ~= 2
                        % Increment abort counter and mark as found
                        match_found = true;
                        % behav_data_shk_to_abort(rows_for_new_table+1, :) = data([i, j], :);
                        behav_data_extracted(rows_for_new_table, :) = data(j, :);
                        rows_for_new_table = rows_for_new_table + 1;
                        % Check if the trial has a bigSmall event (bigSmall == 1)
                    % elseif data.bigSmall(j) == 1.2 || data.bigSmall(j) == 0.3
                    %     % Increment no-abort counter and mark as found
                    %     total_shock_no_abort_trials = total_shock_no_abort_trials + 1;
                    %     match_found = true;
                    %     behav_data_extracted(rows_for_new_table, :) = data(j, :);
                    %     rows_for_new_table = rows_for_new_table + 1;
                    else
                        % Move to the next row if no match is found
                        j = j - 1;
                    end
                end
                
                % Set i to j to continue from the row after the found condition
                i = i + 1;
            else
                % Move to the next row if no shock was found
                i = i + 1;
            end
        end

        % % Display the results
        % disp(['Total number of trials where a shock was followed by an abort: ', num2str(total_shock_abort_trials)]);
        % disp(['Total number of trials where a shock was followed by a bigSmall event (no abort): ', num2str(total_shock_no_abort_trials)]);
        % total_shock_abort_trials_array(ii) = total_shock_abort_trials;
        % total_shock_no_abort_trials_array(ii) = total_shock_no_abort_trials;
        behav_data_extracted_array{ii} = behav_data_extracted;
        clear behav_data_extracted
    end
end

%% Get Ca data
iter = 0;
iter = iter+1;
neuron_num = 0;
animalIDs = fieldnames(final);

for ii = 1:numel(animalIDs)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = behav_data_extracted_array{1, ii};
        ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);

        % Comment out if you don't want to z-score traces prior to analysis
        % ca = zscore(ca, 0, 2);

        if strcmp(ca_data_type, 'S')
            ca = full(ca);
        end

        num_samples = size(ca, 2);
        sampling_frequency = final.(currentanimal).(session_to_analyze).uv.dt * 100;
        time_array = final.(currentanimal).(session_to_analyze).time;

        eTS = BehavData.(epoc_to_align); % Get timestamps
        zb_session = mean(ca, 2);
        zsd_session = std(ca, [], 2);

        % Calculate time windows for each event
        evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
        numMeasurements = round(evtWinSpan / uv.dt); % Rounding due to odd frame rate

        % Extract only odd rows
        odd_indices = 1:2:size(eTS, 1)-1; % Ensuring not to exceed size
        eTS_filtered = eTS(odd_indices);

        % Get the subsequent row's shock values
        subsequent_shock_values = BehavData.shock(odd_indices + 1);  

        % Initialize arrays
        caTraceTrials_shock = {};
        caTraceTrials_no_shock = {};
        zall_mean_all_shock = [];
        zall_mean_all_no_shock = [];

        for u = 1:size(ca, 1)
            neuron_num = neuron_num + 1;
            caTraceTrials = NaN(length(eTS_filtered), numMeasurements);
            unitTrace = ca(u, :);

            if isempty(eTS_filtered) || length(eTS_filtered) == 1
                caTraceTrials(1, 1:size(ts1, 2)) = nan;
                zall(1, 1:size(ts1, 2)) = nan;
                sem_all(neuron_num, :) = nan;
                zall_mean_all(neuron_num, :) = nan;
            else
                % Align and process calcium traces
                [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times] = ...
                    align_and_zscore(BehavData, unitTrace, eTS_filtered, uv, time_array, zb_session, zsd_session, u, use_normalized_time);

                StartChoiceCollect_times_array{neuron_num} = StartChoiceCollect_times;
                caTraceTrials = caTraceTrials(:, 1:size(ts1, 2));
                zall = zall_window(:, 1:size(ts1, 2));

                % Apply Savitzky-Golay filter
                for z = 1:size(zall, 1)
                    zall(z, :) = sgolayfilt(zall(z, :), 9, 21);
                end

                % Store categorized data based on subsequent row's shock value
                shock_mask = subsequent_shock_values == 1; % Logical array where shock = 1
                no_shock_mask = ~shock_mask; % Logical array where shock ≠ 1 (or NaN)

                % Allocate trials to separate arrays
                caTraceTrials_shock{neuron_num} = zall(shock_mask, :);
                caTraceTrials_no_shock{neuron_num} = zall(no_shock_mask, :);

                % Store other necessary variables
                zall_array{neuron_num} = zall;
                zall_mean_all(neuron_num, :) = nanmean(zall, 1);
                sem_all(neuron_num, :) = nanstd(zall, 1) / sqrt(size(zall, 1));
                caTraceTrials_shock_mean(neuron_num, :) = nanmean(zall(shock_mask, :));
                caTraceTrials_no_shock_mean(neuron_num, :) = nanmean(zall(no_shock_mask, :));
                clear zall caTraceTrials zb zsd;
            end
        end
    end
end

% Store outputs in cell arrays
zall_mean_all_array(iter) = {zall_mean_all};
sem_all_array(iter) = {sem_all};
caTraceTrials_shock_all{iter} = caTraceTrials_shock;
caTraceTrials_no_shock_all{iter} = caTraceTrials_no_shock;

%%
figure; plot(ts1, mean(caTraceTrials_shock_mean))
hold on; plot(ts1, mean(caTraceTrials_no_shock_mean))


%% Run eventRelatedActivityAndClassification on choiceTime.postchoice_0to2.AA_1

% update respClass_all_array location depending on how many variables you
% have also analyzed using eventRelatedActivity
figure; plot(ts1, mean(caTraceTrials_shock_mean(respClass_all_array{1,2} == 1, :)));
hold on; plot(ts1, mean(caTraceTrials_no_shock_mean(respClass_all_array{1,2} == 1, :)))

