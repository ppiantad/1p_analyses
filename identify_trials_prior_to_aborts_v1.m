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
                % Set `i` to `j` to continue from the row after the found condition
                i = i + 1;
            else
                % Move to the next row if no shock was found
                i = i + 1;
            end
        end

        % Display the results
        disp(['Total number of trials where a shock was followed by an abort: ', num2str(total_shock_abort_trials)]);
        disp(['Total number of trials where a shock was followed by a bigSmall event (no abort): ', num2str(total_shock_no_abort_trials)]);
        total_shock_abort_trials_array(ii) = total_shock_abort_trials;
        total_shock_no_abort_trials_array(ii) = total_shock_no_abort_trials;
        behav_data_extracted_array{ii} = behav_data_extracted;
        clear behav_data_extracted
    end
end

% shk_abort_to_shk_choice_ratio = total_shock_abort_trials_array./total_shock_no_abort_trials_array