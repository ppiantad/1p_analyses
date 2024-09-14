%% overall, this does not seem to change much


uv.ca_data_type = "C_raw"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate

session_to_analyze = 'Pre_RDT_RM';
uv.yoke_data = 0; % set to 1 if you want to be prompted to yoke the number of trials analyzed, set to 0 otherwise

epoc_to_align = 'choiceTime'; % stTime choiceTime collectionTime
period_of_interest = 'prechoice';

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

ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(uv.ca_data_type);

if strcmp(uv.ca_data_type, 'S')
    ca = full(ca);

end
ca_test = ca(90, :);
[~, num_samples] = size(ca);

currentanimal = char(animalIDs(1));
BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
[BehavData,trials,varargin_identity_class]=TrialFilter_test(BehavData, 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1); %'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1    % 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 2, 'BLOCK', 3
shuffled_data = zeros(uv.resamples, num_samples); % Preallocate matrix for efficiency
eTS = BehavData.(epoc_to_align); %get time stamps
time_array = final.(currentanimal).(session_to_analyze).time;


for hh = 1:uv.resamples %1:uv.resamples*10


    shift_val = randi(num_samples); % Generate a random shift value for each signal RUAIRI RECOMMENDED KEEPING THE SAME SHIFT VAL, rather than randomizing per neuron. this is because then you keep the overall correlation b/w the neurons, but disrupt the relationship to the event timestamps

    shuffled_data(hh,:) = circshift(ca_test, shift_val,2); % Perform the circular shuffle




end

for ff = 1:size(ca_test, 1)
    unitTrace = ca_test; %get trace
    [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);
    % [caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials] = align_only(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);
    caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1

    % if strcmp(uv.zscore_to, 'window')
        zall = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
    % elseif strcmp(uv.zscore_to, 'session')
        % zall = zall_session(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
    % elseif strcmp(uv.zscore_to, 'baseline')
        % zall = zall_baselined(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
        % zall = zscored_caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
    % end
    % for some events, the mice have no trials, therefore there are
    % no traces. this line basically skips those neurons (adding 0
    % to the respClass struct), to maintain the same total # of
    % cells throughout (for easy filtering if wanting to check the
    % mean of the respClass.activated = 1 neurons, for example
    % also made it so that if the mouse only has 1 trial, add 0,
    % because these trials have a SEM of 0 and the shuffling method
    % does not work. potentially possible to address this with
    % another method?
    % zall_shuffled{jj} = zall_shuffled;
end






for jj = 1:size(shuffled_data, 1)
    shuffled_trace = shuffled_data(jj, :);
    [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials] = align_and_zscore(BehavData, shuffled_trace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);
    % [caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials] = align_only(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);
    caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1

    % if strcmp(uv.zscore_to, 'window')
        zall_shuff = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
    % elseif strcmp(uv.zscore_to, 'session')
        % zall_shuff = zall_session(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
    % elseif strcmp(uv.zscore_to, 'baseline')
        % zall_shuff = zall_baselined(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
        % zall = zscored_caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
    % end
    % for some events, the mice have no trials, therefore there are
    % no traces. this line basically skips those neurons (adding 0
    % to the respClass struct), to maintain the same total # of
    % cells throughout (for easy filtering if wanting to check the
    % mean of the respClass.activated = 1 neurons, for example
    % also made it so that if the mouse only has 1 trial, add 0,
    % because these trials have a SEM of 0 and the shuffling method
    % does not work. potentially possible to address this with
    % another method?
    zall_shuffled{jj} = zall_shuff;
end

figure; plot(ts1, mean(zall));


%%
% Step 1: Calculate the mean for the real data
zall_mean = mean(zall, 1);  % mean across rows (averaging neuronal activity for each cell)

% Define the time window and create an index for the window [1 3]
time_window = uv.evtSigWin.outcome;  % [1 3]
time_indices = find(ts1 >= time_window(1) & ts1 <= time_window(2));  % indices for the time window

real_mean = mean(zall_mean(time_indices));

% Number of shuffles (1000) and initialize result array
num_shuffles = numel(zall_shuffled);
shuffled_means = zeros(num_shuffles, 1);  % 1000x1 array to store the means

% Step 2: Calculate the mean within the subwindow for each shuffled dataset
for i = 1:num_shuffles
    % Extract the data for the current shuffle and within the time window
    current_data = zall_shuffled{i};  % 30x160 array for each shuffle
    current_subwindow_data = current_data(:, time_indices);  % Extract columns within time window
    
    % Calculate the mean across rows (neurons) for this subwindow
    shuffled_means(i) = mean(current_subwindow_data(:));  % mean across all values in the subwindow
end


% Step 3: Compute the 99th percentile of the shuffled distribution
shuffled_percentile_99 = prctile(shuffled_means, 99);  % 99th percentile for shuffled means

% Step 4: Identify cells whose real mean activity exceeds the 99th percentile
cell_exceeds_99th = zall_mean > shuffled_percentile_99;



%%

% Step 1: Create a histogram of the shuffled means
figure;
histogram(shuffled_means, 30, 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
hold on;

% Step 2: Plot a vertical line for the real mean (zall_mean)
xline(real_mean, 'r', 'LineWidth', 2, 'DisplayName', 'Real Mean');

% Step 3: Plot a vertical line for the 99th percentile (shuffled_percentile_99)
xline(shuffled_percentile_99, 'b', 'LineWidth', 2, 'DisplayName', '99th Percentile');

% Step 4: Add labels and legend
xlabel('Mean Neuronal Activity');
ylabel('Frequency');
title('Comparison of Real Neuronal Activity to Shuffled Distribution');
legend('Shuffled Means', 'Real Mean', '99th Percentile');

% Step 5: Display the plot
hold off;