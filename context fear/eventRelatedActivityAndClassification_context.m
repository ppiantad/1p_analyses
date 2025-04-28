%% Run me first


iter = 0;


load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
load('acton.mat')



%% Load the session you want to examine


% load('BLA-NAcShell_Risk_2023_09_15.mat')


%% Edit these user variables with what you want to look at
uv.evtWin = [-4 4]; %what time do you want to look at around each event [-2 8] [-10 5] [-10 10]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
uv.evtSigWin.outcome = [0 2]; %for SHK or immediate post-choice [0 2] [0 1]
% uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at

ca_data_type = "C_raw"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate


session_to_analyze = 'D1_Afternoon';

% Parameters
session_duration = 12 * 60; % seconds
sampling_rate = 10; % Hz
total_time_points = session_duration * sampling_rate;

shock_start_time = 4 * 60; % First shock in seconds
shock_interval = 60; % Interval between shocks in seconds
shock_duration = 2; % Duration of each shock in seconds
num_shocks = 6;

% Initialize the footshock array
footshock = zeros(1, total_time_points);

% Calculate the indices for each shock
for i = 0:(num_shocks-1)
    shock_start_idx = (shock_start_time + i * shock_interval) * sampling_rate + 1;
    shock_end_idx = shock_start_idx + shock_duration * sampling_rate - 1;
    footshock(shock_start_idx:shock_end_idx) = 1;
end

% Parameters
shock_start_time = 4 * 60; % First shock in seconds
shock_interval = 60; % Interval between shocks in seconds
shock_duration = 2; % Duration of each shock in seconds
num_shocks = 6;

% Initialize variables
shk_on = zeros(1, num_shocks);
shk_off = zeros(1, num_shocks);

% Calculate on and off times for each shock
for i = 0:(num_shocks-1)
    shk_on(i+1) = shock_start_time + i * shock_interval; % Shock start time in seconds
    shk_off(i+1) = shk_on(i+1) + shock_duration; % Shock end time in seconds
end

yoke_data = 0; % 1, set to 1 if you want to be prompted to yoke the number of trials analyzed, set to 0 otherwise

epoc_to_align = 'choiceTime';
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);

neuron_num = 0;
use_normalized_time = 0;

clear neuron_mean neuron_sem neuron_num zall_mean zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all 

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

%% %user selected variables
clear neuron_mean neuron_sem neuron_num trials
uv.chooseFluoresenceOrRate = 1;                                             %set to 1 to classify fluoresence response; set to 2 to classify firing rate responses
uv.sigma = 1.5;  %1.5                                                             %this parameter controls the number of standard deviations that the response must exceed to be classified as a responder. try 1 as a starting value and increase or decrease as necessary.
% uv.evtWin = [-10 10];                                                       %time window around each event in sec relative to event times (use long windows here to see more data)
% % uv.evtSigWin.outcome = [-3 0]; %for trial start
% uv.evtSigWin.outcome = [-4 0]; %for pre-choice   [-4 0]    [-4 1]                              %period within time window that response is classified on (sec relative to event)
% uv.evtSigWin.outcome = [0 2]; %for SHK or immediate post-choice [0 2]
% uv.evtSigWin.outcome = [1 3]; %for REW collection [1 3]
% 
% % uv.evtSigWin.trial_start = [-3 0]; %for trial start
% uv.evtSigWin.prechoice = [-4 0]; %for pre-choice   [-4 0]    [-4 1]                              %period within time window that response is classified on (sec relative to event)
% uv.evtSigWin.postchoice = [0 2]; %for SHK or immediate post-choice [0 2]
% uv.evtSigWin.collect = [1 3]; %for REW collection [1 3]

uv.resamples = 100                                                         %number of resamples to use in shuffle analysis 1000
uv.zscore_to = 'window'; %
% 'window'
% 'baseline'
% 'session'
uv.smoothing.method = 'savgol';
uv.smoothing.params = [9 21];







iter = iter +1;


evtWinSpan = max(uv.evtWin) - min(uv.evtWin);                               %calculate length of each period to examine neural activity in
numMeasurements = evtWinSpan/uv.dt;                                         %calculate the number of measurements in each time window

animalIDs = (fieldnames(final));
neuron_num = 0;
num_trials = 0;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    session_string{iter} = session_to_analyze;

    if isfield(final.(currentanimal), session_to_analyze)

        if exist('full_filter_string', 'var')
            if uv.yoke_data == 1
                
                for i = 1:size(full_filter_string, 2)
                    fprintf('%d. %s\n', i, full_filter_string{1, i});
                end

                % Prompt the user for input
                user_selection = input('Which data would you like to match trials to?: ');

                % Check if the input is valid
                if user_selection >= 1 && user_selection <= size(full_filter_string, 2)
                    selected_data = full_filter_string{1, user_selection};
                    fprintf('You have selected: %s\n', selected_data);
                else
                    disp('Invalid selection. Please run the script again and enter a valid number.');
                end
                size_to_downsample_to = size(trials_per_mouse{ii, user_selection}, 1);
                if size(BehavData, 1) > size_to_downsample_to
                    % Randomly select rows from BehavData
                    rand_indices = randperm(size(BehavData, 1), size_to_downsample_to);
                    BehavData = BehavData(rand_indices, :);
                    trials = trials(rand_indices, :);
                    trials = sortrows(trials);
                    % Sort the filtered BehavData by the Trial column
                    BehavData = sortrows(BehavData, 'Trial');
                else
                    % If the size is not greater, keep BehavData as it is
                    disp('No downsampling needed.');
                end

            else

            end
        end



        
        ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);
        ca_spikes = full(final.(currentanimal).(session_to_analyze).CNMFe_data.S);
        % ca = zscore(ca, 0, 2);
        if strcmp(ca_data_type, 'S')
            ca = full(ca);

        end
        num_samples = size(ca, 2);
        sampling_frequency = (final.(currentanimal).(session_to_analyze).uv.sampling_rate);

        time_array = (0:(num_samples-1)) / sampling_frequency;
        eTS = shk_on'; %get time stamps
        for e = 1:size(eTS,1)                                                      %for each event
            evtWinIdx = ts1 >= uv.evtSigWin.outcome(1,1) &...          %calculate logical index for each event period
                ts1 <= uv.evtSigWin.outcome(1,2);

        end
        clear evtWinSpan e



        zb_session = mean(ca,2);
        zsd_session = std(ca,[],2);
        % caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


        %calculate time windows for each event
        evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
        numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate

        for u = 1:size(ca,1)
            neuron_num = neuron_num+1;

            caTraceTrials = NaN(size(eTS,1),numMeasurements); %
            unitTrace = ca(u,:); %get trace
            caTraceTrials_spikes = NaN(size(eTS,1),numMeasurements); %
            unitTrace_spikes = ca_spikes(u,:); %get trace
            if isempty(eTS)
                caTraceTrials(1, 1:size(ts1, 2)) = NaN;
                zall(1, 1:size(ts1, 2)) = NaN;
            else
                for t = 1:size(eTS,1)
                    % set each trial's temporal boundaries


                    timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
                    BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
                    unitTrace_zscored = zscore(unitTrace);


                    if min(timeWin) > min(time_array) && max(timeWin) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                        % get unit event counts in trials
                        % get unit ca traces in trials
                        idx = time_array >= min(timeWin) & time_array < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                        bl_idx = time_array > min(BL_win) & time_array < max(BL_win);
                        %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                        caTraceTrials(t,:) = unitTrace(idx);
                        zscored_caTraceTrials(t, :) = unitTrace_zscored(idx);
                        zb(t,:) = nanmean(unitTrace(bl_idx)); %baseline mean
                        zb_window(t,:) = nanmean(caTraceTrials(t,:));
                        zsd(t,:) = nanstd(unitTrace(bl_idx)); %baseline std
                        zsd_window(t,:) = nanstd(caTraceTrials(t,:));
                        tmp = 0;
                        for j = 1:size(caTraceTrials,2)
                            tmp = tmp+1;
                            zall_baselined(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
                            zall_window(t,tmp) = (caTraceTrials(t,j) - zb_window(t))/zsd_window(t);
                            zall_session(t,tmp) = (caTraceTrials(t,j) - zb_session(u))/zsd_session(u);
                        end
                        clear j;

                    end

                end
                caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                if strcmp(uv.zscore_to, 'window') 
                    zall = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                elseif strcmp(uv.zscore_to, 'session') 
                    zall = zall_session(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                elseif strcmp(uv.zscore_to, 'baseline') 
                    zall = zall_baselined(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                % zall = zscored_caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                elseif strcmp(uv.zscore_to, 'luthi')
                    zall = zall_luthi(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                end
                % for some events, the mice have no trials, therefore there are
                % no traces. this line basically skips those neurons (adding 0
                % to the respClass struct), to maintain the same total # of
                % cells throughout (for easy filtering if wanting to check the
                % mean of the respClass.activated = 1 neurons, for example
                % also made it so that if the mouse only has 1 trial, add 0,
                % because these trials have a SEM of 0 and the shuffling method
                % does not work. potentially possible to address this with
                % another method?
                if isempty(caTraceTrials) || size(caTraceTrials, 1) == 1
                    respClass_all(neuron_num) = 0;
                    % respClass.(session_to_analyze).(identity_classification_str).(filter_args).activated(neuron_num,1) = 0;
                    % respClass.(session_to_analyze).(identity_classification_str).(filter_args).inhibited(neuron_num,1) = 0;
                    % respClass.(session_to_analyze).(identity_classification_str).(filter_args).neutral(neuron_num,1) = 0;
                    respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args)(u,1) = 0;
                    % respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).activated(qq,1) = 0;
                    % respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).inhibited(qq,1) = 0;
                    % respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).neutral(qq,1) = 0;
                    neuron_mean(neuron_num,:) = NaN;
                    neuron_sem(neuron_num,:) = NaN;
                elseif ~isempty(caTraceTrials) && size(caTraceTrials, 1) > 1

                    % Loop through each row of zall
                    for z = 1:size(zall, 1)
                        if strcmp(uv.smoothing.method, 'savgol')
                            % Apply Savitzky-Golay filter to each row
                            zall(z, :) = sgolayfilt(zall(z, :), uv.smoothing.params(1), uv.smoothing.params(2));
                        end
                    end
                    mouse_cells(iter, neuron_num) = {currentanimal};
                    zall_array(iter, neuron_num) = {zall};
                    % neuron_mean(neuron_num,:) = sgolayfilt((mean(zall,1)), 9, 21);
                    neuron_mean(neuron_num,:) = (mean(zall,1));
                    if size(zall, 1) == 1
                        neuron_sem(neuron_num,:) = zeros(1, size(zall, 2));
                    else
                        neuron_sem(neuron_num,:) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                    end
                    zall_mouse{ii, iter}(u) = {zall};
                    caTraceTrials_spikes_mouse{ii, iter}(u) = {caTraceTrials_spikes};
                    caTraceTrials_spikes_all(iter, neuron_num) = {caTraceTrials_spikes};
                    caTraceTrials_mouse{ii, iter}(u) = {caTraceTrials};
                    neuron_mean_mouse_unnormalized{ii, iter}(u,: ) = mean(caTraceTrials, 1);
                    neuron_sem_mouse_unnormalized{ii, iter}(u,: ) = nanstd(caTraceTrials,1)/(sqrt(size(caTraceTrials, 1)));
                    neuron_mean_mouse{ii, iter}(u,: ) = mean(zall, 1);
                    neuron_sem_mouse{ii, iter}(u,: ) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                    caTraceTrials_unnormalized_array(iter, neuron_num) = {caTraceTrials};
                    caTraceTrials = zall;
                    
                    clear zall zb zsd zall_baselined zall_window zall_session;

                    trialCt = size(caTraceTrials,1);                                    %number of trials for currently analyzed event

                    
                    for g = 1:size(caTraceTrials_spikes, 1)                                         %for each resampling of the data
                        current_trial_data = caTraceTrials_spikes(g, :);
                        [~, lc] = findpeaks(current_trial_data(:, evtWinIdx));
                        
                        current_trial_peaks(g) = size(lc, 2);

                    end
                    trial_peaks_array(iter, neuron_num) = {current_trial_peaks};
                    trial_peaks_mouse_array{ii, iter}(u) = {current_trial_peaks};
                    clear lc pk current_trial_data current_trial_peaks


                    for g = 1:uv.resamples                                              %for each resampling of the data
                        [~,shuffledIDX] = sort(randi...                                 %generate a matrix of random integers from 1 to the number of measurements in each time window (not each number is generated/some are repeated)
                            (numMeasurements,trialCt,numMeasurements),2);               %sort the data index to create a new list of indices
                        for t = 1:trialCt                                               %for each trial
                            shuffledTrace(t,:) = caTraceTrials(t,shuffledIDX(t,:));     %shuffle the calcium trace
                            %         shuffledEvtRate(t,:) = caEvtRateTrials(t,shuffledIDX(t,:)); %shuffle the event rate
                        end
                        nullDistTrace(g,:) = nanmean(shuffledTrace);                    %calculate the NaN mean of the shuffled traces
                        %     nullDistEvtRate(g,:) = nanmean(shuffledEvtRate);                %calculate the NaN mean of the shuffled event rates
                    end
                    clear shuffled* g t trialCt

                    
                    
                    % % 9/13/2024
                    % % alternate way of shuffling, maybe worth trying
                    % [trial_num, sample_num] = size(caTraceTrials);
                    % % [trial_num, sample_num] = size(caTraceTrials);
                    % % shift_val = randi(sample_num)
                    % for g = 1:uv.resamples                                              %for each resampling of the data
                    %     %sort the data index to create a new list of indices
                    %     for t = 1:trial_num
                    % 
                    %         shift_val = randi(sample_num); %for each trial
                    %         shuffledTrace(t,:) = circshift(caTraceTrials(t,:), shift_val,2);     %shuffle the calcium trace
                    %         %         shuffledEvtRate(t,:) = caEvtRateTrials(t,shuffledIDX(t,:)); %shuffle the event rate
                    %     end
                    %     nullDistTrace(g,:) = nanmean(shuffledTrace);                    %calculate the NaN mean of the shuffled traces
                    %     %     nullDistEvtRate(g,:) = nanmean(shuffledEvtRate);                %calculate the NaN mean of the shuffled event rates
                    % end
                    % clear shuffled* g t trialCt





                    %% choose to classify fluoresence or event rates
                    if uv.chooseFluoresenceOrRate == 1                                  %if user selected to classify the fluoresence
                        nullDist = nullDistTrace;                                       %direct transfer
                        empiricalTrialWin = nanmean...
                            (caTraceTrials(:,evtWinIdx));                   %NaN mean of the fluorescent response across trials, within the time window. this gets the within trial mean.
                        empiricalSEM = nansem...
                            (caTraceTrials(:,:));
                        % otherPeriodWin = nanmean...
                        %     (caTraceTrials(:, ~evtWinIdx));
                        % otherEventWin1 = nanmean...
                        %     (caTraceTrials(:, other_evtWinIdx1{iter,:}));
                        % otherEventWin2 = nanmean...
                        %     (caTraceTrials(:, other_evtWinIdx2{iter,:}));
                        empiricalWinAvg = nanmean(empiricalTrialWin);                   %across trial mean
                        empiricalSEMAvg = nanmean(empiricalSEM);                   %across trial mean
                        % otherPeriodWinAvg = nanmean(otherPeriodWin);
                        % otherEventWinAvg1 = nanmean(otherEventWin1);
                        % otherEventWinAvg2 = nanmean(otherEventWin2);
                    elseif uv.chooseFluoresenceOrRate == 2                              %if user selected to classify the event rates
                        nullDist = nullDistEvtRate;                                     %direct transfer
                        empiricalTrialWin = nanmean...
                            (caEvtRateTrials(:,evtWinIdx.(evts{e})),2);                 %within trial average
                        empiricalWinAvg = nanmean(empiricalTrialWin);                   %across trial mean
                    end
                    clear empiricalTrialWin otherPeriodWin empiricalSEM otherEventWin1 ootherEventWin2
                    %%
                    % sdNull = nanstd(nullDist(:,evtWinIdx));                                          %calculate the standard deviation of the null distribution
                    sdNull = nanstd(nullDist(:));                                          %calculate the standard deviation of the null distribution
                    upperSD = nanmean(nullDist(:)) + (uv.sigma*sdNull);                    %calculate upper limit of interval around the mean
                    lowerSD = nanmean(nullDist(:)) - (uv.sigma*sdNull);                    %calculate lower limit of interval around the mean
                    clear sdNull nullDist
                    % Assuming you have already computed empiricalWinAvg, upperSD, lowerSD, otherPeriodWinAvg, and caTraceTrials

                    % Step 1: Calculate the mean of the maxima of individual traces
                    maxima = max(caTraceTrials, [], 2); % Find the maximum value in each row (trial)
                    mean_maxima = mean(maxima);

                    % Step 2: Calculate the standard error of the mean (SEM) of the maxima of individual traces
                    n_trials = size(caTraceTrials, 1); % Number of trials
                    sem_maxima = std(maxima) / sqrt(n_trials);

                    % Step 3: Use Student's t distribution to determine the critical value
                    alpha = 0.05; % Significance level (95% confidence interval)
                    df = n_trials - 1; % Degrees of freedom
                    critical_value = tinv(1 - alpha/2, df); % Two-tailed critical value

                    % Step 4: Multiply the SEM by the critical value to obtain the margin of error
                    margin_of_error = critical_value * sem_maxima;

                    % Step 5: Construct the confidence interval around zero
                    confidence_interval = [-margin_of_error, margin_of_error];
                    clear margin_of_error critical_value df sem_maxima n_trials maxima mean_maxima


                    %%
                    % respClass.(session_to_analyze).(identity_classification_str).(filter_args).activated(neuron_num,1) = empiricalWinAvg > upperSD;     %classify as activated if empirical response exceeds upper limit
                    % respClass.(session_to_analyze).(identity_classification_str).(filter_args).inhibited(neuron_num,1) = empiricalWinAvg < lowerSD;     %classify as inhibited if empirical response exceeds lower limit
                    % respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).activated(qq,1) = empiricalWinAvg > upperSD;     %classify as activated if empirical response exceeds upper limit
                    % respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).inhibited(qq,1) = empiricalWinAvg < lowerSD;     %classify as inhibited if empirical response exceeds lower limit
                    % respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).neutral(qq,1) = respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).activated(qq,1) == 0 & respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).inhibited(qq,1) == 0;
                    % Check if empiricalWinAvg exceeds the 95% confidence interval boundary
                    if empiricalWinAvg > upperSD  %empiricalWinAvg > upperSD & empiricalWinAvg > otherEventWinAvg1 & empiricalWinAvg > otherEventWinAvg2
                        respClass_all(neuron_num) = 1;
                        respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align)(u,1) = 1;
                        respClass_all_array_mouse{ii, iter}(u) = 1;
                    elseif empiricalWinAvg < lowerSD  % empiricalWinAvg < lowerSD & empiricalWinAvg < otherPeriodWinAvg & empiricalWinAvg < otherEventWinAvg2
                        respClass_all(neuron_num) = 2;
                        respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align)(u,1) = 2;
                        respClass_all_array_mouse{ii, iter}(u) = 2;
                    else
                        respClass_all(neuron_num) = 3;
                        respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align)(u,1) = 3;
                        respClass_all_array_mouse{ii, iter}(u) = 3;
                    end
                    clear upperSD lowerSD empiricalWinAvg otherPeriodWinAvg

                    clear upperSD lowerSD empiricalWinAvg otherPeriodWinAvg upper_limit empiricalSEMAvg confidence_interval otherEventWinAvg1
                    %             %% store trial by trial data
                    %             unitXTrials(u).(evts{e}).caEvtCts = caEvtCtTrials;                  %store evoked calcium event counts over all trials
                    %             unitXTrials(u).(evts{e}).caEvtRate = caEvtRateTrials;               %store evoked calcium event rates over all trials
                    %             unitXTrials(u).(evts{e}).caTraces = caTraceTrials;                  %store evoked calcium traces over all trials
                    %             %% store unit averaged data
                    %             unitAVG.(evts{e}).caEvtCts(u,:) = nanmean(caEvtCtTrials);           %store trial averaged event counts
                    %             unitAVG.(evts{e}).caEvtRates(u,:) = nanmean(caEvtRateTrials);       %store trial averaged event rates
                    %             unitAVG.(evts{e}).caTraces(u,:) = nanmean(caTraceTrials);           %store trial averaged calcium traces
                    %             clear caEvtCtTrials caTraceTrials caEvtRateTrials
                end
                clear caTraceTrials;
            end
        end
    end
    clear ca BehavData

end

%%


for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    session_string{iter} = session_to_analyze;

    if isfield(final.(currentanimal), session_to_analyze)
        current_animal_treatment{ii} = final.(currentanimal).experimental_grp;





    end
end

num_cells_experimental = 0;
num_cells_one_context = 0;
num_cells_no_shock = 0;

total_active_experimental = 0;
total_active_one_context = 0;
total_active_no_shock = 0;

total_inhibited_experimental = 0;
total_inhibited_one_context = 0;
total_inhibited_no_shock = 0;

for qq = 1:size(respClass_all_array_mouse, 1)
    if ~isempty(respClass_all_array_mouse{qq})
        current_treatment = current_animal_treatment{qq};
        if strcmp(current_treatment, 'Experimental')
            num_cells_experimental = num_cells_experimental + size(respClass_all_array_mouse{qq}, 2);
            total_active_experimental = total_active_experimental + sum(respClass_all_array_mouse{qq} == 1);
            total_inhibited_experimental = total_inhibited_experimental + sum(respClass_all_array_mouse{qq} == 2);
        
        elseif strcmp(current_treatment, 'One Context')
            num_cells_one_context = num_cells_one_context + size(respClass_all_array_mouse{qq}, 2);
            total_active_one_context = total_active_one_context + sum(respClass_all_array_mouse{qq} == 1);
            total_inhibited_one_context = total_inhibited_one_context + sum(respClass_all_array_mouse{qq} == 2);
            
        elseif strcmp(current_treatment, 'No Shock')
            num_cells_no_shock = num_cells_no_shock + size(respClass_all_array_mouse{qq}, 2);
            total_active_no_shock = total_active_no_shock + sum(respClass_all_array_mouse{qq} == 1);
            total_inhibited_no_shock = total_inhibited_no_shock + sum(respClass_all_array_mouse{qq} == 2);

        end

    end



end


percent_active_experimental = (total_active_experimental/num_cells_experimental)*100
percent_inhibited_experimental = (total_inhibited_experimental/num_cells_experimental)*100
percent_not_active_experimental = 100-[percent_active_experimental + percent_inhibited_experimental] 

percent_active_one_context = (total_active_one_context/num_cells_one_context)*100
percent_inhibited_one_context = (total_inhibited_one_context/num_cells_one_context)*100
percent_not_active_one_context = 100-[percent_active_one_context + percent_inhibited_one_context] 

percent_active_no_shock = (total_active_no_shock/num_cells_no_shock)*100
percent_inhibited_no_shock = (total_inhibited_no_shock/num_cells_no_shock)*100
percent_not_active_no_shock = 100-[percent_active_no_shock + percent_inhibited_no_shock] 

figure; pie([percent_active_experimental percent_inhibited_experimental percent_not_active_experimental])
figure; pie([percent_active_one_context percent_inhibited_one_context percent_not_active_one_context])
figure; pie([percent_active_no_shock percent_inhibited_no_shock percent_not_active_no_shock])




