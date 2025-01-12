%% Run me first


iter = 0;


load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
load('acton.mat')



%% Load the session you want to examine


% load('BLA-NAcShell_Risk_2023_09_15.mat')

% load('BLA-NAcShell_Risk_2024_04_16.mat')

% load('BLA_panneuronal_Risk_2023_07_06.mat')

% load('BLA_panneuronal_Risk_2024_01_04.mat')

load('BLA_panneuronal_Risk_2024_11_21_just_CNMFe_and_BehavData.mat')

% load('NAcSh_D2_Cre-OFF_GCAMP_all.mat')

% load('BLA_panneuronal_matched_RM_D1_vs_Pre_RDT_RM_01042024.mat')

% load('BLA_panneuronal_matched_Pre_RDT_RM_vs_RDT_D1_01042024.mat')

% load('BLA_panneuronal_Risk_matched_RDT_D1_vs_RDT_D2.mat')

% load('BLA_panneuronal_Risk_matched_RDT_D2_vs_RDT_D3.mat')

% load('BLA_panneuronal_Risk_matched_RDT_D1_vs_SHOCK_TEST.mat')

% load('BLA_NAcSh_Risk_matched_RM_D1_vs_Pre_RDT_RM.mat');

% load('BLA_NAcSh_Risk_matched_Pre_RDT_RM_vs_RDT_D1.mat')

% load('BLA_NAcSh_Risk_matched_RDT_D1_vs_RDT_D2.mat')


%% Edit these user variables with what you want to look at
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5]
% uv.BLper = [-10 -5];
uv.BLper = [uv.evtWin(1) uv.evtWin(1)+3];
uv.dt = 0.1; %what is your frame rate
% uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at

uv.ca_data_type = "C_raw"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate

session_to_analyze = 'RDT_D1';
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
elseif strcmp(epoc_to_align, 'PressTime')
    if strcmp(period_of_interest, 'prechoice')
        uv.evtSigWin.outcome = [-4 0]; %for pre-choice   [-4 0]    [-4 1]
    elseif strcmp(period_of_interest, 'postchoice')
        uv.evtSigWin.outcome = [0 2]; %for SHK or immediate post-choice [0 2]
    end
end

ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
animalIDs = (fieldnames(final));
neuron_num = 0;
use_normalized_time = 0;
clear neuron_mean neuron_sem neuron_num zall_mean zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all 


% these are mice that did not complete the entire session - kinda have to
% toss them to do some comparisons during RDT


if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39', 'BLA_Insc_13'};

    for i = 1:length(fieldsToRemove)
        if isfield(final, fieldsToRemove{i})
            final = rmfield(final, fieldsToRemove{i});
        end
    end
elseif strcmp('RM_D1', session_to_analyze)

    fieldsToRemove = {'BLA_Insc_32'}

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



% sub_window_idx = ts1 >= uv.evtSigWin.outcome(1) & ts1 <= uv.evtSigWin.outcome(2);
% 
% identity_classification_win = 'Outcome';
% identity_classification_str = join(string(uv.evtSigWin.outcome), 'to');
% 
% if contains(identity_classification_str, '-')
%     identity_classification_str = strrep(identity_classification_str, '-', 'Minus_');
% end
% 
% if contains(identity_classification_str, '.')
%     identity_classification_str = strrep(identity_classification_str, '.', 'point');
% end
% 
% identity_classification_str = join([identity_classification_win, identity_classification_str],'_');

sub_window_idx = ts1 >= uv.evtSigWin.outcome(1) & ts1 <= uv.evtSigWin.outcome(2);

identity_classification_win = period_of_interest;
identity_classification_str = join(string(uv.evtSigWin.outcome), 'to');

if contains(identity_classification_str, '-')
    identity_classification_str = strrep(identity_classification_str, '-', 'Minus_');
end

if contains(identity_classification_str, '.')
    identity_classification_str = strrep(identity_classification_str, '.', 'point');
end

identity_classification_str = join([identity_classification_win, identity_classification_str],'_');





iter = iter +1;


evtWinSpan = max(uv.evtWin) - min(uv.evtWin);                               %calculate length of each period to examine neural activity in
numMeasurements = evtWinSpan/uv.dt;                                         %calculate the number of measurements in each time window

animalIDs = (fieldnames(final));
neuron_num = 0;
num_trials = 0;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    session_string{iter} = session_to_analyze;
    event_classification_string{iter} = identity_classification_str;
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData; % final.(currentanimal).(session_to_analyze).uv.BehavData;
        % for labeling shock + 1 trials, for subsequent analysis
        if contains(session_to_analyze, 'RDT')
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
            end
        end
        % [BehavData,trials, varargin_identity_class] = TrialFilter_PR(BehavData, 'ALL', 1);
        [BehavData,trials,varargin_identity_class]=TrialFilter_test(BehavData, 'OMITALL', 0, 'BLANK_TOUCH', 0); %'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1    % 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 2, 'BLOCK', 3
        varargin_strings = string(varargin_identity_class);
        varargin_strings = strrep(varargin_strings, '0.3', 'Small');
        varargin_strings = strrep(varargin_strings, '1.2', 'Large');
        filter_args = strjoin(varargin_strings,'_');
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



        % uncomment if you want to test specifically for particular ranges
        % during shock test
        % shk_intensity_range = [0.42 0.50];
        % BehavData = BehavData(round(BehavData.shockIntensity, 2) >= shk_intensity_range(1) & round(BehavData.shockIntensity, 2) <= shk_intensity_range(2), :);
        % trials = trials(round(BehavData.shockIntensity, 2) >=  shk_intensity_range(1) & round(BehavData.shockIntensity, 2) <= shk_intensity_range(2), :);

        % % Create a logi0.10cal index array based on your conditions
        % logical_index = BehavData.stTime - BehavData.TrialPossible >= 10 & BehavData.stTime - BehavData.TrialPossible <= 50;
        % 
        % % Use the logical index array to subset BehavData
        % BehavData = BehavData(logical_index,: );
        % trials = trials(logical_index);

        num_trials = num_trials+sum(numel(trials));
        if ~contains(session_to_analyze, 'PR')
            if ~strcmp('stTime',BehavData.Properties.VariableNames)
                BehavData.stTime = BehavData.TrialPossible - 5;
            end
            if ~strcmp('collectionTime',BehavData.Properties.VariableNames)
                BehavData.collectionTime = BehavData.choiceTime + 5;
                
            end
        end

        if contains(session_to_analyze, 'PR')
            BehavData.choiceTime = BehavData.PressTime;
            if ~strcmp('stTime',BehavData.Properties.VariableNames)
                BehavData.stTime = BehavData.PressTime - 5;
            end
            if ~strcmp('collectionTime',BehavData.Properties.VariableNames)
                BehavData.collectionTime = BehavData.choiceTime + 5;
                
            end
        end

        % if strcmpi('stTime',epoc_to_align)
        %     time2EPOC = BehavData.stTime(:) - BehavData.choiceTime(:);
        % elseif strcmpi('collectionTime', epoc_to_align)
        %     BL_time = BehavData.collectionTime(:)-BehavData.choiceTime(:);
        %     time2EPOC = BehavData.choiceTime(:)-BehavData.collectionTime(:); %time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
        % elseif strcmpi('choiceTime', epoc_to_align)
        %     BL_time = BehavData.stTime(:)-BehavData.choiceTime(:);
        %     %     time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
        %     time2EPOC = BehavData.stTime(:)-BehavData.choiceTime(:);
        % end
        
        behav_tbl_temp{ii,:} = BehavData;

        
        trials = cell2mat(trials);
        ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(uv.ca_data_type);
        % ca = zscore(ca, 0, 2);
        if strcmp(uv.ca_data_type, 'S')
            ca = full(ca);

        end
        num_samples = size(ca, 2);
        sampling_frequency = (final.(currentanimal).(session_to_analyze).uv.dt)*100;
        time_array = final.(currentanimal).(session_to_analyze).time;
        % time_array = (0:(num_samples-1)) / sampling_frequency;
        
        for e = 1:size(BehavData,1)                                                      %for each event
            evtWinIdx = ts1 >= uv.evtSigWin.outcome(1,1) &...          %calculate logical index for each event period
                ts1 <= uv.evtSigWin.outcome(1,2);

        end
        clear evtWinSpan e
        eTS = BehavData.(epoc_to_align); %get time stamps
        eTS_collection = BehavData.collectionTime; %get time stamps
        if strcmp(epoc_to_align, 'choiceTime') & uv.evtSigWin.outcome == [-4 0] %for REW collection
            other_evtWinIdx1{iter,:} = ts1 >= 0 & ts1 <= 2; %ts1 >= 0 & ts1 <= 2;
            other_evtWinIdx2{iter,:} = ts1 >= 1 & ts1 <= 3;
        elseif strcmp(epoc_to_align, 'choiceTime') & uv.evtSigWin.outcome == [0 2] %for SHK or immediate post-choice [0 2]
            other_evtWinIdx1{iter,:} = ts1 >= -1 & ts1 <= 0;
            other_evtWinIdx2{iter,:} = ts1 >= -1 & ts1 <= 0; %might not need this one? ts1 >= 2 & ts1 <= 4
            % other_evtWinIdx2 = other_evtWinIdx1; 
        elseif strcmp(epoc_to_align, 'collectionTime') & uv.evtSigWin.outcome == [1 3] %for SHK or immediate post-choice [1 3]
            other_evtWinIdx1{iter,:} = ts1 >= -4 & ts1 <= 0;
            other_evtWinIdx2{iter,:} = ts1 >= -2 & ts1 <= 0;
        elseif strcmp(epoc_to_align, 'stTime') & uv.evtSigWin.outcome == [0 2] %for SHK or immediate post-choice [1 3]
            other_evtWinIdx1{iter,:} = ts1 >= 0 & ts1 <= 2;
            other_evtWinIdx2{iter,:} = ts1 >= 0 & ts1 <= 2;

        end

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
            if isempty(eTS)
                caTraceTrials(1, 1:size(ts1, 2)) = NaN;
                zall(1, 1:size(ts1, 2)) = NaN;
            else
                [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);
                % [caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials] = align_only(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);
                caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1

                if strcmp(uv.zscore_to, 'window') 
                    zall = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                elseif strcmp(uv.zscore_to, 'session') 
                    zall = zall_session(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                elseif strcmp(uv.zscore_to, 'baseline') 
                    zall = zall_baselined(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                % zall = zscored_caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
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
                    caTraceTrials_mouse{ii, iter}(u) = {caTraceTrials};
                    neuron_mean_mouse_unnormalized{ii, iter}(u,: ) = mean(caTraceTrials, 1);
                    neuron_sem_mouse_unnormalized{ii, iter}(u,: ) = nanstd(caTraceTrials,1)/(sqrt(size(caTraceTrials, 1)));
                    neuron_mean_mouse{ii, iter}(u,: ) = mean(zall, 1);
                    neuron_sem_mouse{ii, iter}(u,: ) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                    caTraceTrials_unnormalized_array(iter, neuron_num) = {caTraceTrials};
                    trials_per_mouse{ii, iter} = trials;
                    caTraceTrials = zall;
                    clear zall zb zsd zall_baselined zall_window zall_session;

                    trialCt = size(caTraceTrials,1);                                    %number of trials for currently analyzed event



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
                        respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args)(u,1) = 1;
                        respClass_all_array_mouse{ii, iter}(u) = 1;
                    elseif empiricalWinAvg < lowerSD  % empiricalWinAvg < lowerSD & empiricalWinAvg < otherPeriodWinAvg & empiricalWinAvg < otherEventWinAvg2
                        respClass_all(neuron_num) = 2;
                        respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args)(u,1) = 2;
                        respClass_all_array_mouse{ii, iter}(u) = 2;
                    else
                        respClass_all(neuron_num) = 3;
                        respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args)(u,1) = 3;
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
varargin_list{iter,:} = varargin_identity_class;
behav_tbl_iter{iter, :} = behav_tbl_temp;
user_vars{iter} = uv; 

clear behav_tbl_temp

epoc_to_align_all{iter,:} = epoc_to_align;
identity_class_string_all{iter,:} = identity_classification_str;
all_filter_args{iter,:} = filter_args;

full_filter_string{iter} = strcat(epoc_to_align_all{iter,:}, '.', identity_class_string_all{iter,:}, '.', all_filter_args{iter,:});

num_trials_per_event(iter) = num_trials;
sum_activated(iter) = sum(respClass_all == 1);
sum_inhibited(iter) = sum(respClass_all == 2);
sum_neutral(iter) = sum(respClass_all == 3);
sum_inhibited_percent(iter) = (sum_inhibited(iter)/neuron_num)*100;
sum_activated_percent(iter) = (sum_activated(iter)/neuron_num)*100;
sum_neutral_percent(iter) = (sum_neutral(iter)/neuron_num)*100;

figure; pie([sum_activated(iter) sum_inhibited(iter) sum_neutral(iter)])
title({"Event: " + strrep(full_filter_string{iter}, '_', '-'), "Num of events (across mice): " + num_trials_per_event(iter), "Num of neurons (across mice): " + neuron_num}, 'FontSize', 9)
labels = {'activated: ' + string(sum_activated(iter)), 'inhibited: ' + string(sum_inhibited(iter)), 'neutral: ' + string(sum_neutral(iter))};
legend(labels)
neuron_mean_array(iter) = {neuron_mean};
neuron_sem_array(iter) = {neuron_sem};
respClass_all_array(:,iter) = {respClass_all}';
respClass_all_array_mouse_array(:,iter) = {respClass_all_array_mouse}';
clear respClass_all 


%% plot activated neurons

figure; plot(ts1, neuron_mean(respClass_all_array{:,iter} == 1,:))
hold on;
figure;
% figure; plot(ts1, mean(neuron_mean(respClass_all(iter,:) == 1,:)))
shadedErrorBar(ts1, nanmean(neuron_mean(respClass_all_array{:,iter} == 1,:)), nanmean(neuron_sem(respClass_all_array{:,iter} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
hold off; 
% hold on;
% figure; plot(ts1, mean(neuron_mean(respClass_all(iter,:) == 1,:)))
% errorplot3(mean(neuron_mean(respClass_all(iter,:) == 1,:))-mean(neuron_sem(respClass_all(iter,:) == 1,:)),mean(neuron_mean(respClass_all(iter,:) == 1,:))+mean(neuron_sem(respClass_all(iter,:) == 1,:)),[-10 10],batlowW(iter,:),.15);
% hold off;
% plot inhibited neurons
figure; plot(ts1, neuron_mean(respClass_all_array{:,iter} == 2,:))
figure; shadedErrorBar(ts1, nanmean(neuron_mean(respClass_all_array{:,iter} == 2,:)), nanmean(neuron_sem(respClass_all_array{:,iter} == 2,:)), 'lineProps', {'color', batlowW(iter,:)});
hold off; 
% plot neutral neurons
figure; plot(ts1, neuron_mean(respClass_all_array{:,iter} == 3,:))
figure; shadedErrorBar(ts1, nanmean(neuron_mean(respClass_all_array{:,iter} == 3,:)), nanmean(neuron_sem(respClass_all_array{:,iter} == 3,:)), 'lineProps', {'color', batlowW(iter,:)});


%% Use this code to plot heatmaps for each individual cell, across trials for all levels of iter
% **most useful for plotting matched cells within the same experiment, e.g., pan-neuronal matched Pre-RDT RM vs. RDT D1**

for ii = 1:size(zall_array, 2)
    figure;
    % Initialize variables to store global max and min for heatmap and line graph
    globalMaxHeatmap = -inf;
    globalMinHeatmap = inf;
    globalMaxYLine = -inf;
    globalMinYLine = inf;
    for qq = 1:iter
        
        % Create subplot with 1 row and 2 columns
        num_columns_plot = iter;
        subplot(2, num_columns_plot, qq);

        % Generate the heatmap
        imagesc(ts1, 1, zall_array{qq, ii});
        colormap gray;
        hold on;
        title({"Cell from " + strrep(mouse_cells{qq, ii}, '_', '-'), "Classified as " + respClass_all_array{qq}(ii), "(Overall cell number = " + (ii) + ")"}, 'FontSize', 9)
        
        % Update global max and min for heatmap
        localMaxHeatmap = max(zall_array{qq, ii}(:));
        localMinHeatmap = min(zall_array{qq, ii}(:));
        globalMaxHeatmap = max(globalMaxHeatmap, localMaxHeatmap);
        globalMinHeatmap = min(globalMinHeatmap, localMinHeatmap);
  
        % Find the row index in animalIDs that matches mouse_cells{qq, ii}
        isMatch = find(ismember(animalIDs,       mouse_cells{qq, ii}));

        if ~isempty(isMatch)  
            % Access beha     v_tbl_iter using the first index (assuming there's only one match)
            behav_tbl = behav_tbl_iter{qq, 1}{isMatch, 1};
            time2Collect = behav_tbl.collectionTime(:) - behav_tbl.choiceTime(:);
            trialStartTime = behav_tbl.stTime(:) - behav_tbl.choiceTime(:);
            median_trialStartTime = median(trialStartTime)
            median_time2Collect = median(time2Collect)
            xline(median_trialStartTime)
            xline(median_time2Collect)
            [numTrials, ~] = size(behav_tbl.collectionTime(:));
            Tris = [1:numTrials]';
            % scatter(time2Collect, Tris               , 'Marker', 'p', 'MarkerFaceColor', 'w', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.7)
            % scatter(trialStartTime, Tris, 'Marker', 's', 'MarkerFaceColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.7)
            xline(0, 'Color', 'w')
        end     
        colorbar;
        
        hold off;
        clear time2Collect trialStartTime numTrials Tris behav_tbl

        % Create subplot for the mean and raw data
        subplot(2, num_columns_plot, iter + qq);

        % Plot the mean as a thick black line
        meanData = mean(zall_array{qq, ii});
        plot(ts1, meanData, 'r', 'LineWidth', 2);
        hold on;

        % Plot the raw data in grey with transparency
        for trial = 1:size(zall_array{qq, ii}, 1)
            plot(ts1, zall_array{qq, ii}(trial, :), 'Color', [0.1, 0.1, 0.1, 0.1]);
            hold on;
        end

        title("Mean and Raw Data", 'FontSize', 9)
        % Update global max and min for line graph
        localMaxYLine = max(zall_array{qq, ii}, [], "all");
        localMinYLine = min(zall_array{qq, ii}, [], "all");
        globalMaxYLine = max(globalMaxYLine, localMaxYLine);
        globalMinYLine = min(globalMinYLine, localMinYLine);
        
        % hold off;        
    end
    % Set the same colorbar scale for all heatmap subplots
    for qq = 1:iter
        subplot(2, num_columns_plot, qq);
        %set limits for the heatmap, use a flat value to exaggerate diffs,
        %or set based on the max and min of the range
        clim([-1, 1]); %clim([globalMinHeatmap-1, globalMaxHeatmap-1]);
        colorbar;
    end

    % Set the same Y-axis scale for all line graph subplots
    for qq = 1:iter
        subplot(2, num_columns_plot, iter + qq);
        ylim([-4 4]);
        xlim(uv.evtWin);
    end

    pause
    hold
    close
end





%%
excited_to_excited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1;
inhibited_to_inhibited = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 2;
excited_to_excited_sum = sum(excited_to_excited);
indices_excited_to_excited = find(excited_to_excited == 1);
inhibited_to_inhibited_sum = sum(inhibited_to_inhibited);

excited_to_inhibited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 2;
excited_to_inhibited_sum = sum(excited_to_inhibited);
excited_to_neutral = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3;
excited_to_neutral_sum = sum(excited_to_neutral);

inhibited_to_excited = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 1;
inhibited_to_excited_sum = sum(inhibited_to_excited);
inhibited_to_neutral = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 3;
inhibited_to_neutral_sum = sum(inhibited_to_neutral);

% above this line are good for single session stacked plots (e.g., how do
% neurons on session 1 that are activated (or inhibited) change to session
% 2?) 

neutral_to_excited = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1;
neutral_to_excited_sum = sum(neutral_to_excited);

neutral_to_inhibited = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 2;
neutral_to_inhibited_sum = sum(neutral_to_inhibited);

% need to add excited_to_excited twice because they are present in BOTH
% sessions
total_activated_possible = excited_to_excited_sum+ excited_to_excited_sum + excited_to_inhibited_sum+ excited_to_neutral_sum+ neutral_to_excited_sum + inhibited_to_excited_sum;
% need to add inhibited_to_inhibited twice because they are present in BOTH
% sessions
total_inhibited_possible = inhibited_to_inhibited_sum + inhibited_to_inhibited_sum + inhibited_to_excited_sum + inhibited_to_neutral_sum + neutral_to_inhibited_sum + excited_to_inhibited_sum;






neutral_to_neutral = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3;
neutral_to_neutral_sum = sum(neutral_to_neutral);

exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1);
exclusive_activated_session_2 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2);



exclusive_activated_session_1_indices = find(exclusive_activated_session_1(1,:) == 1);


exclusive_inhibited_session_1 = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 3;
exclusive_inhibited_session_1_sum = sum(exclusive_inhibited_session_1);
exclusive_inhibited_session_2 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 2;
exclusive_inhibited_session_2_sum = sum(exclusive_inhibited_session_2);
% exclusive_activated_sum = sum(exclusive_activated);
% exclusive_inhibited_sum = sum(exclusive_inhibited);

data = [exclusive_activated_session_1_sum exclusive_activated_session_2_sum exclusive_inhibited_session_1_sum exclusive_inhibited_session_2_sum ...
    excited_to_excited_sum inhibited_to_inhibited_sum excited_to_inhibited_sum inhibited_to_excited_sum, neutral_to_neutral_sum]


pie_labels = {'large reward consumption activated ONLY', 'small reward consumption activated ONLY', 'large reward consumption inhibited ONLY', ...
    'small reward consumption inhibited ONLY', 'large AND small reward consumption excited', 'large AND small reward consumption inhibited', ...
    'large reward consumption excited AND small reward consumption inhibited', 'large reward consumption inhibited AND small reward consumption excited', 'neutral'};

figure; pie(data);
legend(pie_labels)

% exclusive_modulated = exclusive_activated_sum + exclusive_inhibited_sum;


%%
% Example 2: Nested pie chart with custom colors for each wedge

% Initialize data points
inner_pie = [exclusive_activated_session_1_sum/neuron_num,...
            exclusive_inhibited_session_1_sum/neuron_num,...
            exclusive_activated_session_2_sum/neuron_num,...
            exclusive_inhibited_session_2_sum/neuron_num,...
            excited_to_excited_sum/neuron_num,...
            inhibited_to_inhibited_sum/neuron_num,...
            excited_to_inhibited_sum/neuron_num,...
            inhibited_to_excited_sum/neuron_num,...
            neutral_to_neutral_sum/neuron_num];
outer_pie = [(exclusive_activated_session_1_sum+exclusive_inhibited_session_1_sum)/neuron_num, ...
            (exclusive_activated_session_2_sum+exclusive_inhibited_session_2_sum)/neuron_num, ...
            (excited_to_excited_sum+inhibited_to_inhibited_sum+excited_to_inhibited_sum+inhibited_to_excited_sum)/neuron_num,...
            (neutral_to_neutral_sum/neuron_num)];
C = {...
    outer_pie,... % Inner to outer layer
    inner_pie};

% Custom colors
inner_colors = [...
    0 0.4470 0.7410;...
    0 0.4470 0.7410;...
    1 0 0;...
    1 0 0;...
    0 1 0;...
    0 1 0;...
    0 1 0;...
    0 1 0;...
    0.8 0.8 0.8];
outer_colors = [...
    0 0.4470 0.7410;...
    1 0 0;...
    0 1 0;...
    0.8 0.8 0.8];
wedge_colors = {...
    outer_colors,...
    inner_colors};

% Spider plot
nested_pie(C,...
    'WedgeColors', wedge_colors,...
    'LegendOrder', [1, 2],...
    'RhoLower', 0.5);

% Title
title('Nested Pie Chart');


hold off; 
figure; donutchart(outer_pie, 'InnerRadius', 0.8)
figure; pie(inner_pie)

%%
%CREATE A STACKED BAR PLOT TO SHOW THE PROPORTION OF MODULATED (CROSS
%SESSION) NEURONS THAT CO OR EXCLUSIVELY MODULATED

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the total height of the stacked bar
total_height = sum(total_activated_possible)-excited_to_excited_sum;
% total_height = sum(total_inhibited_possible)-inhibited_to_inhibited_sum;

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Create a matrix for the bar plot data
stacked_plot_data = [excited_to_excited_sum, excited_to_inhibited_sum, excited_to_neutral_sum, neutral_to_excited_sum, inhibited_to_excited_sum];
% stacked_plot_data = [inhibited_to_inhibited_sum, inhibited_to_excited_sum, inhibited_to_neutral_sum, neutral_to_inhibited_sum, excited_to_inhibited_sum];


% Calculate the remaining portion
% remaining_portion = total_height - sum(stacked_plot_data);
figure;
% Create a vector for the x-axis values
x = 1;

% Create a stacked bar plot
bar(x, [stacked_plot_data], 'stacked');

% Set the colors for each section
colormap([0.3, 0.6, 0.9; 0.5, 0.5, 0.5; 0.2, 0.2, 0.2]); % Blue, Gray, and Dark Gray

% % Set the x-axis tick labels
% xticklabels({'Activation Types'});

% Set the y-axis limits to match the total height
ylim([0, total_height]);

% Add labels and a legend
% ylabel('Total Activation');
legend('Co-Excited', 'Excited to Inhibited', 'Excited to Neutral', 'Neutral to Excited', 'Inhibited to Excited');
% legend('Co-Inhibited', 'Inhibited to Excited', 'Inhibited to Neutral', 'Neutral to Inhibited', 'Excited to Inhibited');

%%
%CREATE A STACKED BAR PLOT TO SHOW THE PROPORTION OF MODULATED neurons
%within a single session

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the total height of the stacked bar
total_height = sum(sum_activated(1));
% total_height = sum(sum_inhibited);

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Create a matrix for the bar plot data
stacked_plot_data = [excited_to_excited_sum, excited_to_inhibited_sum, excited_to_neutral_sum];
% stacked_plot_data = [co_inhibited_sum, exclusive_inhibited_sum, inhibited_to_excited_sum, inhibited_to_neutral_sum];


% Calculate the remaining portion
% remaining_portion = total_height - sum(stacked_plot_data);
figure;
% Create a vector for the x-axis values
x = 1;

% Create a stacked bar plot
bar(x, [stacked_plot_data], 'stacked');

% Set the colors for each section
colormap([0.3, 0.6, 0.9; 0.5, 0.5, 0.5; 0.2, 0.2, 0.2]); % Blue, Gray, and Dark Gray

% % Set the x-axis tick labels
% xticklabels({'Activation Types'});

% Set the y-axis limits to match the total height
ylim([0, total_height]);

legend_string = strcat('Co-Excited =', num2str(excited_to_excited_sum));

% Add labels and a legend
% ylabel('Total Activation');
legend(strcat('Co-Excited = ', num2str(excited_to_excited_sum)), strcat('Excited to Inhibited = ', num2str(excited_to_inhibited_sum)), strcat('Excited to Neutral = ', num2str(excited_to_neutral_sum)));
% legend('Co-Inhibited', 'Exclusive Inhibited', 'Inhib to Excite', 'Inhib to Neutral');

%%
%CREATE A STACKED BAR PLOT TO SHOW THE PROPORTION OF MODULATED (CROSS
%SESSION) NEURONS THAT CO OR EXCLUSIVELY MODULATED

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the total height (total number of neurons)
% total_height = sum(sum_activated);
total_height = sum(sum_inhibited);


% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the percentages for each section
% percentage_co_activated = (co_activated_sum / total_height) * 100;
% percentage_exclusive_activated = (exclusive_activated_sum / total_height) * 100;
% percentage_remaining = 100 - (percentage_co_activated + percentage_exclusive_activated);

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the percentages for each section
percentage_co_inhibited = (co_inhibited_sum / total_height) * 100;
percentage_exclusive_inhibited = (exclusive_inhibited_sum / total_height) * 100;
percentage_remaining = 100 - (percentage_co_inhibited + percentage_exclusive_inhibited);

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Create a matrix for the bar plot data
% stacked_plot_data = [percentage_co_activated, percentage_exclusive_activated, percentage_remaining];
stacked_plot_data = [percentage_co_inhibited, percentage_exclusive_inhibited, percentage_remaining];

figure;
% Create a vector for the x-axis values
x = 1;

% Create a stacked bar plot
bar(x, stacked_plot_data, 'stacked');

% Set the colors for each section
colormap([0.3, 0.6, 0.9; 0.5, 0.5, 0.5; 0.2, 0.2, 0.2]); % Blue, Gray, and Dark Gray

% Set the x-axis tick labels
xticklabels({'Activation Types'});

% Set the y-axis limits to be percentages (0 to 100)
ylim([0, 100]);


legend('Co', 'Exclusive', 'Remaining');



%%
Block_1_activated = respClass_all_array{1,1} == 1;
Block_1_activated_sum = sum(Block_1_activated);
Block_2_activated = respClass_all_array{1,2} == 1;
Block_2_activated_sum = sum(Block_2_activated);
Block_2_consistent_activated = Block_1_activated & respClass_all_array{1,2} == 1;
disp(sum(Block_2_consistent_activated));
Block_3_consistent_activated = Block_1_activated & respClass_all_array{1,3} == 1;
disp(sum(Block_3_consistent_activated));

Block_1_responsive = respClass_all_array{1,1} == 1 | respClass_all_array{1,2} == 2;
disp(sum(Block_1_responsive));
Block_2_consistent_responsive = (respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1) | (respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 2);
disp(sum(Block_2_consistent_responsive));
Block_3_consistent_responsive = (respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1) & respClass_all_array{1,3} == 1 | (respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 2) & respClass_all_array{1,3} == 2;
disp(sum(Block_3_consistent_responsive));


Block_2_large_to_small = (respClass_all_array{1,2} ~= Block_2_consistent_responsive) & (respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 1) | (respClass_all_array{1,2} == 2 & respClass_all_array{1,3} == 2);
disp(sum(Block_2_large_to_small));
Block_3_large_to_small = (respClass_all_array{1,2} ~= Block_2_consistent_responsive) & (respClass_all_array{1,1} == 1 & respClass_all_array{1,4} == 1) | (respClass_all_array{1,2} == 2 & respClass_all_array{1,4} == 2);
disp(sum(Block_3_large_to_small));

% this is the start of checking if neurons are MORE active than during
% other events, i.e. if you wanted to check if REW activated neurons are
% significantly differentially activated by SHK. preliminary poking around
% seems to suggest that few large reward active neurons have their activity
% increase in response to SHK
co_activated_indices = find(excited_to_excited(1,:) == 1);
co_activated_indices_sum = numel(co_activated_indices);
for qq = 1:size(co_activated_indices, 2)
    [h(qq),p(qq),ci{qq},stats{qq}] = ttest(neuron_mean_array{1, 1}(co_activated_indices(qq),sub_window_idx),neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx));
    mean_diff(qq) = mean(neuron_mean_array{1, 1}(co_activated_indices(qq),sub_window_idx) - mean(neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx)));
end

sig_increase_shk_from_large = co_activated_indices(p < 0.05 & mean_diff < 0);
sig_increase_shk_from_large_sum = numel(sig_increase_shk_from_large);
sig_increase_shk_from_large_ind = zeros(1, size(respClass_all_array{1,2}, 2));

sig_increase_shk_from_large_ind(:, sig_increase_shk_from_large) = 1;

shk_activated = respClass_all_array{1,2} == exclusive_activated_session_2 |  respClass_all_array{1,2} == sig_increase_shk_from_large_ind;
shk_activated_sum = sum(shk_activated)

figure; plot(ts1, mean(neuron_mean_array{1, 2}(shk_activated,:))); hold on; plot(ts1,  mean(neuron_mean_array{1, 1}(Block_1_activated,:)));



figure;
hold on; 
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(shk_activated,:)), mean(neuron_sem_array{1, 2}(shk_activated,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 1}(Block_1_activated,:)), mean(neuron_sem_array{1, 1}(Block_1_activated,:)), 'lineProps', {'color', batlowW(iter,:)});


% total_modulated = [sum_activated + sum_inhibited];
% total_co_modulated = [co_activated_sum + co_inhibited_sum];


total_modulated = [(Block_1_activated_sum/neuron_num)*100 (shk_activated_sum/neuron_num)*100];

A = total_modulated;
I = (sig_increase_shk_from_large_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end




total_modulated = [(Block_1_activated_sum/neuron_num)*100 (Block_2_activated_sum/neuron_num)*100];

A = total_modulated;
I = (co_activated_indices_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end


figure;
hold on; 
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(exclusive_activated_session_2,:)), mean(neuron_sem_array{1, 2}(exclusive_activated_session_2,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 1}(exclusive_activated_session_1,:)), mean(neuron_sem_array{1, 1}(exclusive_activated_session_1,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 1}(excited_to_excited,:)), mean(neuron_sem_array{1, 1}(excited_to_excited,:)), 'lineProps', {'color', batlowW(iter,:)});



%%
total_modulated = [sum_activated + sum_inhibited];
total_co_modulated = [co_activated_sum + co_inhibited_sum];

A = (total_modulated/neuron_num)*100;
I = (total_co_modulated/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
  %Now label each zone 
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end




A = (sum_activated/neuron_num)*100;
I = (excited_to_excited_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end


A = (sum_activated(iter)/sum_activated(iter-1))*100;
I = (excited_to_excited_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end


A = (sum_activated/neuron_num)*100;
I = (exclusive_activated_sum/co_activated_sum)*100
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end



A = (sum_inhibited/neuron_num)*100;
I = (co_inhibited_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end



large_pre_and_shock_post = (respClass_all_array{1,1} == 1 & shk_activated == 1);
large_pre_and_shock_post_sum = sum(large_pre_and_shock_post)
large_pre_to_small_pre = (respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 1);
large_pre_to_small_pre_sum = sum(large_pre_to_small_pre)

overlap = large_pre_and_shock_post == 1 & large_pre_to_small_pre == 1; 
overlap_sum = sum(overlap)

large_responsive_sum = sum(respClass_all_array{1,1} == 1);


total_modulated = [(large_pre_and_shock_post_sum/large_responsive_sum)*100 (large_pre_to_small_pre_sum/large_responsive_sum)*100];
overlap_sum = (sum(overlap)/large_responsive_sum)*100;

A = total_modulated;
I = overlap_sum;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end



%% calculates RESPONSIVE neuron percents for a series of filtering steps. dependent on what filters are used, so beware! 

large_pre_and_shock_post = (respClass_all_array{1,1} == 1 | rrespClass_all_array{1,1} == 2) & (shk_activated == 1);
large_pre_and_shock_post_sum = sum(large_pre_and_shock_post)
large_pre_to_small_pre = (respClass_all_array{1,1} == 1 | respClass_all_array{1,1} == 2) & (respClass_all_array{1,3} == 1 | respClass_all_array{1,3} == 2);
large_pre_to_small_pre_sum = sum(large_pre_to_small_pre)
large_pre_still_large_pre = (respClass_all_array{1,1} == 1) & (respClass_all_array{1,4} == 1);
large_pre_still_large_sum = sum(large_pre_still_large_pre)


overlap_shock_small = large_pre_and_shock_post == 1 & large_pre_to_small_pre == 1; 
overlap_shock_small_sum = sum(overlap_shock_small)

overlap_shock_large = large_pre_and_shock_post == 1 & large_pre_still_large_pre == 1; 
overlap_shock_large_sum = sum(overlap_shock_large)

overlap_small_large = large_pre_to_small_pre == 1 & large_pre_still_large_pre == 1; 
overlap_small_large_sum = sum(overlap_small_large)

overlap_all = large_pre_to_small_pre == 1 & large_pre_still_large_pre == 1 & large_pre_and_shock_post == 1; 
overlap_all_sum = sum(overlap_all)

large_responsive_sum = sum(respClass_all_array{1,1} == 1 | respClass_all_array{1,1} == 2)


total_modulated = [(large_pre_and_shock_post_sum/large_responsive_sum)*100 (large_pre_to_small_pre_sum/large_responsive_sum)*100 (large_pre_still_large_sum/large_responsive_sum)*100];
overlap_sum = [(overlap_shock_small_sum/large_responsive_sum)*100 (overlap_shock_large_sum/large_responsive_sum)*100 (overlap_small_large_sum/large_responsive_sum)*100 (overlap_all_sum/large_responsive_sum)*100]

A = total_modulated;
I = overlap_sum;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
[H, S] = venn(A,I);
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end



%%
[median_choice_time_block_1, median_choice_time_block_2, median_choice_time_block_3, median_collect_time_block_1, median_collect_time_block_2, median_collect_time_block_3, median_collect_time_from_choice] = get_median_choice_and_collect_fn(behav_tbl_iter);

figure;
hold on; 
shadedErrorBar(ts1, mean(neuron_mean(respClass_all_array{:,iter} == 1,:)), mean(neuron_sem(respClass_all_array{:,iter} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean(respClass_all_array{:,iter} == 2,:)), mean(neuron_sem(respClass_all_array{:,iter} == 2,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean(respClass_all_array{:,iter} == 3,:)), mean(neuron_sem(respClass_all_array{:,iter} == 3,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
%%

pun_responsive_B1_large_to_B2_null = (respClass_all_array{1,1} == 1 & respClass_all_array{1,1} == shk_activated) & respClass_all_array{1,3} == 2;
sum(pun_responsive_B1_large_to_B2_null)

%%
consistent_increase = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 1 & respClass_all_array{1, 3} == 1;
consumption_encoding_lost = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} ~= 1 & respClass_all_array{1, 3} ~= 1;
consumption_encoding_gained = respClass_all_array{1, 1} ~= 1 & respClass_all_array{1, 2} == 1 & respClass_all_array{1, 3} == 1;
