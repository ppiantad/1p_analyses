iter = 0

%%
% load('BLA-NAcShell_Risk_2024_04_16.mat')

% load('BLA_panneuronal_Risk_2024_01_04.mat')

load('BLA_panneuronal_Risk_2024_04_19_just_CNMFe_and_BehavData.mat')

load('BLA_panneuronal_SLEAP_data_2024_02_26_plus_XY_correction.mat')

% load('NAcSh_D2_Cre-OFF_GCAMP_all.mat')

% load('BLA_panneuronal_matched_Pre_RDT_RM_vs_RDT_D1_01042024.mat')

% load('BLA_panneuronal_Risk_matched_RM_D1_vs_Pre_RDT_RM.mat')

% load('BLA_NAcSh_Risk_matched_Pre_RDT_RM_vs_RDT_D1.mat')

%% Edit these uservariables with what you want to look at
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
% uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at

ca_data_type = "C_raw"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate

session_to_analyze = 'RDT_D1';
epoc_to_align = 'choiceTime';
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
animalIDs = (fieldnames(final));
neuron_num = 0;
use_normalized_time = 0;

clear neuron_mean neuron_sem neuron_num zall_mean zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all 


%% FILTER TO GET UN-SHUFFLED DATA
iter = iter+1;
neuron_num = 0;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze) && isfield(final_SLEAP.(currentanimal), session_to_analyze)
        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        % block_1 = [BehavData.stTime(BehavData.Block == 1) BehavData.collectionTime(BehavData.Block == 1)]; 
        % block_1 = [block_1(1, 1) block_1(end, 2)];
        % block_2 = [BehavData.stTime(BehavData.Block == 2) BehavData.collectionTime(BehavData.Block == 2)]; 
        % block_2 = [block_2(1, 1) block_2(end, 2)];
        % block_3 = [BehavData.stTime(BehavData.Block == 3) BehavData.collectionTime(BehavData.Block == 3)];
        % block_3 = [block_3(1, 1) block_3(end, 2)];
        [BehavData,trials,varargin]=TrialFilter(BehavData,'SHK', 1);
        trials = cell2mat(trials);
        behav_tbl_temp{ii,:} = BehavData;
        % % BehavData = BehavData(BehavData.shockIntensity >= 0.08 & BehavData.shockIntensity <= 0.13, :);
        % % trials = trials(BehavData.shockIntensity >= 0.08 & BehavData.shockIntensity <= 0.13, :);
        %
        % %Create a logical index array based on your conditions
        % logical_index = BehavData.stTime - BehavData.TrialPossible >= 10 & BehavData.stTime - BehavData.TrialPossible <= 50;
        %
        % % Use the logical index array to subset BehavData
        % BehavData = BehavData(logical_index,: );
        % trials = trials(logical_index);


        ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);

        if strcmp(ca_data_type, 'S')
            ca = full(ca);

        end

        velocity_data = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.vel_cm_s';
        velocity_data = zscore(velocity_data);
        velocity_data = sgolayfilt(velocity_data, 9, 21);
        zall_motion = velocity_data;
        num_samples = size(ca, 2);
        sampling_frequency = (final.(currentanimal).(session_to_analyze).uv.dt)*100;
        time_array = final.(currentanimal).(session_to_analyze).time;
        if ~strcmp('stTime',BehavData.Properties.VariableNames)
            BehavData.stTime = BehavData.TrialPossible - 5;
        end
        if ~strcmp('collectionTime',BehavData.Properties.VariableNames)
            BehavData.collectionTime = BehavData.choiceTime + 5;
        end



        % time_array = (0:(num_samples-1)) / sampling_frequency;
        eTS = BehavData.(epoc_to_align); %get time stamps

        zb_session = mean(ca,2);
        zsd_session = std(ca,[],2);
        % caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


        %calculate time windows for each event
        evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
        numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate

        % because the first trial possible is ALWAYS 60 seconds after ABET is
        % issued, you can determine what adjustment has been made to this column
        % (adding time to account for calcium recording starting first) by
        % subtracting off 60 from the first element
        if BehavData.TrialPossible(1) < 60
            stTime = BehavData.TrialPossible(1);
        else
            stTime = BehavData.TrialPossible(1)-60;
        end

        
        trim_length = size(velocity_data, 2);

        % time_array = time_array(:, 1:trim_length);

        % it is necessary to add some empty columns because the way it is right
        % now, SLEAP data starts when ABET starts. might make more sense to just
        % trim off the initial samples of the session-long calcium data instead
        % though?
        % num_columns_added_to_SLEAP = round(sampling_frequency * stTime);
        % 
        % buffer_to_SLEAP = zeros(1, num_columns_added_to_SLEAP);
        % 
        % zall_motion = [buffer_to_SLEAP velocity_data];
        % 
        samples_to_remove_from_ca = round(sampling_frequency * stTime);
        
        % get trial x trial velocity data for every mouse
        for t = 1:size(eTS,1)
            time_array_motion = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.idx_time';
            % if strcmp(currentanimal, 'BLA_Insc_09')
            %     time_array_motion = time_array_motion+3.732;
            % end
            timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
            BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
            if min(timeWin) > min(time_array_motion) && max(timeWin) < max(time_array_motion)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                % get unit event counts in trials
                % get unit ca traces in trials
                idx = time_array_motion > min(timeWin) & time_array_motion < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                bl_idx = time_array_motion > min(BL_win) & time_array_motion < max(BL_win);
                %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                motion_TraceTrials(t,1:sum(idx)) = velocity_data(idx);
                % zscored_caTraceTrials(t, 1:sum(idx)) = unitTrace_zscored(idx);
                % zb(t,:) = nanmean(unitTrace(bl_idx)); %baseline mean
                % zb_window(t,:) = nanmean(caTraceTrials(t,:));
                % zsd(t,:) = nanstd(unitTrace(bl_idx)); %baseline std
                % zsd_window(t,:) = nanstd(caTraceTrials(t,:));
                % tmp = 0;
                % for j = 1:size(caTraceTrials,2)
                %     tmp = tmp+1;
                %     zall_baselined(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
                %     zall_window(t,tmp) = (caTraceTrials(t,j) - zb_window(t))/zsd_window(t);
                %     zall_session(t,tmp) = (caTraceTrials(t,j) - zb_session(u))/zsd_session(u);
                % end
                % clear j;

            end
        end
        motion_TraceTrials_mouse{ii} =  motion_TraceTrials;
        clear motion_TraceTrials
        for u = 1:size(ca,1)
            neuron_num = neuron_num+1;
            % initialize trial matrices
            caTraceTrials = NaN(size(eTS,1),numMeasurements); %
            unitTrace = ca(u,:); %get trace
            unitTrace_trimmed = unitTrace(:,samples_to_remove_from_ca:end);
            unitTrace_trimmed = normalize(unitTrace_trimmed); %use normalize instead of zscore in case using spike_prob, because this by default includes NaN values @ the start and the end
            unitTrace_trimmed = sgolayfilt(unitTrace_trimmed, 9, 21);
            if length(unitTrace_trimmed) < length(zall_motion)
                velocity_ca_correlation_neuron = corrcoef(zall_motion(:, 1:size(unitTrace_trimmed, 2)), unitTrace_trimmed, 'rows', 'complete');
                velocity_ca_correlation_neuron_array(neuron_num) = velocity_ca_correlation_neuron(1,2);
                zall_motion = zall_motion(:, 1:size(unitTrace_trimmed, 2));
                ca_and_velocity{neuron_num} = [unitTrace_trimmed; zall_motion];

            else
                velocity_ca_correlation_neuron = corrcoef(zall_motion, unitTrace_trimmed(:, 1:size(zall_motion, 2)), 'rows', 'complete');
                velocity_ca_correlation_neuron_array(neuron_num) = velocity_ca_correlation_neuron(1,2);
                unitTrace_trimmed = unitTrace_trimmed(:, 1:size(zall_motion, 2));
                ca_and_velocity{neuron_num} = [unitTrace_trimmed; zall_motion];

            end
            clear unitTrace_trimmed velocity_ca_correlation_neuron
            if isempty(eTS)
                caTraceTrials(1, 1:size(ts1, 2)) = NaN;
                zall(1, 1:size(ts1, 2)) = NaN;
            else
                [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);
                % [normalized_trial_ca, concatenated_normalized_trial_ca] = normalize_trials_in_time_fn(trial_ca);
                % time_normalized_ca{neuron_num} = concatenated_normalized_trial_ca;
                StartChoiceCollect_times_array{neuron_num} = StartChoiceCollect_times;
                caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                zall = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                for z = 1:size(zall, 1)
                    % Apply Savitzky-Golay filter to each row
                    zall(z, :) = sgolayfilt(zall(z, :), 9, 21);
                end

                zall_array_session{neuron_num} = zall_session(:, 1:size(ts1, 2));
                neuron_mean_unnormalized(neuron_num,:) = nanmean(caTraceTrials,1);
                zall_array{neuron_num} = zall(:, 1:size(ts1, 2));
                zall_mouse{ii, iter}(u) = {zall(:, 1:size(ts1, 2))};
                sem_mouse{ii, iter}(u) = {nanstd(zall,1)/(sqrt(size(zall, 1)))};
                caTraceTrials_mouse{ii, iter}(u) = {caTraceTrials(:, 1:size(ts1, 2))};
                zall_mean_all(neuron_num,:) = nanmean(zall(:, 1:size(ts1, 2)));

                if size(zall, 1) == 1
                    sem_all(neuron_num, size(zall, 2)) = NaN;

                else
                    sem_temp = nanstd(zall,1)/(sqrt(size(zall, 1)));
                    sem_all(neuron_num,:) = sem_temp(:, 1:size(ts1, 2));
                end


                trials_per_mouse{ii, iter} = trials;
                clear zall caTraceTrials zb zsd sem_temp;
            end
        end
    end
end
zall_mean_all_array(iter) = {zall_mean_all};
neuron_mean_all_unnormalized(iter) = {neuron_mean_unnormalized};
sem_all_array(iter) = {sem_all};
varargin_list{iter,:} = varargin;
behav_tbl_iter{iter, :} = behav_tbl_temp;
velocity_ca_correlation_neuron_array_iter(iter) = {velocity_ca_correlation_neuron_array}; 



% Plot histogram
figure;
histogram(velocity_ca_correlation_neuron_array, 'Normalization', 'probability');
title('Distribution of Correlation Coefficients for Session Long Ca and Motion');
xlabel('Correlation Coefficient');
ylabel('Probability');



clear behav_tbl_temp

%% FILTER TO GET, ALIGNED SHUFFLED CA DATA

neuron_num = 0;
iter = iter+1;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
    [BehavData,trials,varargin]=TrialFilter(BehavData,'REW',1.2);
    trials = cell2mat(trials);
    ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);

    % shuffle data for comparison
    [num_cells, num_samples] = size(ca);
    shuffled_data = zeros(num_cells, num_samples); % Preallocate matrix for efficiency
    shift_val = randi(num_samples); % Generate a random shift value for each signal RUAIRI RECOMMENDED KEEPING THE SAME SHIFT VAL, rather than randomizing per neuron. this is because then you keep the overall correlation b/w the neurons, but disrupt the relationship to the event timestamps
    for i = 1:num_cells
        % shift_val = randi(num_signals); % Generate a random shift value for each signal
        shuffled_data(i,:) = circshift(ca(i,:), shift_val,2); % Perform the circular shuffle
    end
    ca = shuffled_data;
    if strcmp(ca_data_type, 'S')
        ca = full(ca);

    end
    num_samples = size(ca, 2);
    sampling_frequency = (final.(currentanimal).(session_to_analyze).uv.dt)*100;
    time_array = (0:(num_samples-1)) / sampling_frequency;
    eTS = BehavData.(epoc_to_align); %get time stamps

    zb_session = mean(ca,2);
    zsd_session = std(ca,[],2);
    % caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


    %calculate time windows for each event
    evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
    numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
    for u = 1:size(ca,1)
        neuron_num = neuron_num+1; 
        % initialize trial matrices
        caTraceTrials = NaN(size(eTS,1),numMeasurements); %
        unitTrace = ca(u,:); %get trace
        [zall_baselined, zall_window, zall_session, caTraceTrials] = align_and_zscore(unitTrace, eTS, uv, time_array, zb_session, zsd_session, u);

        caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
        zall = zall_session(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
        zall_array_session{neuron_num} = zall_session(:, 1:size(ts1, 2));
        zall_array{neuron_num} = zall(:, 1:size(ts1, 2));
        zall_mouse{ii, iter}(u) = {zall(:, 1:size(ts1, 2))};
        caTraceTrials_mouse{ii, iter}(u) = {caTraceTrials(:, 1:size(ts1, 2))};
        zall_mean(neuron_num,:) = mean(zall(:, 1:size(ts1, 2)));
        trials_per_mouse{ii, iter} = trials;
        clear zall caTraceTrials zb zsd;

    end
end
