iter = 0

%% Edit these uservariables with what you want to look at
uv.evtWin = [-2 8]; %what time do you want to look at around each event
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
% uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at

session_to_analyze = 'RDT_D1';
epoc_to_align = 'collectionTime';
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
animalIDs = (fieldnames(final));
neuron_num = 0;

clear neuron_mean neuron_sem neuron_num zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized


%%

neuron_num = 0;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    


    
    BehavData = final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.BehavData;
    [BehavData,trials,varargin]=TrialFilter(BehavData,'REW',1.2, 'BLOCK', 1);
    trials = cell2mat(trials);
    ca = final.(currentanimal).(session_to_analyze).CNMFe_data.C_raw;


    num_samples = size(ca, 2);
    sampling_frequency = (final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.dt)*100;
    time_array = (0:(num_samples-1)) / sampling_frequency;
    eTS = BehavData.(epoc_to_align); %get time stamps

    zb_session = mean(ca,2);
    zsd_session = std(ca,[],2);
    caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


    %calculate time windows for each event
    evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
    numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
    for u = 1:size(ca,1)
        neuron_num = neuron_num+1; 
        % initialize trial matrices
        caTraceTrials = NaN(size(eTS,1),numMeasurements); %
        unitTrace = ca(u,:); %get trace
        %         current_animal = char(upper(mat_strings(1)));
        % current_animal = matlab.lang.makeValidName(current_animal);
        % current_session = char(folder_strings(end));
        % current_session = regexprep(current_session,{' ', '-'}, '_');
        %         current_session = matlab.lang.makeValidName(current_session);

        %             %%
        for t = 1:size(eTS,1)
            % set each trial's temporal boundaries
            timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
            BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
            if min(timeWin) > min(time_array) & max(timeWin) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                % get unit event counts in trials
                % get unit ca traces in trials
                idx = time_array > min(timeWin) & time_array< max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                bl_idx = time_array > min(BL_win) & time_array < max(BL_win);
                %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                caTraceTrials(t,1:sum(idx)) = unitTrace(idx);
                zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
                zb_window(t,:) = mean(caTraceTrials(t,:));
                zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
                zsd_window(t,:) = std(caTraceTrials(t,:));
                tmp = 0;
                for j = 1:size(caTraceTrials,2)
                    tmp = tmp+1;
                    zall(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
                    zall_window(t,tmp) = (caTraceTrials(t,j) - zb_window(t))/zsd_window(t);
                    zall_session(t,tmp) = (caTraceTrials(t,j) - zb_session(u))/zsd_session(u);
                end
                clear j;



            end
        end
        zall_array{neuron_num} = zall;
        zall_mouse{ii, iter+1}(u) = {zall};
        caTraceTrials_mouse{ii, iter+1}(u) = {caTraceTrials};
        zall_mean(neuron_num,:) = mean(zall);
        trials_per_mouse{ii, iter+1} = trials;
        clear zall caTraceTrials zb zsd; 

    end
end

iter = iter+1;

%% SHUFFLE CA DATA

neuron_num = 0;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    


    
    BehavData = final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.BehavData;
    [BehavData,trials,varargin]=TrialFilter(BehavData,'REW',1.2, 'BLOCK', 1);
    trials = cell2mat(trials);
    ca = final.(currentanimal).(session_to_analyze).CNMFe_data.C_raw;

    for g = 1:100 %for each resampling of the data g = 1:uv.resamples
        [num_timepoints, num_signals] = size(ca);
        shuffled_data = zeros(num_timepoints, num_signals); % Preallocate matrix for efficiency

        for i = 1:num_timepoints
            shift_val = randi(num_signals); % Generate a random shift value for each signal
            shuffled_data(i,:) = circshift(ca(i,:), shift_val,2); % Perform the circular shuffle
        end
        % nullDistTrace(g,:) = nanmean(X); %calculate the NaN mean of the shuffled traces
        %calculate the NaN mean of the shuffled event rates
    end

    ca = shuffled_data;
    
    num_samples = size(ca, 2);
    sampling_frequency = (final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.dt)*100;
    time_array = (0:(num_samples-1)) / sampling_frequency;
    eTS = BehavData.(epoc_to_align); %get time stamps

    zb_session = mean(ca,2);
    zsd_session = std(ca,[],2);
    caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


    %calculate time windows for each event
    evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
    numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
    for u = 1:size(ca,1)
        neuron_num = neuron_num+1; 
        % initialize trial matrices
        caTraceTrials = NaN(size(eTS,1),numMeasurements); %
        unitTrace = ca(u,:); %get trace
        %         current_animal = char(upper(mat_strings(1)));
        % current_animal = matlab.lang.makeValidName(current_animal);
        % current_session = char(folder_strings(end));
        % current_session = regexprep(current_session,{' ', '-'}, '_');
        %         current_session = matlab.lang.makeValidName(current_session);

        %             %%
        for t = 1:size(eTS,1)
            % set each trial's temporal boundaries
            timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
            BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
            if min(timeWin) > min(time_array) & max(timeWin) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                % get unit event counts in trials
                % get unit ca traces in trials
                idx = time_array > min(timeWin) & time_array< max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                bl_idx = time_array > min(BL_win) & time_array < max(BL_win);
                %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                caTraceTrials(t,1:sum(idx)) = unitTrace(idx);
                zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
                zb_window(t,:) = mean(caTraceTrials(t,:));
                zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
                zsd_window(t,:) = std(caTraceTrials(t,:));
                tmp = 0;
                for j = 1:size(caTraceTrials,2)
                    tmp = tmp+1;
                    zall(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
                    zall_window(t,tmp) = (caTraceTrials(t,j) - zb_window(t))/zsd_window(t);
                    zall_session(t,tmp) = (caTraceTrials(t,j) - zb_session(u))/zsd_session(u);
                end
                clear j;



            end
        end
        zall_array{neuron_num} = zall;
        zall_mouse{ii, iter+1}(u) = {zall};
        caTraceTrials_mouse{ii, iter+1}(u) = {caTraceTrials};
        zall_mean(neuron_num,:) = mean(zall);
        trials_per_mouse{ii, iter+1} = trials;
        clear zall caTraceTrials zb zsd; 

    end
end

iter = iter+1;