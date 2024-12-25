
% 
% num_samples = size(final.b78764.D1_Afternoon.CNMFe_data.C_raw, 2);
% 
% time_step = final.b78764.D1_Afternoon.uv.sampling_rate/100;
% 
% time_series = (1:num_samples) * time_step;
% 
% figure; plot(time_series, final.b78764.D1_Afternoon.CNMFe_data.C_raw(1, :))
%%
iter = 0

%%
% load('BLA-NAcShell_Risk_2024_01_04.mat')

%% Edit these uservariables with what you want to look at
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5] [-10 10]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
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
%% FILTER TO GET UN-SHUFFLED DATA
iter = iter+1;
neuron_num = 0;
animalIDs = (fieldnames(final));


for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        if exist('full_filter_string', 'var')
            if yoke_data == 1

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
        % comment out below if you don't want to zscore traces prior to the
        % rest of the analysis
        % ca = zscore(ca, 0, 2);

        if strcmp(ca_data_type, 'S')
            ca = full(ca);

        end
        num_samples = size(ca, 2);
        sampling_frequency = (final.(currentanimal).(session_to_analyze).uv.sampling_rate)/100;
        time_series = (0:num_samples) * sampling_frequency;


        time_array = (0:(num_samples-1)) / sampling_frequency;
        eTS = shk_on'; %get time stamps
        zb_session = [];
        zsd_session = [];
        zb_session = mean(ca,2);
        zsd_session = std(ca,[],2);
        % caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


        %calculate time windows for each event
        evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
        numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
        for u = 1:size(ca,1)
            neuron_num = neuron_num+1;
            % initialize trial matrices
            % caTraceTrials = NaN(size(eTS,1),numMeasurements); %
            unitTrace = ca(u,:); %get trace
            if isempty(eTS) || size(eTS, 1) == 1
                caTraceTrials(1, 1:size(ts1, 2)) = nan;
                zall(1, 1:size(ts1, 2)) = nan;
                sem_all(neuron_num, size(zall, 2)) = nan;
                zall_mean_all(neuron_num,:) = nan;
            else
                
                for t = 1:size(eTS,1)
                    % set each trial's temporal boundaries


                    timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
                    BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
                    unitTrace_zscored = zscore(unitTrace);


                    if min(timeWin) > min(time_series) && max(timeWin) < max(time_series)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                        % get unit event counts in trials
                        % get unit ca traces in trials
                        idx = time_series >= min(timeWin) & time_series < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                        bl_idx = time_series > min(BL_win) & time_series < max(BL_win);
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
                zall = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                for z = 1:size(zall, 1)
                    % Apply Savitzky-Golay filter to each row
                    zall(z, :) = sgolayfilt(zall(z, :), 9, 21);
                end

                zall_array_session{neuron_num} = zall_session(:, 1:size(ts1, 2));
                neuron_mean_unnormalized(neuron_num,:) = nanmean(caTraceTrials,1);
                zall_array{neuron_num} = zall;
                zall_mouse{ii, iter}(u) = {zall};
                sem_mouse{ii, iter}(u) = {nanstd(zall,1)/(sqrt(size(zall, 1)))};
                caTraceTrials_mouse{ii, iter}(u) = {caTraceTrials};
                neuron_mean_mouse_unnormalized{ii, iter}(u,: ) = mean(caTraceTrials, 1);
                unnormalized_by_mouse{ii, iter}(u) = {caTraceTrials(:, 1:size(ts1, 2))};
                neuron_sem_mouse_unnormalized{ii, iter}(u,: ) = nanstd(caTraceTrials,1)/(sqrt(size(caTraceTrials, 1)));
                neuron_mean_mouse{ii, iter}(u,: ) = mean(zall, 1);
                neuron_sem_mouse{ii, iter}(u,: ) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                zall_mean_all(neuron_num,:) = nanmean(zall);
                mouse_cells(iter, neuron_num) = {currentanimal};
                if size(zall, 1) == 1
                   

                else
                    sem_temp = nanstd(zall,1)/(sqrt(size(zall, 1)));
                    sem_all(neuron_num,:) = sem_temp;
                end



                clear zall caTraceTrials zb zsd sem_temp;
            end
        end
    end
end
zall_mean_all_array(iter) = {zall_mean_all};
neuron_mean_all_unnormalized(iter) = {neuron_mean_unnormalized};
sem_all_array(iter) = {sem_all};

