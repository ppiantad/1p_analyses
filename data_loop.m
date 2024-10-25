iter = 0

%%
% load('BLA-NAcShell_Risk_2024_01_04.mat')

% load('BLA_panneuronal_Risk_2024_01_04.mat')

load('BLA_panneuronal_Risk_2024_03_07_just_CNMFe_and_BehavData.mat')

% load('NAcSh_D2_Cre-OFF_GCAMP_all.mat')

% load('BLA_panneuronal_matched_Pre_RDT_RM_vs_RDT_D1_01042024.mat')

% load('BLA_panneuronal_Risk_matched_RM_D1_vs_Pre_RDT_RM.mat')

% load('BLA_NAcSh_Risk_matched_Pre_RDT_RM_vs_RDT_D1.mat')

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


session_to_analyze = 'Pre_RDT_RM';

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
        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        % for labeling shock + 1 trials, for subsequent analysis
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
        % block_1 = [BehavData.stTime(BehavData.Block == 1) BehavData.collectionTime(BehavData.Block == 1)]; 
        % block_1 = [block_1(1, 1) block_1(end, 2)];
        % block_2 = [BehavData.stTime(BehavData.Block == 2) BehavData.collectionTime(BehavData.Block == 2)]; 
        % block_2 = [block_2(1, 1) block_2(end, 2)];
        % block_3 = [BehavData.stTime(BehavData.Block == 3) BehavData.collectionTime(BehavData.Block == 3)];
        % block_3 = [block_3(1, 1) block_3(end, 2)];
        [BehavData,trials, varargin_identity_class]=TrialFilter_test(BehavData, 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 3);
        % BehavData = BehavData(1:2, :);
        % trials = trials(1:2, :);
        varargin_strings = string(varargin_identity_class);
        varargin_strings = strrep(varargin_strings, '0.3', 'Small');
        varargin_strings = strrep(varargin_strings, '1.2', 'Large');
        filter_args = strjoin(varargin_strings,'_');

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
        % comment out below if you don't want to zscore traces prior to the
        % rest of the analysis
        % ca = zscore(ca, 0, 2);

        if strcmp(ca_data_type, 'S')
            ca = full(ca);

        end
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
            caTraceTrials = NaN(size(eTS,1),numMeasurements); %
            unitTrace = ca(u,:); %get trace
            if isempty(eTS) || size(eTS, 1) == 1
                caTraceTrials(1, 1:size(ts1, 2)) = nan;
                zall(1, 1:size(ts1, 2)) = nan;
                sem_all(neuron_num, size(zall, 2)) = nan;
                zall_mean_all(neuron_num,:) = nan;
            else
                [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);
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
                neuron_mean_mouse_unnormalized{ii, iter}(u,: ) = mean(caTraceTrials, 1);
                unnormalized_by_mouse{ii, iter}(u) = {caTraceTrials(:, 1:size(ts1, 2))};
                neuron_sem_mouse_unnormalized{ii, iter}(u,: ) = nanstd(caTraceTrials,1)/(sqrt(size(caTraceTrials, 1)));
                neuron_mean_mouse{ii, iter}(u,: ) = mean(zall, 1);
                neuron_sem_mouse{ii, iter}(u,: ) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                zall_mean_all(neuron_num,:) = nanmean(zall(:, 1:size(ts1, 2)));

                if size(zall, 1) == 1
                   

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
varargin_list{iter,:} = varargin_identity_class;
behav_tbl_iter{iter, :} = behav_tbl_temp;
epoc_to_align_all{iter,:} = epoc_to_align;
all_filter_args{iter,:} = filter_args;

full_filter_string{iter} = strcat(epoc_to_align_all{iter,:}, '.', all_filter_args{iter,:});

clear behav_tbl_temp

%% FILTER TO GET, ALIGNED SHUFFLED CA DATA

neuron_num = 0;
iter = iter+1;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
    [BehavData,trials,varargin_identity_class]=TrialFilter(BehavData,'REW', 1.2, 'BLOCK', 1);
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
