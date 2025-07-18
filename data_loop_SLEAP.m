iter = 0

%%


%% Edit these uservariables with what you want to look at
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5] [-10 10]
uv.BLper = [-10 -5];
uv.dt = 1/10; %what is your frame rate
% uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at

animalIDs = (fieldnames(final_SLEAP));

session_to_analyze = 'RDT_D1';

yoke_data = 0; % 1, set to 1 if you want to be prompted to yoke the number of trials analyzed, set to 0 otherwise

epoc_to_align = 'collectionTime';
ts1 = (uv.evtWin(1):uv.dt:uv.evtWin(2)-uv.dt);

% neuron_num = 0;
use_normalized_time = 0;

clear neuron_mean neuron_sem zall_mean zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all 

if any(startsWith(string(animalIDs), 'RDT-F'))
    if strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
        fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39', 'BLA_Insc_41'};

        for i = 1:length(fieldsToRemove)
            if isfield(final_SLEAP, fieldsToRemove{i})
                final_SLEAP = rmfield(final_SLEAP, fieldsToRemove{i});
            end
        end
    elseif strcmp('RDT_D2', session_to_analyze)

        fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_39'};

        for i = 1:length(fieldsToRemove)
            if isfield(final_SLEAP, fieldsToRemove{i})
                final_SLEAP = rmfield(final_SLEAP, fieldsToRemove{i});
            end
        end
    end

else
    if strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
        fieldsToRemove = {'BLA_Insc_35', 'BLA_Insc_38', 'BLA_Insc_41'};

        for i = 1:length(fieldsToRemove)
            if isfield(final_SLEAP, fieldsToRemove{i})
                final_SLEAP = rmfield(final_SLEAP, fieldsToRemove{i});
            end
        end

    end

end


%% FILTER TO GET UN-SHUFFLED DATA
iter = iter+1;
animalIDs = (fieldnames(final_SLEAP));
for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    if ~isfield(final_SLEAP.(currentanimal), session_to_analyze)
        caTraceTrials(1, 1:size(ts1, 2)) = nan;
        zall(1, 1:size(ts1, 2)) = nan;
        sem_all(ii, 1:size(ts1, 2)) = nan;
        zall_mean_all(ii,1:size(ts1, 2)) = nan;
    elseif isempty(fieldnames(final_SLEAP.(currentanimal).(session_to_analyze)))
        caTraceTrials(1, 1:size(ts1, 2)) = NaN;
        zall(1, 1:size(ts1, 2)) = NaN;
        sem_all(ii, 1:size(ts1, 2)) = NaN;
        zall_mean_all(ii,1:size(ts1, 2)) = NaN;
    elseif isfield(final_SLEAP.(currentanimal), session_to_analyze)
        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        % for labeling shock + 1 trials, for subsequent analysis
        if ~strcmp(session_to_analyze, 'SHOCK_TEST')
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

        [BehavData,trials, varargin_identity_class]=TrialFilter_test(BehavData, 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 3, 'SHK', 0);

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


        velocity_data = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.vel_cm_s'; %velocity_data = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.vel_filtered_2';
        % comment out below if you don't want to zscore traces prior to the
        % rest of the analysis
        % velocity_data = zscore(velocity_data, 0, 2);

        num_samples = size(velocity_data, 2);

        time_array_SLEAP = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data.idx_time;
        if ~strcmp('stTime',BehavData.Properties.VariableNames)
            BehavData.stTime = BehavData.TrialPossible - 5;
        end
        if ~strcmp('collectionTime',BehavData.Properties.VariableNames)
            BehavData.collectionTime = BehavData.choiceTime + 5;
        end



        % time_array = (0:(num_samples-1)) / sampling_frequency;
        eTS = BehavData.(epoc_to_align); %get time stamps

        zb_session = mean(velocity_data,2);
        zsd_session = std(velocity_data,[],2);
        % caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


        %calculate time windows for each event
        evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
        numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
        for u = 1:size(velocity_data,1)

            % initialize trial matrices
            caTraceTrials = NaN(size(eTS,1),numMeasurements); %
            unitTrace = velocity_data(u,:); %get trace
            if isempty(eTS) || size(eTS, 1) == 1
                caTraceTrials(1, 1:size(ts1, 2)) = nan;
                zall(1, 1:size(ts1, 2)) = nan;
                sem_all(ii, size(zall, 2)) = nan;
                zall_mean_all(ii,:) = nan;
            else
                [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array_SLEAP, zb_session, zsd_session, u, use_normalized_time);
                % [normalized_trial_ca, concatenated_normalized_trial_ca] = normalize_trials_in_time_fn(trial_ca);
                % time_normalized_ca{neuron_num} = concatenated_normalized_trial_ca;
             
                caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                zall = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                for z = 1:size(zall, 1)
                    % Apply Savitzky-Golay filter to each row
                    zall(z, :) = sgolayfilt(zall(z, :), 9, 21);
                end

                zall_array_session{ii} = zall_session(:, 1:size(ts1, 2));
                neuron_mean_unnormalized(ii,:) = nanmean(caTraceTrials,1);
                neuron_sem_unnormalized(ii,:) = nanstd(caTraceTrials,1)/(sqrt(size(caTraceTrials, 1)));
                zall_array{ii} = zall(:, 1:size(ts1, 2));
                zall_mouse{ii, iter}(u) = {zall(:, 1:size(ts1, 2))};
                unnormalized_mouse{ii} = caTraceTrials(:, 1:size(ts1, 2));
                unnormalized_sem_mouse{ii} = nanstd(caTraceTrials,1)/(sqrt(size(caTraceTrials, 1)));
                sem_mouse{ii, iter}(u) = {nanstd(zall,1)/(sqrt(size(zall, 1)))};
                caTraceTrials_mouse{ii, iter}(u) = {caTraceTrials(:, 1:size(ts1, 2))};
                neuron_mean_mouse_unnormalized{ii, iter}(u,: ) = mean(caTraceTrials, 1);
                neuron_sem_mouse_unnormalized{ii, iter}(u,: ) = nanstd(caTraceTrials,1)/(sqrt(size(caTraceTrials, 1)));
                neuron_mean_mouse{ii, iter}(u,: ) = mean(zall, 1);
                neuron_sem_mouse{ii, iter}(u,: ) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                zall_mean_all(ii,:) = nanmean(zall(:, 1:size(ts1, 2)));

                if size(zall, 1) == 1
                   

                else
                    sem_temp = nanstd(zall,1)/(sqrt(size(zall, 1)));
                    sem_all(ii,:) = sem_temp(:, 1:size(ts1, 2));
                end


                trials_per_mouse{ii, iter} = trials;
                clear zall caTraceTrials zb zsd sem_temp;
            end
        end
    end
end
zall_mean_all_array(iter) = {zall_mean_all};
neuron_mean_all_unnormalized(iter) = {neuron_mean_unnormalized};
neuron_sem_all_unnormalized(iter) = {neuron_sem_unnormalized};
sem_all_array(iter) = {sem_all};
varargin_list{iter,:} = varargin_identity_class;
behav_tbl_iter{iter, :} = behav_tbl_temp;
epoc_to_align_all{iter,:} = epoc_to_align;
all_filter_args{iter,:} = filter_args;

full_filter_string{iter} = strcat(epoc_to_align_all{iter,:}, '.', all_filter_args{iter,:});

clear behav_tbl_temp


%% 
% Get AUCs for the relevant periods for the 3 defined events
% Define time windows
pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [2 4];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

% Find indices corresponding to each time window
pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);


% Initialize arrays to store AUCs
% auc_pre_choice = zeros(size(neuron_mean_array));
% auc_post_choice = zeros(size(neuron_mean_array));
% auc_consumption = zeros(size(neuron_mean_array));
mouse_count = 0;
% Iterate over each element of neuron_mean_array
for j = 1:size(unnormalized_mouse, 2)
    prechoice_mean(:, j) = mean(neuron_mean_all_unnormalized{1, j}(:, pre_choice_indices), 2);
    postchoice_mean(:, j) = mean(neuron_mean_all_unnormalized{1, j}(:, post_choice_indices), 2);
    consumption_mean(:, j) = mean(neuron_mean_all_unnormalized{1, j}(:, consumption_indices), 2);
    for i = 1:size(unnormalized_mouse, 1)
        % Select data where exclusive_activated_session_1 is 1


        selected_data = unnormalized_mouse{i, j}{1, 1};
        % % Extract time variable (assuming it's named 'ts1')
        % ts1_data = ts1{i}(exclusive_activated_session_1{i} == 1);



        % Compute AUC for each time window
        % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
        auc_pre_choice(i, j) = mean(trapz(selected_data(:, pre_choice_indices), 2));
        auc_post_choice(i, j) = mean(trapz(selected_data(:, post_choice_indices), 2));
        auc_consumption(i, j) = mean(trapz(selected_data(:, consumption_indices), 2));

        mean_pre_choice(i, j) = mean(mean(selected_data(:, pre_choice_indices), 2));
        mean_post_choice(i, j) = mean(mean(selected_data(:, post_choice_indices), 2));
        mean_consumption(i, j) = mean(mean(selected_data(:, consumption_indices), 2));


    end
end



%% FILTER TO GET, ALIGNED SHUFFLED CA DATA

neuron_num = 0;
iter = iter+1;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
    [BehavData,trials,varargin_identity_class]=TrialFilter(BehavData,'REW',1.2);
    trials = cell2mat(trials);
    velocity_data = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);

    % shuffle data for comparison
    [num_cells, num_samples] = size(velocity_data);
    shuffled_data = zeros(num_cells, num_samples); % Preallocate matrix for efficiency
    shift_val = randi(num_samples); % Generate a random shift value for each signal RUAIRI RECOMMENDED KEEPING THE SAME SHIFT VAL, rather than randomizing per neuron. this is because then you keep the overall correlation b/w the neurons, but disrupt the relationship to the event timestamps
    for i = 1:num_cells
        % shift_val = randi(num_signals); % Generate a random shift value for each signal
        shuffled_data(i,:) = circshift(velocity_data(i,:), shift_val,2); % Perform the circular shuffle
    end
    velocity_data = shuffled_data;
    if strcmp(ca_data_type, 'S')
        velocity_data = full(velocity_data);

    end
    num_samples = size(velocity_data, 2);
    sampling_frequency = (final.(currentanimal).(session_to_analyze).uv.dt)*100;
    time_array_SLEAP = (0:(num_samples-1)) / sampling_frequency;
    eTS = BehavData.(epoc_to_align); %get time stamps

    zb_session = mean(velocity_data,2);
    zsd_session = std(velocity_data,[],2);
    % caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


    %calculate time windows for each event
    evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
    numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
    for u = 1:size(velocity_data,1)
        neuron_num = neuron_num+1; 
        % initialize trial matrices
        caTraceTrials = NaN(size(eTS,1),numMeasurements); %
        unitTrace = velocity_data(u,:); %get trace
        [zall_baselined, zall_window, zall_session, caTraceTrials] = align_and_zscore(unitTrace, eTS, uv, time_array_SLEAP, zb_session, zsd_session, u);

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


%%
%These data can be used to plot the median or mean choice
% time on a PCA graph, for example

behav_tbl_iter_single = behav_tbl_iter(1);

% Initialize the concatenated table
concatenatedTable = table();

% Iterate through the 3x1 cell array
for i = 1:numel(behav_tbl_iter_single)
    % Assuming each cell contains a 12x1 cell array of tables
    twelveByOneCellArray = behav_tbl_iter_single{i};
    
    % Initialize a temporary table to store the concatenated tables for this cell
    tempTable = table();
    
    % Iterate through the 12x1 cell array
    for j = 1:numel(twelveByOneCellArray)
        % Assuming each cell in the 12x1 cell array contains a table
        currentTable = twelveByOneCellArray{j};
        
        % Concatenate the current table to the temporary table vertically
        tempTable = vertcat(tempTable, currentTable);
    end
    
    % Concatenate the temporary table to the overall concatenated table vertically
    concatenatedTable = vertcat(concatenatedTable, tempTable);
end

median_choice_time_block_1 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_choice_time_block_2 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_choice_time_block_3 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));

median_collect_time_block_1 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_collect_time_block_2 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_collect_time_block_3 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));





median_start_time_from_choice = median(concatenatedTable.stTime - concatenatedTable.choiceTime);
median_collect_time_from_choice = median(concatenatedTable.collectionTime - concatenatedTable.choiceTime);

median_start_time_from_choice_large = median(concatenatedTable.stTime(concatenatedTable.bigSmall == 1.2) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 1.2));
median_start_time_from_choice_small = median(concatenatedTable.stTime(concatenatedTable.bigSmall == 0.3) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 0.3));

median_collect_time_from_choice_large = median(concatenatedTable.collectionTime(concatenatedTable.bigSmall == 1.2) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 1.2));
median_collect_time_from_choice_small = median(concatenatedTable.collectionTime(concatenatedTable.bigSmall == 0.3) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 0.3));




%% for hM4Di vs mCherry





% Find indices where Animals match valid_animalIDs
valid_idx = ismember(animalIDs, hM4Di_IDs);

% Filter path_length_table to only include valid animals
valid_animalIDs = animalIDs(valid_idx, :);
filtered_neuron_mean_unnormalized = neuron_mean_unnormalized(valid_idx, :);
filtered_neuron_sem_unnormalized = neuron_sem_unnormalized(valid_idx, :);

% Extract TreatmentCondition for the matched valid_animalIDs
[~, loc] = ismember(valid_animalIDs, hM4Di_IDs);

% Get the corresponding TreatmentCondition from risk_table
filtered_treatment_conditions = hM4Di_treatment_groups(loc);

% Find indices where TreatmentCondition is 'mCherry'
mCherry_idx = strcmp(filtered_treatment_conditions, 'mCherry');
hM4Di_idx = strcmp(filtered_treatment_conditions, 'hM4Di');


filtered_neuron_mean_unnormalized_mCherry = filtered_neuron_mean_unnormalized(mCherry_idx == 1, :);
filtered_neuron_mean_unnormalized_hM4Di = filtered_neuron_mean_unnormalized(hM4Di_idx == 1, :);
filtered_neuron_sem_unnormalized_mCherry = filtered_neuron_sem_unnormalized(mCherry_idx == 1, :);
filtered_neuron_sem_unnormalized_hM4Di = filtered_neuron_sem_unnormalized(hM4Di_idx == 1, :);


mean_large = nanmean(filtered_neuron_mean_unnormalized_mCherry, 1);
mean_small = nanmean(filtered_neuron_mean_unnormalized_hM4Di, 1);
sem_large = nanstd(filtered_neuron_mean_unnormalized_mCherry, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_mCherry, 1));
sem_small = nanstd(filtered_neuron_mean_unnormalized_hM4Di, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_hM4Di, 1));


behav_data_hM4Di = behav_tbl_iter{1, 1}(hM4Di_idx);
behav_data_mCherry = behav_tbl_iter{1, 1}(mCherry_idx);


concatenatedTable_mCherry = table();

for i = 1:numel(behav_data_mCherry)
    currentTable = behav_data_mCherry{i};
    concatenatedTable_mCherry = vertcat(concatenatedTable_mCherry, currentTable);
end

median_start_time_from_choice_mCherry = median(concatenatedTable_mCherry.stTime - concatenatedTable_mCherry.choiceTime);
median_collect_time_from_choice_mCherry = median(concatenatedTable_mCherry.collectionTime - concatenatedTable_mCherry.choiceTime);

concatenatedTable_hM4Di = table();

for i = 1:numel(behav_data_hM4Di)
    currentTable = behav_data_hM4Di{i};
    concatenatedTable_hM4Di = vertcat(concatenatedTable_hM4Di, currentTable);
end

median_start_time_from_choice_hM4Di = median(concatenatedTable_hM4Di.stTime - concatenatedTable_hM4Di.choiceTime);
median_collect_time_from_choice_hM4Di = median(concatenatedTable_hM4Di.collectionTime - concatenatedTable_hM4Di.choiceTime);



 mean_data_array = {filtered_neuron_mean_unnormalized_mCherry, filtered_neuron_mean_unnormalized_hM4Di};
 sem_data_array = {filtered_neuron_sem_unnormalized_mCherry, filtered_neuron_sem_unnormalized_hM4Di};

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-6 6], [0 5], 3);




figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 450; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xline(0);
xline(median_start_time_from_choice_mCherry, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_mCherry, 'r', {'Median', 'collect', 'latency'})
xline(median_start_time_from_choice_hM4Di, '--g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_hM4Di, '--r', {'Median', 'collect', 'latency'})
xlim([-6 6]);
ylim([0 5]);
% Set X-axis ticks
set(gca, 'XTick', [-6, 0, 6], 'YTick', [0 2.5 5]);
shadedErrorBar(ts1, mean_large, sem_large, 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, mean_small, sem_small, 'lineProps', {'color', 'k'});
plot(ts1, perm_p_sig+1)


%% for stGtACR vs mCherry





% Find indices where Animals match valid_animalIDs
valid_idx = ismember(animalIDs, stGtACR_IDs);

% Filter path_length_table to only include valid animals
valid_animalIDs = animalIDs(valid_idx, :);
filtered_neuron_mean_unnormalized = neuron_mean_unnormalized(valid_idx, :);
filtered_neuron_sem_unnormalized = neuron_sem_unnormalized(valid_idx, :);

% Extract TreatmentCondition for the matched valid_animalIDs
[~, loc] = ismember(valid_animalIDs, stGtACR_IDs);

% Get the corresponding TreatmentCondition from risk_table
filtered_treatment_conditions = stGtACR_treatment_groups(loc);

% Find indices where TreatmentCondition is 'mCherry'
mCherry_idx = strcmp(filtered_treatment_conditions, 'mCherry');
hM4Di_idx = strcmp(filtered_treatment_conditions, 'stGtACR');


filtered_neuron_mean_unnormalized_mCherry = filtered_neuron_mean_unnormalized(mCherry_idx == 1, :);
filtered_neuron_mean_unnormalized_hM4Di = filtered_neuron_mean_unnormalized(hM4Di_idx == 1, :);
filtered_neuron_sem_unnormalized_mCherry = filtered_neuron_sem_unnormalized(mCherry_idx == 1, :);
filtered_neuron_sem_unnormalized_hM4Di = filtered_neuron_sem_unnormalized(hM4Di_idx == 1, :);


mean_large = nanmean(filtered_neuron_mean_unnormalized_mCherry, 1);
mean_small = nanmean(filtered_neuron_mean_unnormalized_hM4Di, 1);
sem_large = nanstd(filtered_neuron_mean_unnormalized_mCherry, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_mCherry, 1));
sem_small = nanstd(filtered_neuron_mean_unnormalized_hM4Di, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_hM4Di, 1));


behav_data_hM4Di = behav_tbl_iter{1, 1}(hM4Di_idx);
behav_data_mCherry = behav_tbl_iter{1, 1}(mCherry_idx);


concatenatedTable_mCherry = table();

for i = 1:numel(behav_data_mCherry)
    currentTable = behav_data_mCherry{i};
    concatenatedTable_mCherry = vertcat(concatenatedTable_mCherry, currentTable);
end

median_start_time_from_choice_mCherry = median(concatenatedTable_mCherry.stTime - concatenatedTable_mCherry.choiceTime);
median_collect_time_from_choice_mCherry = median(concatenatedTable_mCherry.collectionTime - concatenatedTable_mCherry.choiceTime);

concatenatedTable_hM4Di = table();

for i = 1:numel(behav_data_hM4Di)
    currentTable = behav_data_hM4Di{i};
    concatenatedTable_hM4Di = vertcat(concatenatedTable_hM4Di, currentTable);
end

median_start_time_from_choice_hM4Di = median(concatenatedTable_hM4Di.stTime - concatenatedTable_hM4Di.choiceTime);
median_collect_time_from_choice_hM4Di = median(concatenatedTable_hM4Di.collectionTime - concatenatedTable_hM4Di.choiceTime);



 mean_data_array = {filtered_neuron_mean_unnormalized_mCherry, filtered_neuron_mean_unnormalized_hM4Di};
 sem_data_array = {filtered_neuron_sem_unnormalized_mCherry, filtered_neuron_sem_unnormalized_hM4Di};

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-6 6], [0 5], 3);




figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 450; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xline(0);
xline(median_start_time_from_choice_mCherry, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_mCherry, 'r', {'Median', 'collect', 'latency'})
xline(median_start_time_from_choice_hM4Di, '--g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_hM4Di, '--r', {'Median', 'collect', 'latency'})
xlim([-6 6]);
ylim([0 5]);
% Set X-axis ticks
set(gca, 'XTick', [-6, 0, 6], 'YTick', [0 2.5 5]);
shadedErrorBar(ts1, mean_large, sem_large, 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, mean_small, sem_small, 'lineProps', {'color', 'k'});
plot(ts1, perm_p_sig+1)

%% for PdCO vs mCherry





% Find indices where Animals match valid_animalIDs
valid_idx = ismember(animalIDs, PdCO_IDs);

% Filter path_length_table to only include valid animals
valid_animalIDs = animalIDs(valid_idx, :);
filtered_neuron_mean_unnormalized = neuron_mean_unnormalized(valid_idx, :);
filtered_neuron_sem_unnormalized = neuron_sem_unnormalized(valid_idx, :);

% Extract TreatmentCondition for the matched valid_animalIDs
[~, loc] = ismember(valid_animalIDs, PdCO_IDs);

% Get the corresponding TreatmentCondition from risk_table
filtered_treatment_conditions = PdCO_treatment_groups(loc);

% Find indices where TreatmentCondition is 'mCherry'
mCherry_idx = strcmp(filtered_treatment_conditions, 'mCherry');
hM4Di_idx = strcmp(filtered_treatment_conditions, 'PdCO');


filtered_neuron_mean_unnormalized_mCherry = filtered_neuron_mean_unnormalized(mCherry_idx == 1, :);
filtered_neuron_mean_unnormalized_hM4Di = filtered_neuron_mean_unnormalized(hM4Di_idx == 1, :);
filtered_neuron_sem_unnormalized_mCherry = filtered_neuron_sem_unnormalized(mCherry_idx == 1, :);
filtered_neuron_sem_unnormalized_hM4Di = filtered_neuron_sem_unnormalized(hM4Di_idx == 1, :);


mean_large = nanmean(filtered_neuron_mean_unnormalized_mCherry, 1);
mean_small = nanmean(filtered_neuron_mean_unnormalized_hM4Di, 1);
sem_large = nanstd(filtered_neuron_mean_unnormalized_mCherry, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_mCherry, 1));
sem_small = nanstd(filtered_neuron_mean_unnormalized_hM4Di, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_hM4Di, 1));


behav_data_hM4Di = behav_tbl_iter{1, 1}(hM4Di_idx);
behav_data_mCherry = behav_tbl_iter{1, 1}(mCherry_idx);


concatenatedTable_mCherry = table();

for i = 1:numel(behav_data_mCherry)
    currentTable = behav_data_mCherry{i};
    concatenatedTable_mCherry = vertcat(concatenatedTable_mCherry, currentTable);
end

median_start_time_from_choice_mCherry = median(concatenatedTable_mCherry.stTime - concatenatedTable_mCherry.choiceTime);
median_collect_time_from_choice_mCherry = median(concatenatedTable_mCherry.collectionTime - concatenatedTable_mCherry.choiceTime);

concatenatedTable_hM4Di = table();

for i = 1:numel(behav_data_hM4Di)
    currentTable = behav_data_hM4Di{i};
    concatenatedTable_hM4Di = vertcat(concatenatedTable_hM4Di, currentTable);
end

median_start_time_from_choice_hM4Di = median(concatenatedTable_hM4Di.stTime - concatenatedTable_hM4Di.choiceTime);
median_collect_time_from_choice_hM4Di = median(concatenatedTable_hM4Di.collectionTime - concatenatedTable_hM4Di.choiceTime);



 mean_data_array = {filtered_neuron_mean_unnormalized_mCherry, filtered_neuron_mean_unnormalized_hM4Di};
 sem_data_array = {filtered_neuron_sem_unnormalized_mCherry, filtered_neuron_sem_unnormalized_hM4Di};

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-6 6], [0 5], 3);




figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 450; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xline(0);
xline(median_start_time_from_choice_mCherry, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_mCherry, 'r', {'Median', 'collect', 'latency'})
xline(median_start_time_from_choice_hM4Di, '--g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_hM4Di, '--r', {'Median', 'collect', 'latency'})
xlim([-6 6]);
ylim([0 5]);
% Set X-axis ticks
set(gca, 'XTick', [-6, 0, 6], 'YTick', [0 2.5 5]);
shadedErrorBar(ts1, mean_large, sem_large, 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, mean_small, sem_small, 'lineProps', {'color', 'k'});
plot(ts1, perm_p_sig+1)

%%
%% for Females vs Males





% Find indices where Animals match valid_animalIDs
valid_idx = ismember(animalIDs, ChrimsonR_IDs);

% Filter path_length_table to only include valid animals
valid_animalIDs = animalIDs(valid_idx, :);
filtered_neuron_mean_unnormalized = neuron_mean_unnormalized(valid_idx, :);
filtered_neuron_sem_unnormalized = neuron_sem_unnormalized(valid_idx, :);

% Extract TreatmentCondition for the matched valid_animalIDs
[~, loc] = ismember(valid_animalIDs, ChrimsonR_IDs);

% Get the corresponding TreatmentCondition from risk_table
filtered_treatment_conditions = ChrimsonR_treatment_groups(loc);

% Find indices where TreatmentCondition is 'mCherry'
mCherry_idx = strcmp(filtered_treatment_conditions, 'Female');
hM4Di_idx = strcmp(filtered_treatment_conditions, 'Male');


filtered_neuron_mean_unnormalized_mCherry = filtered_neuron_mean_unnormalized(mCherry_idx == 1, :);
filtered_neuron_mean_unnormalized_hM4Di = filtered_neuron_mean_unnormalized(hM4Di_idx == 1, :);
filtered_neuron_sem_unnormalized_mCherry = filtered_neuron_sem_unnormalized(mCherry_idx == 1, :);
filtered_neuron_sem_unnormalized_hM4Di = filtered_neuron_sem_unnormalized(hM4Di_idx == 1, :);


mean_large = nanmean(filtered_neuron_mean_unnormalized_mCherry, 1);
mean_small = nanmean(filtered_neuron_mean_unnormalized_hM4Di, 1);
sem_large = nanstd(filtered_neuron_mean_unnormalized_mCherry, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_mCherry, 1));
sem_small = nanstd(filtered_neuron_mean_unnormalized_hM4Di, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_hM4Di, 1));


behav_data_hM4Di = behav_tbl_iter{iter, 1}(hM4Di_idx);
behav_data_mCherry = behav_tbl_iter{iter, 1}(mCherry_idx);


concatenatedTable_mCherry = table();

for i = 1:numel(behav_data_mCherry)
    currentTable = behav_data_mCherry{i};
    concatenatedTable_mCherry = vertcat(concatenatedTable_mCherry, currentTable);
end

median_start_time_from_choice_mCherry = median(concatenatedTable_mCherry.stTime - concatenatedTable_mCherry.choiceTime);
median_collect_time_from_choice_mCherry = median(concatenatedTable_mCherry.collectionTime - concatenatedTable_mCherry.choiceTime);

concatenatedTable_hM4Di = table();

for i = 1:numel(behav_data_hM4Di)
    currentTable = behav_data_hM4Di{i};
    concatenatedTable_hM4Di = vertcat(concatenatedTable_hM4Di, currentTable);
end

median_start_time_from_choice_hM4Di = median(concatenatedTable_hM4Di.stTime - concatenatedTable_hM4Di.choiceTime);
median_collect_time_from_choice_hM4Di = median(concatenatedTable_hM4Di.collectionTime - concatenatedTable_hM4Di.choiceTime);



 mean_data_array = {filtered_neuron_mean_unnormalized_mCherry, filtered_neuron_mean_unnormalized_hM4Di};
 sem_data_array = {filtered_neuron_sem_unnormalized_mCherry, filtered_neuron_sem_unnormalized_hM4Di};

 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-4 4], [0 15], 3);




figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 450; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xline(0);
xline(median_start_time_from_choice_mCherry, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_mCherry, 'g', {'Median', 'collect', 'latency'})
xline(median_start_time_from_choice_hM4Di, 'k', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice_hM4Di, 'k', {'Median', 'collect', 'latency'})
xlim([-8 8]);
ylim([0 15]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8], 'YTick', [0:5:15]);
shadedErrorBar(ts1, mean_large, sem_large, 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, mean_small, sem_small, 'lineProps', {'color', 'k'});
plot(ts1, perm_p_sig+1)


pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 1];     % Post-choice period: 0 to 2 s
consumption_window = [1 3];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

% Find indices corresponding to each time window
pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);


mouse_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(filtered_neuron_mean_unnormalized_mCherry, 1)
    mouse_count = mouse_count+1;
    % Compute AUC for each time window
    % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
    action_auc_pre_choice_females(iter, mouse_count) = trapz(filtered_neuron_mean_unnormalized_mCherry(i, pre_choice_indices));
    action_auc_post_choice_females(iter, mouse_count) = trapz(filtered_neuron_mean_unnormalized_mCherry(i, post_choice_indices));
    action_auc_consumption_females(iter, mouse_count) = trapz(filtered_neuron_mean_unnormalized_mCherry(i, consumption_indices));
end

mouse_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(filtered_neuron_mean_unnormalized_hM4Di, 1)
    mouse_count = mouse_count+1;
    % Compute AUC for each time window
    % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
    action_auc_pre_choice_males(iter, mouse_count) = trapz(filtered_neuron_mean_unnormalized_hM4Di(i, pre_choice_indices));
    action_auc_post_choice_males(iter, mouse_count) = trapz(filtered_neuron_mean_unnormalized_hM4Di(i, post_choice_indices));
    action_auc_consumption_males(iter, mouse_count) = trapz(filtered_neuron_mean_unnormalized_hM4Di(i, consumption_indices));
end


%%
%%
% Create sample data for testing (replace with your actual data)
% mCherry_data_means = [action_auc_pre_choice_females];
% hM4Di_data_mean = [action_auc_pre_choice_males];

% mCherry_data_means = [action_auc_post_choice_females];
% hM4Di_data_mean = [action_auc_post_choice_males];

mCherry_data_means = [action_auc_consumption_females];
hM4Di_data_mean = [action_auc_consumption_males];

% Combine the datasets for easier handling
all_data = {mCherry_data_means, hM4Di_data_mean};

% Calculate means for bar heights (mean of each row)
means = [mean(mCherry_data_means, 2), mean(hM4Di_data_mean, 2)];

% Calculate SEM for error bars (SEM = std/sqrt(n))
sem = [std(mCherry_data_means, 0, 2) ./ sqrt(size(mCherry_data_means, 2)), ...
       std(hM4Di_data_mean, 0, 2) ./ sqrt(size(hM4Di_data_mean, 2))];

% Grouped positions for the bars
x = [1, 2, 3; 5, 6, 7]; % First group at positions 1,2,3; second group at 5,6,7

% Create figure
figure('Position', [100, 100, 300, 450]); % Made wider to accommodate 6 bars
hold on;

% Loop through each dataset
for i = 1:size(all_data, 2)
    data = all_data{i}; % Current dataset
    
    % Plot bars for each row (3 bars per dataset)
    for row = 1:3
        bar_x = x(i, row); % Position for the current bar
        bar(bar_x, means(row, i), 0.6); % Plot bar
        
        % Add error bars
        errorbar(bar_x, means(row, i), sem(row, i), 'k', 'LineWidth', 1.5, 'CapSize', 8);
    end
    
    % Overlay scatter points and connect with lines
    jittered_x = zeros(size(data)); % To store jittered x-coordinates
    
    if i <= 1
        % Add jitter to x-coordinates and plot scatter points
        for col = 1:size(data, 2) % For each column (data point)
            for row = 1:3 % For each row
                scatter_x = x(i, row) + (rand(1, 1) - 0.5) * 0.3; % Add jitter
                jittered_x(row, col) = scatter_x; % Store jittered x-coordinate
                scatter(scatter_x, data(row, col), 50, 'green', 'filled', 'square');
            end
        end
        % Connect scatter points with lines (each column connected across rows)
        for col = 1:size(data, 2)
            plot(jittered_x(:, col), data(:, col), 'green-', 'LineWidth', 0.5);
        end
    else
        % Add jitter to x-coordinates and plot scatter points
        for col = 1:size(data, 2) % For each column (data point)
            for row = 1:3 % For each row
                scatter_x = x(i, row) + (rand(1, 1) - 0.5) * 0.3; % Add jitter
                jittered_x(row, col) = scatter_x; % Store jittered x-coordinate
                scatter(scatter_x, data(row, col), 40, 'k', 'filled');
            end
        end
        % Connect scatter points with lines (each column connected across rows)
        for col = 1:size(data, 2)
            plot(jittered_x(:, col), data(:, col), 'k-', 'LineWidth', 0.5);
        end
    end
end

% Adjustments for aesthetics
set(gca, 'XTick', [2, 6], 'XTickLabel', {'F', 'M'});
% xlabel('Groups');
% ylabel('Values');
ylim([0 40]);
yticks([0:10:40]);
xlim([0 8]);
hold off;

%% for analyzing shock test, grouping ranges of intensities together

% ranges = [2:7; 8:13; 14:19; 20:25]

ranges = [2:13; 14:25]

for hh = 1:size(unnormalized_mouse, 2)


    velocity_data_extracted = unnormalized_mouse{1, hh};
    if isempty(velocity_data_extracted)

        mouse_velocity_in_ranges_all{hh} = nan;
        mouse_sem_in_ranges{hh} = nan;
    else
        for gg = 1:size(velocity_data_extracted, 2)
            for mm = 1:size(ranges, 1)

                mouse_velocity_in_ranges(mm, :) = mean(velocity_data_extracted(ranges(mm, :), :));
                mouse_sem_in_ranges(mm, :) = std(velocity_data_extracted(ranges(mm, :), :)/sqrt(size(velocity_data_extracted(ranges(mm, :), :), 1)));

            end


        end

        mouse_velocity_in_ranges_all{hh} = mouse_velocity_in_ranges;
        mouse_sem_in_ranges_all{hh} = mouse_sem_in_ranges;
        clear mouse_velocity_in_ranges mouse_sem_in_ranges;
    end

end

% Find indices where Animals match valid_animalIDs
% valid_idx = ismember(animalIDs, hM4Di_IDs);
valid_idx = ismember(animalIDs, stGtACR_IDs);

% Filter path_length_table to only include valid animals
valid_animalIDs = animalIDs(valid_idx, :);
filtered_mouse_velocity_in_ranges_all = mouse_velocity_in_ranges_all(1, valid_idx);
filtered_mouse_sem_in_ranges_all = mouse_sem_in_ranges_all(1, valid_idx);

% Initialize the output matrices
row1 = zeros(size(valid_animalIDs, 1), size(ts1, 2));
row2 = zeros(size(valid_animalIDs, 1), size(ts1, 2));
row3 = zeros(size(valid_animalIDs, 1), size(ts1, 2));
row4 = zeros(size(valid_animalIDs, 1), size(ts1, 2));


sem_data_row1 = zeros(size(valid_animalIDs, 1), size(ts1, 2));
sem_data_row2 = zeros(size(valid_animalIDs, 1), size(ts1, 2));
sem_data_row3 = zeros(size(valid_animalIDs, 1), size(ts1, 2));
sem_data_row4 = zeros(size(valid_animalIDs, 1), size(ts1, 2));

% Loop through each cell and extract the rows
for i = 1:size(valid_animalIDs, 1)
    data = filtered_mouse_velocity_in_ranges_all{i}; % Extract the 4x60 double array
    row1(i, :) = data(1, :);
    row2(i, :) = data(2, :);
    % row3(i, :) = data(3, :);
    % row4(i, :) = data(4, :);
    sem_data = filtered_mouse_sem_in_ranges_all{i}; % Extract the 4x60 double array
    sem_data_row1(i, :) = sem_data(1, :);
    sem_data_row2(i, :) = sem_data(2, :);
    % sem_data_row3(i, :) = sem_data(3, :);
    % sem_data_row4(i, :) = sem_data(4, :);



end


%% for PLOTTING SHOCK TEST DATA, BROKEN INTO SHOCK RANGES


data_to_plot = row1; 
sem_to_plot = sem_data_row1;

% Extract TreatmentCondition for the matched valid_animalIDs
[~, loc] = ismember(valid_animalIDs, hM4Di_IDs);

% Get the corresponding TreatmentCondition from risk_table
filtered_treatment_conditions = hM4Di_treatment_groups(loc);

% Find indices where TreatmentCondition is 'mCherry'
mCherry_idx = strcmp(filtered_treatment_conditions, 'mCherry');
hM4Di_idx = strcmp(filtered_treatment_conditions, 'hM4Di');


filtered_neuron_mean_unnormalized_mCherry = data_to_plot(mCherry_idx == 1, :);
filtered_neuron_mean_unnormalized_hM4Di = data_to_plot(hM4Di_idx == 1, :);

filtered_neuron_sem_unnormalized_mCherry = sem_to_plot(mCherry_idx == 1, :);
filtered_neuron_sem_unnormalized_hM4Di = sem_to_plot(hM4Di_idx == 1, :);

mean_large = nanmean(filtered_neuron_mean_unnormalized_mCherry, 1);
mean_small = nanmean(filtered_neuron_mean_unnormalized_hM4Di, 1);
sem_large = nanstd(filtered_neuron_mean_unnormalized_mCherry, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_mCherry, 1));
sem_small = nanstd(filtered_neuron_mean_unnormalized_hM4Di, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_hM4Di, 1));



figure;
hold on
% Create a histogram for allCorrelations

width = 150; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-2 6]);
ylim([0 10]);
% Set X-axis ticks
set(gca, 'XTick', [-2, 0, 2, 4, 6], 'YTick', [0 5 10]);
shadedErrorBar(ts1, mean_large, sem_large, 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, mean_small, sem_small, 'lineProps', {'color', 'k'});

mean_data_array = {filtered_neuron_mean_unnormalized_mCherry, filtered_neuron_mean_unnormalized_hM4Di};
sem_data_array = {filtered_neuron_sem_unnormalized_mCherry, filtered_neuron_sem_unnormalized_hM4Di};

[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-6 6], [0 10], 3);



%% for PLOTTING SHOCK TEST DATA, BROKEN INTO SHOCK RANGES


data_to_plot = row2; 
sem_to_plot = sem_data_row2;


% Extract TreatmentCondition for the matched valid_animalIDs
[~, loc] = ismember(valid_animalIDs, stGtACR_IDs);

% Get the corresponding TreatmentCondition from risk_table
filtered_treatment_conditions = stGtACR_treatment_groups(loc);

% Find indices where TreatmentCondition is 'mCherry'
mCherry_idx = strcmp(filtered_treatment_conditions, 'mCherry');
hM4Di_idx = strcmp(filtered_treatment_conditions, 'stGtACR');


filtered_neuron_mean_unnormalized_mCherry = data_to_plot(mCherry_idx == 1, :);
filtered_neuron_mean_unnormalized_hM4Di = data_to_plot(hM4Di_idx == 1, :);

filtered_neuron_sem_unnormalized_mCherry = sem_to_plot(mCherry_idx == 1, :);
filtered_neuron_sem_unnormalized_hM4Di = sem_to_plot(hM4Di_idx == 1, :);

mean_large = nanmean(filtered_neuron_mean_unnormalized_mCherry, 1);
mean_small = nanmean(filtered_neuron_mean_unnormalized_hM4Di, 1);
sem_large = nanstd(filtered_neuron_mean_unnormalized_mCherry, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_mCherry, 1));
sem_small = nanstd(filtered_neuron_mean_unnormalized_hM4Di, 0, 1) ./ sqrt(size(filtered_neuron_mean_unnormalized_hM4Di, 1));



figure;
hold on
% Create a histogram for allCorrelations

width = 150; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-2 6]);
ylim([0 10]);
% Set X-axis ticks
set(gca, 'XTick', [-2, 0, 2, 4, 6], 'YTick', [0 5 10]);
shadedErrorBar(ts1, mean_large, sem_large, 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, mean_small, sem_small, 'lineProps', {'color', 'k'});

mean_data_array = {filtered_neuron_mean_unnormalized_mCherry, filtered_neuron_mean_unnormalized_hM4Di};
sem_data_array = {filtered_neuron_sem_unnormalized_mCherry, filtered_neuron_sem_unnormalized_hM4Di};

[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-6 6], [0 10], 3);
