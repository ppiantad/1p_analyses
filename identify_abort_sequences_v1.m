%% This code produces thesaame output as Tommy Gunawan's R code, so good to go

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

%% get trials where shocks were followed by abort (get first abort) or choice (any choice, large/small)

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
            if data.shock(i) == 1
                % Initialize a flag to determine if we found a matching condition
                match_found = false;
                behav_data_extracted(rows_for_new_table, :) = data(i, :);
                rows_for_new_table = rows_for_new_table + 1;
                % Look for a matching condition in the subsequent rows
                j = i + 1;
                while j <= height(data) && ~match_found
                    % Check if the trial is an abort (type_binary == 1)
                    if data.type_binary(j) == 1 || data.type_binary(j) == 2
                        % Increment abort counter and mark as found
                        total_shock_abort_trials = total_shock_abort_trials + 1;
                        match_found = true;
                        % behav_data_shk_to_abort(rows_for_new_table+1, :) = data([i, j], :);
                        behav_data_extracted(rows_for_new_table, :) = data(j, :);
                        rows_for_new_table = rows_for_new_table + 1;
                        % Check if the trial has a bigSmall event (bigSmall == 1)
                    elseif data.bigSmall(j) == 1.2 || data.bigSmall(j) == 0.3
                        % Increment no-abort counter and mark as found
                        total_shock_no_abort_trials = total_shock_no_abort_trials + 1;
                        match_found = true;
                        behav_data_extracted(rows_for_new_table, :) = data(j, :);
                        rows_for_new_table = rows_for_new_table + 1;
                    else
                        % Move to the next row if no match is found
                        j = j + 1;
                    end
                end
                % Set `i` to `j` to continue from the row after the found condition
                i = j;
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

shk_abort_to_shk_choice_ratio = total_shock_abort_trials_array./total_shock_no_abort_trials_array

%% get abort sequences following shocks
shock_abort_sequence = 0;
neuron_num = 0;
animalIDs = (fieldnames(final));
for ii = 1:size(animalIDs, 1)
    sequence = 0;
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
            if data.shock(i) == 1
                % Initialize a flag to determine if we found a matching condition
                match_found = false;
                behav_data_extracted(rows_for_new_table, :) = data(i, :);
                rows_for_new_table = rows_for_new_table + 1;
                % Look for a matching condition in the subsequent rows
                shock_abort_sequence = 0;
                sequence = sequence + 1;
                j = i + 1;
                while j <= height(data) && ~match_found
                    % Check if the trial is an abort (type_binary == 1)
                    if data.type_binary(j) == 1 || data.type_binary(j) == 2
                        % Increment abort counter and mark as found

                        shock_abort_sequence = shock_abort_sequence + 1;
                        % behav_data_shk_to_abort(rows_for_new_table+1, :) = data([i, j], :);

                        j = j + 1;
                        % Check if the trial has a bigSmall event (bigSmall == 1)
                    elseif data.bigSmall(j) == 1.2 || data.bigSmall(j) == 0.3
                        % Increment no-abort counter and mark as found

                        match_found = true;

                    else
                        j = j + 1;
                    end
                end
                shock_abort_sequence_all(sequence) = shock_abort_sequence;
                % Set `i` to `j` to continue from the row after the found condition
                i = j;
            else
                % Move to the next row if no shock was found
                i = i + 1;
            end
        end
        behav_data_extracted_array{ii} = behav_data_extracted;
        shock_abort_sequence_all_array{ii} = shock_abort_sequence_all;
        clear behav_data_extracted shock_abort_sequence_all
    end
end

%% Plot https://stackoverflow.com/questions/54528239/boxplot-for-paired-observations

combined_data = [total_shock_abort_trials_array ; total_shock_no_abort_trials_array]';
figure();
coordLineStyle = 'k.';
boxplot(combined_data, 'Symbol', coordLineStyle); hold on;
parallelcoords(combined_data, 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
  'Marker', '.', 'MarkerSize', 10);

TF = isoutlier(combined_data, 'grubbs');


%% get Ca data
iter = 0;
iter = iter+1;
neuron_num = 0;
animalIDs = (fieldnames(final));


for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = behav_data_extracted_array{1, ii};
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



                clear zall caTraceTrials zb zsd sem_temp;
            end
        end
    end
end
zall_mean_all_array(iter) = {zall_mean_all};
neuron_mean_all_unnormalized(iter) = {neuron_mean_unnormalized};
sem_all_array(iter) = {sem_all};
% varargin_list{iter,:} = varargin_identity_class;
% behav_tbl_iter{iter, :} = behav_tbl_temp;
% epoc_to_align_all{iter,:} = epoc_to_align;
% all_filter_args{iter,:} = filter_args;

% full_filter_string{iter} = strcat(epoc_to_align_all{iter,:}, '.', all_filter_args{iter,:});

% clear behav_tbl_temp

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


%%
array_for_means = 1; 

% Initialize the new cell array to store the mean values
meanZallMouse = cell(size(zall_mouse, 2), 1);

% Define the time range for 0 to 2 seconds
% timeRange = (ts1 >= -4) & (ts1 <= 0);
timeRange = (ts1 >= 0) & (ts1 <= 2);
% timeRange = (ts1 >= 1) & (ts1 <= 3);


% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_mouse)
    % Get the current nested cell array
    nestedCellArray_1 = zall_mouse{i, array_for_means};

    % Initialize the nested cell array for storing mean values
    meanNestedCellArray = cell(size(nestedCellArray_1));
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(nestedCellArray_1)
        % Get the current double array
        currentArray = nestedCellArray_1{j};

        meanValues = mean(currentArray(:, timeRange), 2);
        % meanValues = max(currentArray(:, timeRange), [], 2);
        
        % Store the mean values in the corresponding cell of the nested cell array
        meanNestedCellArray{j} = meanValues;
    end
    
    % Store the nested cell array of mean values in the corresponding cell of the main cell array
    meanZallMouse{i} = meanNestedCellArray;
end

% Now, meanZallMouse contains the mean activity for each row in the time period 0 to 2 seconds
% Each cell in meanZallMouse contains a nested cell array with the


%%
% Initialize arrays to store the data

shockResponses_all_mice = [];
trialChoices_all_mice = [];
animal_IDs_all_mice = [];  

% Iterate through each level of meanZallMouse
for kk = 1:length(meanZallMouse) %1:length(meanZallMouse)
    shockResponses = [];
    trialChoices = [];
    animal_IDs = [];  
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{kk};
    
    % comment out if necessary
    % limit to just SHOCK activated neurons
    % get this by loading the 10x variable dataset, or similar, where array
    % 4 is data from SHK == 1
    % meanNestedCellArray = meanNestedCellArray(respClass_all_array_mouse{kk, 4} == 1);

    % Extract the data for easy reference
    data = behav_data_extracted_array{1, kk};

    % Find indices of even-numbered rows
    % because of the above code, the even rows contain data from shk + 1
    % trials, while the odd rows contain data from shk trials themselves
    even_rows = 2:2:height(data);
    odd_rows = 1:2:height(data);
    % Initialize the categorical_outcome variable with NaNs
    categorical_outcome = nan(height(data), 1);

    % Loop through each even-numbered row
    % determine whether a trial a choice (bigSmall > 0) or an abort
    % (type_binary > 0)
    for i = even_rows
        % Check if bigSmall > 0
        if data.bigSmall(i) > 0
            categorical_outcome(i) = 1;
            % Check if type_binary > 0
        elseif data.type_binary(i) > 0
            categorical_outcome(i) = 0;
        end
    end
    categorical_outcome_mouse{kk} = categorical_outcome;
    categorical_outcome = categorical_outcome(~isnan(categorical_outcome));

    % Get the trial choices for the current mouse
    currentTrialChoices = categorical_outcome;
    % currentTrialChoices = currentTrialChoices(currentTrialChoices(:,1) ~= 1, 2);
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        
        meanValues = meanNestedCellArray{j};
        % Get the current mean values array
        % only use data from shk trials, these will be the predictors
        meanValues = meanValues(odd_rows);
        
        % Here we use the meanValues as is, no averaging across trials
        % Flatten the meanValues to a single row vector, if needed
        shockResponses = [shockResponses; meanValues];
        
        % Append the corresponding trial choice
        trialChoices = [trialChoices; currentTrialChoices];
        animal_IDs = [animal_IDs; repmat(kk, length(meanValues), 1)];
    end
    shockResponses_all_mice = [shockResponses_all_mice; shockResponses];
    trialChoices_all_mice = [trialChoices_all_mice; trialChoices];
    animal_IDs_all_mice = [animal_IDs_all_mice; animal_IDs];
end

%% save data if necessary (sending to Tommy to see if he can help w/ regression

% % Save animal_IDs_all_mice to a CSV file
% writematrix(animal_IDs_all_mice, 'animal_IDs_all_mice.csv');
% 
% % Save shockResponses_all_mice to a CSV file
% writematrix(shockResponses_all_mice, 'shockResponses_all_mice.csv');
% 
% % Save trialChoices_all_mice to a CSV file
% writematrix(trialChoices_all_mice, 'trialChoices_all_mice.csv');



%%
% Step 1: Encode the choices
% Convert trialChoices_all_mice to a binary variable (1 for 1.2, 0 for 0.3)
trialChoices_binary = trialChoices_all_mice;

% Step 2: Fit a logistic regression model
% Use the calcium imaging data in shockResponses_all_mice as the predictor
X = shockResponses_all_mice;  % Predictor (calcium response)
y = trialChoices_binary;       % Response (binary choices)

% Fit logistic regression model
% Use MATLAB's fitglm function with 'binomial' distribution for logistic regression
mdl = fitglm(X, y, 'Distribution', 'binomial');

% Display model summary
disp(mdl)

% Step 3: Calculate Odds Ratios
% Extract coefficients
coefficients = mdl.Coefficients.Estimate;

% Calculate odds ratios by exponentiating the coefficients
odds_ratios = exp(coefficients);

% Display odds ratios
fprintf('Odds Ratio for Intercept: %.4f\n', odds_ratios(1));
fprintf('Odds Ratio for Predictor (shockResponses_all_mice): %.4f\n', odds_ratios(2));

% Step 4: Make predictions
% Get predicted probabilities for the choices
predictedProbabilities = predict(mdl, X);

% Optional: Classify based on probability threshold (e.g., 0.5)
predictedChoices = predictedProbabilities > 0.5;

% Step 5: Assess model accuracy (if desired)
% Calculate accuracy as the proportion of correctly predicted choices
accuracy = mean(predictedChoices == y);
fprintf('Prediction Accuracy: %.2f%%\n', accuracy * 100);

%%

% Step 1: Prepare data for plotting
% Sort shockResponses_all_mice for a smoother curve
[X_sorted, sortIdx] = sort(shockResponses_all_mice);
predictedProbabilities_sorted = predict(mdl, X_sorted);

% Step 2: Plot the actual data points
figure;
hold on;
% Plot raw data: shockResponses vs actual choices
scatter(shockResponses_all_mice, trialChoices_binary, 'b', 'filled', 'DisplayName', 'Actual Choices');
ylabel('Choice (1 = Large, 0 = Small)');
xlabel('Shock Response');
title('Logistic Regression: Shock Response vs Trial Choice');

% Step 3: Plot the fitted logistic regression curve
plot(X_sorted, predictedProbabilities_sorted, 'r-', 'LineWidth', 2, 'DisplayName', 'Predicted Probability');

% Add legend
legend('Location', 'best');

% Step 4: Improve plot aesthetics (optional)
ylim([-0.1, 1.1]); % Set y-axis limits to make the plot clearer
grid on;
hold off;
%%
% Assume shockResponses_all_mice and trialChoices_all_mice are defined

% Filter shock responses based on trial choices
shock_responses_choice_0 = shockResponses_all_mice(trialChoices_all_mice == 0);
shock_responses_choice_1 = shockResponses_all_mice(trialChoices_all_mice == 1);

% Plot the histograms in a single figure with two subplots
figure;

% Plot histogram for trials where trialChoices_all_mice == 0
subplot(2, 1, 1);
histogram(shock_responses_choice_0);
title('Histogram of Shock Responses (trialChoices == Abort)');
xlabel('Shock Response');
ylabel('Frequency');

% Plot histogram for trials where trialChoices_all_mice == 1
subplot(2, 1, 2);
histogram(shock_responses_choice_1);
title('Histogram of Shock Responses (trialChoices == Choice)');
xlabel('Shock Response');
ylabel('Frequency');

% Adjust layout for clarity
% sgtitle('Shock Responses for Different Trial Choices');

%%
neuron_num = 0;

for zz = 1:size(zall_mouse, 1)
    first_level_data = zall_mouse{zz, 1};
    get_trial_data = categorical_outcome_mouse{zz}; 
    % Find indices of rows where the value is 0, and get the index of the previous row
    rows_preceding_0 = find(get_trial_data(2:end) == 0);

    % Find indices of rows where the value is 1, and get the index of the previous row
    rows_preceding_1 = find(get_trial_data(2:end) == 1);
    for xx = 1:size(first_level_data, 2)
        second_level_data = first_level_data{xx};
        neuron_num=neuron_num+1;
        mean_0_trials(neuron_num, :) = mean(second_level_data(rows_preceding_0, :));
        mean_1_trials(neuron_num, :) = mean(second_level_data(rows_preceding_1, :));
        



    end


end

figure; plot(ts1, mean(mean_0_trials)); hold on; plot(ts1, mean(mean_1_trials))

%%
%%
% Initialize arrays to store the data

shockResponses_all_mice = [];
trialChoices_all_mice = [];
animal_IDs_all_mice = [];  

% Iterate through each level of meanZallMouse
for kk = 1:length(meanZallMouse) %1:length(meanZallMouse)
    shockResponses = [];
    trialChoices = [];
    animal_IDs = [];  
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{kk};
    
    % comment out if necessary
    % limit to just SHOCK activated neurons
    % get this by loading the 10x variable dataset, or similar, where array
    % 4 is data from SHK == 1
    % meanNestedCellArray = meanNestedCellArray(respClass_all_array_mouse{kk, 4} == 1);
    
    % Extract the data for easy reference
    data = behav_data_extracted_array{1, kk};


    % Get the trial choices for the current mouse
    currentTrialChoices = shock_abort_sequence_all_array{kk}';

    % uncomment to only use trials where an abort occurred
    % currentTrialChoices = currentTrialChoices(currentTrialChoices(:,1) ~= 0, 1);
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        
        meanValues = meanNestedCellArray{j};
        
        % uncomment to only use trials where an abort occurred
        % meanValues = meanValues(currentTrialChoices(:,1) ~= 0, 1);
        
        % Here we use the meanValues as is, no averaging across trials
        % Flatten the meanValues to a single row vector, if needed
        shockResponses = [shockResponses; meanValues];
        
        % Append the corresponding trial choice
        trialChoices = [trialChoices; currentTrialChoices];
        animal_IDs = [animal_IDs; repmat(kk, length(meanValues), 1)];
    end
    shockResponses_all_mice = [shockResponses_all_mice; shockResponses];
    trialChoices_all_mice = [trialChoices_all_mice; trialChoices];
    animal_IDs_all_mice = [animal_IDs_all_mice; animal_IDs];
end

%%
% Step 1: Define the predictor and response variables
X = shockResponses_all_mice;    % Predictor (calcium response)
y = trialChoices_all_mice;       % Response (count data for choices)

% Step 2: Fit a Poisson regression model
% Use MATLAB's fitglm function with 'poisson' distribution for count data
mdl = fitglm(X, y, 'Distribution', 'poisson');

% Display model summary
disp(mdl);

% Step 3: Calculate the Exponentiated Coefficients (Incidence Rate Ratios)
% Extract coefficients
coefficients = mdl.Coefficients.Estimate;

% Calculate incidence rate ratios (IRRs) by exponentiating the coefficients
% These are analogous to odds ratios in logistic regression
incidence_rate_ratios = exp(coefficients);

% Display incidence rate ratios
fprintf('Incidence Rate Ratio for Intercept: %.4f\n', incidence_rate_ratios(1));
fprintf('Incidence Rate Ratio for Predictor (shockResponses_all_mice): %.4f\n', incidence_rate_ratios(2));

% Step 4: Make predictions
% Get predicted counts for each trial based on the model
predictedCounts = predict(mdl, X);

% Optional: Assess model accuracy with a measure of goodness-of-fit
% Calculate residuals to evaluate how well the model fits
residuals = y - predictedCounts;
sse = sum(residuals.^2); % Sum of squared errors
fprintf('Sum of Squared Errors (SSE): %.2f\n', sse);



