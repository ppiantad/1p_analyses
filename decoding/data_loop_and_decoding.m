iter = 0;


%%
num_iterations = 1; 
caTraceTrials_mouse_iterations = cell(1, num_iterations);
% iter = 0;
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
counts = 0;
num_comparisons = 2; 
epoc_to_align = 'choiceTime';


ca_data_type = "C_raw"; % C % C_raw
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes

use_normalized_time = 0;
shuffle_confirm = 1; %1 if you want shuffle, 0 if you don't

session_to_analyze = 'Pre_RDT_RM';

if strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
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


clear neuron_mean neuron_sem neuron_num zall_mean zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all caTraceTrials_mouse caTraceTrials_current caTraceTrials_mouse_decoding
%%

animalIDs = (fieldnames(final));

iter = iter + 1;

for num_iteration = 1:num_iterations
    fprintf('The current iteration is: %d\n', num_iteration);


    event_to_analyze = {'BLOCK',1,'REW',1.2};

    if exist('iter', 'var') == 1

    elseif exist('iter', 'var') == 0
        iter = 0;
    end
    % disp(['iter: ' string(iter)]);
    clear neuron_mean neuron_sem neuron_num zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized

    


    neuron_num = 0;
    

    sum_trials_per_iter = 0;
    filter_names_idx = cellfun(@ischar,event_to_analyze);
    filter_strings = string(event_to_analyze(filter_names_idx));
    for num_comparison = 1:num_comparisons
        if num_comparison == 1 %num_comparison == 3 num_comparison == 1 % if you want to force shuffle, swap to num_comparisons 3 (if doing 2 comparisons) and change the shuffle below. prob should make this a little more intuitive in the future
            neuron_num = 0;
            % neuron_sem = zeros(1, size(ts1, 2));
            for ii = 1:size(fieldnames(final),1)
                currentanimal = char(animalIDs(ii));
                if isfield(final.(currentanimal), session_to_analyze)
                    BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
                    [BehavData,trials,varargin]=TrialFilter_test(BehavData,'REW', 1.2); %'OMITALL', 0, 'BLANK_TOUCH', 0
                    trials = cell2mat(trials);
                    ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);
                    
                    % use if you want to check an event without copious
                    % other filters
                    % ca = ca(respClass_all_array_mouse{ii, 3} == 1, :);
                    
                    % uncomment & edit me if you want to examine a subset
                    % of neurons! make sure to do in both parts of else /
                    % elseif statement!!
                    % ca = ca(respClass_all_array_mouse_pre_choice_active{ii, 1} == 1, :);
                    % ca = ca(respClass_all_array_mouse_post_choice_reward{ii, 1} == 1, :);
                    % ca = ca(respClass_all_array_mouse_consumption{ii, 1} == 1, :);

                    % uncomment below if you want to examine a subset of
                    % neurons that were not responsive to any of the
                    % Pre_RDT_RM events
                    % ca = ca(respClass_all_array_mouse_true_neutral{ii, 1} == 1, :);

                    num_samples = size(ca, 2);
                    sampling_frequency = (final.(currentanimal).(session_to_analyze).uv.dt)*100;
                    time_array = final.(currentanimal).(session_to_analyze).time;
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
                        [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);

                        caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                        zall = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                        zall_array{neuron_num} = zall;
                        zall_mouse{ii, num_comparison}(u) = {zall};
                        caTraceTrials_mouse{ii, num_comparison}(u) = {caTraceTrials};
                        zall_mean(neuron_num,:) = mean(zall);
                        trials_per_mouse{ii, num_comparison} = trials;
                        clear zall caTraceTrials zb zsd;
                    end

                end
            end


            % disp(['iter = ' string(iter)])
        elseif num_comparison == 2 %num_comparison == 2 || num_comparison == 1 %use this one on the left if you want to do shuffle vs shuffle
            neuron_num = 0;
            for ii = 1:size(fieldnames(final),1)
                currentanimal = char(animalIDs(ii));
                if isfield(final.(currentanimal), session_to_analyze)
                    BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
                    [BehavData,trials,varargin]=TrialFilter_test(BehavData,'REW', 1.2);
                    trials = cell2mat(trials);

                    ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);
                    
                    % use if you want to check an event without copious
                    % other filters
                    % ca = ca(respClass_all_array_mouse{ii, 3} == 1, :);                   
                   
                    % uncomment & edit me if you want to examine a subset
                    % of neurons! make sure to do in both parts of else /
                    % elseif statement!!
                    % ca = ca(respClass_all_array_mouse_pre_choice_active{ii, 1} == 1, :);
                    % ca = ca(respClass_all_array_mouse_post_choice_reward{ii, 1} == 1, :);
                    % ca = ca(respClass_all_array_mouse_consumption{ii, 1} == 1, :);

                    % uncomment below if you want to examine a subset of
                    % neurons that were not responsive to any of the
                    % Pre_RDT_RM events
                    % ca = ca(respClass_all_array_mouse_true_neutral{ii, 1} == 1, :);

                    % shuffle data for comparison
                    [num_cells, num_samples] = size(ca);
                    shuffled_data = zeros(num_cells, num_samples); % Preallocate matrix for efficiency
                    shift_val = randi(num_samples); % Generate a random shift value for each signal RUAIRI RECOMMENDED KEEPING THE SAME SHIFT VAL, rather than randomizing per neuron. this is because then you keep the overall correlation b/w the neurons, but disrupt the relationship to the event timestamps
                    for i = 1:num_cells
                        % shift_val = randi(num_signals); % Generate a random shift value for each signal
                        shuffled_data(i,:) = circshift(ca(i,:), shift_val,2); % Perform the circular shuffle
                    end
                    ca = shuffled_data;

                    num_samples = size(ca, 2);
                    sampling_frequency = (final.(currentanimal).(session_to_analyze).uv.dt)*100;
                    time_array = final.(currentanimal).(session_to_analyze).time;
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
                        [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time);

                        caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                        zall = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
                        zall_array{neuron_num} = zall;
                        zall_mouse{ii, num_comparison}(u) = {zall};
                        caTraceTrials_mouse{ii, num_comparison}(u) = {caTraceTrials};
                        zall_mean(neuron_num,:) = mean(zall);
                        trials_per_mouse{ii, num_comparison} = trials;
                        clear zall caTraceTrials zb zsd;

                    end
                end
            
            clear zall caTraceTrials zb zsd;
            end
        end
    end

    caTraceTrials_mouse_iterations(1, num_iteration) = {caTraceTrials_mouse};
    zall_mouse_iterations(1, num_iteration) = {zall_mouse};
end

data_for_decoding = caTraceTrials_mouse_iterations;

%% only run this section if you DO NOT WANT TO decode across time!  
% attempting to focus decoding on particular epoch. this section takes the mean activity in a "relevant_period" (pre-choice, etc) for each neuron to be decoded, plus its shuffle. 
% The final data structure is the exact same as
% "caTraceTrials_mouse_iterations", except each final level contains means,
% not all data from the time series. one critical thing is that, since the
% code below expects timeseries (not mean), the "ts1" variable which
% represents the time window must be set to "1", meaning you are
% effectively decoding 1 "sample" (the means)

ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
relevant_period = [-4 0]
sub_window_idx = ts1 >= relevant_period(1) & ts1 <= relevant_period(2);

% Preallocate memory for caTraceTrials_mouse_iterations_means
caTraceTrials_mouse_iterations_means = cell(1, size(caTraceTrials_mouse_iterations, 2));
for ff = 1:size(caTraceTrials_mouse_iterations, 2)
    for gg = 1:size(caTraceTrials_mouse_iterations{1, ff}, 1)
        current_level = caTraceTrials_mouse_iterations{1, ff}(gg,:);
        for hh = 1:size(current_level, 2)

            % Check if current_level{1, hh} contains a cell array
            if iscell(current_level{1, hh})
                for jj = 1:size(current_level{1, hh}, 2)
                    cell_level = current_level{1, hh}{jj};

                    cell_mean = mean(cell_level(:, sub_window_idx), 2);

                    % Initialize if necessary
                    if isempty(caTraceTrials_mouse_iterations_means{1, 1})
                        caTraceTrials_mouse_iterations_means{1, 1} = cell(size(caTraceTrials_mouse_iterations{1, ff}));
                    end

                    % Initialize if necessary at the gg, hh level
                    if isempty(caTraceTrials_mouse_iterations_means{1, 1}{gg, hh})
                        caTraceTrials_mouse_iterations_means{1, 1}{gg, hh} = cell(size(current_level{1, hh}));
                    end

                    caTraceTrials_mouse_iterations_means{1, ff}{gg, hh}{jj} = cell_mean;
                end
            end
        end
    end
end

iter = iter+1;
ts1 = 1;

data_for_decoding = caTraceTrials_mouse_iterations_means;



%%

caTraceTrials_current = []
empty_rows_indices = []


for uu = 1:size(data_for_decoding, 2)
    caTraceTrials_current = data_for_decoding{:,uu};
    % caTraceTrials_current(all(cellfun(@isempty,caTraceTrials_current),2), : ) = [];
    % check if there are mice with no data - delete these & the
    % corresponding related arrays
    empty_rows_indices = find(cellfun(@isempty, caTraceTrials_current(:,1)));
    caTraceTrials_current(empty_rows_indices, :) = [];
    % animalIDs(empty_rows_indices,:) = [];
    % trials_per_mouse(empty_rows_indices, :) = [];
    % zall_mouse(empty_rows_indices, :) = [];
    fprintf('The current iteration is: %d\n', uu);
    for bb = 1:size(caTraceTrials_current, 1)
        caTraceTrials_mouse_decoding = caTraceTrials_current(bb,:);
        % disp('The current mouse being decoded is:' + string(animalIDs(bb)))
        fprintf('The current mouse being decoded is: %d\n', bb);
        [trimmed_concatenatedColumns_offsets,...
            trimmed_concatenatedColumns_time_offsets,...
            trimmed_concatenatedColumns_trials_offsets,...
            trimmed_concatenatedEvents_offsets]...
            = flatten_data_for_offset_decoding_fn(caTraceTrials_mouse_decoding, ts1, num_comparisons);

        % additions from Ruairi 01/12/2024
        k = 5; %number of cross-validation folds 10
        accuracy_by_offset = zeros(size(trimmed_concatenatedColumns_offsets, 1), 1);
        numTrees = 100; % Number of decision trees in the forest
        for p = 1:size(trimmed_concatenatedColumns_offsets, 1)
            offset_1_GCAMP = cell2mat(trimmed_concatenatedColumns_offsets(p,:));
            offset_1_events_offset = cell2mat(trimmed_concatenatedEvents_offsets(p,:));
            offset_1_trials_offset = cell2mat(trimmed_concatenatedColumns_trials_offsets(p,:));
            offset_1_time_offset = cell2mat(trimmed_concatenatedColumns_time_offsets(p,:));

            y = offset_1_events_offset(:,1);
            y = y -1; %values need to be 0 or 1 for fitcnb
            X = offset_1_GCAMP;
            X = zscore(X);
            idx = randperm(size(X, 1));
            X = X(idx, :);
            y = y(idx, :);
            cv = cvpartition(size(X, 1),"KFold", k);


            for i = 1:k
                xTrain = X(cv.training(i),:);
                yTrain = y(cv.training(i),:);
                xTest = X(cv.test(i), :);
                yTest = y(cv.test(i), :);
                % model = TreeBagger(numTrees, xTrain, yTrain, 'Method', 'classification');
                % model = fitglm(xTrain, yTrain, 'Distribution', 'binomial' , 'Link', 'logit');
                % model = fitcsvm(xTrain, yTrain);
                model = fitcnb(xTrain, yTrain);
                yPred = predict(model,xTest);
                accuracy(i) = sum(yPred == yTest)/numel(yTest);
            end
            accuracy_by_offset(p) = mean(accuracy);
        end

        % figure; plot(ts1, accuracy_by_offset);
        accuracy_at_loop(:, bb) = accuracy_by_offset;
    end
    accuracy_per_iteration{iter}(uu) = {accuracy_at_loop};
    cross_mouse_accuracy_per_iteration{iter}(:, uu) = mean(accuracy_at_loop, 2);
    sem_accuracy_per_iteration{iter}(:, uu) = std(accuracy_at_loop,[],2)/sqrt(size(accuracy_at_loop, 2));

    clear accuracy_at_loop
end




%% plot figure - specifics need to be tailored to decoding above
for zz = 1:size(accuracy_per_iteration, 2)
    current_iter = accuracy_per_iteration{1, zz}
    for qq = 1:size(current_iter, 2)
        sem_data{zz}(:, qq) = std(current_iter{1, qq}, [], 2)/sqrt(size(current_iter{1, qq},2));
    end
end


[concatenatedTable_all] = get_median_choice_and_collect_fn(behav_tbl_iter)
for zz = 1:size(concatenatedTable_all, 2)
    temp_table = concatenatedTable_all{zz}; 
    % median_choice_time_all(zz) = median(temp_table.choiceTime - temp_table.stTime);
    median_collect_time_all(zz) = median(temp_table.collectionTime - temp_table.choiceTime);
    median_start_time_all(zz) = median(temp_table.stTime - temp_table.choiceTime);
    % median_choice_time_block_2_all(zz) = median(temp_table.choiceTime(temp_table.Block == 2) - temp_table.stTime(temp_table.Block == 2));
    % median_choice_time_block_3_all(zz) = median(temp_table.choiceTime(temp_table.Block == 3) - temp_table.stTime(temp_table.Block == 3));
    clear temp_table
end




figure;
shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 1}, 2), nanmean(sem_data{1, 1}, 2), 'lineProps', {'color', 'k'});
hold on;shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 2}, 2), nanmean(sem_data{1, 2}, 2), 'lineProps', {'color', 'g'});
hold on;shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 3}, 2), nanmean(sem_data{1, 3}, 2), 'lineProps', {'color', 'c'});
hold on;shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 4}, 2), nanmean(sem_data{1, 4}, 2), 'lineProps', {'color', 'b'});
hold on;shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 5}, 2), nanmean(sem_data{1, 5}, 2), 'lineProps', {'color', [.7 .7 .7]});


xline(0);
xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'all neurons', 'pre-choice active', 'post-choice reward active', 'consumption active', 'shuffle (consumption shuff vs. consumption shuff)'}, 'Location','northwest')