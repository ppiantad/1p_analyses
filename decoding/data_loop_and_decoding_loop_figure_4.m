% to use this code, use decoding_variables_to_use to load data
% plus need to load relevant 10x dataset, run block_wise_changes_v1 etc
iter = 0;
num_iterations = 5;
caTraceTrials_mouse_iterations = cell(1, num_iterations);
uv.evtWin = [-5 1]; %what time do you want to look at around each event [-2 8] [-10 5]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
ca_data_type = "C_raw"; % C % C_raw
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes

session_to_analyze = 'RDT_D1';

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





use_normalized_time = 0;
% shuffle_confirm = 1; %1 if you want shuffle, 0 if you don't
counts = 0;
num_comparisons = 2;
animalIDs = (fieldnames(final));

%%
for dd = 1:size(neuron_subgroup, 1)

    %%

    % iter = 0;

    ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);

    epoc_to_align = current_epoc{dd};
    to_shuffle = shuffle_confirm(dd);

    caTraceTrials_mouse_iterations = cell(1, num_iterations);





    clear neuron_mean neuron_sem neuron_num zall_mean zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all caTraceTrials_mouse caTraceTrials_current caTraceTrials_mouse_decoding
    %%



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
        if shuffle_confirm(dd) == 0
            for num_comparison = 1:num_comparisons
                if num_comparison == 1 %num_comparison == 3 num_comparison == 1 % if you want to force shuffle, swap to num_comparisons 3 (if doing 2 comparisons) and change the shuffle below. prob should make this a little more intuitive in the future
                    neuron_num = 0;
                    % neuron_sem = zeros(1, size(ts1, 2));
                    for ii = 1:size(fieldnames(final),1)
                        currentanimal = char(animalIDs(ii));
                        if isfield(final.(currentanimal), session_to_analyze)
                            BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
                            [BehavData,trials,varargin]=TrialFilter_test(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 2, 'BLOCK', 3); %'OMITALL', 0, 'BLANK_TOUCH', 0
                            trials = cell2mat(trials);
                            ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);

                            % use if you want to check an event without copious
                            % other filters
                            % ca = ca(respClass_all_array_mouse{ii, 3} == 1, :);

                            % uncomment & edit me if you want to examine a subset
                            % of neurons! make sure to do in both parts of else /
                            % elseif statement!!
                            if strcmp(neuron_subgroup{dd}, 'all')
                                % ca = ca;
                            elseif strcmp(neuron_subgroup{dd}, 'pre-choice')
                                ca = ca(prechoice_remapped_mouse{ii, 1} ~= 1, :);
                            elseif strcmp(neuron_subgroup{dd}, 'post-choice reward')
                                ca = ca(postchoice_reward_block_1_mouse{ii, 1} == 1, :);
                            elseif strcmp(neuron_subgroup{dd}, 'consumption')
                                ca = ca(collect_block_1_mouse{ii, 1} == 1, :);
                            elseif strcmp(neuron_subgroup{dd}, 'true neutral')
                                ca = ca(true_neutral_block_1_mouse{ii, 1} == 1, :);
                            end
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
                elseif num_comparison == 2  %num_comparison == 2 || num_comparison == 1 %use this one on the left if you want to do shuffle vs shuffle
                    neuron_num = 0;
                    for ii = 1:size(fieldnames(final),1)
                        currentanimal = char(animalIDs(ii));
                        if isfield(final.(currentanimal), session_to_analyze)
                            BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
                            [BehavData,trials,varargin]=TrialFilter_test(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 2, 'BLOCK', 3);
                            trials = cell2mat(trials);

                            ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);

                            % use if you want to check an event without copious
                            % other filters
                            % ca = ca(respClass_all_array_mouse{ii, 3} == 1, :);

                            % uncomment & edit me if you want to examine a subset
                            % of neurons! make sure to do in both parts of else /
                            % elseif statement!!
                            if strcmp(neuron_subgroup{dd}, 'all')
                                % ca = ca;
                            elseif strcmp(neuron_subgroup{dd}, 'pre-choice')
                                ca = ca(prechoice_remapped_mouse{ii, 1} ~= 1, :);
                            elseif strcmp(neuron_subgroup{dd}, 'post-choice reward')
                                ca = ca(postchoice_reward_block_1_mouse{ii, 1} == 1, :);
                            elseif strcmp(neuron_subgroup{dd}, 'consumption')
                                ca = ca(collect_block_1_mouse{ii, 1} == 1, :);
                            elseif strcmp(neuron_subgroup{dd}, 'true neutral')
                                ca = ca(true_neutral_block_1_mouse{ii, 1} == 1, :);
                            end
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
        elseif shuffle_confirm(dd) == 1
            for num_comparison = 1:num_comparisons
                if num_comparison == 1 || num_comparison == 2  %num_comparison == 2 || num_comparison == 1 %use this one on the left if you want to do shuffle vs shuffle
                    neuron_num = 0;
                    for ii = 1:size(fieldnames(final),1)
                        currentanimal = char(animalIDs(ii));
                        if isfield(final.(currentanimal), session_to_analyze)
                            BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
                            [BehavData,trials,varargin]=TrialFilter_test(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 2, 'BLOCK', 3);
                            trials = cell2mat(trials);

                            ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);

                            % use if you want to check an event without copious
                            % other filters
                            % ca = ca(respClass_all_array_mouse{ii, 3} == 1, :);

                            % uncomment & edit me if you want to examine a subset
                            % of neurons! make sure to do in both parts of else /
                            % elseif statement!!
                            if strcmp(neuron_subgroup{dd}, 'all')
                                % ca = ca;
                            elseif strcmp(neuron_subgroup{dd}, 'pre-choice')
                                ca = ca(prechoice_remapped_mouse{ii, 1} ~= 1, :);
                            elseif strcmp(neuron_subgroup{dd}, 'post-choice reward')
                                ca = ca(postchoice_reward_block_1_mouse{ii, 1} == 1, :);
                            elseif strcmp(neuron_subgroup{dd}, 'consumption')
                                ca = ca(collect_block_1_mouse{ii, 1} == 1, :);
                            elseif strcmp(neuron_subgroup{dd}, 'true neutral')
                                ca = ca(true_neutral_block_1_mouse{ii, 1} == 1, :);
                            end
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
    end
    data_for_decoding = zall_mouse_iterations;

    %% only run this section if you want to not decode across time!
    % attempting to focus decoding on particular epoch. this section takes the mean activity in a "relevant_period" (pre-choice, etc) for each neuron to be decoded, plus its shuffle.
    % The final data structure is the exact same as
    % "caTraceTrials_mouse_iterations", except each final level contains means,
    % not all data from the time series. one critical thing is that, since the
    % code below expects timeseries (not mean), the "ts1" variable which
    % represents the time window must be set to "1", meaning you are
    % effectively decoding 1 "sample" (the means)

    ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
    relevant_period = ca_period_to_use(dd, :);
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

    % iter = iter+1;
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

                accuracy = zeros(k,1);
                precision = zeros(k,1);
                recall = zeros(k,1);
                f1_score = zeros(k,1);
                
                for i = 1:k
                    xTrain = X(cv.training(i),:);
                    yTrain = y(cv.training(i),:);
                    xTest = X(cv.test(i), :);
                    yTest = y(cv.test(i), :);
                    % model = TreeBagger(numTrees, xTrain, yTrain, 'Method', 'classification');
                    % model = fitglm(xTrain, yTrain, 'Distribution', 'binomial' , 'Link', 'logit');
                    model = fitcsvm(xTrain, yTrain);
                    % model = fitcnb(xTrain, yTrain);
                    yPred = predict(model,xTest);
                    % Compute confusion matrix components
                    TP = sum((yPred == 1) & (yTest == 1));
                    TN = sum((yPred == 0) & (yTest == 0));
                    FP = sum((yPred == 1) & (yTest == 0));
                    FN = sum((yPred == 0) & (yTest == 1));

                    % Calculate metrics
                    accuracy(i) = (TP + TN) / (TP + TN + FP + FN);
                    if (TP + FP) > 0
                        precision(i) = TP / (TP + FP);
                    else
                        precision(i) = NaN;
                    end
                    if (TP + FN) > 0
                        recall(i) = TP / (TP + FN);
                    else
                        recall(i) = NaN;
                    end
                    if (precision(i) + recall(i)) > 0
                        f1_score(i) = 2 * (precision(i) * recall(i)) / (precision(i) + recall(i));
                    else
                        f1_score(i) = NaN;
                    end
                end

                % Store averaged results across folds
                accuracy_by_offset(p) = mean(accuracy, 'omitnan');
                precision_by_offset(p) = mean(precision, 'omitnan');
                recall_by_offset(p) = mean(recall, 'omitnan');
                f1_score_by_offset(p) = mean(f1_score, 'omitnan');
            end

            accuracy_at_loop(:, bb) = accuracy_by_offset;
            precision_at_loop(:, bb) = precision_by_offset;
            recall_at_loop(:, bb) = recall_by_offset;
            f1_score_at_loop(:, bb) = f1_score_by_offset;
        end

        accuracy_per_iteration{iter}(uu) = {accuracy_at_loop};
        precision_per_iteration{iter}(uu) = {precision_at_loop};
        recall_per_iteration{iter}(uu) = {recall_at_loop};
        f1_score_per_iteration{iter}(uu) = {f1_score_at_loop};

        cross_mouse_accuracy_per_iteration{iter}(:, uu) = mean(accuracy_at_loop, 2);
        cross_mouse_precision_per_iteration{iter}(:, uu) = mean(precision_at_loop, 2);
        cross_mouse_recall_per_iteration{iter}(:, uu) = mean(recall_at_loop, 2);
        cross_mouse_f1_per_iteration{iter}(:, uu) = mean(f1_score_at_loop, 2);

        sem_accuracy_per_iteration{iter}(:, uu) = std(accuracy_at_loop,[],2)/sqrt(size(accuracy_at_loop, 2));
        sem_precision_per_iteration{iter}(:, uu) = std(precision_at_loop,[],2)/sqrt(size(precision_at_loop, 2));
        sem_recall_per_iteration{iter}(:, uu) = std(recall_at_loop,[],2)/sqrt(size(recall_at_loop, 2));
        sem_f1_per_iteration{iter}(:, uu) = std(f1_score_at_loop,[],2)/sqrt(size(f1_score_at_loop, 2));

        clear accuracy_at_loop precision_at_loop recall_at_loop f1_score_at_loop
    end

end


%% plot figure - specifics need to be tailored to decoding above
% for zz = 1:size(accuracy_per_iteration, 2)
%     current_iter = accuracy_per_iteration{1, zz}
%     for qq = 1:size(current_iter, 2)
%         sem_data{zz}(:, qq) = std(current_iter{1, qq}, [], 2)/sqrt(size(current_iter{1, qq},2));
%     end
% end
% 
% 
% [concatenatedTable_all] = get_median_choice_and_collect_fn(behav_tbl_iter)
% for zz = 1:size(concatenatedTable_all, 2)
%     temp_table = concatenatedTable_all{zz}; 
%     % median_choice_time_all(zz) = median(temp_table.choiceTime - temp_table.stTime);
%     median_collect_time_all(zz) = median(temp_table.collectionTime - temp_table.choiceTime);
%     median_start_time_all(zz) = median(temp_table.stTime - temp_table.choiceTime);
%     % median_choice_time_block_2_all(zz) = median(temp_table.choiceTime(temp_table.Block == 2) - temp_table.stTime(temp_table.Block == 2));
%     % median_choice_time_block_3_all(zz) = median(temp_table.choiceTime(temp_table.Block == 3) - temp_table.stTime(temp_table.Block == 3));
%     clear temp_table
% end
% 
% 
% 
% 
% figure;
% shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 1}, 2), nanmean(sem_data{1, 1}, 2), 'lineProps', {'color', 'k'});
% hold on;shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 2}, 2), nanmean(sem_data{1, 2}, 2), 'lineProps', {'color', 'g'});
% hold on;shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 3}, 2), nanmean(sem_data{1, 3}, 2), 'lineProps', {'color', 'c'});
% hold on;shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 4}, 2), nanmean(sem_data{1, 4}, 2), 'lineProps', {'color', 'b'});
% hold on;shadedErrorBar(ts1, nanmean(cross_mouse_accuracy_per_iteration{1, 5}, 2), nanmean(sem_data{1, 5}, 2), 'lineProps', {'color', [.7 .7 .7]});
% 
% 
% xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_all, 'r', {'Median', 'collect', 'latency'})
% xlabel('Time from Large Rew Choice (s)');
% legend({'all neurons', 'pre-choice active', 'post-choice reward active', 'consumption active', 'shuffle (consumption shuff vs. consumption shuff)'}, 'Location','northwest')

%%
% Initialize variables to store means and standard deviations for each group
group_means = zeros(4, 2); % 5 groups, 4 bars per group
group_stddevs = zeros(4, 2);

% Calculate means and standard deviations for each group
for i = 1:4
    group_data = cell2mat(cross_mouse_accuracy_per_iteration(:, (i-1)*2 + 1:i*2));
    for j = 1:2
        set_data = group_data(:, (j-1)*5 + 1:j*5);
        group_means(i, j) = mean(set_data, 'all');
        group_stddevs(i, j) = std(set_data, 0, 'all');
    end
end

% Define proper Y positions for horizontal bars
num_groups = size(group_means, 1); % Number of groups (5 in this case)
num_items_per_group = size(group_means, 2); % Number of items per group (4 in this case)
group_spacing = 1.5; % Space between groups
item_spacing = 0.3; % Space between items within a group

% Compute Y positions for all items in all groups
barPositions = zeros(num_items_per_group, num_groups);
for group_idx = 1:num_groups
    start_pos = (group_idx - 1) * group_spacing;
    for item_idx = 1:num_items_per_group
        barPositions(item_idx, group_idx) = start_pos + (item_idx - 1) * item_spacing;
    end
end

% Reverse the order of the groups in barPositions
barPositions = fliplr(barPositions); % Reverse column order for groups
barPositions = flip(barPositions);
% group_means = flip(group_means, 2);


% Reverse group labels for proper Y-axis
group_labels = {'All', 'Pre-choice', 'Post-choice reward', 'Consum.', 'Non-responsive'};

% Plot horizontal bars using errorbar
figure;
hold on;
for i = 1:num_items_per_group
    % Horizontal error bars
    errorbar(group_means(:, i), barPositions(i, :), group_stddevs(:, i), 'horizontal', '.', 'LineWidth', 1, 'DisplayName', sprintf('Epoch %d', i));
end
hold off;

% Customize the plot
% xlabel('Mean Accuracy');
% ylabel('Group');
% title('Mean Accuracy by Group');

% Compute Y-tick positions in ascending order for the reversed group order
yticks(mean(flip(barPositions, 2), 1)); % Flip the columns back to get ascending order for ticks
yticklabels(flip(group_labels)); % Correct order of labels
xlim([0.5 1])
xticks([0.5 0.75 1])
legend('Pre-choice epoch', 'Post-choice epoch', 'Consumption epoch', 'Shuffle (Consumption)', 'Location', 'eastoutside');


%%

% Initialize vectors to store results
num_pairs = 2; % Since you have 8 cells, there are 4 pairs
h_values = zeros(1, num_pairs);
p_values = zeros(1, num_pairs);

for i = 1:num_pairs
    % Get the index for the paired cells
    idx1 = (i - 1) * 2 + 1;
    idx2 = idx1 + 1;
    
    % Convert each cell to a matrix and compute mean accuracy
    data_matrix_1 = cell2mat(accuracy_per_iteration{1, idx1}');
    mean_accuracy_1 = mean(data_matrix_1, 1);
    
    data_matrix_2 = cell2mat(accuracy_per_iteration{1, idx2}');
    mean_accuracy_2 = mean(data_matrix_2, 1);
    
    % Perform paired-samples t-test
    [h, p] = ttest(mean_accuracy_1, mean_accuracy_2);
    
    % Store results
    h_values(i) = h;
    p_values(i) = p;
end

% Display results
disp('Paired t-test results:');
disp(table((1:num_pairs)', h_values', p_values', 'VariableNames', {'Pair', 'h', 'p'}));


% for testing "real" vs "real" - check the indices to make sure it still
% matches
% Convert each cell to a matrix and compute mean accuracy
data_matrix_1 = cell2mat(f1_score_per_iteration{1, 1}');
mean_accuracy_1 = mean(data_matrix_1, 1);

data_matrix_2 = cell2mat(f1_score_per_iteration{1, 3}');
mean_accuracy_2 = mean(data_matrix_2, 1);

% Perform paired-samples t-test
[h, p] = ttest(mean_accuracy_1, mean_accuracy_2);

% Store results
h_values(i) = h;
p_values(i) = p;

%%
% Define figure width and height (adjust as needed)
fig_width = 300; % Width in pixels
fig_height = 600; % Height in pixels

% Create figure with custom size
figure('Position', [100, 100, fig_width, fig_height]); 
hold on;

% Number of pairs
num_pairs = 2; 
mean_values = zeros(num_pairs, 2); % Store means for each pair
sem_values = zeros(num_pairs, 2); % Store SEM for each pair
all_data = cell(num_pairs, 2); % Store individual data for scatter

% Process data for plotting
for i = 1:num_pairs
    idx1 = (i - 1) * 2 + 1;
    idx2 = idx1 + 1;
    
    % Convert cell contents to matrices and compute means
    data_matrix_1 = cell2mat(f1_score_per_iteration{1, idx1}');
    data_matrix_2 = cell2mat(f1_score_per_iteration{1, idx2}');
    
    % Mean values
    mean_values(i, 1) = mean(mean(data_matrix_1, 1));
    mean_values(i, 2) = mean(mean(data_matrix_2, 1));
    
    % Standard Error of the Mean (SEM)
    sem_values(i, 1) = std(mean(data_matrix_1, 1)) / sqrt(size(data_matrix_1, 1));
    sem_values(i, 2) = std(mean(data_matrix_2, 1)) / sqrt(size(data_matrix_2, 1));
    
    % Store data for scatter plot
    all_data{i, 1} = mean(data_matrix_1, 1);
    all_data{i, 2} = mean(data_matrix_2, 1);
end

% Reverse order of bars
mean_values = flipud(mean_values);
sem_values = flipud(sem_values);
all_data = flipud(all_data);

% Create horizontal bar plot
y_positions = 1:num_pairs; % Set positions for the pairs
bar_width = 0.6; % Bar width to allow space for scatter points

% Plot bars
bar_handles = barh(y_positions, mean_values, bar_width, 'grouped');
bar_handles(1).FaceColor = [0.2, 0.4, 0.8]; % Color for first condition
bar_handles(2).FaceColor = [0.8, 0.4, 0.2]; % Color for second condition

% Get the actual bar positions
x_offsets = [bar_handles(1).XEndPoints; bar_handles(2).XEndPoints];

% Add error bars centered on bars
errorbar(mean_values(:,1), x_offsets(1,:), sem_values(:,1), 'horizontal', 'k', 'linestyle', 'none', 'linewidth', 1.5);
errorbar(mean_values(:,2), x_offsets(2,:), sem_values(:,2), 'horizontal', 'k', 'linestyle', 'none', 'linewidth', 1.5);

% Add scatter points & connect them
for i = 1:num_pairs
    % Get scatter y-positions directly over bars
    y1 = x_offsets(1, i) * ones(size(all_data{i, 1}));
    y2 = x_offsets(2, i) * ones(size(all_data{i, 2}));
    
    % Scatter individual points
    scatter(all_data{i, 1}, y1, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(all_data{i, 2}, y2, 50, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Connect points with lines
    plot([all_data{i, 1}; all_data{i, 2}], [y1; y2], 'k-', 'LineWidth', 1);
end

% Labels & Aesthetics
yticks(y_positions);
yticklabels({'Non-resp', 'Consum.', 'Post-choiceRew', 'Pre-choice'}); % Reverse order
% xlabel('Mean Accuracy');
% title('Paired Accuracy Comparisons');
% legend({'Condition 1', 'Condition 2'}, 'Location', 'SouthEast');
set(gca, 'FontSize', 12);
% Set X-axis limits and define ticks explicitly
xlim([0.35, 1]); 
xticks(0.4:0.1:1); % Tick marks at each value from 0.4 to 0.9
box off;
hold off;
