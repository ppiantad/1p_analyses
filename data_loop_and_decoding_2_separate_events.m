num_iterations = 1; 
caTraceTrials_mouse_iterations = cell(1, num_iterations);
iter = 0;
uv.evtWin = [-10 5]; %what time do you want to look at around each event [-2 8] [-10 5]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
counts = 0;
num_comparisons = 2; 

ca_data_type = "C_raw"; % C % C_raw
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes


shuffle_confirm = 1; %1 if you want shuffle, 0 if you don't

clear neuron_mean neuron_sem neuron_num zall_mean zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all caTraceTrials_mouse caTraceTrials_current caTraceTrials_mouse_decoding
%%
for num_iteration = 1:num_iterations
    fprintf('The current iteration is: %d\n', num_iteration);

    session_to_analyze = 'Pre_RDT_RM';
    epoc_to_align = 'choiceTime';
    event_to_analyze = {'BLOCK',1,'REW',1.2};

    if exist('iter', 'var') == 1

    elseif exist('iter', 'var') == 0
        iter = 0;
    end
    % disp(['iter: ' string(iter)]);
    clear neuron_mean neuron_sem neuron_num zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized

    

    animalIDs = (fieldnames(final));
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
                    [BehavData,trials,varargin]=TrialFilter(BehavData,'REW', 1.2, 'BLOCK', 1);
                    trials = cell2mat(trials);
                    ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);


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
                        [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u);

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

            iter = iter+1;
            % disp(['iter = ' string(iter)])
        elseif num_comparison == 2  %num_comparison == 2 %num_comparison == 2 || num_comparison == 1
            neuron_num = 0;
            for ii = 1:size(fieldnames(final),1)
                currentanimal = char(animalIDs(ii));
                if isfield(final.(currentanimal), session_to_analyze)
                    BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
                    [BehavData,trials,varargin]=TrialFilter(BehavData,'REW', 1.2, 'BLOCK', 2);
                    trials = cell2mat(trials);

                    ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);
                    % % shuffle data for comparison
                    % [num_cells, num_samples] = size(ca);
                    % shuffled_data = zeros(num_cells, num_samples); % Preallocate matrix for efficiency
                    % shift_val = randi(num_samples); % Generate a random shift value for each signal RUAIRI RECOMMENDED KEEPING THE SAME SHIFT VAL, rather than randomizing per neuron. this is because then you keep the overall correlation b/w the neurons, but disrupt the relationship to the event timestamps
                    % for i = 1:num_cells
                    %     % shift_val = randi(num_signals); % Generate a random shift value for each signal
                    %     shuffled_data(i,:) = circshift(ca(i,:), shift_val,2); % Perform the circular shuffle
                    % end
                    % ca = shuffled_data;

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
                        [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u);

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
    clear iter
    caTraceTrials_mouse_iterations(1, num_iteration) = {caTraceTrials_mouse};
    zall_mouse_iterations(1, num_iteration) = {zall_mouse};

end



%%

for uu = 1:size(caTraceTrials_mouse_iterations, 2)
    caTraceTrials_current = caTraceTrials_mouse_iterations{:,uu};
    % caTraceTrials_current(all(cellfun(@isempty,caTraceTrials_current),2), : ) = [];
    % check if there are mice with no data - delete these & the
    % corresponding related arrays
    empty_rows_indices = find(cellfun(@isempty, caTraceTrials_current(:,1)));
    caTraceTrials_current(empty_rows_indices, :) = [];
    % animalIDs(empty_rows_indices,:) = [];
    % trials_per_mouse(empty_rows_indices, :) = [];
    zall_mouse(empty_rows_indices, :) = [];
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
        k = 5; %set number of cross validation folks. 10
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
    accuracy_per_iteration(uu) = {accuracy_at_loop};
    cross_mouse_accuracy_per_iteration(:, uu) = mean(accuracy_at_loop, 2);
    clear accuracy_at_loop
end




%% old - for one mouse. likely to delete once decoding_all_beta is finished
select_mouse = 'BLA_Insc_28';

select_mouse_index = find(strcmp(animalIDs, select_mouse));




trials_per_mouse_decoding = trials_per_mouse(select_mouse_index,:);



for uu = 1:size(caTraceTrials_mouse_iterations, 2)
    caTraceTrials_mouse_current = caTraceTrials_mouse_iterations{:,uu};
    caTraceTrials_mouse_decoding = caTraceTrials_mouse_current(select_mouse_index,:);
    
    [trimmed_concatenatedColumns_offsets,...
        trimmed_concatenatedColumns_time_offsets,...
        trimmed_concatenatedColumns_trials_offsets,...
        trimmed_concatenatedEvents_offsets]...
        = flatten_data_for_offset_decoding_fn(caTraceTrials_mouse_decoding, ts1, num_comparisons);

    % additions from Ruairi 01/12/2024
    k = 10;
    accuracy_by_offset = zeros(size(trimmed_concatenatedColumns_offsets, 1), 1);
    numTrees = 100; % Number of decision trees in the forest
    for p = 1:size(trimmed_concatenatedColumns_offsets, 1)
        offset_1_GCAMP = cell2mat(trimmed_concatenatedColumns_offsets(p,:));
        offset_1_events_offset = cell2mat(trimmed_concatenatedEvents_offsets(p,:));
        offset_1_trials_offset = cell2mat(trimmed_concatenatedColumns_trials_offsets(p,:));
        offset_1_time_offset = cell2mat(trimmed_concatenatedColumns_time_offsets(p,:));

        y = offset_1_events_offset(:,1);
        y = y -1;
        X = offset_1_GCAMP;
        X = zscore(X);
        idx = randperm(size(X, 1));
        X = X(idx, :);
        y = y(idx,:);
        cv = cvpartition(size(X, 1),"KFold", k);


        for i = 1:k
            xTrain = X(cv.training(i),:);
            yTrain = y(cv.training(i),:);
            xTest = X(cv.test(i), :);
            yTest = y(cv.test(i), :);
            % model = TreeBagger(numTrees, xTrain, yTrain, 'Method', 'classification');
            % model = fitglm(xTrain, yTrain, 'Distribution', 'binomial' , 'Link', 'logit');
            model = fitcnb(xTrain, yTrain);
            yPred = predict(model,xTest);
            accuracy(i) = sum(yPred == yTest)/numel(yTest);
        end
        accuracy_by_offset(p) = mean(accuracy);
    end

    figure; plot(ts1, accuracy_by_offset);
    accuracy_at_loop(:, uu) = accuracy_by_offset;

end