num_iterations = 1; 
caTraceTrials_mouse_iterations = cell(1, num_iterations);
iter = 0;
uv.evtWin = [-2 8]; %what time do you want to look at around each event [-2 8]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
counts = 0;
num_comparisons = 2; 

ca_data_type = "C_raw"; % C % C_raw
% CNMFe_data.S = inferred spikes
% CNMFe_data.C = denoised traces
% CNMFe_data.C_raw = CNMFe data
%%
for num_iteration = 1:num_iterations
    fprintf('The current iteration is: %d\n', num_iteration);

    session_to_analyze = 'Pre_RDT_RM';
    epoc_to_align = 'collectionTime';
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
        if num_comparison == 1 %num_comparison == 3
            neuron_num = 0;
            % neuron_sem = zeros(1, size(ts1, 2));
            for ii = 1:size(fieldnames(final),1)
                currentanimal = char(animalIDs(ii));
                if isfield(final.(currentanimal), session_to_analyze)
                    BehavData = final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.BehavData;
                    [BehavData,trials,varargin]=TrialFilter(BehavData,'REW',1.2);
                    trials = cell2mat(trials);
                    ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);


                    num_samples = size(ca, 2);
                    sampling_frequency = (final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.dt)*100;
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
                                idx = time_array >= min(timeWin) & time_array< max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
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
        elseif num_comparison == 2 %num_comparison == 2 %num_comparison == 2 || num_comparison == 1
            neuron_num = 0;
            for ii = 1:size(fieldnames(final),1)
                currentanimal = char(animalIDs(ii));
                if isfield(final.(currentanimal), session_to_analyze)
                    BehavData = final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.BehavData;
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

                    num_samples = size(ca, 2);
                    sampling_frequency = (final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.dt)*100;
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
                                idx = time_array >= min(timeWin) & time_array< max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
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
select_mouse = 'BLA_Insc_28';

select_mouse_index = find(strcmp(animalIDs, select_mouse));




trials_per_mouse_decoding = trials_per_mouse(select_mouse_index,:);

%% old - for one mouse. likely to delete once decoding_all_beta is finished

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