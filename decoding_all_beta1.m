


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
        k = 10;
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