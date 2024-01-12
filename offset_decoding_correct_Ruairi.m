%% first flatten data using flatten_data_for_offset_decoding, ensure that resulting arrays are correct! 

% offset_to_decode = 2; 

% Set the number of folds
numFolds = 1;

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

%% shuffled
% first flatten data using flatten_data_for_offset_decoding, ensure that resulting arrays are correct! 

% offset_to_decode = 2; 

% Set the number of folds
numFolds = 1;

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
    for g = 1:1000 %for each resampling of the data g = 1:uv.resamples
        [num_timepoints, num_signals] = size(X);
        shuffled_data = zeros(num_timepoints, num_signals); % Preallocate matrix for efficiency

        for i = 1:num_timepoints
            shift_val = randi(num_signals); % Generate a random shift value for each signal
            shuffled_data(i,:) = circshift(X(i,:), shift_val,2); % Perform the circular shuffle
        end
        nullDistTrace(g,:) = nanmean(X); %calculate the NaN mean of the shuffled traces
        %calculate the NaN mean of the shuffled event rates
    end
    X = shuffled_data;
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