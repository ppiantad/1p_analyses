%% first flatten data using flatten_data_for_offset_decoding, ensure that resulting arrays are correct! 

% offset_to_decode = 2; 

% Set the number of folds
numFolds = 5;


for p = 1:size(trimmed_concatenatedColumns_offsets, 1)

    for fold = 1:numFolds

        offset_1_GCAMP = cell2mat(trimmed_concatenatedColumns_offsets(p,:));
        offset_1_events_offset = cell2mat(trimmed_concatenatedEvents_offsets(p,:));
        offset_1_trials_offset = cell2mat(trimmed_concatenatedColumns_trials_offsets(p,:));

        y = offset_1_events_offset;
        X = offset_1_GCAMP;



        [trainIdx, testIdx] = crossvalind('HoldOut', size(X, 1), 0.3);

        XTrain = X(trainIdx, 1:end);
        yTrain = y(trainIdx);
        XTest = X(testIdx, 1:end);
        yTest = y(testIdx);

        % Calculate the mean and standard deviation for each column
        column_means = mean(XTrain);
        column_sds = std(XTrain);

        % Initialize a matrix to store the z-scored values
        zscored_XTrain = zeros(size(XTrain));

        % Loop through each column
        for col = 1:size(X, 2)
            % Calculate z-score for each element in the column
            zscored_XTrain(:, col) = (XTrain(:, col) - column_means(col)) / column_sds(col);
        end

        % X_normalized = zscore(X, 0, 1);
        %
        % X_column_means(1, :) = mean(X, 1);
        % X_column_SD(1,:) = std(X, 1);


        % Train a Random Forest classification model
        numTrees = 100; % Number of decision trees in the forest
        model = TreeBagger(numTrees, zscored_XTrain, yTrain, 'Method', 'classification');

        %z-score test data based on XTrain values
        
        for col = 1:size(XTest, 2)
            zscored_XTest(:,col) = (XTest(:,col) - column_means(col)) / column_sds(col);
        end

        % Predict trial block labels using the trained model
        predictions = predict(model, zscored_XTest);

        % Convert cell array of strings to numeric labels
        predictions = str2double(predictions);

        % Evaluate the model's performance
        correctPredictions = sum(predictions == yTest);
        totalSamples = numel(yTest);
        accuracy(p, fold) = correctPredictions / totalSamples;
        predictions_array(:,p) = predictions;
    end
    
    final_accuracy(p) = mean(accuracy(p,: ));
    

end

%% Attempt @ shuffling data
numFolds = 5;


for p = 1:size(trimmed_concatenatedColumns_offsets, 1)

    for fold = 1:numFolds

        offset_1_GCAMP = cell2mat(trimmed_concatenatedColumns_offsets(p,:));
        offset_1_events_offset = cell2mat(trimmed_concatenatedEvents_offsets(p,:));
        offset_1_trials_offset = cell2mat(trimmed_concatenatedColumns_trials_offsets(p,:));

        y = offset_1_events_offset;
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

        [trainIdx, testIdx] = crossvalind('HoldOut', size(X, 1), 0.3);

        XTrain = X(trainIdx, 1:end);
        yTrain = y(trainIdx);
        XTest = X(testIdx, 1:end);
        yTest = y(testIdx);

        % Calculate the mean and standard deviation for each column
        column_means = mean(XTrain);
        column_sds = std(XTrain);

        % Initialize a matrix to store the z-scored values
        zscored_XTrain = zeros(size(XTrain));

        % Loop through each column
        for col = 1:size(X, 2)
            % Calculate z-score for each element in the column
            zscored_XTrain(:, col) = (XTrain(:, col) - column_means(col)) / column_sds(col);
        end

        % X_normalized = zscore(X, 0, 1);
        %
        % X_column_means(1, :) = mean(X, 1);
        % X_column_SD(1,:) = std(X, 1);


        % Train a Random Forest classification model
        numTrees = 100; % Number of decision trees in the forest
        model = TreeBagger(numTrees, zscored_XTrain, yTrain, 'Method', 'classification');

        %z-score test data based on XTrain values
        
        for col = 1:size(XTest, 2)
            zscored_XTest(:,col) = (XTest(:,col) - column_means(col)) / column_sds(col);
        end

        % Predict trial block labels using the trained model
        predictions = predict(model, zscored_XTest);

        % Convert cell array of strings to numeric labels
        predictions = str2double(predictions);

        % Evaluate the model's performance
        correctPredictions = sum(predictions == yTest);
        totalSamples = numel(yTest);
        accuracy(p, fold) = correctPredictions / totalSamples;
        predictions_array(:,p) = predictions;
    end
    
    final_accuracy(p) = mean(accuracy(p,: ));
    

end




%%
confusionMat = confusionmat(yTest, predictions);

% Visualize the confusion matrix
heatmap(confusionMat, 'XLabel', 'Predicted', 'YLabel', 'Actual');

disp(['Accuracy: ' num2str(accuracy)]);
