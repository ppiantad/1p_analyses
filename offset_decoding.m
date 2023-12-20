%% first flatten data using flatten_data_for_offset_decoding, ensure that resulting arrays are correct! 

% offset_to_decode = 2; 

% Set the number of folds
numFolds = 5;


for p = 1:size(trimmed_concatenatedColumns_offsets, 1)

    for fold = 1:numFolds

        offset_1_GCAMP = cell2mat(trimmed_concatenatedColumns_offsets(p,:));
        offset_1_events_offset = cell2mat(trimmed_concatenatedEvents_offsets(p,:));
        offset_1_trials_offset = cell2mat(trimmed_concatenatedColumns_trials_offsets(p,:));

        y = offset_1_events_offset(:, p);
        X = offset_1_GCAMP;

        % Calculate the mean and standard deviation for each column
        column_means = mean(X);
        column_sds = std(X);


        % Initialize a matrix to store the z-scored values
        zscored_X = zeros(size(X));

        % Loop through each column
        for col = 1:size(X, 2)
            % Calculate z-score for each element in the column
            zscored_X(:, col) = (X(:, col) - column_means(col)) / column_sds(col);
        end




        % X_normalized = zscore(X, 0, 1);
        %
        % X_column_means(1, :) = mean(X, 1);
        % X_column_SD(1,:) = std(X, 1);

        

        [trainIdx, testIdx] = crossvalind('HoldOut', size(zscored_X, 1), 0.3);

        XTrain = zscored_X(trainIdx, 1:end);
        yTrain = y(trainIdx);
        XTest = zscored_X(testIdx, 1:end);
        yTest = y(testIdx);


        % Train a Random Forest classification model
        numTrees = 100; % Number of decision trees in the forest
        model = TreeBagger(numTrees, XTrain, yTrain, 'Method', 'classification');

        % model = fitcsvm(XTrain, yTrain);

        % Predict trial block labels using the trained model
        predictions = predict(model, XTest);

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
