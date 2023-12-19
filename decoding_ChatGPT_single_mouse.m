% Initialize a cell array to store the concatenated data
concatenated_mean_data = cell(12, 1);
% FINISH EDITING THIS TO ALSO GET THE SEM SO THAT YOU CAN PLOT INDIVIDUAL
% MOUSE DATA MORE EASILY
concatenated_sem_data = cell(12, 1);

% Loop through each row of the neuron_mean_mouse cell array
for i = 1:size(neuron_mean_mouse, 1)
    % Extract the data from the current row for columns {1,1}, {1,2}, and {1,3}
    current_row_data = neuron_mean_mouse{i, 1};
    for q = 1:iter-1
        % current_row_data = neuron_mean_mouse{i, q};
        current_row_data = vertcat(current_row_data, neuron_mean_mouse{i, q+1});
        

        % Store the concatenated data in the result cell array
        concatenated_data{i} = current_row_data;
        neurons_per_mouse(i) = size(neuron_mean_mouse{i, 1}, 1);
    end
end
%%
% select the # animal to decode
animalNum_to_decode = 1;

%get the number of neurons that that mouse has (for calculations later)
numNeurons_mouse = neurons_per_mouse(animalNum_to_decode);

neuron_mean_concat = concatenated_data{animalNum_to_decode, 1};




% Create label vector y (corresponding to trial blocks)
numNeuronsPerCondition = numNeurons_mouse;
% change depending on the number of behaviors to decode!
numConditions = size(neuron_mean_concat, 1)/numNeuronsPerCondition;


%% FINISH CODE ABOVE SO THAT YOU CAN USE SINGLE MOUSE DATA, NEED TO FINISH THE SEM CODE ABOVE

load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
colormap(batlowW); % c1 = colorbar; 
% Create label vector y (corresponding to trial blocks)



y = repelem(1:numConditions, numNeuronsPerCondition);



figure;
hold on; 
ylim([-0.8 0.8])
xticks([-10 -5 0 5 10])
shadedErrorBar(ts1, mean(neuron_mean_concat(:,y == 1)), mean(neuron_sem_concat(:,y == 1)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(neuron_mean_concat(:,y == 2)), mean(neuron_sem_concat(:,y == 2)), 'lineProps', {'color', batlowW(1,:)});
shadedErrorBar(ts1, mean(neuron_mean_concat(:,y == 3)), mean(neuron_sem_concat(:,y == 3)), 'lineProps', {'color', batlowW(1,:)});
hold off

%% SVM
% code created with the aid of ChatGPT on 8/28/2023

neuron_mean_concat = concatenated_data{animalNum_to_decode, 1};




% Create label vector y (corresponding to trial blocks)
numNeuronsPerCondition = numNeurons_mouse;
% change depending on the number of behaviors to decode!
numConditions = size(neuron_mean_concat, 1)/numNeuronsPerCondition;




numSamples = size(neuron_mean_concat, 2);
y = repelem(1:numConditions, numNeuronsPerCondition);

% Split data into training and testing sets
[trainIdx, testIdx] = crossvalind('HoldOut', size(neuron_mean_concat, 1), 0.2);
XTrain = neuron_mean_concat(trainIdx, :);
yTrain = y(trainIdx);
XTest = neuron_mean_concat(testIdx, :);
yTest = y(testIdx);

% % Train a classification model (SVM in this example)
% model = fitcsvm(XTrain, yTrain);

% Train a multiclass classification model using ECOC
model = fitcecoc(XTrain, yTrain);

% Predict trial block labels using the trained model
predictions = predict(model, XTest);

% Evaluate the model's performance
correctPredictions = sum(predictions == yTest');
totalSamples = numel(yTest);
accuracy = correctPredictions / totalSamples;
confusionMat = confusionmat(yTest, predictions);

% Visualize the confusion matrix
heatmap(confusionMat, 'XLabel', 'Predicted', 'YLabel', 'Actual');

disp(['Accuracy: ' num2str(accuracy)]);


% Iterate through time samples and compute accuracy
for t = 1:numSamples
    XTestSample = XTest(:, t);
    predictions = predict(model, XTestSample);
    accuracy = sum(predictions == yTest) / numel(yTest);
    accuracies(t) = accuracy;
end

% Plot the accuracies over time
figure;
plot(1:numSamples, accuracies);
xlabel('Time Sample');
ylabel('Accuracy');
title('Accuracy Over Time');






figure; 
hold on;
for ii = 1:numConditions
    plot(ts1, mean(neuron_mean_concat(y==ii, :)));
end




%% Random Forest (this seemed to work best to decode Blocks of Large Rew Trials (B1 Large Choice Aligned, B2 Large Choice Aligned, B3 Large Choice Aligned)
% Assuming you have already loaded and preprocessed your data
% neuron_mean_concats is a 300x200 array with 300 rows (neurons) and 200 columns (samples)







y = repelem(1:numConditions, numNeuronsPerCondition);

% Split data into training and testing sets
[trainIdx, testIdx] = crossvalind('HoldOut', size(neuron_mean_concat, 1), 0.2);
XTrain = neuron_mean_concat(trainIdx, 1:end);
yTrain = y(trainIdx);
XTest = neuron_mean_concat(testIdx, 1:end);
yTest = y(testIdx);

% Train a Random Forest classification model
numTrees = 100; % Number of decision trees in the forest
model = TreeBagger(numTrees, XTrain, yTrain, 'Method', 'classification');

% Predict trial block labels using the trained model
predictions = predict(model, XTest);

% Convert cell array of strings to numeric labels
predictions = str2double(predictions);

% Evaluate the model's performance
correctPredictions = sum(predictions == yTest');
totalSamples = numel(yTest);
accuracy = correctPredictions / totalSamples;
confusionMat = confusionmat(yTest, predictions);

% Visualize the confusion matrix
heatmap(confusionMat, 'XLabel', 'Predicted', 'YLabel', 'Actual');

disp(['Accuracy: ' num2str(accuracy)]);


figure; 
hold on;
for ii = 1:numConditions
    plot(ts1, mean(neuron_mean_concat(y==ii, 1:end)));
end
hold off;




%% Random Forest data organized with separate conditions concatenated in columns, rows are neurons (this seemed to work best to decode Blocks of Large Rew Trials (B1 Large Choice Aligned, B2 Large Choice Aligned, B3 Large Choice Aligned)
% Assuming you have already loaded and preprocessed your data
% neuron_mean_concats is a 300x200 array with 300 rows (neurons) and 200 columns (samples)




% Create label vector y (corresponding to trial blocks)
numNeuronsPerCondition = neuron_num;
% change depending on the number of behaviors to decode!
numConditions = size(neuron_mean_concat, 1)/numMeasurements;




y = repelem(1:numConditions, numMeasurements);

% Split data into training and testing sets
[trainIdx, testIdx] = crossvalind('HoldOut', size(neuron_mean_concat, 2), 0.2);
XTrain = neuron_mean_concat(1:end, trainIdx)';
yTrain = y(trainIdx)';
XTest = neuron_mean_concat(1:end, testIdx)';
yTest = y(testIdx)';

% Train a Random Forest classification model
numTrees = 100; % Number of decision trees in the forest
model = TreeBagger(numTrees, XTrain, yTrain, 'Method', 'classification');

% Predict trial block labels using the trained model
predictions = predict(model, XTest);

% Convert cell array of strings to numeric labels
predictions = str2double(predictions);

% Evaluate the model's performance
correctPredictions = sum(predictions == yTest);
totalSamples = numel(yTest);
accuracy = correctPredictions / totalSamples;
confusionMat = confusionmat(yTest, predictions);

% Visualize the confusion matrix
heatmap(confusionMat, 'XLabel', 'Predicted', 'YLabel', 'Actual');

disp(['Accuracy: ' num2str(accuracy)]);


figure; 
hold on;
for ii = 1:numConditions
    plot(ts1, mean(neuron_mean_concat(1:end, y==ii)));
end
hold off;


%% Random Forest decoder accuracy over time - this takes a while to run but seems to work. requires numTrees from previous step 
% Initialize variables to store accuracy at each time point

numTimePoints = size(XTest, 2);  % Number of time points
accuracyAtTime = zeros(1, numTimePoints);

%Train a Random Forest classification model for each time point
for t = 1:numTimePoints
    XTrainTime = XTrain(:, 1:t);  % Use data up to the current time point
    XTestTime = XTest(:, 1:t);    % Use data up to the current time point for testing

    % Train a Random Forest classification model for this time point
    modelTime = TreeBagger(numTrees, XTrainTime, yTrain, 'Method', 'classification');

    % Predict trial block labels using the trained model for this time point
    predictionsTime = predict(modelTime, XTestTime);

    % Convert cell array of strings to numeric labels
    predictionsTime = str2double(predictionsTime);

    % Calculate accuracy for this time point
    correctPredictionsTime = sum(predictionsTime == yTest');
    totalSamplesTime = numel(yTest);
    accuracyAtTime(t) = correctPredictionsTime / totalSamplesTime;
end

%Plot accuracy at each time point
figure;
plot(1:numTimePoints, accuracyAtTime);
xlabel('Time Point');
ylabel('Accuracy');
title('Accuracy at Each Time Point');

%% Random forest accuracy over time, but not cumulative. This results in much lower accuracy since it is working with only one sample at a time, and may not be that useful
% Initialize variables to store accuracy at each time point
numTimePoints = size(XTest, 2);  % Number of time points
accuracyAtTime = zeros(1, numTimePoints);

% Train a Random Forest classification model for each time point
for t = 1:numTimePoints
    XTrainTime = XTrain(:, t);  % Use data at the current time point
    XTestTime = XTest(:, t);    % Use data at the current time point for testing

    % Train a Random Forest classification model for this time point
    modelTime = TreeBagger(numTrees, XTrainTime, yTrain, 'Method', 'classification');

    % Predict trial block labels using the trained model for this time point
    predictionsTime = predict(modelTime, XTestTime);

    % Convert cell array of strings to numeric labels
    predictionsTime = str2double(predictionsTime);

    % Calculate accuracy for this time point
    correctPredictionsTime = sum(predictionsTime == yTest');
    totalSamplesTime = numel(yTest);
    accuracyAtTime(t) = correctPredictionsTime / totalSamplesTime;
end

% Plot accuracy at each time point
figure;
plot(1:numTimePoints, accuracyAtTime);
xlabel('Time Point');
ylabel('Accuracy');
title('Accuracy at Each Time Point');


%% shuffle and resample the data to generate data for a circularly shuffled Random Forest decoding
trialCt = size(neuron_mean_concat,1); %number of trials for currently analyzed event
for g = 1:1000 %for each resampling of the data g = 1:uv.resamples
    [num_timepoints, num_signals] = size(neuron_mean_concat);
    shuffled_data = zeros(num_timepoints, num_signals); % Preallocate matrix for efficiency

    for i = 1:num_timepoints
        shift_val = randi(num_signals); % Generate a random shift value for each signal
        shuffled_data(i,:) = circshift(neuron_mean_concat(i,:), shift_val,2); % Perform the circular shuffle
    end
    nullDistTrace(g,:) = nanmean(shuffled_data); %calculate the NaN mean of the shuffled traces
    %calculate the NaN mean of the shuffled event rates
end
% clear shuffled* g t trialCt





y = repelem(1:numConditions, numNeuronsPerCondition);

% Split data into training and testing sets
[trainIdx, testIdx] = crossvalind('HoldOut', size(shuffled_data, 1), 0.2);
XTrain = shuffled_data(trainIdx, 1:end-1);
yTrain = y(trainIdx);
XTest = shuffled_data(testIdx, 1:end-1);
yTest = y(testIdx);

% Train a Random Forest classification model
numTrees = 100; % Number of decision trees in the forest
model = TreeBagger(numTrees, XTrain, yTrain, 'Method', 'classification');

% Predict trial block labels using the trained model
predictions = predict(model, XTest);

% Convert cell array of strings to numeric labels
predictions = str2double(predictions);

% Evaluate the model's performance
correctPredictions = sum(predictions == yTest');
totalSamples = numel(yTest);
accuracy = correctPredictions / totalSamples;
confusionMat = confusionmat(yTest, predictions);

% Visualize the confusion matrix
heatmap(confusionMat, 'XLabel', 'Predicted', 'YLabel', 'Actual');

disp(['Accuracy: ' num2str(accuracy)]);


figure; 
hold on;
for ii = 1:numConditions
    plot(ts1, mean(shuffled_data(y==ii, 1:end)));
end
hold off;




%% Naive Bayes

y = repelem(1:numConditions, numNeuronsPerCondition);

% Split data into training and testing sets
[trainIdx, testIdx] = crossvalind('HoldOut', size(neuron_mean_concat, 1), 0.2);
XTrain = neuron_mean_concat(trainIdx, :);
yTrain = y(trainIdx);
XTest = neuron_mean_concat(testIdx, :);
yTest = y(testIdx);

% Train a Gaussian Naive Bayes classifier
model = fitcnb(XTrain, yTrain);

% Predict trial block labels using the trained model
predictions = predict(model, XTest);

% Evaluate the model's performance
correctPredictions = sum(predictions == yTest');
totalSamples = numel(yTest);
accuracy = correctPredictions / totalSamples;
confusionMat = confusionmat(yTest, predictions);
% Visualize the confusion matrix
heatmap(confusionMat, 'XLabel', 'Predicted', 'YLabel', 'Actual');

disp(['Accuracy: ' num2str(accuracy)]);


figure; 
hold on;
for ii = 1:numConditions
    plot(ts1, mean(neuron_mean_concat(y==ii, :)));
end


% Define the number of classes
numClasses = numConditions;

% Initialize arrays to store precision, recall, and F1-score for each class
precisionArray = zeros(1, numClasses);
recallArray = zeros(1, numClasses);
f1ScoreArray = zeros(1, numClasses);

% Initialize arrays to store TP, FP, and FN for micro-average
TP_micro = zeros(1, numClasses);
FP_micro = zeros(1, numClasses);
FN_micro = zeros(1, numClasses);

for classIdx = 1:numClasses
    % True positives, false positives, false negatives, and true negatives for the current class
    TP = sum(predictions' == classIdx & yTest == classIdx); % True Positives
    FP = sum(predictions' == classIdx & yTest ~= classIdx); % False Positives
    FN = sum(predictions' ~= classIdx & yTest == classIdx); % False Negatives
    
    % Precision, Recall, and F1-Score calculation for the current class
    precision = TP / (TP + FP);
    recall = TP / (TP + FN);
    f1Score = 2 * (precision * recall) / (precision + recall);

    % Store the metrics for the current class
    precisionArray(classIdx) = precision;
    recallArray(classIdx) = recall;
    f1ScoreArray(classIdx) = f1Score;
    
    % Store TP, FP, FN for micro-average
    TP_micro(classIdx) = TP;
    FP_micro(classIdx) = FP;
    FN_micro(classIdx) = FN;

    % Display precision, recall, and F1-score for the current class
    disp(['Class ' num2str(classIdx) ' Precision: ' num2str(precision)]);
    disp(['Class ' num2str(classIdx) ' Recall: ' num2str(recall)]);
    disp(['Class ' num2str(classIdx) ' F1-Score: ' num2str(f1Score)]);
end

% Calculate macro-average precision, recall, and F1-Score
macroPrecision = mean(precisionArray);
macroRecall = mean(recallArray);
macroF1Score = mean(f1ScoreArray);

% Calculate micro-average precision, recall, and F1-Score
microPrecision = sum(TP_micro) / (sum(TP_micro) + sum(FP_micro));
microRecall = sum(TP_micro) / (sum(TP_micro) + sum(FN_micro));
microF1Score = 2 * (microPrecision * microRecall) / (microPrecision + microRecall);

% Display macro and micro averages
disp(['Macro-Average Precision: ' num2str(macroPrecision)]);
disp(['Macro-Average Recall: ' num2str(macroRecall)]);
disp(['Macro-Average F1-Score: ' num2str(macroF1Score)]);
disp(['Micro-Average Precision: ' num2str(microPrecision)]);
disp(['Micro-Average Recall: ' num2str(microRecall)]);
disp(['Micro-Average F1-Score: ' num2str(microF1Score)]);




%% CNN - does not work for some reason
% Assuming you have already loaded and preprocessed your data
% neuron_mean_concats is a 4059x201 array with 4059 rows (neurons) and 201 columns (samples)

numNeuronsPerCondition = neuron_num;
numConditions = size(neuron_mean_concat, 1)/numNeuronsPerCondition;
numSamples = size(neuron_mean_concat, 2);

% Reshape data to [height x width x channels x samples]
height = numNeuronsPerCondition;
width = 1;
channels = numConditions;
inputData = reshape(neuron_mean_concat, height, width, channels, numSamples);

% Create label vector y (corresponding to trial blocks)
y = repelem(1:numConditions, numNeuronsPerCondition);

% Create indices for training and testing
totalNeurons = numNeuronsPerCondition * numConditions;
[trainIdx, testIdx] = crossvalind('HoldOut', totalNeurons, 0.2);

% Split data into training and testing sets
XTrain = inputData(:, :, :, :);
yTrain = categorical(y(trainIdx));
XTest = inputData(:, :, :, testIdx);
yTest = categorical(y(testIdx));

% Define and build a simple CNN model
layers = [
    imageInputLayer([height width channels])
    convolution2dLayer(3, 16)
    reluLayer
    fullyConnectedLayer(numConditions)
    softmaxLayer
    classificationLayer
];

% Compile and train the model
options = trainingOptions('sgdm', 'MaxEpochs', 10, 'Verbose', false);
model = trainNetwork(XTrain, yTrain, layers, options);

% Predict trial block labels using the trained model
predictions = classify(model, XTest);

% Evaluate the model's performance
accuracy = sum(predictions == yTest) / numel(yTest);
confusionMat = confusionmat(yTest, predictions);

% Visualize the confusion matrix
heatmap(confusionMat, 'XLabel', 'Predicted', 'YLabel', 'Actual');

disp(['Accuracy: ' num2str(accuracy)]);


figure; 
hold on;
for ii = 1:numConditions
    plot(ts1, mean(neuron_mean_concat(y==ii, :)));
end
hold off;

%% k Nearest Neighbors

numNeighbors = 5; % Number of neighbors for k-NN

% Create label vector y (corresponding to trial blocks)
y = repelem(1:numConditions, numNeuronsPerCondition);

% Split data into training and testing sets
[trainIdx, testIdx] = crossvalind('HoldOut', size(neuron_mean_concat, 1), 0.2);
XTrain = neuron_mean_concat(trainIdx, :);
yTrain = y(trainIdx);
XTest = neuron_mean_concat(testIdx, :);
yTest = y(testIdx);

% Train a k-Nearest Neighbors classifier
knnModel = fitcknn(XTrain, yTrain, 'NumNeighbors', numNeighbors);

% Predict trial block labels using the trained k-NN model
predictions = predict(knnModel, XTest);

% Evaluate the model's performance
correctPredictions = sum(predictions == yTest');
totalSamples = numel(yTest);
accuracy = correctPredictions / totalSamples;
confusionMat = confusionmat(yTest, predictions);
% Visualize the confusion matrix
heatmap(confusionMat, 'XLabel', 'Predicted', 'YLabel', 'Actual');

disp(['Accuracy: ' num2str(accuracy)]);


figure; 
hold on;
for ii = 1:numConditions
    plot(ts1, mean(neuron_mean_concat(y==ii, :)));
end
hold off;

%% Linear Discriminant Analysis
% Assuming you have already loaded and preprocessed your data
% neuron_mean_concats is a 300x200 array with 300 rows (neurons) and 200 columns (samples)

numNeighbors = 5; % Number of neighbors for k-NN

% Create label vector y (corresponding to trial blocks)
y = repelem(1:numConditions, numNeuronsPerCondition);

% Split data into training and testing sets
[trainIdx, testIdx] = crossvalind('HoldOut', size(neuron_mean_concat, 1), 0.2);
XTrain = neuron_mean_concat(trainIdx, :);
yTrain = y(trainIdx);
XTest = neuron_mean_concat(testIdx, :);
yTest = y(testIdx);

% Train a Linear Discriminant Analysis (LDA) classifier
ldaModel = fitcdiscr(XTrain, yTrain);

% Predict trial block labels using the trained LDA model
predictions = predict(ldaModel, XTest);

% Evaluate the model's performance
correctPredictions = sum(predictions == yTest');
totalSamples = numel(yTest);
accuracy = correctPredictions / totalSamples;
confusionMat = confusionmat(yTest, predictions);
% Visualize the confusion matrix
heatmap(confusionMat, 'XLabel', 'Predicted', 'YLabel', 'Actual');

disp(['Accuracy: ' num2str(accuracy)]);


figure; 
hold on;
for ii = 1:numConditions
    plot(ts1, mean(neuron_mean_concat(y==ii, :)));
end
hold off;
