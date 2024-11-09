%% Run data_loop
% [BehavData,trials, varargin_identity_class]=TrialFilter_test(BehavData, 'SHK', 1);
% [BehavData,trials, varargin_identity_class]=TrialFilter_test(BehavData, 'LOSS_PLUS_ONE', 1);
% 


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

        % uncomment below if you want to mean center
        % currentArray_mean = mean(currentArray, 2);
        % currentArray = currentArray-currentArray_mean;
        % Compute the mean activity for each row in the time range 0 to 2 seconds
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



%% attempting logistic regression to predict large vs. small vs. omit from magnitude of SHK response


% Initialize arrays to store the data

shockResponses_all_mice = [];
trialChoices_all_mice = [];

% Iterate through each level of meanZallMouse
for i = 1:length(meanZallMouse) %1:length(meanZallMouse)
    shockResponses = [];
    trialChoices = [];
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{i};
    
    % Get the trial choices for the current mouse
    currentTrialChoices = [behav_tbl_iter{2, 1}{i}.ForceFree behav_tbl_iter{2, 1}{i}.bigSmall];
    currentTrialChoices = currentTrialChoices(currentTrialChoices(:,1) ~= 1, 2);
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        meanValues = meanNestedCellArray{j};
        % Get the current mean values array
        meanValues = meanNestedCellArray{j};
        meanValues = meanValues(currentTrialChoices(:,1) ~= 1);
        
        % Here we use the meanValues as is, no averaging across trials
        % Flatten the meanValues to a single row vector, if needed
        shockResponses = [shockResponses; meanValues];
        
        % Append the corresponding trial choice
        trialChoices = [trialChoices; currentTrialChoices];
    end
    shockResponses_all_mice = [shockResponses_all_mice; shockResponses];
    trialChoices_all_mice = [trialChoices_all_mice; trialChoices];
end

% Convert trialChoices to categorical if not already
% trialChoices = categorical(trialChoices, {'LargeReward', 'SmallReward', 'NoChoice'});

y = trialChoices_all_mice;
mdl = fitmnr(shockResponses_all_mice, y);

% Display the model coefficients
disp('Model Coefficients:');
disp(mdl);

%%
% Step 1: Encode the choices
% Convert trialChoices_all_mice to a binary variable (1 for 1.2, 0 for 0.3)
trialChoices_binary = trialChoices_all_mice == 1.2;

% Step 2: Fit a logistic regression model
% Use the calcium imaging data in shockResponses_all_mice as the predictor
X = shockResponses_all_mice;  % Predictor (calcium response)
y = trialChoices_binary;       % Response (binary choices)

% Fit logistic regression model
% Use MATLAB's fitglm function with 'binomial' distribution for logistic regression
mdl = fitglm(X, y, 'Distribution', 'binomial');

% Display model summary
disp(mdl)

% Step 3: Make predictions
% Get predicted probabilities for the choices
predictedProbabilities = predict(mdl, X);

% Optional: Classify based on probability threshold (e.g., 0.5)
predictedChoices = predictedProbabilities > 0.5;

% Step 4: Assess model accuracy (if desired)
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
