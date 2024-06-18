% % Define the time range for 0 to 2 seconds
% timeRange = (ts1 >= 0) & (ts1 <= 2);
% 
% % Iterate through each cell in the cell array
% for i = 1:length(zall_array)
%     % Get the current double array
%     currentArray = zall_array{i};
% 
%     % Compute the mean activity for each row in the time range 0 to 2 seconds
%     meanValues = mean(currentArray(:, timeRange), 2);
% 
%     % Store the mean values in the corresponding cell of the new cell array
%     meanCellArray{i} = meanValues;
% end

%% use these data for mouse x mouse, which is likely better
% run eventRelatedActivityAndClassification with:
%   choiceTime.Outcome_0to2.SHK_1
%   choiceTime.Outcome_0to2.LOSS_PLUS_ONE_1

% Initialize the new cell array to store the mean values
meanZallMouse = cell(length(zall_mouse), 1);

% Define the time range for 0 to 2 seconds
timeRange = (ts1 >= 0) & (ts1 <= 2);
% timeRange = (ts1 >= 1) & (ts1 <= 3);


% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_mouse)
    % Get the current nested cell array
    nestedCellArray_2 = zall_mouse{i};
    
    % Initialize the nested cell array for storing mean values
    meanNestedCellArray = cell(size(nestedCellArray_2));
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(nestedCellArray_2)
        % Get the current double array
        currentArray = nestedCellArray_2{j};
        
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

%%
for q = 1:length (behav_tbl_iter{1, 1})
    nestedCellArray_1 = behav_tbl_iter{1, 1}{q};
    nestedCellArray_2 = behav_tbl_iter{2, 1}{q};
    trial_choice_times = nestedCellArray_2.choiceTime - nestedCellArray_2.stTime;
    delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
    delay_to_collect_post_shk = nestedCellArray_1.collectionTime - nestedCellArray_1.choiceTime;
    trial_choice_times_by_mouse{q} = trial_choice_times;
    delay_to_initiation_by_mouse{q} = delay_to_initiation;
    delay_to_collect_post_shk_by_mouse{q} = delay_to_collect_post_shk;
    clear trial_choice_times delay_to_initiation delay_to_collect_post_shk



end

%%

% Assuming the following variables are defined:
% meanZallMouse: 14x1 cell array where each cell contains another cell array with mean values
% trial_choice_times_by_mouse: 1x11 cell array containing values to correlate with

% Initialize the new cell array to store the correlation results
correlationResults = cell(size(meanZallMouse));

variable_to_correlate = delay_to_initiation_by_mouse;

% Iterate through each level of meanZallMouse
for i = 1:length(meanZallMouse)
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{i};
    
    % Initialize the nested cell array for storing correlation results
    correlationNestedArray = zeros(size(meanNestedCellArray));
    
    % Determine the corresponding index in trial_choice_times_by_mouse
    % Adjust this logic based on how the indices are mapped
    trialIndex = mod(i-1, length(variable_to_correlate)) + 1;
    
    % Get the corresponding trial choice times
    trialChoiceTimes = variable_to_correlate{trialIndex};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        % Get the current mean values array
        meanValues = meanNestedCellArray{j};
        
        % Check if trialChoiceTimes has the same length as meanValues
        if length(trialChoiceTimes) == length(meanValues)
            % Compute the correlation
            correlationCoeff = corr(meanValues, trialChoiceTimes(:));
        else
            % If lengths do not match, handle the mismatch (e.g., set correlation to NaN)
            correlationCoeff = NaN;
        end
        
        % Store the correlation coefficient in the nested cell array
        correlationNestedArray(j) = correlationCoeff;
    end
    
    % Store the nested cell array of correlation coefficients in the main cell array
    correlationResults{i} = correlationNestedArray;
end

% Now, correlationResults contains the correlation coefficients for each nested structure in meanZallMouse
%%
% Assuming correlationResults is defined and contains the correlation coefficients

% Initialize an empty array to collect all correlation coefficients
allCorrelations = [];

% Iterate through each level of correlationResults
for i = 1:length(correlationResults)
    % Get the current nested cell array of correlation coefficients
    correlationNestedArray = correlationResults{i};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(correlationNestedArray)
        % Get the current correlation coefficient
        correlationCoeff = correlationNestedArray(j);
        
        % Check if the coefficient is not NaN (if applicable)
        if ~isnan(correlationCoeff)
            % Append the coefficient to the allCorrelations array
            allCorrelations = [allCorrelations; correlationCoeff];
        end
    end
end

% Now, allCorrelations contains all the correlation coefficients
% Create a histogram of the correlation coefficients
figure;
histogram(allCorrelations);
xlabel('Correlation Coefficient');
ylabel('Frequency');
title('Histogram of Correlation Coefficients');

% Optionally, you can add a vertical line at 0 for reference
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;

%% SHK responsive neurons assumed to be stored in respClass_all_array{1, 1} for this purpose - change as necessary
only_shk_responsive_corrs = allCorrelations(respClass_all_array{1, 1}==1);
not_shk_responsive_corrs = allCorrelations(respClass_all_array{1, 1}==3);
% Now, allCorrelations contains all the correlation coefficients
% Create a histogram of the correlation coefficients
figure;
histogram(only_shk_responsive_corrs);
xlabel('Correlation Coefficient');
ylabel('Frequency');
title('Histogram of Correlation Coefficients');

% Optionally, you can add a vertical line at 0 for reference
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;

%%
% Assuming the following variables are defined:
% allCorrelations: array containing correlation coefficients from correlationResults
% only_shk_responsive_corrs: array containing correlation coefficients from a different variable

% Create a histogram for allCorrelations
figure;
histogram(not_shk_responsive_corrs , 'Normalization', 'probability', 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'BinWidth', 0.05);
hold on;

% Create a histogram for only_shk_responsive_corrs on the same figure
histogram(only_shk_responsive_corrs, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.5, 'BinWidth', 0.05);

% Add labels and title
xlabel('Correlation Coefficient');
ylabel('Probability');
title('Histograms of Correlation Coefficients');
legend('All Correlations', 'Only SHK Responsive Correlations');

% Optionally, you can add a vertical line at 0 for reference
yLimits = ylim;
plot([0 0], yLimits, 'k--', 'LineWidth', 2);
xtickformat('%.2f');
ytickformat('%.2f');
hold off;

% Perform a Kolmogorov-Smirnov test to compare the two distributions
[h, p] = kstest2(not_shk_responsive_corrs , only_shk_responsive_corrs);

% Display the results of the statistical test
fprintf('Kolmogorov-Smirnov test result:\n');
fprintf('h = %d (0 means the null hypothesis cannot be rejected, 1 means it can be rejected)\n', h);
fprintf('p-value = %.4f\n', p);

[h,p,ci,stats] = ttest2(not_shk_responsive_corrs , only_shk_responsive_corrs)

%%

% Calculate means
mean_only_shk = mean(only_shk_responsive_corrs);
mean_not_shk = mean(not_shk_responsive_corrs);



% % Plot means as bars
% figure;
% hold on;
% bar(1, mean_only_shk, 'FaceColor', 'r'); % Red bar for 'only_shk_responsive_corrs'
% bar(2, mean_not_shk, 'FaceColor', 'b'); % Blue bar for 'not_shk_responsive_corrs'
% 
% swarmchart(1, only_shk_responsive_corrs, 5)
% 
% % Scatter individual data points
% scatter(ones(size(only_shk_responsive_corrs)), only_shk_responsive_corrs, 'r', 'filled', 'jitter', 'on', 'jitterAmount', 0.15); % Red points with jitter
% scatter(2 * ones(size(not_shk_responsive_corrs)), not_shk_responsive_corrs, 'b', 'filled', 'jitter', 'on', 'jitterAmount', 0.15); % Blue points with jitter
% 
% % Customize plot
% xlim([0.5, 2.5]);
% xticks([1 2]);
% xticklabels({'Only Shk Responsive', 'Not Shk Responsive'});
% ylabel('Correlation Values');
% title('Correlation Values and Means');
% grid on;
% 
% % Add legend
% legend({'Mean Only Shk', 'Mean Not Shk', 'Individual Only Shk', 'Individual Not Shk'}, 'Location', 'best');
% 
% hold off;

figure;
width = 500; % Width of the figure
height = 1000; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
swarmchart(ones(1, length(only_shk_responsive_corrs)), only_shk_responsive_corrs)
hold on
swarmchart(ones(1, length(not_shk_responsive_corrs))*2, not_shk_responsive_corrs)

% yline(mean(only_shk_responsive_corrs), ones(length(only_shk_responsive_corrs)))
plot([0.5; 1.5], [mean(only_shk_responsive_corrs); mean(only_shk_responsive_corrs)], 'LineWidth',3)
plot([1.5; 2.5], [mean(not_shk_responsive_corrs); mean(not_shk_responsive_corrs)], 'LineWidth',3)
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off


%% attempting logistic regression to predict large vs. small vs. omit from magnitude of SHK response


% Initialize arrays to store the data
shockResponses = [];
trialChoices = [];

% Iterate through each level of meanZallMouse
for i = 1:1 %1:length(meanZallMouse)
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
end

% Convert trialChoices to categorical if not already
trialChoices = categorical(trialChoices, {'LargeReward', 'SmallReward', 'NoChoice'});

y = trialChoices;
mdl = fitmnr(shockResponses, y);

% Display the model coefficients
disp('Model Coefficients:');
disp(mdl);

%% minor attempts to correlate with MOTION (velocity after shock)

meanNestedCellArray = meanZallMouse{6};

for j = 1:length(meanNestedCellArray)
    % Get the current mean values array
    meanValues = meanNestedCellArray{j};

    % Check if trialChoiceTimes has the same length as meanValues
    if length(meanValues_motion) == length(meanValues)
        % Compute the correlation
        correlationCoeff = corr(meanValues, meanValues_motion(:));
    else
        % If lengths do not match, handle the mismatch (e.g., set correlation to NaN)
        correlationCoeff = NaN;
    end

    % Store the correlation coefficient in the nested cell array
    correlationNestedArray(j) = correlationCoeff;
end
