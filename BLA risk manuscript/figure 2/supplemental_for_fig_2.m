% to check whether correlations are strong across blocks load Pre_RDT_RM 10
% variable data, run data_loop.m with [BehavData,trials,varargin_identity_class]=TrialFilter(BehavData,'REW', 1.2, 'BLOCK', 1);
% run code below as necessary

%%
for q = 1:length (behav_tbl_iter{11, 1}  )
    nestedCellArray_1 = behav_tbl_iter{11, 1}  {q};
    % nestedCellArray_2 = behav_tbl_iter{2, 1}{q};
    % if size(nestedCellArray_1, 1) > size(nestedCellArray_2, 1)
    %     delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime(1:end-1,:);
    % else 
    %     delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
    % end
    % 

    for zz = 1:size(nestedCellArray_1, 1)
        valid_start_times = nestedCellArray_1.stTime(2:end);
        valid_choice_times = nestedCellArray_1.choiceTime(1:end-1);
        delay_to_initiation = valid_start_times - valid_choice_times;
    end

    trial_choice_times = nestedCellArray_1.choiceTime - nestedCellArray_1.stTime;
    % delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
    delay_to_collect_post_shk = nestedCellArray_1.collectionTime - nestedCellArray_1.choiceTime;
    trial_choice_times_by_mouse{q} = trial_choice_times;
    delay_to_initiation_by_mouse{q} = delay_to_initiation;
    delay_to_collect_post_shk_by_mouse{q} = delay_to_collect_post_shk;
    clear trial_choice_times delay_to_initiation delay_to_collect_post_shk



end


variable_to_correlate = trial_choice_times_by_mouse;

%%
%%
array_for_means = 11; 

% Initialize the new cell array to store the mean values
meanZallMouse = cell(size(zall_mouse, 2), 1);

% Define the time range for 0 to 2 seconds
timeRange = (ts1 >= -4) & (ts1 <= 0);
% timeRange = (ts1 >= 0) & (ts1 <= 2);
% timeRange = (ts1 >= 1) & (ts1 <= 3);


% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_mouse)
    % Get the current nested cell array
    nestedCellArray_1 = zall_mouse{i, array_for_means};
    nestedCellArray_2 = zall_mouse{i, 1};
    % Initialize the nested cell array for storing mean values
    meanNestedCellArray = cell(size(nestedCellArray_1));
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(nestedCellArray_1)
        % Get the current double array
        currentArray = nestedCellArray_1{j};
        comparisonArray_for_size = nestedCellArray_2{j};
        
        if isequal(variable_to_correlate, delay_to_collect_post_shk_by_mouse)
            currentArray = nestedCellArray_1{j};
        else
            if size(currentArray, 1) > size(comparisonArray_for_size, 1)
                currentArray = currentArray(1:end-1,:);
            else

            end
        end
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

%%
correlationResults = cell(size(meanZallMouse));



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
    trialChoiceTimes = variable_to_correlate{i};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        % Get the current mean values array
        meanValues = meanNestedCellArray{j};
        
        % Check if trialChoiceTimes has the same length as meanValues
        if length(trialChoiceTimes) == length(meanValues)
            % Compute the correlation
            correlationCoeff = corr(meanValues, trialChoiceTimes(:));
        elseif length(trialChoiceTimes) < length(meanValues)
            correlationCoeff = corr(meanValues(1:end-1), trialChoiceTimes(:));
        else
            % If lengths do not match, handle the mismatch (e.g., set correlation to NaN)
            correlationCoeff = NaN;
        end
        
        % Store the correlation coefficient in the nested cell array
        correlationNestedArray(j) = correlationCoeff;
    end
    clear meanValues
    % Store the nested cell array of correlation coefficients in the main cell array
    correlationResults{i} = correlationNestedArray;
end
%%
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

%%
only_shk_responsive_corrs = allCorrelations(prechoice_blocks_2_and_3==1);
not_shk_responsive_corrs = allCorrelations(prechoice_blocks_2_and_3~=1);
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
%%
% Assuming the following variables are defined:
% allCorrelations: array containing correlation coefficients from correlationResults
% only_shk_responsive_corrs: array containing correlation coefficients from a different variable

% Calculate means
mean_only_shk = mean(only_shk_responsive_corrs);
mean_not_shk = mean(not_shk_responsive_corrs);

% Create a histogram for allCorrelations
figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
histogram(not_shk_responsive_corrs , 'Normalization', 'probability', 'FaceColor', 'blue','BinWidth', 0.05,'LineStyle','none');
hold on;

% Create a histogram for only_shk_responsive_corrs on the same figure
histogram(only_shk_responsive_corrs, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05, 'LineStyle','none');
xline(mean_only_shk, 'r')
xline(mean_not_shk, 'g')
% Add labels and title
xlabel('Correlation Coefficient');
ylabel('Probability');
% title('Histograms of Correlation Coefficients');
% legend('All Correlations', 'Only SHK Responsive Correlations');

% Optionally, you can add a vertical line at 0 for reference
yLimits = ylim;
plot([0 0], yLimits, 'k', 'LineWidth', 2);
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
bar_separation_value = 3;

figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
swarmchart(ones(1, length(only_shk_responsive_corrs)), only_shk_responsive_corrs)
hold on
swarmchart(ones(1, length(not_shk_responsive_corrs))*bar_separation_value, not_shk_responsive_corrs)

% yline(mean(only_shk_responsive_corrs), ones(length(only_shk_responsive_corrs)))
plot([0.5; 1.5], [mean(only_shk_responsive_corrs); mean(only_shk_responsive_corrs)], 'LineWidth',3)
plot([bar_separation_value-.5; bar_separation_value+.5], [mean(not_shk_responsive_corrs); mean(not_shk_responsive_corrs)], 'LineWidth',3)
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off



%% FOR BELOW - RUN data_loop_SLEAP RDT_D1

% [BehavData,trials, varargin_identity_class]=TrialFilter_test(BehavData, 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1);

% load the 10 variable RDT_D1 .mat file, run above code or
% block_wise_changes_v1.m
%%
array_for_means = 11; 

% Initialize the new cell array to store the mean values
meanZallMouse_motion = cell(size(zall_mouse, 2), 1);

% Define the time range for 0 to 2 seconds
timeRange = (ts1 >= -4) & (ts1 <= 0);
% timeRange = (ts1 >= 0) & (ts1 <= 2);
% timeRange = (ts1 >= 1) & (ts1 <= 3);


% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_mouse)
    % Get the current nested cell array
    nestedCellArray_1 = zall_mouse{i, array_for_means};

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
    meanZallMouse_motion{i} = meanNestedCellArray;
end

% Now, meanZallMouse contains the mean activity for each row in the time period 0 to 2 seconds
% Each cell in meanZallMouse contains a nested cell array with the
%%
array_for_means = 1; 

% Initialize the new cell array to store the mean values
meanZallMouse = cell(size(zall_mouse, 2), 1);

% Define the time range for 0 to 2 seconds
timeRange = (ts1 >= -4) & (ts1 <= 0);
% timeRange = (ts1 >= 0) & (ts1 <= 2);
% timeRange = (ts1 >= 1) & (ts1 <= 3);


% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_mouse)
    % Get the current nested cell array
    nestedCellArray_1 = zall_mouse{i, array_for_means};
    nestedCellArray_2 = zall_mouse{i, 1};
    % Initialize the nested cell array for storing mean values
    meanNestedCellArray = cell(size(nestedCellArray_1));
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(nestedCellArray_1)
        % Get the current double array
        currentArray = nestedCellArray_1{j};
        comparisonArray_for_size = nestedCellArray_2{j};
        

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

%%
correlationResults = cell(size(meanZallMouse));



% Iterate through each level of meanZallMouse
for i = 1:length(meanZallMouse)
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{i};
    
    % Initialize the nested cell array for storing correlation results
    correlationNestedArray = zeros(size(meanNestedCellArray));
    

    % Get the corresponding trial choice times
    trialChoiceTimes = meanZallMouse_motion{i, 1}{1, 1};
    % trialChoiceTimes = b1_path_length_mouse{1, i};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        % Get the current mean values array
        meanValues = meanNestedCellArray{j};
        
        % Check if trialChoiceTimes has the same length as meanValues
        if length(trialChoiceTimes) == length(meanValues)
            % Compute the correlation
            correlationCoeff = corr(meanValues, trialChoiceTimes(:));
        elseif length(trialChoiceTimes) < length(meanValues)
            correlationCoeff = corr(meanValues(1:end-1), trialChoiceTimes(:));
        else
            % If lengths do not match, handle the mismatch (e.g., set correlation to NaN)
            correlationCoeff = NaN;
        end
        
        % Store the correlation coefficient in the nested cell array
        correlationNestedArray(j) = correlationCoeff;
    end
    clear meanValues
    % Store the nested cell array of correlation coefficients in the main cell array
    correlationResults{i} = correlationNestedArray;
end

%%
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

%%
only_shk_responsive_corrs = allCorrelations(prechoice_block_1==1);
not_shk_responsive_corrs = allCorrelations(prechoice_block_1~=1);
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
%%
% Assuming the following variables are defined:
% allCorrelations: array containing correlation coefficients from correlationResults
% only_shk_responsive_corrs: array containing correlation coefficients from a different variable

% Calculate means
mean_only_shk = mean(only_shk_responsive_corrs);
mean_not_shk = mean(not_shk_responsive_corrs);

% Create a histogram for allCorrelations
figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
histogram(not_shk_responsive_corrs , 'Normalization', 'probability', 'FaceColor', 'blue','BinWidth', 0.05,'LineStyle','none');
hold on;

% Create a histogram for only_shk_responsive_corrs on the same figure
histogram(only_shk_responsive_corrs, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05, 'LineStyle','none');
xline(mean_only_shk, 'r')
xline(mean_not_shk, 'g')
% Add labels and title
xlabel('Correlation Coefficient');
ylabel('Probability');
% title('Histograms of Correlation Coefficients');
% legend('All Correlations', 'Only SHK Responsive Correlations');

% Optionally, you can add a vertical line at 0 for reference
yLimits = ylim;
plot([0 0], yLimits, 'k', 'LineWidth', 2);
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
bar_separation_value = 3;

figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
swarmchart(ones(1, length(only_shk_responsive_corrs)), only_shk_responsive_corrs)
hold on
swarmchart(ones(1, length(not_shk_responsive_corrs))*bar_separation_value, not_shk_responsive_corrs)

% yline(mean(only_shk_responsive_corrs), ones(length(only_shk_responsive_corrs)))
plot([0.5; 1.5], [mean(only_shk_responsive_corrs); mean(only_shk_responsive_corrs)], 'LineWidth',3)
plot([bar_separation_value-.5; bar_separation_value+.5], [mean(not_shk_responsive_corrs); mean(not_shk_responsive_corrs)], 'LineWidth',3)
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off

