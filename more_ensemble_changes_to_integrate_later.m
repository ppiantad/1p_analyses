% load dataset where you have done prechoice, post-choice, consumption (all block 1), shk
% and also prechoice, post-choice, consumption (blocks 2 and 3)

exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} ~= 1 & respClass_all_array{1,4} ~= 1;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1)
exclusive_activated_session_2 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} ~= 1 & respClass_all_array{1,4} ~= 1;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2)
exclusive_activated_session_3 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} == 1 & respClass_all_array{1,4} ~= 1;
exclusive_activated_session_3_sum = sum(exclusive_activated_session_3)
exclusive_activated_session_4 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} ~= 1 & respClass_all_array{1,4} == 1;
exclusive_activated_session_4_sum = sum(exclusive_activated_session_4)




exclusive_activated_session_prechoice_risky = respClass_all_array{1,5} == 1 & respClass_all_array{1,6} ~= 1 & respClass_all_array{1,7} ~= 1 & respClass_all_array{1,4} ~= 1;
exclusive_activated_session_prechoice_risky_sum = sum(exclusive_activated_session_prechoice_risky)

sum(exclusive_activated_session_prechoice_risky == 1 & exclusive_activated_session_1 == 1)


exclusive_activated_session_post_choice_rew_risky = respClass_all_array{1,5} ~= 1 & respClass_all_array{1,6} == 1 & respClass_all_array{1,7} ~= 1 & respClass_all_array{1,4} ~= 1;
exclusive_activated_session_post_choice_rew_risky_sum = sum(exclusive_activated_session_post_choice_rew_risky)


exclusive_activated_session_consumption_risky = respClass_all_array{1,5} ~= 1 & respClass_all_array{1,6} ~= 1 & respClass_all_array{1,7} == 1 & respClass_all_array{1,4} ~= 1;
exclusive_activated_session_consumption_risky_sum = sum(exclusive_activated_session_consumption_risky)

sum(exclusive_activated_session_consumption_risky == 1 & exclusive_activated_session_3 == 1)

figure; plot(ts1, mean(neuron_mean_array{1,3}(exclusive_activated_session_3 == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1,7}(exclusive_activated_session_consumption_risky == 1, :)))
figure; plot(ts1, mean(neuron_mean_array{1,1}(exclusive_activated_session_1 == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1,5}(exclusive_activated_session_prechoice_risky == 1, :)))
figure; plot(ts1, mean(neuron_mean_array{1,2}(exclusive_activated_session_2 == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1,6}(exclusive_activated_session_post_choice_rew_risky == 1, :)))


%% trying to create code to make a modified version of Courtin's correlation matrices from the outcome changes
% this is a little tricky with the RDT, because mice have diff #s of
% trials. so, it might be best to take the first, median, and last trial
% (assuming mice have 3 or > trials, this should work)
pre_choice_only = zall_array(:, exclusive_activated_session_1);

consumption_neurons_only = zall_array(:, exclusive_activated_session_3);
initial_consumption_neurons_during_risk = zall_array(:, exclusive_activated_session_3);

consumption_neurons_only_means = neuron_mean_array{1, 3}(exclusive_activated_session_3, :);

initial_consumption_neurons_only_means_during_risk = neuron_mean_array{1, 7}(exclusive_activated_session_3, :);

new_consumption_neurons_only_means_risky = neuron_mean_array{1, 7}(exclusive_activated_session_consumption_risky, :);

new_consumption_neurons_only_means_block_1 = neuron_mean_array{1, 3}(exclusive_activated_session_consumption_risky, :);

array_to_use = consumption_neurons_only;
% use if you want to calc correlations to entire window 
time_index = ts1 == ts1;
% use if you want ot calc correlatiosn to some subwindow
% time_index= ts1 >= 1 & ts1 <= 3;
% Assuming 'consumption_neurons_only' is your 7x77 cell array
numColumns = size(array_to_use, 2);  % Number of columns
firstRows = {};  % To store the first row from each column
medianRows = {};  % To store the median row from each column
finalRows = {};  % To store the final row from each column

for col = 1:numColumns
    % rows 1 to 3 are from block 1, rows 5 to 7 are from RDT
    columnData = array_to_use{3, col};  % Get the data in the first row of the current column
    
    if ~isempty(columnData)
        % Get the first row
        firstRows{col} = columnData(1, :);
        
        % Get the median row
        medianRowIndex = round(size(columnData, 1) / 2);
        medianRows{col} = columnData(medianRowIndex, :);
        
        % Get the final row
        finalRows{col} = columnData(end, :);
    end
end

% Convert cell arrays to matrices (if all rows have the same length)
firstRows = cell2mat(firstRows.');
medianRows = cell2mat(medianRows.');
finalRows = cell2mat(finalRows.');

% Display the results
% disp('First Rows:');
% disp(firstRows);
% 
% disp('Median Rows:');
% disp(medianRows);
% 
% disp('Final Rows:');
% disp(finalRows);

%%

% Initialize correlation matrices
correlationsFirstRows = zeros(size(firstRows, 1));
correlationsMedianRows = zeros(size(medianRows, 1));
correlationsFinalRows = zeros(size(finalRows, 1));

% Calculate correlations for firstRows
for i = 1:size(firstRows, 1)
    for j = 1:size(firstRows, 1)
        correlationsFirstRows(i, j) = corr(firstRows(i, time_index)', firstRows(j, time_index)');
    end
end

% Calculate correlations for medianRows
for i = 1:size(medianRows, 1)
    for j = 1:size(medianRows, 1)
        correlationsMedianRows(i, j) = corr(medianRows(i, time_index)', medianRows(j, time_index)');
    end
end

% Calculate correlations for finalRows
for i = 1:size(finalRows, 1)
    for j = 1:size(finalRows, 1)
        correlationsFinalRows(i, j) = corr(finalRows(i, time_index)', finalRows(j, time_index)');
    end
end

% Calculate correlations for firstRows
for i = 1:size(consumption_neurons_only_means, 1)
    for j = 1:size(consumption_neurons_only_means, 1)
        correlations_consumption_neurons_only_means(i, j) = corr(consumption_neurons_only_means(i, time_index)', consumption_neurons_only_means(j, time_index)');
    end
end

for i = 1:size(initial_consumption_neurons_only_means_during_risk, 1)
    for j = 1:size(initial_consumption_neurons_only_means_during_risk, 1)
        correlations_initial_consumption_neurons_only_means_during_risk(i, j) = corr(initial_consumption_neurons_only_means_during_risk(i, time_index)', initial_consumption_neurons_only_means_during_risk(j, time_index)');
    end
end


for i = 1:size(new_consumption_neurons_only_means_risky, 1)
    for j = 1:size(new_consumption_neurons_only_means_risky, 1)
        correlations_new_consumption_neurons_only_means_risky(i, j) = corr(new_consumption_neurons_only_means_risky(i, time_index)', new_consumption_neurons_only_means_risky(j, time_index)');
    end
end


for i = 1:size(new_consumption_neurons_only_means_block_1, 1)
    for j = 1:size(new_consumption_neurons_only_means_block_1, 1)
        correlations_new_consumption_neurons_only_means_block_1(i, j) = corr(new_consumption_neurons_only_means_block_1(i, time_index)', new_consumption_neurons_only_means_block_1(j, time_index)');
    end
end

% Calculate mean correlations, including the diagonal
meanCorrFirstRows = mean(correlationsFirstRows(:));
meanCorrMedianRows = mean(correlationsMedianRows(:));
meanCorrFinalRows = mean(correlationsFinalRows(:));
meanCorr_consumption_neurons_only_means = mean(correlations_consumption_neurons_only_means(:));
meanCorr_correlations_initial_consumption_neurons_only_means_during_risk = mean(correlations_initial_consumption_neurons_only_means_during_risk(:));
meanCorr_correlations_new_consumption_neurons_only_means_risky = mean(correlations_new_consumption_neurons_only_means_risky(:));
meanCorr_correlations_new_consumption_neurons_only_means_block_1 = mean(correlations_new_consumption_neurons_only_means_block_1(:));
% Display the results
disp('Mean pairwise correlation for First Rows:');
disp(meanCorrFirstRows);

disp('Mean pairwise correlation for Median Rows:');
disp(meanCorrMedianRows);

disp('Mean pairwise correlation for Final Rows:');
disp(meanCorrFinalRows);
