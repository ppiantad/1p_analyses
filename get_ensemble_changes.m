% Initialize empty arrays to store the indices
indices_2_1 = [];
indices_3_1 = [];

% Iterate over the cells in the 2nd column ({1,2})
for i = 1:numel(respClass_all_array{1,2})
    % Check if the value is 1 in the 2nd column, but not in the 1st column
    if respClass_all_array{1,2}(i) == 1 && respClass_all_array{1,1}(i) ~= 1
        indices_2_1 = [indices_2_1, i]; % Store the index
    end
end

% Iterate over the cells in the 3rd column ({1,3})
for i = 1:numel(respClass_all_array{1,3})
    % Check if the value is 1 in the 3rd column, but not in the 1st column
    if respClass_all_array{1,3}(i) == 1 && respClass_all_array{1,1}(i) ~= 1
        indices_3_1 = [indices_3_1, i]; % Store the index
    end
end

% % Display the indices for the comparison between {1,1} and {1,2}
% disp("Indices for comparison between {1,1} and {1,2}:");
% disp(indices_2_1);
% 
% % Display the indices for the comparison between {1,1} and {1,3}
% disp("Indices for comparison between {1,1} and {1,3}:");
% disp(indices_3_1);


%%
%classic loss of ensemble strength
figure; 
plot(nanmean(neuron_mean_array{1,1}(respClass_all_array{1,1} == 1,:)))
hold on; plot(nanmean(neuron_mean_array{1,2}(respClass_all_array{1,1} == 1,:)))
hold on; plot(nanmean(neuron_mean_array{1,3}(respClass_all_array{1,1} == 1,:)))



figure;
shadedErrorBar(ts1, nanmean(neuron_mean_array{1,1}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 1}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 2}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 2}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')



%uut actually ensemble strength is maintained, just by different neurons
figure; 
plot(nanmean(neuron_mean_array{1,1}(respClass_all_array{1,1} == 1,:)))
hold on; plot(nanmean(neuron_mean_array{1,2}(respClass_all_array{1,2} == 1,:)))
hold on; plot(nanmean(neuron_mean_array{1,3}(respClass_all_array{1,3} == 1,:)))

%here are the different neurons
figure; 
% plot(nanmean(neuron_mean_array{1,1}(intersect(indices_2_1, indices_3_1),:)))
plot(nanmean(neuron_mean_array{1,1}(union(indices_2_1, indices_3_1),:)))
hold on; plot(nanmean(neuron_mean_array{1,2}(indices_2_1,:)))
hold on; plot(nanmean(neuron_mean_array{1,3}(indices_3_1,:)))

%%
% Get neurons from zall array that correspond to each event

zall_array_filtered = cell(size(zall_array));
for hh = 1:size(zall_array, 1)
    zall_original_ensemble_neurons_x_trial{hh} = zall_array(hh, respClass_all_array{1, 1}==1);
    zall_block_2_ensemble_neurons_x_trial{hh} = zall_array(hh, indices_2_1);
    zall_block_3_ensemble_neurons_x_trial{hh} = zall_array(hh, indices_3_1);
end

%% THIS DOESN'T WORK BECAUSE MICE HAVE DIFFERENT #S OF TRIALS. MAYBE IT MAKES MORE SENSE TO JUST ESCHEW THIS DETAILED LOOK & INSTEAD CALCULATE CORRELATIONS FROM THE MEAN DATA FOR EACH NEURON THAT IS PRESENT IN AN ENSEMBLE, ACROSS BLOCKS
% Initialize cell arrays to store the means
zall_original_ensemble_neurons_x_trial_means = cell(size(zall_original_ensemble_neurons_x_trial));

for i = 1:numel(zall_original_ensemble_neurons_x_trial)
    level_cell_array = zall_original_ensemble_neurons_x_trial{i};
    % Assuming your cell array is named 'level_cell_array'
    num_rows = size(level_cell_array{1}, 1); % Number of rows in each 24x160 matrix
    num_cols = numel(level_cell_array); % Number of columns in the cell array
    result = cell(num_rows, 1);

    for row = 1:num_rows
        row_values = [];
        for col = 1:num_cols
            row_values = [row_values; level_cell_array{col}(row, :)];
        end
        result{row} = row_values;
        level_means(row, :) = mean(row_values);
    end
    results{i} = result;
    zall_original_ensemble_neurons_x_trial_means{i} = level_means;
    clear level_means
end


figure; plot(ts1, mean(original_ensemble_neurons_x_trial_means{1, 1}))
hold on; plot(ts1, mean(original_ensemble_neurons_x_trial_means{1, 2}))
hold on; plot(ts1, mean(original_ensemble_neurons_x_trial_means{1, 3}))

%%
block_1_og_ensemble = neuron_mean_array{1, 1}(respClass_all_array{1,1} == 1, :);

% Assuming your array is named results{1, 1}{1, 1}
array = block_1_og_ensemble;

% Get the number of rows
numRows = size(array, 1);

% Preallocate a matrix to store the correlation coefficients
correlationMatrix = zeros(numRows, numRows);

% Calculate the correlation coefficients for each pair of rows
for i = 1:numRows
    for j = 1:numRows
        correlationMatrix(i, j) = corr(array(i, :)', array(j, :)');
    end
end

% Display the correlation matrix
% disp(correlationMatrix);

% If you only want unique correlations, you can use the triu function
uniqueCorrelationMatrix = triu(correlationMatrix, 1);
uniqueCorrelations = unique(uniqueCorrelationMatrix);

mean(uniqueCorrelations)

%%
block_2_og_ensemble = neuron_mean_array{1, 2}(respClass_all_array{1,1} == 1, :);

% Assuming your array is named results{1, 1}{1, 1}
array = block_2_og_ensemble;

% Get the number of rows
numRows = size(array, 1);

% Preallocate a matrix to store the correlation coefficients
correlationMatrix = zeros(numRows, numRows);

% Calculate the correlation coefficients for each pair of rows
for i = 1:numRows
    for j = 1:numRows
        correlationMatrix(i, j) = corr(array(i, :)', array(j, :)');
    end
end

% Display the correlation matrix
% disp(correlationMatrix);

% If you only want unique correlations, you can use the triu function
uniqueCorrelationMatrix = triu(correlationMatrix, 1);
uniqueCorrelations = unique(uniqueCorrelationMatrix);

mean(uniqueCorrelations)

%%
block_3_og_ensemble = neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1, :);

% Assuming your array is named results{1, 1}{1, 1}
array = block_3_og_ensemble;

% Get the number of rows
numRows = size(array, 1);

% Preallocate a matrix to store the correlation coefficients
correlationMatrix = zeros(numRows, numRows);

% Calculate the correlation coefficients for each pair of rows
for i = 1:numRows
    for j = 1:numRows
        correlationMatrix(i, j) = corr(array(i, :)', array(j, :)');
    end
end

% Display the correlation matrix
% disp(correlationMatrix);

% If you only want unique correlations, you can use the triu function
uniqueCorrelationMatrix = triu(correlationMatrix, 1);
uniqueCorrelations = unique(uniqueCorrelationMatrix);

nanmean(uniqueCorrelations)

%%
block_1_new_ensemble = neuron_mean_array{1, 1}(indices_2_1, :);

% Assuming your array is named results{1, 1}{1, 1}
array = block_1_new_ensemble;

% Get the number of rows
numRows = size(array, 1);

% Preallocate a matrix to store the correlation coefficients
correlationMatrix = zeros(numRows, numRows);

% Calculate the correlation coefficients for each pair of rows
for i = 1:numRows
    for j = 1:numRows
        correlationMatrix(i, j) = corr(array(i, :)', array(j, :)');
    end
end

% Display the correlation matrix
% disp(correlationMatrix);

% If you only want unique correlations, you can use the triu function
uniqueCorrelationMatrix = triu(correlationMatrix, 1);
uniqueCorrelations = unique(uniqueCorrelationMatrix);

mean(uniqueCorrelations)

%%
block_2_new_ensemble = neuron_mean_array{1, 2}(indices_2_1, :);

% Assuming your array is named results{1, 1}{1, 1}
array = block_2_new_ensemble;

% Get the number of rows
numRows = size(array, 1);

% Preallocate a matrix to store the correlation coefficients
correlationMatrix = zeros(numRows, numRows);

% Calculate the correlation coefficients for each pair of rows
for i = 1:numRows
    for j = 1:numRows
        correlationMatrix(i, j) = corr(array(i, :)', array(j, :)');
    end
end

% Display the correlation matrix
% disp(correlationMatrix);

% If you only want unique correlations, you can use the triu function
uniqueCorrelationMatrix = triu(correlationMatrix, 1);
uniqueCorrelations = unique(uniqueCorrelationMatrix);

mean(uniqueCorrelations)

%%
block_3_new_ensemble = neuron_mean_array{1, 3}(indices_2_1, :);

% Assuming your array is named results{1, 1}{1, 1}
array = block_3_new_ensemble;

% Get the number of rows
numRows = size(array, 1);

% Preallocate a matrix to store the correlation coefficients
correlationMatrix = zeros(numRows, numRows);

% Calculate the correlation coefficients for each pair of rows
for i = 1:numRows
    for j = 1:numRows
        correlationMatrix(i, j) = corr(array(i, :)', array(j, :)');
    end
end

% Display the correlation matrix
% disp(correlationMatrix);

% If you only want unique correlations, you can use the triu function
uniqueCorrelationMatrix = triu(correlationMatrix, 1);
uniqueCorrelations = unique(uniqueCorrelationMatrix);

nanmean(uniqueCorrelations)










%%
%BELOW USES THE ENSEMBLE MEANS FOR CORRELATIONS, WHICH MIGHT BE BETTTER
%THAN WHAT IS DONE ABOVE?
% Initialize empty arrays to store the indices
% need to make sure you have run eventRelatedActivity with mice excluded
% who did not finish the RDT session

indices_2_1 = [];
indices_3_1 = [];

for zz = 1:size(respClass_all_array_mouse, 1)
    step_1_data =respClass_all_array_mouse(zz, :);

    indices_2_1{zz} = find(step_1_data{1,2} == 1 & step_1_data{1,1} ~= 1);
    indices_3_1{zz} = find(step_1_data{1,3} == 1 & step_1_data{1,1} ~= 1);



end






%%

for zz = 1:size(neuron_mean_mouse, 1)
    step_1_data = neuron_mean_mouse(zz, :);
    for qq = 1:size(step_1_data, 2)
        ensemble_1_mean{qq}(zz, :) = mean(step_1_data{1, qq}(respClass_all_array_mouse{zz,1}==1, :));
        ensemble_2_mean{qq}(zz, :) = mean(step_1_data{1, qq}(indices_2_1{1, zz}, :));
        ensemble_3_mean{qq}(zz, :) = mean(step_1_data{1, qq}(indices_3_1{1, zz}, :));
    end



end


%%
% block_1_og_ensemble = neuron_mean_array{1, 1}(respClass_all_array{1,1} == 1, :);

% Assuming your array is named results{1, 1}{1, 1}
array = ensemble_1_mean{1,1};

% Get the number of rows
numRows = size(array, 1);

% Preallocate a matrix to store the correlation coefficients
correlationMatrix = zeros(numRows, numRows);

% Calculate the correlation coefficients for each pair of rows
for i = 1:numRows
    for j = 1:numRows
        correlationMatrix(i, j) = corr(array(i, :)', array(j, :)');
    end
end

% Display the correlation matrix
% disp(correlationMatrix);

% If you only want unique correlations, you can use the triu function
uniqueCorrelationMatrix = triu(correlationMatrix, 1);
uniqueCorrelations = unique(uniqueCorrelationMatrix);

mean(uniqueCorrelations)



%%
% block_1_og_ensemble = neuron_mean_array{1, 1}(respClass_all_array{1,1} == 1, :);

% Assuming your array is named results{1, 1}{1, 1}
array = ensemble_2_mean{1,1};

% Get the number of rows
numRows = size(array, 1);

% Preallocate a matrix to store the correlation coefficients
correlationMatrix = zeros(numRows, numRows);

% Calculate the correlation coefficients for each pair of rows
for i = 1:numRows
    for j = 1:numRows
        correlationMatrix(i, j) = corr(array(i, :)', array(j, :)');
    end
end

% Display the correlation matrix
% disp(correlationMatrix);

% If you only want unique correlations, you can use the triu function
uniqueCorrelationMatrix = triu(correlationMatrix, 1);
uniqueCorrelations = unique(uniqueCorrelationMatrix);

nanmean(uniqueCorrelations)


%%
% block_1_og_ensemble = neuron_mean_array{1, 1}(respClass_all_array{1,1} == 1, :);

% Assuming your array is named results{1, 1}{1, 1}
array = ensemble_3_mean{1,1};

% Get the number of rows
numRows = size(array, 1);

% Preallocate a matrix to store the correlation coefficients
correlationMatrix = zeros(numRows, numRows);

% Calculate the correlation coefficients for each pair of rows
for i = 1:numRows
    for j = 1:numRows
        correlationMatrix(i, j) = corr(array(i, :)', array(j, :)');
    end
end

% Display the correlation matrix
% disp(correlationMatrix);

% If you only want unique correlations, you can use the triu function
uniqueCorrelationMatrix = triu(correlationMatrix, 1);
uniqueCorrelations = unique(uniqueCorrelationMatrix);

nanmean(uniqueCorrelations)
