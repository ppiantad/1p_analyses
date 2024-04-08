% Assuming zall_array is your cell array

% Create a function to extract the first row of each item
extract_first_row = @(x) x(1,:);

% Use cellfun to apply the function to each cell
first_rows = cellfun(extract_first_row, zall_array, 'UniformOutput', false);

% Concatenate the first rows vertically
result = vertcat(first_rows{:});

%% get rid of data where there are < 90 trials (because the next part requires all trials)
% Initialize logical indices to keep track of columns to keep
columns_to_keep = false(1, size(zall_array, 2));

for col = 1:size(zall_array, 2)
    % Check the size of the first element in the column
    if size(zall_array{1, col}, 1) == 90
        % If the size is 90, mark the column to keep
        columns_to_keep(col) = true;
    end
end


% Create a new cell array to store the filtered columns
filtered_zall_array = zall_array(:, columns_to_keep);

%%
% Initialize a cell array to store result vectors
result_vectors = cell(size(filtered_zall_array));

for i = 1:size(filtered_zall_array, 2)
    % Initialize a counter for row indexing
    row_index = 1;
    
    % Initialize a counter for result_vectors indexing
    result_index = 1;
    
    % Iterate over the rows in steps of 5 until the end of the matrix
    while row_index <= size(filtered_zall_array{1, i}, 1)
        % Extract the next 5 rows of the matrix
        rows = filtered_zall_array{1, i}(row_index:min(row_index+4, end), :);
        
        % Calculate the average of the rows
        avg_rows = mean(rows, 1);
        
        % Store the result vector
        result_vectors{result_index, i} = avg_rows;
        
        % Update row index for the next iteration
        row_index = row_index + 5;
        
        % Update result index for the next iteration
        result_index = result_index + 1;
    end
end
%%
% Determine the number of columns
num_columns = size(result_vectors, 2);

% Initialize a variable to store concatenated matrices
concatenated_matrices = [];

% Loop over each column
for i = 1:num_columns
    % Access the matrix stored in result_vectors{1, i}
    matrix = result_vectors{1, i};
    
    % Vertically concatenate the matrix to the existing concatenated matrices
    concatenated_matrices = vertcat(concatenated_matrices, matrix);
end

% respClass_all = respClass_all_array{1, 8}';
[peak_values, time_of_peak_activity] = max(concatenated_matrices, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
neuron_mean_sorted = concatenated_matrices(sort_indices, :);


% Sort the rows of activated_neuron_mean based on peak_times.
% [~, sort_indices] = sort(time_of_peak_activity);
% activated_neuron_mean_sorted = activated_rows(sort_indices, :);

% Now, activated_neuron_mean_sorted contains the rows of neuron_mean filtered by respClass_all == 1
% and sorted by the time of peak activity.

figure;
% Generate the heatmap
imagesc(ts1, 1, neuron_mean_sorted);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');
clim([-1 1])


%%
% Determine the number of columns
num_columns = size(result_vectors, 2);

% Initialize a variable to store concatenated matrices
concatenated_matrices = [];

% Loop over each column
for i = 1:num_columns
    % Access the matrix stored in result_vectors{1, i}
    matrix = result_vectors{end, i};
    
    % Vertically concatenate the matrix to the existing concatenated matrices
    concatenated_matrices = vertcat(concatenated_matrices, matrix);
end

% respClass_all = respClass_all_array{1, 8}';
[peak_values, time_of_peak_activity] = max(concatenated_matrices, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
neuron_mean_sorted = concatenated_matrices(sort_indices, :);


% Sort the rows of activated_neuron_mean based on peak_times.
% [~, sort_indices] = sort(time_of_peak_activity);
% activated_neuron_mean_sorted = activated_rows(sort_indices, :);

% Now, activated_neuron_mean_sorted contains the rows of neuron_mean filtered by respClass_all == 1
% and sorted by the time of peak activity.

figure;
% Generate the heatmap
imagesc(ts1, 1, neuron_mean_sorted);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');
clim([-1 1])