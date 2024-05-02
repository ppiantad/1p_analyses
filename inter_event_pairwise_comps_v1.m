select_mouse = 'BLA_Insc_27';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';

select_mouse_cells_indices = find(strcmp(mouse_cells(1,:), select_mouse));
% Assuming respClass_all_array_filtered is your 1x3 cell array
respClass_all_array_filtered = cell(1, 3); % Initialize filtered array


for i = 1:3
    respClass_all_array_filtered{i} = respClass_all_array{i}(select_mouse_cells_indices);

end



% Initialize empty arrays to store the indices
indices_2_1 = [];
indices_3_1 = [];

% Iterate over the cells in the 2nd column ({1,2})
for i = 1:numel(respClass_all_array_filtered{1,2})
    % Check if the value is 1 in the 2nd column, but not in the 1st column
    if respClass_all_array_filtered{1,2}(i) == 1 && respClass_all_array_filtered{1,1}(i) ~= 1
        indices_2_1 = [indices_2_1, i]; % Store the index
    end
end

% Iterate over the cells in the 3rd column ({1,3})
for i = 1:numel(respClass_all_array_filtered{1,3})
    % Check if the value is 1 in the 3rd column, but not in the 1st column
    if respClass_all_array_filtered{1,3}(i) == 1 && respClass_all_array_filtered{1,1}(i) ~= 1
        indices_3_1 = [indices_3_1, i]; % Store the index
    end
end


%%

% original_consumption_ensemble = [neuron_mean_mouse{select_mouse_index, 1}(respClass_all_array_filtered{1, 1}==1, :)];



zall_mouse_filtered = zall_mouse(select_mouse_index,:);

for hh = 1:size(zall_mouse_filtered, 2)
    original_ensemble_neurons_x_trial{hh} = zall_mouse_filtered{1, hh}(respClass_all_array_filtered{1, 1}==1);
    block_2_ensemble_neurons_x_trial{hh} = zall_mouse_filtered{1, hh}(indices_2_1);
    block_3_ensemble_neurons_x_trial{hh} = zall_mouse_filtered{1, hh}(indices_3_1);
end




% Initialize cell arrays to store the means
original_ensemble_neurons_x_trial_means = cell(size(original_ensemble_neurons_x_trial));

% % Loop through the first level of cell array
% for i = 1:numel(original_ensemble_neurons_x_trial)
%     % Get the cell array at the current level
%     level_cell_array = original_ensemble_neurons_x_trial{i};
% 
%     % Initialize array to store means for the current level
%     level_means = zeros(size(level_cell_array{1}, 1), size(level_cell_array{1}, 2));
% 
%     % Loop through each row of the cell arrays in the current level
%     for row = 1:size(level_cell_array{1}, 1)
%         % Initialize array to store values for the current row
%         % row_values = zeros(1, numel(level_cell_array));
% 
%         % Loop through each cell array in the current level
%         for j = 1:numel(level_cell_array)
%             % Extract the current row from the current cell array
%             current_row = level_cell_array{j}(row,:);
% 
%             % Store the values for the current row
%             row_values(j,:) = current_row;
%         end
% 
%         % Store the mean of the current row
%         level_means(row, :) = mean(row_values);
%     end
% 
%     % Store the means for the current level
%     original_ensemble_neurons_x_trial_means{i} = level_means;
% end


%  THIS WORKS - BUILD OFF OF THIS CODE TOMORROW

for i = 1:numel(original_ensemble_neurons_x_trial)
    level_cell_array = original_ensemble_neurons_x_trial{i};
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
    original_ensemble_neurons_x_trial_means{i} = level_means;
    clear level_means
end


figure; plot(ts1, mean(original_ensemble_neurons_x_trial_means{1, 1}))
hold on; plot(ts1, mean(original_ensemble_neurons_x_trial_means{1, 2}))
hold on; plot(ts1, mean(original_ensemble_neurons_x_trial_means{1, 3}))


%%

% original_plus_block_2_emergence = [test; neuron_mean_mouse{select_mouse_index, 1}(respClass_all_array_filtered{1, 1}~=1 & respClass_all_array_filtered{1, 2}==1 & respClass_all_array_filtered{1, 3}~=1, :)];
% 
% 
% test = [neuron_mean_mouse{select_mouse_index, 1}(respClass_all_array_filtered{1, 1}==1 & respClass_all_array_filtered{1, 2}~=1 & respClass_all_array_filtered{1, 3}~=1, :)];
% test = [test; neuron_mean_mouse{select_mouse_index, 2}(indices_2_1,:)];
% test = [test; neuron_mean_mouse{select_mouse_index, 3}(indices_3_1,:)];
% % test = [test; neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 2}~=1 & respClass_all_array_filtered{1, 1}~=1 & respClass_all_array_filtered{1,3}~=1,:)]
% 
% 
% test = [neuron_mean_mouse{select_mouse_index, 1}(indices_2_1,:)];
% test = [test; neuron_mean_mouse{select_mouse_index, 2}(indices_2_1,:)];
% test = [test; neuron_mean_mouse{select_mouse_index, 3}(indices_2_1,:)];



% Initialize an empty array to store concatenated means
concatenated_means = [];

% Loop through each level of original_ensemble_neurons_x_trial_means
for i = 1:numel(original_ensemble_neurons_x_trial_means)
    % Vertically concatenate the means from the current level
    concatenated_means = vertcat(concatenated_means, original_ensemble_neurons_x_trial_means{i});
end


% data = test;


data = concatenated_means;

alpha = 0.0001;

% Initialize matrices to store correlation coefficients and p-values
correlation_matrix = zeros(size(data, 1));
p_value_matrix = zeros(size(data, 1));

% Calculate correlation coefficients and p-values between rows
for i = 1:size(data, 1)
    for j = 1:size(data, 1)
        [corr_coeff, p_value] = corrcoef(data(i, :)', data(j, :)');
        correlation_matrix(i, j) = corr_coeff(1, 2); % Store correlation coefficient
        p_value_matrix(i, j) = p_value(1, 2); % Store p-value
    end
end

% Plot the correlation matrix
figure;
imagesc(correlation_matrix);
colorbar; % Add a colorbar to the plot
axis square; % Make the plot square for better visualization
title('Correlation Matrix');
xlabel('Trial Number');
ylabel('Trial Number');

% Show row and column indices on the plot
xticks(0:5:size(data, 1));
yticks(0:5:size(data, 1));

size_start = 0;
for ii = 1:size(original_ensemble_neurons_x_trial_means, 2)
    size_start = size_start + size(original_ensemble_neurons_x_trial_means{1, ii}, 1);
    xline(size_start, 'LineWidth', 2)
    yline(size_start, 'LineWidth', 2)
    
end


% If you want to customize the color map, you can use colormap function
% For example, using a blue-white-red colormap:
colormap(bluewhitered);

% If you want to limit the color scale to the range [0, 1]
caxis([-1 1]); % Assuming correlations range from -1 to 1

% % Display p-values as text on the plot
% for i = 1:size(data, 1)
%     for j = 1:size(data, 1)
%         text(j, i, sprintf('p = %.4f', p_value_matrix(i, j)), ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     end
% end

og_ensemble = zall_array(:, respClass_all_array{1, 1} == 1);


for zz = 1:size(og_ensemble, 2)
    og_ensemble_data = og_ensemble{1, zz};
   


end



% Assuming your array is named results{1, 1}{1, 1}
array = results{1, 1}{2, 1}  ;

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

% Display the unique correlations
% disp(uniqueCorrelations);
% mean(uniqueCorrelations)


%%

% Preallocate a cell array to store the mean correlation for each level
meanCorrelationArray = cell(size(results));

% Loop through each column of results
for col = 1:numel(results)
    % Get the subarray from the current column
    subarray = results{col}; % Assuming the subarray is 1xN cell array
    
    % Preallocate a cell array to store mean correlations for each level
    meanCorrelationLevels = [];
    
    % Loop through each level of subarray
    for level = 1:numel(subarray)
        array = subarray{level};
        
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

        % If you only want unique correlations, you can use the triu function
        uniqueCorrelationMatrix = triu(correlationMatrix, 1);
        uniqueCorrelations = unique(uniqueCorrelationMatrix);

        % Calculate the mean of the unique correlation coefficients
        meanUniqueCorrelation = mean(uniqueCorrelations);

        % Store the meanUniqueCorrelation in the cell array
        meanCorrelationLevels(level) = meanUniqueCorrelation;
    end
    
    % Store the meanCorrelationLevels in the meanCorrelationArray
    meanCorrelationArray{col} = meanCorrelationLevels;
    disp(mean(meanCorrelationLevels));
end

% Display the meanCorrelationArray
% disp(meanCorrelationArray);



%%
% original_consumption_ensemble = [neuron_mean_mouse{select_mouse_index, 1}(respClass_all_array_filtered{1, 1}==1, :)];




% Initialize cell arrays to store the means
block_2_ensemble_neurons_x_trial_means = cell(size(block_2_ensemble_neurons_x_trial));

% % Loop through the first level of cell array
% for i = 1:numel(original_ensemble_neurons_x_trial)
%     % Get the cell array at the current level
%     level_cell_array = original_ensemble_neurons_x_trial{i};
% 
%     % Initialize array to store means for the current level
%     level_means = zeros(size(level_cell_array{1}, 1), size(level_cell_array{1}, 2));
% 
%     % Loop through each row of the cell arrays in the current level
%     for row = 1:size(level_cell_array{1}, 1)
%         % Initialize array to store values for the current row
%         % row_values = zeros(1, numel(level_cell_array));
% 
%         % Loop through each cell array in the current level
%         for j = 1:numel(level_cell_array)
%             % Extract the current row from the current cell array
%             current_row = level_cell_array{j}(row,:);
% 
%             % Store the values for the current row
%             row_values(j,:) = current_row;
%         end
% 
%         % Store the mean of the current row
%         level_means(row, :) = mean(row_values);
%     end
% 
%     % Store the means for the current level
%     original_ensemble_neurons_x_trial_means{i} = level_means;
% end


%  THIS WORKS - BUILD OFF OF THIS CODE TOMORROW

for i = 1:numel(block_2_ensemble_neurons_x_trial)
    level_cell_array = block_2_ensemble_neurons_x_trial{i};
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
    block_2_ensemble_neurons_x_trial_means{i} = level_means;
    clear level_means
end


figure; plot(ts1, mean(block_2_ensemble_neurons_x_trial_means{1, 1}))
hold on; plot(ts1, mean(block_2_ensemble_neurons_x_trial_means{1, 2}))
hold on; plot(ts1, mean(block_2_ensemble_neurons_x_trial_means{1, 3}))


%%

% original_plus_block_2_emergence = [test; neuron_mean_mouse{select_mouse_index, 1}(respClass_all_array_filtered{1, 1}~=1 & respClass_all_array_filtered{1, 2}==1 & respClass_all_array_filtered{1, 3}~=1, :)];
% 
% 
% test = [neuron_mean_mouse{select_mouse_index, 1}(respClass_all_array_filtered{1, 1}==1 & respClass_all_array_filtered{1, 2}~=1 & respClass_all_array_filtered{1, 3}~=1, :)];
% test = [test; neuron_mean_mouse{select_mouse_index, 2}(indices_2_1,:)];
% test = [test; neuron_mean_mouse{select_mouse_index, 3}(indices_3_1,:)];
% % test = [test; neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 2}~=1 & respClass_all_array_filtered{1, 1}~=1 & respClass_all_array_filtered{1,3}~=1,:)]
% 
% 
% test = [neuron_mean_mouse{select_mouse_index, 1}(indices_2_1,:)];
% test = [test; neuron_mean_mouse{select_mouse_index, 2}(indices_2_1,:)];
% test = [test; neuron_mean_mouse{select_mouse_index, 3}(indices_2_1,:)];



% Initialize an empty array to store concatenated means
concatenated_means = [];

% Loop through each level of original_ensemble_neurons_x_trial_means
for i = 1:numel(block_2_ensemble_neurons_x_trial_means)
    % Vertically concatenate the means from the current level
    concatenated_means = vertcat(concatenated_means, block_2_ensemble_neurons_x_trial_means{i});
end

%%
% data = test;


data = concatenated_means;

alpha = 0.0001;

% Initialize matrices to store correlation coefficients and p-values
correlation_matrix = zeros(size(data, 1));
p_value_matrix = zeros(size(data, 1));

% Calculate correlation coefficients and p-values between rows
for i = 1:size(data, 1)
    for j = 1:size(data, 1)
        [corr_coeff, p_value] = corrcoef(data(i, :)', data(j, :)');
        correlation_matrix(i, j) = corr_coeff(1, 2); % Store correlation coefficient
        p_value_matrix(i, j) = p_value(1, 2); % Store p-value
    end
end

% Plot the correlation matrix
figure;
imagesc(correlation_matrix);
colorbar; % Add a colorbar to the plot
axis square; % Make the plot square for better visualization
title('Correlation Matrix');
xlabel('Trial Number');
ylabel('Trial Number');

% Show row and column indices on the plot
xticks(0:5:size(data, 1));
yticks(0:5:size(data, 1));

size_start = 0;
for ii = 1:size(block_2_ensemble_neurons_x_trial_means, 2)
    size_start = size_start + size(block_2_ensemble_neurons_x_trial_means{1, ii}, 1);
    xline(size_start, 'LineWidth', 2)
    yline(size_start, 'LineWidth', 2)
    
end

% If you want to customize the color map, you can use colormap function
% For example, using a blue-white-red colormap:
colormap(bluewhitered);

% If you want to limit the color scale to the range [0, 1]
caxis([-1 1]); % Assuming correlations range from -1 to 1

% % Display p-values as text on the plot
% for i = 1:size(data, 1)
%     for j = 1:size(data, 1)
%         text(j, i, sprintf('p = %.4f', p_value_matrix(i, j)), ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     end
% end



%%

% Preallocate a cell array to store the mean correlation for each level
meanCorrelationArray = cell(size(results));

% Loop through each column of results
for col = 1:numel(results)
    % Get the subarray from the current column
    subarray = results{col}; % Assuming the subarray is 1xN cell array
    
    % Preallocate a cell array to store mean correlations for each level
    meanCorrelationLevels = [];
    
    % Loop through each level of subarray
    for level = 1:numel(subarray)
        array = subarray{level};
        
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

        % If you only want unique correlations, you can use the triu function
        uniqueCorrelationMatrix = triu(correlationMatrix, 1);
        uniqueCorrelations = unique(uniqueCorrelationMatrix);

        % Calculate the mean of the unique correlation coefficients
        meanUniqueCorrelation = mean(uniqueCorrelations);

        % Store the meanUniqueCorrelation in the cell array
        meanCorrelationLevels(level) = meanUniqueCorrelation;
    end
    
    % Store the meanCorrelationLevels in the meanCorrelationArray
    meanCorrelationArray{col} = meanCorrelationLevels;
    disp(mean(meanCorrelationLevels));
end

% Display the meanCorrelationArray
% disp(meanCorrelationArray);
%%
% original_consumption_ensemble = [neuron_mean_mouse{select_mouse_index, 1}(respClass_all_array_filtered{1, 1}==1, :)];




% Initialize cell arrays to store the means
block_3_ensemble_neurons_x_trial_means = cell(size(block_3_ensemble_neurons_x_trial));

% % Loop through the first level of cell array
% for i = 1:numel(original_ensemble_neurons_x_trial)
%     % Get the cell array at the current level
%     level_cell_array = original_ensemble_neurons_x_trial{i};
% 
%     % Initialize array to store means for the current level
%     level_means = zeros(size(level_cell_array{1}, 1), size(level_cell_array{1}, 2));
% 
%     % Loop through each row of the cell arrays in the current level
%     for row = 1:size(level_cell_array{1}, 1)
%         % Initialize array to store values for the current row
%         % row_values = zeros(1, numel(level_cell_array));
% 
%         % Loop through each cell array in the current level
%         for j = 1:numel(level_cell_array)
%             % Extract the current row from the current cell array
%             current_row = level_cell_array{j}(row,:);
% 
%             % Store the values for the current row
%             row_values(j,:) = current_row;
%         end
% 
%         % Store the mean of the current row
%         level_means(row, :) = mean(row_values);
%     end
% 
%     % Store the means for the current level
%     original_ensemble_neurons_x_trial_means{i} = level_means;
% end


%  THIS WORKS - BUILD OFF OF THIS CODE TOMORROW

for i = 1:numel(block_3_ensemble_neurons_x_trial)
    level_cell_array = block_3_ensemble_neurons_x_trial{i};
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
    block_3_ensemble_neurons_x_trial_means{i} = level_means;
    clear level_means
end


figure; plot(ts1, mean(block_3_ensemble_neurons_x_trial_means{1, 1}))
hold on; plot(ts1, mean(block_3_ensemble_neurons_x_trial_means{1, 2}))
hold on; plot(ts1, mean(block_3_ensemble_neurons_x_trial_means{1, 3}))


%%

% original_plus_block_2_emergence = [test; neuron_mean_mouse{select_mouse_index, 1}(respClass_all_array_filtered{1, 1}~=1 & respClass_all_array_filtered{1, 2}==1 & respClass_all_array_filtered{1, 3}~=1, :)];
% 
% 
% test = [neuron_mean_mouse{select_mouse_index, 1}(respClass_all_array_filtered{1, 1}==1 & respClass_all_array_filtered{1, 2}~=1 & respClass_all_array_filtered{1, 3}~=1, :)];
% test = [test; neuron_mean_mouse{select_mouse_index, 2}(indices_2_1,:)];
% test = [test; neuron_mean_mouse{select_mouse_index, 3}(indices_3_1,:)];
% % test = [test; neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 2}~=1 & respClass_all_array_filtered{1, 1}~=1 & respClass_all_array_filtered{1,3}~=1,:)]
% 
% 
% test = [neuron_mean_mouse{select_mouse_index, 1}(indices_2_1,:)];
% test = [test; neuron_mean_mouse{select_mouse_index, 2}(indices_2_1,:)];
% test = [test; neuron_mean_mouse{select_mouse_index, 3}(indices_2_1,:)];



% Initialize an empty array to store concatenated means
concatenated_means = [];

% Loop through each level of original_ensemble_neurons_x_trial_means
for i = 1:numel(block_3_ensemble_neurons_x_trial_means)
    % Vertically concatenate the means from the current level
    concatenated_means = vertcat(concatenated_means, block_3_ensemble_neurons_x_trial_means{i});
end

%%
% data = test;


data = concatenated_means;

alpha = 0.0001;

% Initialize matrices to store correlation coefficients and p-values
correlation_matrix = zeros(size(data, 1));
p_value_matrix = zeros(size(data, 1));

% Calculate correlation coefficients and p-values between rows
for i = 1:size(data, 1)
    for j = 1:size(data, 1)
        [corr_coeff, p_value] = corrcoef(data(i, :)', data(j, :)');
        correlation_matrix(i, j) = corr_coeff(1, 2); % Store correlation coefficient
        p_value_matrix(i, j) = p_value(1, 2); % Store p-value
    end
end

% Plot the correlation matrix
figure;
imagesc(correlation_matrix);
colorbar; % Add a colorbar to the plot
axis square; % Make the plot square for better visualization
title('Correlation Matrix');
xlabel('Trial Number');
ylabel('Trial Number');

% Show row and column indices on the plot
xticks(0:5:size(data, 1));
yticks(0:5:size(data, 1));

size_start = 0;
for ii = 1:size(block_3_ensemble_neurons_x_trial_means, 2)
    size_start = size_start + size(block_3_ensemble_neurons_x_trial_means{1, ii}, 1);
    xline(size_start, 'LineWidth', 2)
    yline(size_start, 'LineWidth', 2)
    
end

% If you want to customize the color map, you can use colormap function
% For example, using a blue-white-red colormap:
colormap(bluewhitered);

% If you want to limit the color scale to the range [0, 1]
caxis([-1 1]); % Assuming correlations range from -1 to 1

% % Display p-values as text on the plot
% for i = 1:size(data, 1)
%     for j = 1:size(data, 1)
%         text(j, i, sprintf('p = %.4f', p_value_matrix(i, j)), ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     end
% end



%%

% Preallocate a cell array to store the mean correlation for each level
meanCorrelationArray = cell(size(results));

% Loop through each column of results
for col = 1:numel(results)
    % Get the subarray from the current column
    subarray = results{col}; % Assuming the subarray is 1xN cell array
    
    % Preallocate a cell array to store mean correlations for each level
    meanCorrelationLevels = [];
    
    % Loop through each level of subarray
    for level = 1:numel(subarray)
        array = subarray{level};
        
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

        % If you only want unique correlations, you can use the triu function
        uniqueCorrelationMatrix = triu(correlationMatrix, 1);
        uniqueCorrelations = unique(uniqueCorrelationMatrix);

        % Calculate the mean of the unique correlation coefficients
        meanUniqueCorrelation = mean(uniqueCorrelations);

        % Store the meanUniqueCorrelation in the cell array
        meanCorrelationLevels(level) = meanUniqueCorrelation;
    end
    
    % Store the meanCorrelationLevels in the meanCorrelationArray
    meanCorrelationArray{col} = meanCorrelationLevels;
    disp(mean(meanCorrelationLevels));
end

% Display the meanCorrelationArray
% disp(meanCorrelationArray);


