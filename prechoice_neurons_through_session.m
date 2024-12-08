

mouse_to_check = 10; 

prechoice_block_2_3_column_1_data = zall_mouse{mouse_to_check, 1}(prechoice_blocks_2_and_3_mouse{mouse_to_check, 1} == 1);
prechoice_block_2_3_column_2_data = zall_mouse{mouse_to_check, 5} (prechoice_blocks_2_and_3_mouse{mouse_to_check, 1} == 1);

behav_part_1 = behav_tbl_iter{1, 1}{mouse_to_check, 1}  
behav_part_2 = behav_tbl_iter{5, 1}{mouse_to_check, 1} 
concat_behav = vertcat(behav_part_1, behav_part_2 )

% Initialize the concatenated cell array
concatenated_columns = cell(1, size(prechoice_block_2_3_column_1_data, 2));


% Iterate through each column and concatenate the data
for i = 1:size(prechoice_block_2_3_column_1_data, 2)
    concatenated_columns{i} = vertcat(prechoice_block_2_3_column_1_data{i}, ...
                                      prechoice_block_2_3_column_2_data{i});
end

% Initialize an array to store the resulting mean values
[numRows, numCols] = size(concatenated_columns{1}); % Assume all cells have the same dimensions
mean_array = zeros(numRows, numCols); % Resulting double array

% Loop through each cell array
for i = 1:length(concatenated_columns)
    % Add the values from the current cell array to the mean calculation
    mean_array = mean_array + concatenated_columns{i};
end

% Divide by the number of cell arrays to calculate the mean
mean_array = mean_array / length(concatenated_columns);

figure; imagesc(ts1, [], mean_array)

figure; plot(concat_behav.bigSmall)