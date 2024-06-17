% Assuming zall_mouse is your 11x3 cell array

% Initialize the new cell array to store the reorganized data
reorganized_data = cell(11, 3);

for col = 1:3
    for row = 1:11
        % Extract the subarray from the current cell
        subarray = zall_mouse{row, col};
        
        % Determine the number of rows in the subarray
        num_subarray_rows = size(subarray{1}, 1);
        num_subarray_cols = length(subarray);
        
        % Initialize a cell array to store the extracted rows
        extracted_rows = cell(1, num_subarray_rows);
        
        for r = 1:num_subarray_rows
            % Initialize a matrix to store the rows from each element in subarray
            temp_matrix = [];
            
            for c = 1:num_subarray_cols
                % Extract the r-th row from the c-th cell array in subarray
                temp_matrix = [temp_matrix; subarray{c}(r, :)];
            end
            
            % Store the concatenated rows into the extracted_rows cell array
            extracted_rows{r} = temp_matrix;
        end
        
        % Save the extracted_rows into the reorganized_data cell array
        reorganized_data{row, col} = extracted_rows;
    end
end

% 'reorganized_data' now contains the rearranged data

%%
% Define the time ranges
% ts1 = linspace(-2, 1, 160); % Assuming ts1 spans from -2 to 1 and has 160 points
timeRange_1 = (ts1 >= 0) & (ts1 <= 2);
timeRange_2 = (ts1 >= -1) & (ts1 <= 0);

% Initialize the final data structure
final_data = cell(size(reorganized_data));

for col = 1:3
    for row = 1:11
        % Extract the subarray from reorganized_data
        subarray = reorganized_data{row, col};
        
        % Determine the number of rows in the subarray
        num_rows = length(subarray);
        num_cells = size(subarray{1}, 1); % Number of cells in each subarray
        
        % Initialize a cell array to store the results for each row in the subarray
        result_array = cell(1, num_rows);
        
        for r = 1:num_rows
            % Extract the current row data from all cells in the subarray
            data = subarray{r};
            
            % Initialize an array to store the results for each cell in the current row
            row_results = zeros(num_cells, 1);
            
            for cell_idx = 1:num_cells
                % Extract the data for the current cell
                cell_data = data(cell_idx, :);
                
                % Calculate the mean activity in the two time ranges
                mean_activity_1 = mean(cell_data(timeRange_1));
                mean_activity_2 = mean(cell_data(timeRange_2));
                
                % Compare the mean activities
                if mean_activity_1 > mean_activity_2
                    row_results(cell_idx) = 1;
                elseif mean_activity_1 < mean_activity_2
                    row_results(cell_idx) = 2;
                else
                    row_results(cell_idx) = 0;
                end
            end
            
            % Store the row results in the result_array
            result_array{r} = row_results;
        end
        
        % Store the result_array into the final data structure
        final_data{row, col} = result_array;
    end
end

% 'final_data' now contains the comparison results

%%
% Assuming respClass_all_array_mouse is an 11x3 cell array similar to zall_mouse

% Initialize the filtered data structure
filtered_data = cell(size(final_data));
percent_data_1 = cell(size(final_data));
percent_data_all = cell(size(final_data));

for row = 1:11
    % Extract the neuron identities from the first column of respClass_all_array_mouse
    neuron_identities = respClass_all_array_mouse{row, 1};
    
    % Extract the corresponding subarray from the first column of final_data
    final_subarray = final_data{row, 1};
    
    % Initialize a cell array to store the filtered neurons for each row
    filtered_subarray = cell(size(final_subarray));
    
    for subarray_row = 1:length(final_subarray)
        % Get the current row data from final_subarray
        current_data = final_subarray{subarray_row};
        
        % Get the neuron identities for the current row
        % current_neuron_ids = neuron_identities{subarray_row};
        
        % Find the indices of neurons with an identity of 1
        selected_indices = find(neuron_identities == 1);
        
        % Extract the corresponding data from current_data
        filtered_data_row = current_data(selected_indices, :);
        
        % Store the filtered data in the filtered_subarray
        filtered_subarray{subarray_row} = filtered_data_row;

        % Calculate the percentage of 1's in the filtered data row
        if ~isempty(filtered_data_row)
            percent_1 = sum(filtered_data_row == 1) / length(filtered_data_row) * 100;
        else
            percent_1 = 0;
        end
        
        if ~isempty(current_data )
            percent_all = sum(current_data == 1) / length(current_data) * 100;
        else
            percent_all = 0;
        end
        

        % Store the percentage data in the percent_subarray
        percent_subarray(subarray_row) = percent_1;
        percent_all_subarray(subarray_row) = percent_all;
        clear percent_1 percent_all
    end
    
    % Store the filtered subarray in the filtered_data structure
    filtered_data{row, 1} = filtered_subarray;
    percent_data_1{row, 1} = percent_subarray;
    percent_data_all{row, 1} = percent_all_subarray;
    clear percent_subarray percent_all_subarray
end

% 'filtered_data' now contains the filtered neurons from the first column of final_data

%% this is the percentage of the SHK ensemble that is "activated" on each trial

percent_ensemble_activated = percent_data_1(:,1)'; 
