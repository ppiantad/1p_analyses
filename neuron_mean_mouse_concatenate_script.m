% Initialize a cell array to store the concatenated data
concatenated_data = cell(12, 1);

% Loop through each row of the neuron_mean_mouse cell array
for i = 1:size(neuron_mean_mouse, 1)
    % Extract the data from the current row for columns {1,1}, {1,2}, and {1,3}
    current_row_data = neuron_mean_mouse{i, 1};
    current_row_data = vertcat(current_row_data, neuron_mean_mouse{i, 2});
    current_row_data = vertcat(current_row_data, neuron_mean_mouse{i, 3});
    
    % Store the concatenated data in the result cell array
    concatenated_data{i} = current_row_data;
end