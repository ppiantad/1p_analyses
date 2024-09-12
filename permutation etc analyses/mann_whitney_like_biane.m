% Initialize arrays to store p-values
p_values = zeros(1, 160);

% Loop through all 160 columns
for col = 1:160
    % Extract the data for the current column from both cells
    data1 = mean_data_array{1}(:, col);
    data2 = mean_data_array{2}(:, col);
    
    % Perform a two-sided Mann–Whitney U-test
    p = ranksum(data1, data2);
    
    % Store the p-value
    p_values(col) = p;
end

% Display p-values
disp('P-values for each column:');
disp(p_values);
