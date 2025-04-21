% load('BLA_C_raw_no_additional_filtering_RDT_D1_only_completed_sessions_zall_window_base_workspace_10_categories.mat')
% then run data_loop_SLEAP to get the motion data corresponding to whatever
% event you want to correlate with (let's say something like 'OMITALL', 0,
% 'BLANK_TOUCH', 0, 'BLOCK', 1, 'SHK', 0)
% then run block_wise_changes_v1.m

% Define which column of zall_mouse to use
mouse_column = 1;

% Initialize a cell array to store correlation results
% Assuming zall_mouse has the same number of columns as zall_array has cells
num_arrays = size(zall_array, 2);
num_columns = size(zall_mouse{1, mouse_column}, 2);
correlation_results = cell(num_arrays, num_columns);

% Loop through each cell in zall_array
for i = 1:size(zall_array, 2)
    % Extract the current 30x160 array from zall_array
    current_array = zall_array{1, i};
    
    % Access the corresponding row in zall_mouse
    current_mouse_row = zall_mouse{i, mouse_column};
    
    % Loop through all columns in the current cell of zall_mouse
    for j = 1:size(current_mouse_row, 2)
        % Extract the current column data
        current_mouse_data = current_mouse_row{:, j};
        
        % Calculate correlation between current_array and current_mouse_data
        % Using corrcoef which returns correlation coefficient matrix
        % We take the (1,2) element which is the correlation between the two variables
        [R, P] = corrcoef(current_array(:), current_mouse_data(:));
        
        % Store the correlation coefficient and p-value
        correlation_results{i, j} = struct('R', R(1,2), 'P', P(1,2));
    end
end

% Display results
disp('Correlation Results:');
for i = 1:size(correlation_results, 1)
    for j = 1:size(correlation_results, 2)
        fprintf('Array %d, Column %d: R = %.4f, P = %.4f\n', ...
            i, j, correlation_results{i,j}.R, correlation_results{i,j}.P);
    end
end

%%
A = zall_mouse{1, 1}{1, 1}  ;
B = zall_array{1, 1}  ;

r = zeros(30,1); % preallocate for speed and clarity
for i = 1:30
    temp = corrcoef(A(i,:), B(i,:)); % 2x2 matrix
    r(i) = temp(1,2); % extract the actual correlation coefficient
end
mean_r = mean(r);