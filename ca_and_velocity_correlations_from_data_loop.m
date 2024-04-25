% run ca_correlate_session_long_velocity_v1 first! 

% for a single mouse
% Assuming zall_array{1,1} and motion_TraceTrials_mouse{1,1} are your matrices
num_rows = size(zall_array{1,1}, 1); % Assuming both matrices have the same number of rows


correlation_coefficients = zeros(1, num_rows);

for i = 1:num_rows
    % Get the i-th row from both matrices
    row_zall = zall_array{1,1}(i, :);
    row_motion = motion_TraceTrials_mouse{1,1}(i, :);
    
    % Calculate the correlation coefficient between the two rows
    correlation_coefficients(i) = corr(row_zall', row_motion');
end


% Assuming zall_array{1,1} and motion_TraceTrials_mouse{1,1} are your matrices
num_rows = size(zall_array{1,1}, 1); % Assuming both matrices have the same number of rows

% Preallocate correlation matrix
correlation_matrix = zeros(num_rows);

for i = 1:num_rows
    % Get the i-th row from both matrices
    row_zall = zall_array{1,1}(i, :);
    row_motion = motion_TraceTrials_mouse{1,1}(i, :);
    
    % Calculate the correlation coefficient between the two rows
    correlation_matrix(i, :) = corr(row_zall', row_motion');
end

% % Plot heatmap
% figure;
% heatmap(correlation_matrix);
% title('Correlation Matrix Heatmap');
% xlabel('Neuron');
% ylabel('Neuron');

% Assuming correlation_coefficients contains the correlation coefficients

% Plot histogram
figure;
histogram(correlation_coefficients, 'Normalization', 'probability');
title('Distribution of Correlation Coefficients for Single Mouse - Trial x Trial');
xlabel('Correlation Coefficient');
ylabel('Probability');


%%
% Assuming zall_mouse is your cell array containing data for each mouse
num_mice = size(zall_mouse, 1); % Assuming zall_mouse is a column cell array

% calculate correlations for trial-based activity

% Initialize a cell array to store correlation coefficients
correlation_coefficients_mouse = cell(num_mice, 1);

for mouse_idx = 1:num_mice
    % Get the cell array containing data for the current mouse
    mouse_data = zall_mouse{mouse_idx};
    num_cells = size(mouse_data, 2);
    
    % Initialize a cell array to store correlation coefficients for the current mouse
    correlation_coefficients_cell = cell(num_cells, 1);
    
    % Loop through each cell's activity data
    for cell_idx = 1:num_cells
        % Get the activity data for the current cell
        cell_activity = mouse_data{cell_idx};
        
        % Correlate the activity of the current cell with motion data for the current mouse
        correlation_coefficients_cell{cell_idx} = zeros(size(cell_activity, 1), 1); % Initialize
        
        for trial_idx = 1:size(cell_activity, 1)
            % Get the activity data for the current trial
            cell_activity_trial = cell_activity(trial_idx, :);
            motion_trial = motion_TraceTrials_mouse{mouse_idx}(trial_idx, :);
            
            % Calculate the correlation coefficient
            correlation_coefficients_cell{cell_idx}(trial_idx) = corr(cell_activity_trial', motion_trial');
        end
    end
    
    % Store correlation coefficients for the current mouse
    correlation_coefficients_mouse{mouse_idx} = correlation_coefficients_cell;
end

% Now correlation_coefficients_mouse contains the correlation coefficients for each cell of each mouse

%%
% Initialize an array to store mean correlations for each cell
mean_correlation_per_cell = [];

for mouse_idx = 1:num_mice
    correlation_coefficients_cell = correlation_coefficients_mouse{mouse_idx};
    num_cells = numel(correlation_coefficients_cell);
    
    for cell_idx = 1:num_cells
        % Get correlation coefficients for the current cell
        cell_correlation_coefficients = correlation_coefficients_cell{cell_idx};
        
        % Calculate the mean correlation for the current cell
        mean_correlation = mean(cell_correlation_coefficients);
        
        % Store the mean correlation for the current cell
        mean_correlation_per_cell(end + 1) = mean_correlation;
    end
end

% Plot histogram
figure;
histogram(mean_correlation_per_cell, 'Normalization', 'probability');
title('Distribution of Correlation Coefficients Across Mice');
xlabel('Correlation Coefficient');
ylabel('Probability');
