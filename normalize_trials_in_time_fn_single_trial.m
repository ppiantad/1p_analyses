function [normalized_trial_ca, concatenated_normalized_trial_ca] = normalize_trials_in_time_fn_single_trial(trial_ca)


% Initialize the new cell array to store normalized arrays
normalized_trial_ca = cell(size(trial_ca));

% Define the time vector for the normalized data
normalized_time = linspace(0, 10, 100);

% Iterate through each cell of trial_ca
for i = 1:numel(trial_ca)
    % Define the time vector for the original data
    original_time = linspace(0, (length(trial_ca{i}) - 1) * 0.1, length(trial_ca{i}));
    
    % Determine the number of points for interpolation
    num_points = numel(normalized_time);
    
    % Interpolate each array in the current cell to match the normalized_time
    normalized_trial_ca{i} = interp1(original_time, trial_ca{i}, linspace(original_time(1), original_time(end), num_points), 'linear');
end
concatenated_normalized_trial_ca = vertcat(normalized_trial_ca{:});