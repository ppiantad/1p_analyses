% https://www.pnas.org/doi/10.1073/pnas.1901712116
% Goal is to reproduce something like Fig. 3 from this paper

% Assuming neuron_mean is a 600x200 matrix, and respClass_all is a 1x600 vector.
% Transpose respClass_all to match the shape of neuron_mean.
% respClass_all = respClass_all';
% 
% respClass_all = respClass_all_array{1, 4}';


% respClass_all_2 = respClass.Pre_RDT_RM.Outcome_0to2.REW_Large.activated';


% Find the rows in neuron_mean where respClass_all is equal to 1.
% activated_rows = neuron_mean(respClass_all_array{1, 1}  == 1, :);

% Calculate the time of peak activity within each row.
% [peak_values, time_of_peak_activity] = max(activated_rows, [], 2);

% respClass_all = respClass_all_array{1, 8}';
[peak_values, time_of_peak_activity] = max(zall_mean_all_array{1, 1}  , [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
neuron_mean_sorted =zall_mean_all_array{1, 2}  (sort_indices, :);


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

[r, lags] = corrcoef(zall_mean_all_array{1, 1}, zall_mean_all_array{1, 2});



[num_cells, num_samples] = size(zall_mean_all_array{1, 1});
shuffled_data = zeros(num_cells, num_samples); % Preallocate matrix for efficiency
% shift_val = randi(num_samples); % Generate a random shift value for each signal RUAIRI RECOMMENDED KEEPING THE SAME SHIFT VAL, rather than randomizing per neuron. this is because then you keep the overall correlation b/w the neurons, but disrupt the relationship to the event timestamps
for i = 1:num_cells
    shift_val = randi(num_samples); % Generate a random shift value for each signal
    shuffled_data(i,:) = circshift(zall_mean_all_array{1, 1}(i,:), shift_val,2); % Perform the circular shuffle
end
shuffled_neuron_mean_array = shuffled_data;


[peak_values, time_of_peak_activity] = max(shuffled_neuron_mean_array, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
neuron_mean_sorted = shuffled_neuron_mean_array(sort_indices, :);


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
% Assuming zall_array is your cell array

% Create a function to extract the first row of each item
extract_first_row = @(x) x(1,:);

% Use cellfun to apply the function to each cell
first_rows = cellfun(extract_first_row, zall_array, 'UniformOutput', false);

% Concatenate the first rows vertically
result = vertcat(first_rows{:});

result_sorted = result(sort_indices, :);


figure;
% Generate the heatmap
imagesc(ts1, 1, result_sorted);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');
clim([-1 1])

%%
%These data can be used to plot the median or mean choice
% time on a PCA graph, for example

behav_tbl_iter_single = behav_tbl_iter(1);

% Initialize the concatenated table
concatenatedTable = table();

% Iterate through the 3x1 cell array
for i = 1:numel(behav_tbl_iter_single)
    % Assuming each cell contains a 12x1 cell array of tables
    twelveByOneCellArray = behav_tbl_iter_single{i};
    
    % Initialize a temporary table to store the concatenated tables for this cell
    tempTable = table();
    
    % Iterate through the 12x1 cell array
    for j = 1:numel(twelveByOneCellArray)
        % Assuming each cell in the 12x1 cell array contains a table
        currentTable = twelveByOneCellArray{j};
        
        % Concatenate the current table to the temporary table vertically
        tempTable = vertcat(tempTable, currentTable);
    end
    
    % Concatenate the temporary table to the overall concatenated table vertically
    concatenatedTable = vertcat(concatenatedTable, tempTable);
end

median_choice_time_block_1 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_choice_time_block_2 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_choice_time_block_3 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));

median_collect_time_block_1 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_collect_time_block_2 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_collect_time_block_3 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));





median_start_time_from_choice = median(concatenatedTable.stTime - concatenatedTable.choiceTime);
median_collect_time_from_choice = median(concatenatedTable.collectionTime - concatenatedTable.choiceTime);



%% load Pre_RDT_RM base workspace (for particular experiment) and run initial steps of working_with_Pre_RDT_RM_data_v2
only_active_array = neuron_mean_array{1, 1}(exclusive_activated_session_1 | exclusive_activated_session_2 | exclusive_activated_session_3, :);




% Calculate the time of peak activity within each row.
% [peak_values, time_of_peak_activity] = max(activated_rows, [], 2);

% respClass_all = respClass_all_array{1, 8}';
[peak_values, time_of_peak_activity] = max(only_active_array, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
neuron_mean_sorted = only_active_array(sort_indices, :);


% Sort the rows of activated_neuron_mean based on peak_times.
% [~, sort_indices] = sort(time_of_peak_activity);
% activated_neuron_mean_sorted = activated_rows(sort_indices, :);

% Now, activated_neuron_mean_sorted contains the rows of neuron_mean filtered by respClass_all == 1
% and sorted by the time of peak activity.

figure;
% Generate the heatmap
imagesc(ts1, 1, only_active_array);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');
clim([-1 1])
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})

%% load Pre_RDT_RM base workspace (for particular experiment) and run initial steps of working_with_Pre_RDT_RM_data_v2
pre_choice_neurons = neuron_mean_array{1, 1}(exclusive_activated_session_1, :);
post_choice_reward_neurons = neuron_mean_array{1, 1}(exclusive_activated_session_2, :);
consumption_neurons = neuron_mean_array{1, 1}(exclusive_activated_session_3, :);

only_active_array_stacked = [pre_choice_neurons; post_choice_reward_neurons; consumption_neurons];

% Sort the rows of activated_neuron_mean based on peak_times.
% [~, sort_indices] = sort(time_of_peak_activity);
% activated_neuron_mean_sorted = activated_rows(sort_indices, :);

% Now, activated_neuron_mean_sorted contains the rows of neuron_mean filtered by respClass_all == 1
% and sorted by the time of peak activity.

figure;
% Generate the heatmap
imagesc(ts1, 1, only_active_array_stacked);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');
clim([-1 1])
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})


%% load Pre_RDT_RM base workspace (for particular experiment) and run initial steps of working_with_Pre_RDT_RM_data_v2
pre_choice_neurons = neuron_mean_array{1, 1}(exclusive_activated_session_1, :);
post_choice_reward_neurons = neuron_mean_array{1, 1}(exclusive_activated_session_2, :);
consumption_neurons = neuron_mean_array{1, 1}(exclusive_activated_session_3, :);

only_active_array_stacked = [pre_choice_neurons; post_choice_reward_neurons; consumption_neurons];

% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(pre_choice_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
pre_choice_neurons_sorted = pre_choice_neurons(sort_indices, :);



% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(post_choice_reward_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
post_choice_reward_neurons_sorted = post_choice_reward_neurons(sort_indices, :);


% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(consumption_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
consumption_neurons_sorted = consumption_neurons(sort_indices, :);

sorted_only_active_array_stacked = [pre_choice_neurons_sorted; post_choice_reward_neurons_sorted; consumption_neurons_sorted];

% Now, activated_neuron_mean_sorted contains the rows of neuron_mean filtered by respClass_all == 1
% and sorted by the time of peak activity.

figure;
% Generate the heatmap
imagesc(ts1, 1, sorted_only_active_array_stacked);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');
clim([-1 1])
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})