select_mouse = 'BLA_Insc_30';

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
for ii = 1:size(original_ensemble_neurons_x_trial_means, 2)
    size_start = size_start + size(original_ensemble_neurons_x_trial_means{1, i}, 1  );
    xline(size(original_ensemble_neurons_x_trial_means{1, i}, 1  ))
    
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
%             row_values(j
        uu = uu+1;
        % Check if p-value is less than 0.01
        if post_choice_p_value_matrix(i, j) < alpha
            if post_choice_correl_matrix(i, j) > 0
                post_choice_positive_count = post_choice_positive_count + 1;
            elseif post_choice_correl_matrix(i, j) < 0
                post_choice_negative_count = post_choice_negative_count + 1;
            end
        else
            post_choice_no_correlation_count = post_choice_no_correlation_count + 1;
        end
    end
end

post_choice_comparisons_possible = [post_choice_positive_count + post_choice_negative_count + post_choice_no_correlation_count];


disp(['Number of positive correlations: ', num2str(post_choice_positive_count)]);
disp(['Number of negative correlations: ', num2str(post_choice_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(post_choice_no_correlation_count)]);

% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
post_choice_data = [(post_choice_positive_count/post_choice_comparisons_possible)*100, (post_choice_negative_count/post_choice_comparisons_possible)*100, (post_choice_no_correlation_count/post_choice_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, post_choice_data, 'stacked');

% Add labels and title
xlabel('Counts');
ylabel('Correlation Type');
title('Correlation Counts');

% Add legend
legend(labels);

% Adjust x-axis limits
xlim([0.5, 1.5]); % since we only have one set of data, we set the limits to center the bars

% Adjust y-axis limits if needed
% ylim([0, max(data) + 10]); % adjust ylim if needed for better visualization

% Optionally, you can add data labels on each bar
% text(1:length(data), data, num2str(data'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Optionally, you can rotate x-axis labels if needed
% xticklabels(labels);

% Optionally, you can change bar colors
% colormap([0.8 0.2 0.2; 0.2 0.8 0.2; 0.2 0.2 0.8]); % customize colors as needed





%%
consumption_p_value_matrix = p_value_matrix(consumption_index, consumption_index);
consumption_correl_matrix = correlation_matrix(consumption_index, consumption_index);

n = size(consumption_index, 2); % Total number of neurons
k = 2;   % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n, k);
disp(['Number of distinct pairwise combinations: ', num2str(num_combinations)]);


% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
consumption_positive_count = 0;
consumption_negative_count = 0;
consumption_no_correlation_count = 0;

% Get the size of the correlation matrix
matrix_size = size(consumption_correl_matrix, 1);
uu = 1
% Loop over the upper triangular part of the correlation matrix
for i = 1:matrix_size
    for j = i+1:matrix_size % Start from i+1 to exclude the diagonal
        consumption_ensemble_corr_overall(uu) = consumption_correl_matrix(i, j);
        uu = uu+1;
        % Check if p-value is less than 0.01
        if consumption_p_value_matrix(i, j) < alpha
            if consumption_correl_matrix(i, j) > 0
                consumption_positive_count = consumption_positive_count + 1;
            elseif consumption_correl_matrix(i, j) < 0
                consumption_negative_count = consumption_negative_count + 1;
            end
        else
            consumption_no_correlation_count = consumption_no_correlation_count + 1;
        end
    end
end

consumption_comparisons_possible = [consumption_positive_count + consumption_negative_count + consumption_no_correlation_count];


disp(['Number of positive correlations: ', num2str(consumption_positive_count)]);
disp(['Number of negative correlations: ', num2str(consumption_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(consumption_no_correlation_count)]);

% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
consumption_data = [(consumption_positive_count/consumption_comparisons_possible)*100, (consumption_negative_count/consumption_comparisons_possible)*100, (consumption_no_correlation_count/consumption_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, consumption_data, 'stacked');

% Add labels and title
xlabel('Counts');
ylabel('Correlation Type');
title('Correlation Counts');

% Add legend
legend(labels);

% Adjust x-axis limits
xlim([0.5, 1.5]); % since we only have one set of data, we set the limits to center the bars

% Adjust y-axis limits if needed
% ylim([0, max(data) + 10]); % adjust ylim if needed for better visualization

% Optionally, you can add data labels on each bar
% text(1:length(data), data, num2str(data'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Optionally, you can rotate x-axis labels if needed
% xticklabels(labels);

% Optionally, you can change bar colors
% colormap([0.8 0.2 0.2; 0.2 0.8 0.2; 0.2 0.2 0.8]); % customize colors as needed



%%
action_post_choice_p_value_matrix = p_value_matrix(pre_choice_index, post_choice_index);
action_post_choice_correl_matrix = correlation_matrix(pre_choice_index, post_choice_index);

n1 = size(action_post_choice_p_value_matrix, 1); % Number of neurons in the first set
n2 = size(action_post_choice_p_value_matrix, 2); % Number of neurons in the second set
k = 2;    % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n1, k) * nchoosek(n2, k);
disp(['Number of unique pairwise combinations: ', num2str(num_combinations)]);



% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
action_post_choice_positive_count = 0;
action_post_choice_negative_count = 0;
action_post_choice_no_correlation_count = 0;

% % Initialize counters
% positive_count = 0;
% negative_count