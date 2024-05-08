select_mouse = 'BLA_Insc_24';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';

select_mouse_cells_indices = find(strcmp(mouse_cells(1,:), select_mouse));
% Assuming respClass_all_array is your 1x3 cell array
% respClass_all_array = cell(1, 3); % Initialize filtered array
% 
% for i = 1:3
%     respClass_all_array{i} = respClass_all_array{i}(select_mouse_cells_indices);
% end


test = [neuron_mean_array{1, 1}(respClass_all_array{1, 1}==1 & respClass_all_array{1, 2}~=1 & respClass_all_array{1, 3}~=1, :)];
test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 1}~=1 & respClass_all_array{1, 2}==1 & respClass_all_array{1, 3}~=1, :)];
test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 1}~=1 & respClass_all_array{1, 2}~=1 & respClass_all_array{1, 3}==1, :)];
% test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 2}~=1 & respClass_all_array{1, 1}~=1 & respClass_all_array{1,3}~=1,:)]

pre_choice_index = [1:sum(respClass_all_array{1, 1}==1 & respClass_all_array{1, 2}~=1 & respClass_all_array{1, 3}~=1)];
post_choice_index = [pre_choice_index(end)+1:pre_choice_index(end)+sum(respClass_all_array{1, 1}~=1 & respClass_all_array{1, 2}==1 & respClass_all_array{1, 3}~=1)];
consumption_index = [post_choice_index(end)+1:post_choice_index(end)+sum(respClass_all_array{1, 1}~=1 & respClass_all_array{1, 2}~=1 & respClass_all_array{1, 3}==1)];
neutral_index = [consumption_index(end)+1:consumption_index(end)+sum(respClass_all_array{1, 2}~=1 & respClass_all_array{1, 1}~=1 & respClass_all_array{1,3}~=1)];

% tabulate how neurons assigned to neuron_mean_array for the 1st event
% % change across subsequent events
% test = [neuron_mean_array{1, 1}(respClass_all_array{1, 1}==1, :)];
% test = [test; neuron_mean_array{1, 2}(respClass_all_array{1, 1}==1, :)];
% test = [test; neuron_mean_array{1, 3}(respClass_all_array{1, 1}==1, :)];


piechart_data = [sum_activated_percent 100-sum(sum_activated_percent)];
figure; piechart(piechart_data)

%%
data = test;

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
xlabel('Neuron Number');
ylabel('Neuron Number');

% Show row and column indices on the plot
xticks(0:50:size(data, 1));
yticks(0:50:size(data, 1));

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
action_p_value_matrix = p_value_matrix(pre_choice_index, pre_choice_index);
action_correl_matrix = correlation_matrix(pre_choice_index, pre_choice_index);

n = size(pre_choice_index, 2); % Total number of neurons
k = 2;   % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n, k);
disp(['Number of distinct pairwise combinations: ', num2str(num_combinations)]);


% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
action_positive_count = 0;
action_negative_count = 0;
action_no_correlation_count = 0;

% Get the size of the correlation matrix
matrix_size = size(action_correl_matrix, 1);
uu = 1;
% Loop over the upper triangular part of the correlation matrix
for i = 1:matrix_size
    for j = i+1:matrix_size % Start from i+1 to exclude the diagonal
        % get total pairwise correlations between all possible combos of
        % neurons
        action_ensemble_corr_overall(uu) = action_correl_matrix(i, j);
        uu = uu+1;
        % Check if p-value is less than 0.01
        if action_p_value_matrix(i, j) < alpha
            if action_correl_matrix(i, j) > 0
                action_positive_count = action_positive_count + 1;
            elseif action_correl_matrix(i, j) < 0
                action_negative_count = action_negative_count + 1;
            end
        else
            action_no_correlation_count = action_no_correlation_count + 1;
        end
    end
end

action_comparisons_possible = [action_positive_count + action_negative_count + action_no_correlation_count];


disp(['Number of positive correlations: ', num2str(action_positive_count)]);
disp(['Number of negative correlations: ', num2str(action_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(action_no_correlation_count)]);

% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
action_data = [(action_positive_count/action_comparisons_possible)*100, (action_negative_count/action_comparisons_possible)*100, (action_no_correlation_count/action_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, action_data, 'stacked');

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


% uu = 1
% for i = 1:matrix_size
%     for j = i+1:matrix_size % Start from i+1 to exclude the diagonal
%         % Check if p-value is less than 0.01
%         action_ensemble_corr_overall(uu) = action_correl_matrix(i, j);
%         uu = uu+1;
%     end
% end


%%
post_choice_p_value_matrix = p_value_matrix(post_choice_index, post_choice_index);
post_choice_correl_matrix = correlation_matrix(post_choice_index, post_choice_index);

n = size(post_choice_index, 2); % Total number of neurons
k = 2;   % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n, k);
disp(['Number of distinct pairwise combinations: ', num2str(num_combinations)]);


% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
post_choice_positive_count = 0;
post_choice_negative_count = 0;
post_choice_no_correlation_count = 0;

% Get the size of the correlation matrix
matrix_size = size(post_choice_correl_matrix, 1);
uu = 1
% Loop over the upper triangular part of the correlation matrix
for i = 1:matrix_size
    for j = i+1:matrix_size % Start from i+1 to exclude the diagonal
        post_choice_ensemble_corr_overall(uu) = post_choice_correl_matrix(i, j);
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
% negative_count = 0;
% no_correlation_count = 0;

% Get the size of the correlation matrix
[num_neurons_1, num_neurons_2] = size(action_post_choice_correl_matrix);

% Loop through the matrices to count correlations
for i = 1:num_neurons_1
    for j = 1:num_neurons_2
        correlation = action_post_choice_correl_matrix(i, j);
        p_value = action_post_choice_p_value_matrix(i, j);
        if p_value < alpha
            if correlation > 0
                action_post_choice_positive_count = action_post_choice_positive_count + 1;
            elseif correlation < 0
                action_post_choice_negative_count = action_post_choice_negative_count + 1;
            end
        else
            action_post_choice_no_correlation_count = action_post_choice_no_correlation_count + 1;
        end
    end
end


disp(['Number of positive correlations: ', num2str(action_post_choice_positive_count)]);
disp(['Number of negative correlations: ', num2str(action_post_choice_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(action_post_choice_no_correlation_count)]);

action_post_choice_comparisons_possible = [action_post_choice_positive_count+action_post_choice_negative_count+action_post_choice_no_correlation_count];


% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
action_post_choice_data = [(action_post_choice_positive_count/action_post_choice_comparisons_possible)*100, (action_post_choice_negative_count/action_post_choice_comparisons_possible)*100, (action_post_choice_no_correlation_count/action_post_choice_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, action_post_choice_data, 'stacked');

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
action_consumption_p_value_matrix = p_value_matrix(pre_choice_index, consumption_index);
action_consumption_correl_matrix = correlation_matrix(pre_choice_index, consumption_index);

n1 = size(action_consumption_p_value_matrix, 1); % Number of neurons in the first set
n2 = size(action_consumption_p_value_matrix, 2); % Number of neurons in the second set
k = 2;    % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n1, k) * nchoosek(n2, k);
disp(['Number of unique pairwise combinations: ', num2str(num_combinations)]);



% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
action_consumption_positive_count = 0;
action_consumption_negative_count = 0;
action_consumption_no_correlation_count = 0;

% % Initialize counters
% positive_count = 0;
% negative_count = 0;
% no_correlation_count = 0;

% Get the size of the correlation matrix
[num_neurons_1, num_neurons_2] = size(action_consumption_correl_matrix);

% Loop through the matrices to count correlations
for i = 1:num_neurons_1
    for j = 1:num_neurons_2
        correlation = action_consumption_correl_matrix(i, j);
        p_value = action_consumption_p_value_matrix(i, j);
        if p_value < alpha
            if correlation > 0
                action_consumption_positive_count = action_consumption_positive_count + 1;
            elseif correlation < 0
                action_consumption_negative_count = action_consumption_negative_count + 1;
            end
        else
            action_consumption_no_correlation_count = action_consumption_no_correlation_count + 1;
        end
    end
end


disp(['Number of positive correlations: ', num2str(action_consumption_positive_count)]);
disp(['Number of negative correlations: ', num2str(action_consumption_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(action_consumption_no_correlation_count)]);

action_consumption_comparisons_possible = [action_consumption_positive_count+action_consumption_negative_count+action_consumption_no_correlation_count];


% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
action_consumption_data = [(action_consumption_positive_count/action_consumption_comparisons_possible)*100, (action_consumption_negative_count/action_consumption_comparisons_possible)*100, (action_consumption_no_correlation_count/action_consumption_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, action_consumption_data, 'stacked');

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
post_choice_consumption_p_value_matrix = p_value_matrix(post_choice_index, consumption_index);
post_choice_consumption_correl_matrix = correlation_matrix(post_choice_index, consumption_index);

n1 = size(post_choice_consumption_p_value_matrix, 1); % Number of neurons in the first set
n2 = size(post_choice_consumption_p_value_matrix, 2); % Number of neurons in the second set
k = 2;    % Number of neurons chosen for pairwise combinations

num_combinations = nchoosek(n1, k) * nchoosek(n2, k);
disp(['Number of unique pairwise combinations: ', num2str(num_combinations)]);



% Assuming action_correl_matrix and action_p_value_matrix are your matrices

% Initialize counters
post_choice_consumption_positive_count = 0;
post_choice_consumption_negative_count = 0;
post_choice_consumption_no_correlation_count = 0;

% % Initialize counters
% positive_count = 0;
% negative_count = 0;
% no_correlation_count = 0;

% Get the size of the correlation matrix
[num_neurons_1, num_neurons_2] = size(post_choice_consumption_correl_matrix);

% Loop through the matrices to count correlations
for i = 1:num_neurons_1
    for j = 1:num_neurons_2
        correlation = post_choice_consumption_correl_matrix(i, j);
        p_value = post_choice_consumption_p_value_matrix(i, j);
        if p_value < alpha
            if correlation > 0
                post_choice_consumption_positive_count = post_choice_consumption_positive_count + 1;
            elseif correlation < 0
                post_choice_consumption_negative_count = post_choice_consumption_negative_count + 1;
            end
        else
            post_choice_consumption_no_correlation_count = post_choice_consumption_no_correlation_count + 1;
        end
    end
end


disp(['Number of positive correlations: ', num2str(post_choice_consumption_positive_count)]);
disp(['Number of negative correlations: ', num2str(post_choice_consumption_negative_count)]);
disp(['Number of no correlations (p-value > ', num2str(alpha), '): ', num2str(post_choice_consumption_no_correlation_count)]);

post_choice_consumption_comparisons_possible = [post_choice_consumption_positive_count+post_choice_consumption_negative_count+post_choice_consumption_no_correlation_count];


% Assuming you have positive_count, negative_count, and no_correlation_count variables

% Define data for the stacked bar plot
post_choice_consumption_data = [(post_choice_consumption_positive_count/post_choice_consumption_comparisons_possible)*100, (post_choice_consumption_negative_count/post_choice_consumption_comparisons_possible)*100, (post_choice_consumption_no_correlation_count/post_choice_consumption_comparisons_possible)*100];
figure;
% Define labels for the bars
labels = {'Positive Correlation', 'Negative Correlation', 'No Correlation'};

% Create the stacked bar plot
bar(1, post_choice_consumption_data, 'stacked');

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

% Assuming you have action_data, consumption_data, and action_consumption_data matrices

% Define x-values for each set of bars
x = 1:6;

% Create a new figure
figure;

% Plot the stacked bar plots for each dataset
% bar(x, [action_data; post_choice_data; consumption_data; action_post_choice_data; action_consumption_data; post_choice_consumption_data], 'stacked');

barh(x, [action_data; post_choice_data; consumption_data; action_post_choice_data; action_consumption_data; post_choice_consumption_data], 'stacked');
% % Set x-axis tick locations and labels
% xticks(x);
% xticklabels({'aa', 'cc', 'ac'});
% ylabel('% Neuron Pairs');
% title('Title');

% Add legend
legend('Positive correlation', 'Negative correlation', 'No sig correlation');


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

%%

figure; plot(ts1, nanmean(neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 1}==1, :)), 'color', "#D95319")
hold on; plot(ts1, nanmean(neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 2}==1, :)), 'color',  "blue")
hold on; plot(ts1, nanmean(neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 2}==3 & respClass_all_array_filtered{1, 1}==3, :)), 'color',  "black")
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')

%%
figure;
shadedErrorBar(ts1, nanmean(neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 1}==1, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array_filtered{1, 1}==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 2}==1, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array_filtered{1, 2}==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 3}==1, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array_filtered{1, 3}==1, :)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')


%% 
action_matrix = neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 1}==1, :);

consumption_matrix = neuron_mean_mouse{1, 1}(respClass_all_array_filtered{1, 2}==1, :);


% Assuming you have action_matrix and consumption_matrix

% Get the number of columns in each matrix
num_columns_action = size(action_matrix, 1);
num_columns_consumption = size(consumption_matrix, 1);


action_matrix = action_matrix(1:num_columns_consumption, :);

for corr = 1:size(action_matrix, 2)
    
    [corr_coeff, p_value] = corrcoef(action_matrix(:, corr), consumption_matrix(:, corr));
    time_corr(corr) = corr_coeff(1, 2);


end



figure; plot(ts1, time_corr);
hold on; 
plot(ts1, mean(action_matrix))
plot(ts1, mean(consumption_matrix))