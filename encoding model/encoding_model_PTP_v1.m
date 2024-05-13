 ', num2str(alpha), '): ', num2str(consumption_no_correlation_count)]);

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

median_collect_time_block_1 = median(concatenatedTable.col