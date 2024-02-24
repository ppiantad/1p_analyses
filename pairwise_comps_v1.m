test = [neuron_mean_array{1, 1}(respClass_all_array{1, 1}==1, :)];
test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 2}==1, :)];
test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 2}==3 & respClass_all_array{1, 1}==3, :)]


% % tabulate how neurons assigned to neuron_mean_array for the 1st event
% % change across subsequent events
% test = [neuron_mean_array{1, 1}(respClass_all_array{1, 1}==1, :)];
% test = [test; neuron_mean_array{1, 2}(respClass_all_array{1, 1}==1, :)];
% test = [test; neuron_mean_array{1, 3}(respClass_all_array{1, 1}==1, :)];


piechart_data = [sum_activated_percent 100-sum(sum_activated_percent)];
figure; piechart(piechart_data)

%%
data = test;

% Initialize a matrix to store correlation coefficients
correlation_matrix = zeros(size(data, 1));

% Calculate correlation coefficients between rows
for i = 1:size(data, 1)
    for j = 1:size(data, 1)
        correlation_matrix(i, j) = corr(data(i, :)', data(j, :)');
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
caxis([-1 1]); % Assu

%%
%These data can be used to plot the median or mean choice
% time on a PCA graph, for example

behav_tbl_iter = behav_tbl_iter(1);

% Initialize the concatenated table
concatenatedTable = table();

% Iterate through the 3x1 cell array
for i = 1:numel(behav_tbl_iter)
    % Assuming each cell contains a 12x1 cell array of tables
    twelveByOneCellArray = behav_tbl_iter{i};
    
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

median_collect_time_from_choice = median(concatenatedTable.collectionTime - concatenatedTable.choiceTime);

%%

figure; plot(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 1}==1, :)), 'color', "#D95319")
hold on; plot(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 2}==1, :)), 'color',  "blue")
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'consumption active'}, 'Location','northwest')

%%
figure;
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 1}==1, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array{1, 1}==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 2}==1, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array{1, 2}==1, :)), 'lineProps', {'color', batlowW(iter,:)});
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'consumption active'}, 'Location','northwest')