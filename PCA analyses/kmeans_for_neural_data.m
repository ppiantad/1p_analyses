

% data_for_clustering = zall_mean_all_array{1,1};

% data_for_clustering = neuron_mean_array{1, 1}  ;

% data_for_clustering = zscore(data_for_clustering, 0, 2);

data_for_clustering = zall_mean_all_array{1, 1};  

% if you want to pre-specify subwindows to do the clustering on, as is done
% here: https://www.biorxiv.org/content/10.1101/2024.06.26.600895v4
% run data_loop on the relevant data (e.g., Pre_RDT_RM, 'OMIT_ALL', 0,
% 'BLANK_TOUCH', 0)
prechoice_means = mean(zall_mean_all_array{1, 1}(:, ts1 >= -4 & ts1 <= 0), 2);
postchoice_means = mean(zall_mean_all_array{1, 1}(:, ts1 >= 0 & ts1 <= 2), 2);
% also run data_loop for collectionTime
consumption_means = mean(zall_mean_all_array{1, 2}(:, ts1 >= 1 & ts1 <= 3), 2);

all_means = [prechoice_means postchoice_means consumption_means];

data_for_clustering = all_means;




% if you want to pre-specify subwindows to do the clustering on, as is done
% here: https://www.biorxiv.org/content/10.1101/2024.06.26.600895v4
% run data_loop on the relevant data (e.g., Pre_RDT_RM, 'OMIT_ALL', 0,
% 'BLANK_TOUCH', 0)
prechoice_safe_means = mean(zall_mean_all_array{1, 1}(:, ts1 >= -4 & ts1 <= 0), 2);
postchoice_safe_means = mean(zall_mean_all_array{1, 1}(:, ts1 >= 0 & ts1 <= 2), 2);
% also run data_loop for collectionTime
consumption_safe_means = mean(zall_mean_all_array{1, 2}(:, ts1 >= 1 & ts1 <= 3), 2);

prechoice_risky_means = mean(zall_mean_all_array{1, 3}(:, ts1 >= -4 & ts1 <= 0), 2);
postchoice_risky_means = mean(zall_mean_all_array{1, 3}(:, ts1 >= 0 & ts1 <= 2), 2);
% also run data_loop for collectionTime
consumption_risky_means = mean(zall_mean_all_array{1, 4}(:, ts1 >= 1 & ts1 <= 3), 2);

shk_risky_means = mean(zall_mean_all_array{1, 5}(:, ts1 >= 0 & ts1 <= 2), 2);





all_means = [prechoice_safe_means postchoice_safe_means consumption_safe_means prechoice_risky_means postchoice_risky_means consumption_risky_means shk_risky_means];

data_for_clustering = all_means;

data_for_display = horzcat(zall_mean_all_array{1, 1}, zall_mean_all_array{1, 3}, zall_mean_all_array{1, 5});



all_means = [prechoice_safe_means postchoice_safe_means consumption_safe_means prechoice_risky_means postchoice_risky_means consumption_risky_means];
data_for_clustering = all_means;

data_for_display = horzcat(zall_mean_all_array{1, 1}, zall_mean_all_array{1, 3});


%% from Sean to do kmeans on traces

clusters_desired = 9

[kmeans_idx,C,sumdist3] = kmeans(data_for_clustering,clusters_desired,'Distance','correlation','Display','final', 'Replicates', 200,'Start','uniform');


% for ii = 1:clusters_desired
%     figure;
%     plot(ts1, data_for_clustering(kmeans_idx == ii, :));
% 
% end



figure;
for ii = 1:clusters_desired

    plot(ts1, mean(zall_mean_all_array{1, 1}(kmeans_idx == ii, :)));
    hold on;
    legend;
end




figure;
for ii = 1:clusters_desired

    plot(ts1, mean(data_for_clustering(kmeans_idx == ii, :)));
    hold on;
    legend;
end



%% big heatmap for categories



cluster_neurons_neurons_sorted = [];

for ii = 1:clusters_desired
    clustered_neurons = zall_mean_all_array{1, 1}(kmeans_idx == ii, :);
    % Sort the rows of activated_neuron_mean based on peak_times.
    [peak_values, time_of_peak_activity] = max(clustered_neurons, [], 2);
    [~, sort_indices] = sort(time_of_peak_activity);
    neurons_by_cluster{ii} = clustered_neurons;
    cluster_neurons_neurons_sorted = [cluster_neurons_neurons_sorted; clustered_neurons(sort_indices, :)];
end


neurons_sorted_by_cluster = [];

for ii = 1:clusters_desired
    clustered_neurons = data_for_display(kmeans_idx == ii, :);
    % Sort the rows of activated_neuron_mean based on peak_times.
    % [peak_values, time_of_peak_activity] = max(clustered_neurons, [], 2);
    % [~, sort_indices] = sort(time_of_peak_activity);
    % neurons_by_cluster{ii} = clustered_neurons;
    neurons_by_cluster{ii} = clustered_neurons;
    neurons_sorted_by_cluster = [neurons_sorted_by_cluster; clustered_neurons];
end




% 
% cluster_1_neurons = data_for_clustering(kmeans_idx == 1, :);
% cluster_2_neurons = data_for_clustering(kmeans_idx == 2, :);
% cluster_3_neurons = data_for_clustering(kmeans_idx == 3, :);
% cluster_4_neurons = data_for_clustering(kmeans_idx == 4, :);
% 
% 
% 
% % Sort the rows of activated_neuron_mean based on peak_times.
% [peak_values, time_of_peak_activity] = max(cluster_1_neurons, [], 2);
% [~, sort_indices] = sort(time_of_peak_activity);
% cluster_1_neurons_neurons_sorted = cluster_1_neurons(sort_indices, :);
% 
% 
% 
% % Sort the rows of activated_neuron_mean based on peak_times.
% [peak_values, time_of_peak_activity] = max(cluster_2_neurons, [], 2);
% [~, sort_indices] = sort(time_of_peak_activity);
% cluster_2_neurons_neurons_sorted = cluster_2_neurons(sort_indices, :);
% 
% 
% % Sort the rows of activated_neuron_mean based on peak_times.
% [peak_values, time_of_peak_activity] = max(cluster_3_neurons, [], 2);
% [~, sort_indices] = sort(time_of_peak_activity);
% cluster_3_neurons_neurons_sorted = cluster_3_neurons(sort_indices, :);
% 
% 
% % Sort the rows of activated_neuron_mean based on peak_times.
% [peak_values, time_of_peak_activity] = max(cluster_4_neurons, [], 2);
% [~, sort_indices] = sort(time_of_peak_activity);
% cluster_4_neurons_neurons_sorted = cluster_4_neurons(sort_indices, :);
% 
% sorted_only_active_array_stacked = [cluster_1_neurons_neurons_sorted; cluster_2_neurons_neurons_sorted; cluster_3_neurons_neurons_sorted; cluster_4_neurons_neurons_sorted];
% 


%%

figure;
% Generate the heatmap
imagesc(ts1, 1, neurons_sorted_by_cluster);

initial_row = 0

for dd = 1:clusters_desired
    current_cluster = neurons_by_cluster{dd};
    
    row_size = size(current_cluster, 1)

    [xl,xt] = xlin(-1,'Cluster', initial_row,  row_size +initial_row, row_size +initial_row);

    initial_row = row_size +initial_row + 1;

end


% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');
clim([-1 1])
xline(0);


% If you want to customize the color map, you can use colormap function
% For example, using a blue-white-red colormap:
colormap(gray);

% If you want to limit the color scale to the range [0, 1]
caxis([-1 1]); % Assuming correlations range from -1 to 1
% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

%%
figure;

% Create a tiled layout with enough rows for all clusters
tiledlayout(clusters_desired, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% === Heatmap (First Column, Spanning All Rows) ===
ax1 = nexttile([clusters_desired, 1]); % Span all rows
imagesc(ts1, 1, neurons_sorted_by_cluster);
hold on;

initial_row = 0;
for dd = 1:clusters_desired
    current_cluster = neurons_by_cluster{dd};
    row_size = size(current_cluster, 1);
    
    % Draw separator lines
    xlin(-1, 'Cluster', initial_row, row_size + initial_row, row_size + initial_row);
    
    initial_row = row_size + initial_row + 1;
end

colorbar;
xlabel('Time (s)');
ylabel('Neuron');
set(gca, 'YDir', 'reverse');
clim([-1 1]);
xline(0);
colormap(gray);
caxis([-1 1]);

% === Mean Lines (Second Column, One Per Cluster) ===
for dd = 1:clusters_desired
    nexttile;
    
    % Compute mean trace for the current cluster
    current_cluster = neurons_by_cluster{dd};
    mean_trace = mean(current_cluster, 1);
    
    % Normalize values for visibility
    min_val = min(mean_trace);
    max_val = max(mean_trace);
    if max_val > min_val
        normalized_mean_trace = 2 * (mean_trace - min_val) / (max_val - min_val) - 1;
    else
        normalized_mean_trace = zeros(size(mean_trace)); % Avoid divide by zero
    end
    
    % Plot the mean trace
    % plot(1:480, normalized_mean_trace, 'k', 'LineWidth', 1.5);
    plot(1:320, normalized_mean_trace, 'k', 'LineWidth', 1.5);
    % Format subplot
    ylim([-1, 1]); % Keep the scale consistent
    % xlim([1, 480]);
    xlim([1, 320]);
    title(['Cluster ' num2str(dd)]);
    set(gca, 'XTick', [], 'YTick', []); % Hide ticks for clean look
end

