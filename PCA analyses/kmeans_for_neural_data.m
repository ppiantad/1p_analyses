

% data_for_clustering = zall_mean_all_array{1,1};

% data_for_clustering = neuron_mean_array{1, 1}  ;

% data_for_clustering = zscore(data_for_clustering, 0, 2);

data_for_clustering = zall_mean_all_array{1, 1};  

%% from Sean to do kmeans on traces

clusters_desired = 5

[kmeans_idx,C,sumdist3] = kmeans(data_for_clustering,clusters_desired,'Distance','correlation','Display','final', 'Replicates', 200,'Start','uniform');


% for ii = 1:clusters_desired
%     figure;
%     plot(ts1, data_for_clustering(kmeans_idx == ii, :));
% 
% end


figure;
for ii = 1:clusters_desired

    plot(ts1, mean(data_for_clustering(kmeans_idx == ii, :)));
    hold on;
    legend;
end



%% big heatmap for categories


cluster_1_neurons = data_for_clustering(kmeans_idx == 1, :);
cluster_2_neurons = data_for_clustering(kmeans_idx == 2, :);
cluster_3_neurons = data_for_clustering(kmeans_idx == 3, :);
cluster_4_neurons = data_for_clustering(kmeans_idx == 4, :);

cluster_neurons_neurons_sorted = [];

for ii = 1:clusters_desired
    clustered_neurons = data_for_clustering(kmeans_idx == ii, :);
    % Sort the rows of activated_neuron_mean based on peak_times.
    [peak_values, time_of_peak_activity] = max(clustered_neurons, [], 2);
    [~, sort_indices] = sort(time_of_peak_activity);
    neurons_by_cluster{ii} = clustered_neurons;
    cluster_neurons_neurons_sorted = [cluster_neurons_neurons_sorted; clustered_neurons(sort_indices, :)];
end


% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(cluster_1_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
cluster_1_neurons_neurons_sorted = cluster_1_neurons(sort_indices, :);



% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(cluster_2_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
cluster_2_neurons_neurons_sorted = cluster_2_neurons(sort_indices, :);


% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(cluster_3_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
cluster_3_neurons_neurons_sorted = cluster_3_neurons(sort_indices, :);


% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(cluster_4_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
cluster_4_neurons_neurons_sorted = cluster_4_neurons(sort_indices, :);

sorted_only_active_array_stacked = [cluster_1_neurons_neurons_sorted; cluster_2_neurons_neurons_sorted; cluster_3_neurons_neurons_sorted; cluster_4_neurons_neurons_sorted];

% Now, activated_neuron_mean_sorted contains the rows of neuron_mean filtered by respClass_all == 1
% and sorted by the time of peak activity.

figure;
% Generate the heatmap
imagesc(ts1, 1, cluster_neurons_neurons_sorted);

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
