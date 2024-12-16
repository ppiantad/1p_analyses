% use matched dataset, run eventRelatedActivityAndClassification_PTP_v4 on
% something like collectionTime for RM_D1 and Pre_RDT_RM

consum_neurons = neuron_mean_array{1, 1}(respClass_all_array{1, 1}==1, :);

% Sort the rows of activated_neuron_mean based on peak_times.
[peak_values, time_of_peak_activity] = max(consum_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
consum_neurons_sorted = consum_neurons(sort_indices, :);

figure;
% Generate the heatmap
imagesc(ts1, 1, consum_neurons_sorted);

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

late_rm_consum_neurons_on_early_rm = neuron_mean_array{1, 2}(respClass_all_array{1, 1}==1, :);
consum_neurons_sorted_early_rm_sorted = late_rm_consum_neurons_on_early_rm(sort_indices, :);

figure;
% Generate the heatmap
imagesc(ts1, 1, consum_neurons_sorted_early_rm_sorted);

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