% https://www.pnas.org/doi/10.1073/pnas.1901712116
% Goal is to reproduce something like Fig. 3 from this paper

% Assuming neuron_mean is a 600x200 matrix, and respClass_all is a 1x600 vector.
% Transpose respClass_all to match the shape of neuron_mean.
respClass_all = respClass_all';

respClass_all = respClass_all_array{1, 4}';


% respClass_all_2 = respClass.Pre_RDT_RM.Outcome_0to2.REW_Large.activated';


% Find the rows in neuron_mean where respClass_all is equal to 1.
activated_rows = neuron_mean(respClass_all_2 == 1, :);

% Calculate the time of peak activity within each row.
% [peak_values, time_of_peak_activity] = max(activated_rows, [], 2);

respClass_all = respClass_all_array{1, 8}';
[peak_values, time_of_peak_activity] = max(neuron_mean, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
neuron_mean_sorted = neuron_mean(sort_indices, :);


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