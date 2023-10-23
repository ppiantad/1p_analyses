window_sz = (0:.1:20);
ts1 = (-10:.1:10);

% Define the categories
category_1 = 1;
category_2 = 2;
category_3 = 3;

% Define the time window of interest

start_time = 0; % sub-window start time
end_time = 3; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Calculate the mean activity of each neuron within the time window
mean_activity = mean(neuron_mean(:,sub_window_idx), 2);


% Find the indices of the neurons in each category
indices_1 = find(respClass_all == category_1);
indices_2 = find(respClass_all == category_2);
indices_3 = find(respClass_all == category_3);

% Select the three most representative neurons from each category
[~, sorted_indices_1] = sort(mean_activity(indices_1), 'descend');
[~, sorted_indices_2] = sort(mean_activity(indices_2), 'descend');
[~, sorted_indices_3] = sort(mean_activity(indices_3), 'descend');


midway_for_neutral = size(sorted_indices_3, 1)/2;

representative_neurons_1 = indices_1(sorted_indices_1(1:3));
representative_neurons_2 = indices_2(sorted_indices_2(end-3:end));
representative_neurons_3 = indices_3(sorted_indices_3(midway_for_neutral-3:midway_for_neutral));

representative_neurons_all = [neuron_mean(representative_neurons_1,:); neuron_mean(representative_neurons_2,:); neuron_mean(representative_neurons_3,:)];

% Create a figure
figure;

% Plot the activity of each representative neuron as a separate line
hold on;
plot(sgolayfilt((representative_neurons_all(:,:))', 9, 21));
plot(sgolayfilt(neuron_mean(representative_neurons_1,:)', 9, 21), 'LineWidth', 2)';
plot(sgolayfilt(neuron_mean(representative_neurons_2,:)', 9, 21), 'LineWidth', 2)';
plot(sgolayfilt(neuron_mean(representative_neurons_3,:)', 9, 21), 'LineWidth', 2)';

% Add labels and a legend
xlabel('Time (s)');
ylabel('Activity');
legend('Category 1', 'Category 2', 'Category 3');

%%
% Create a figure
figure;

% Define an offset for each category
offset_1 = max(neuron_mean(representative_neurons_1,:),[],2);
offset_2 = max(neuron_mean(representative_neurons_2,:),[],2) + 0.5;
offset_3 = max(neuron_mean(representative_neurons_3,:),[],2) + 1;

offset_all = [offset_1; offset_2; offset_3];

% Plot the activity of each representative neuron as a separate line with offset
hold on;
for i = 1:size(representative_neurons_all, 1)
    plot(ts1, representative_neurons_all(i, :) - offset_all(i)+ sum(offset_all(1:i)), 'LineWidth', 2);
end



%%
figure;
hold on;
for i = 1:length(representative_neurons_1)
    plot(ts1, neuron_mean(representative_neurons_1(i),:), 'LineWidth', 2);
end
hold off;
figure;
hold on;
for i = 1:length(representative_neurons_2)
    plot(ts1, neuron_mean(representative_neurons_2(i),:), 'LineWidth', 2);
end
figure;
hold on;
for i = 1:length(representative_neurons_3)
    plot(ts1, neuron_mean(representative_neurons_3(i),:), 'LineWidth', 2);
end

% Add labels and a legend
xlabel('Time (s)');
ylabel('Activity');
legend('Category 1', 'Category 2', 'Category 3');
