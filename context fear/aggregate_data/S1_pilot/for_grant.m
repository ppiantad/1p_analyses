
load('S1_for_plotting.mat')
load('B6_for_plottingmat')

num_bins = 24;
interval = 2;
total_time = 12;
bins_per_interval = 4; % Number of bins per interval (4 columns per 2 minutes)

figure('Position', [100, 100, 700, 450]); % [left, bottom, width, height]
hold on;

h(1) = shadedErrorBar(1:num_bins, mean(experimental_data), experimental_sem, 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(1:num_bins, mean(experimental_data_B6), experimental_sem_B6, 'lineProps', {'color', 'k'});
% h(2) = shadedErrorBar(1:num_bins, mean(no_shock_data), no_shock_sem, 'lineProps', {'color', 'b'});

% Calculate x-tick positions and labels
x_tick_positions = bins_per_interval:bins_per_interval:num_bins; % [4, 8, 12, ...]
% time_points = 0:interval:total_time - interval; % Time intervals [0, 2, 4, ...]
time_points = interval:interval:total_time; % Time intervals [0, 2, 4, ...]
% Set the x-ticks and labels
xticks(x_tick_positions); % Exact positions on the axis
xticklabels(arrayfun(@num2str, time_points, 'UniformOutput', false)); % Labels for time
xlabel('Time (min)');


% Set axis limits
xlim([1 num_bins]);
ylim([0 0.9]); % Set y-axis limits


%%

load('S1_for_bar_plot.mat')
load('B6_for_bar_plot.mat')

experimental_data_means = experimental_data_means;
EXPERIMENTAL_DATA_S1_data_means = experimental_data_means_S1;
% no_shock_data_means = [mean(no_shock_data_aversive, 2) mean(no_shock_data_safe, 2)];



% Combine the datasets for easier handling
all_data = {experimental_data_means, EXPERIMENTAL_DATA_S1_data_means};

% Calculate means for bar heights
means = [mean(experimental_data_means); 
         mean(EXPERIMENTAL_DATA_S1_data_means);
        ];

% Grouped positions for the bars
x = [1, 2; 3.5, 4.5; 6, 7]; % Adjust spacing as needed

% Bar plot
figure('Position', [100, 100, 300, 450]); % [left, bottom, width, height]
hold on;

% Loop through each group to plot bars, scatter points, and lines
for i = 1:size(all_data, 2)
    % Bar plot for each group
    for col = 1:2
        bar_x = x(i, col); % Position for the current bar
        bar(bar_x, means(i, col), 0.4, 'FaceAlpha', 0.7); % Plot each bar
    end

    % Overlay scatter points and connect with lines for the current variable
    data = all_data{i}; % Current variable's data
    jittered_x = zeros(size(data)); % To store jittered x-coordinates
    for j = 1:size(data, 1)
        % Scatter points for the current row
        scatter_x = x(i, :) + (rand(1, 2) - 0.5) * 0.2; % Add jitter
        jittered_x(j, :) = scatter_x; % Store jittered x-coordinates
        scatter(scatter_x, data(j, :), 40, 'k', 'filled');
    end

    % Connect scatter points with a line using jittered x-coordinates
    for j = 1:size(data, 1)
        plot(jittered_x(j, :), data(j, :), 'k-', 'LineWidth', 0.5);
    end
end

% Adjustments for aesthetics
set(gca, 'XTick', mean(x, 2), 'XTickLabel', {'Experimental', 'One Context', 'No Shock'});
ylim([0 0.8])
hold off;