% Create a combined array for LargeRew and SmallRew
% Add a flag: 1 for LargeRew, 2 for SmallRew
combined_data = [LargeRew', ones(size(LargeRew')); SmallRew', ones(size(SmallRew')) * 2];

% Sort by time (first column)
combined_data = sortrows(combined_data, 1);

% Extract sorted times and labels
sorted_times = combined_data(:, 1);
sorted_labels = combined_data(:, 2);

% Plot the data
figure;
hold on;

% Define colors
blue = [0 0 0.753]; % LargeRew
green = [0 0.353 0]; % SmallRew
red = [0.8 0.2 0.2]; % AA_large

% Loop through sorted times to plot reward markers and AA_large lines
for i = 1:length(sorted_times)
    % Plot LargeRew or SmallRew marker
    if sorted_labels(i) == 1
        % LargeRew marker
        plot([0, 0], [i - 0.4, i + 0.4], 'Color', blue, 'LineWidth', 2);
        
        % Identify AA_large events after this LargeRew and before the next event
        current_time = sorted_times(i);
        if i < length(sorted_times)
            next_time = sorted_times(i + 1);
        else
            next_time = Inf; % No upper limit for the last reward
        end
        
        % Find relevant AA_large events
        relevant_AA_large = AA_large(AA_large >= current_time & AA_large < next_time);
        
        % Plot each AA_large event as a separate line
        for j = 1:length(relevant_AA_large)
            % Adjust X-coordinates for the line
            x_start = 0.2 + 0.6 * (j - 1); % Slightly offset lines horizontally
            x_end = x_start + 0.5;         % Line width
            plot([x_start, x_end], [i, i], 'Color', red, 'LineWidth', 2);
        end
    elseif sorted_labels(i) == 2
        % SmallRew marker
        plot([0, 0], [i - 0.4, i + 0.4], 'Color', green, 'LineWidth', 2);
    end
end

% Customize axis
xlabel('Event Type');
ylabel('Trial Order');
set(gca, 'YDir', 'reverse'); % Flip Y-axis to have the first trial at the top
yticks(1:length(sorted_times));
yticklabels(1:length(sorted_times));
xlim([-1, max(length(AA_large) * 0.6, 1)]); % Adjust X-axis for markers and lines
ylim([0.5, length(sorted_times) + 0.5]);

% Add legend
legend({'Large Reward', 'Small Reward', 'AA Large'}, 'Location', 'best');

hold off;
