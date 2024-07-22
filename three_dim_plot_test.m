% Create the 3D plot


% Find the indices in ts1 closest to median_start_time
[~, idx_start1] = min(abs(ts1 - median_start_time(1)));
[~, idx_start2] = min(abs(ts1 - median_start_time(2)));
[~, idx_collect1] = min(abs(ts1 - median_collect_time(1)));
[~, idx_collect2] = min(abs(ts1 - median_collect_time(2)));
[~, idx_zero1] = min(abs(ts1 - 0)); % Index for ts1 = 0 for the first line
[~, idx_zero2] = min(abs(ts1 - 0)); % Index for ts1 = 0 for the second line


figure;
plot3(PCScore{1, 1}(1, :), PCScore{1, 1}(2, :), PCScore{1, 1}(3, :)); % '-o' adds circle markers to the line
title(['trajectories for ', num2str(uv.evtWin), ' aligned to ', epoc_to_align]);
xlabel(['PC1 (Variance explained: ' num2str(explained(1), '%.2f') '%)']);
ylabel(['PC2 (Variance explained: ' num2str(explained(2), '%.2f') '%)']);
zlabel(['PC3 (Variance explained: ' num2str(explained(3), '%.2f') '%)']);
% grid on; % Optional: adds a grid to the plot

hold on; plot3(PCScore{1, 2}(1, :), PCScore{1, 2}(2, :), PCScore{1, 2}(3, :)); % '-o' adds circle markers to the line
% hold on; plot3(PCScore{1, 3}(1, :), PCScore{1, 3}(2, :), PCScore{1, 3}(3, :), '-o'); % '-o' adds circle markers to the line


% Plot circles at the median start times
plot3(PCScore{1, 1}(1, idx_start1), PCScore{1, 1}(2, idx_start1), PCScore{1, 1}(3, idx_start1), 'r^', 'MarkerSize', 10, 'LineWidth', 2); % Red circle for the first line
plot3(PCScore{1, 2}(1, idx_start2), PCScore{1, 2}(2, idx_start2), PCScore{1, 2}(3, idx_start2), 'g^', 'MarkerSize', 10, 'LineWidth', 2); % Green circle for the second line

% Plot squares at the median collect times
plot3(PCScore{1, 1}(1, idx_collect1), PCScore{1, 1}(2, idx_collect1), PCScore{1, 1}(3, idx_collect1), 'rs', 'MarkerSize', 10, 'LineWidth', 2); % Red square for the first line
plot3(PCScore{1, 2}(1, idx_collect2), PCScore{1, 2}(2, idx_collect2), PCScore{1, 2}(3, idx_collect2), 'gs', 'MarkerSize', 10, 'LineWidth', 2); % Green square for the second line

% Plot circles at ts1 = 0
plot3(PCScore{1, 1}(1, idx_zero1), PCScore{1, 1}(2, idx_zero1), PCScore{1, 1}(3, idx_zero1), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Blue circle for ts1 = 0 on the first line
plot3(PCScore{1, 2}(1, idx_zero2), PCScore{1, 2}(2, idx_zero2), PCScore{1, 2}(3, idx_zero2), 'go', 'MarkerSize', 10, 'LineWidth', 2); % Blue circle for ts1 = 0 on the second line
% Optionally, add a legend
legend('Block 1 - no shock', 'Block 2 and 3 - no shock', 'Median Start Time 1', 'Median Start Time 2', 'Median Collect Time 1', 'Median Collect Time 2', 'ts1 = 0');
grid("off")
%%
% Create the 3D plot
figure;
plot(PCScore{1, 1}(1, :), PCScore{1, 1}(2, :), '-o'); % '-o' adds circle markers to the line
title('3D Plot of PCScore');
xlabel(['PC1 (Variance explained: ' num2str(explained(1), '%.2f') '%)']);
ylabel(['PC2 (Variance explained: ' num2str(explained(2), '%.2f') '%)']);
grid on; % Optional: adds a grid to the plot

hold on; plot(PCScore{1, 2}(1, :), PCScore{1, 2}(2, :), '-o'); % '-o' adds circle markers to the line
% hold on; plot3(PCScore{1, 3}(1, :), PCScore{1, 3}(2, :), PCScore{1, 3}(3, :), '-o'); % '-o' adds circle markers to the line

%%

% Create the 3D plot
figure;
plot(ts1, PCScore{1, 1}(1, :), '-o'); % '-o' adds circle markers to the line
title(['trajectories for ', num2str(uv.evtWin), ' aligned to ', epoc_to_align]);
xlabel('X values (Row 1)');
ylabel('Y values (Row 2)');
zlabel('Z values (Row 3)');
grid on; % Optional: adds a grid to the plot

hold on; plot(ts1, PCScore{1, 2}(1, :), '-o'); % '-o' adds circle markers to the line
% hold on; plot3(PCScore{1, 3}(1, :), PCScore{1, 3}(2, :), PCScore{1, 3}(3, :), '-o'); % '-o' adds circle markers to the line
%%

% Create the 3D plot
figure;
plot(ts1, PCScore{1, 1}(2, :), '-o'); % '-o' adds circle markers to the line
title(['trajectories for ', num2str(uv.evtWin), ' aligned to ', epoc_to_align]);
xlabel('X values (Row 1)');
ylabel('Y values (Row 2)');
zlabel('Z values (Row 3)');
grid on; % Optional: adds a grid to the plot

hold on; plot(ts1, PCScore{1, 2}(2, :), '-o'); % '-o' adds circle markers to the line
% hold on; plot3(PCScore{1, 3}(1, :), PCScore{1, 3}(2, :), PCScore{1, 3}(3, :), '-o'); % '-o' adds circle markers to the line

%%

% Create the 3D plot
figure;
plot(ts1, PCScore{1, 1}(3, :), '-o'); % '-o' adds circle markers to the line
title(['trajectories for ', num2str(uv.evtWin), ' aligned to ', epoc_to_align]);
xlabel('X values (Row 1)');
ylabel('Y values (Row 2)');
zlabel('Z values (Row 3)');
grid on; % Optional: adds a grid to the plot

hold on; plot(ts1, PCScore{1, 2}(3, :), '-o'); % '-o' adds circle markers to the line
% hold on; plot3(PCScore{1, 3}(1, :), PCScore{1, 3}(2, :), PCScore{1, 3}(3, :), '-o'); % '-o' adds circle markers to the line
