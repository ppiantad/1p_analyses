% Define the custom colormap from white to orange
% custom_colormap = [
%     1, 1, 1; % white
%     1, 0.9, 0.8;
%     1, 0.8, 0.6;
%     1, 0.7, 0.4;
%     1, 0.6, 0.2;
%     1, 0.5, 0; % orange
% ];

plot_num = 29; 


% custom_colormap = [
%     1, 1, 1;       % white
%     0.9, 0.95, 0.9;
%     0.8, 0.9, 0.8;
%     0.6, 0.8, 0.6;
%     0.4, 0.7, 0.4;
%     0.2, 0.6, 0.2;
%     0.13, 0.55, 0.13; % forest green
% ];

% custom_colormap = [
%     1, 1, 1;         % white
%     0.9, 0.95, 0.95;
%     0.8, 0.9, 0.9;
%     0.6, 0.85, 0.85;
%     0.4, 0.8, 0.8;
%     0.2, 0.8, 0.8;
%     0.0, 0.8, 0.8;   % robin's egg blue
% ];
% 
custom_colormap = [
    1, 1, 1;         % white
    0.9, 0.9, 0.95;
    0.8, 0.8, 0.9;
    0.6, 0.6, 0.8;
    0.4, 0.4, 0.7;
    0.2, 0.2, 0.6;
    0.0, 0.0, 0.55;   % dark blue
];

% Generate more intermediate colors for a smoother transition
n = 256; % Number of colors
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
ax1 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, 0, zall_array{1, plot_num});

% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar(ax1, 'eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(zall_array{1, plot_num}, 1)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(zall_array{1, plot_num}, 1)]);
xline(0)
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;

% Plot the mean as a thick black line
meanData = mean(zall_array{1, plot_num});
plot(ts1, meanData, 'r', 'LineWidth', 2, 'Color', custom_colormap(end, :));

% Plot the raw data in grey with transparency
for trial = 1:size(zall_array{1, plot_num}, 1)
    plot(ts1, zall_array{1, plot_num}(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
end
ylim([-4 4]);
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
xline(0)
yline(0)
hold off;
