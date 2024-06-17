% Define the custom colormap from white to orange
% custom_colormap = [
%     1, 1, 1; % white
%     1, 0.9, 0.8;
%     1, 0.8, 0.6;
%     1, 0.7, 0.4;
%     1, 0.6, 0.2;
%     1, 0.5, 0; % orange
% ];

% green
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
figure;
subplot(2, 1, 1);
hold on;
% 949, 1001, 946
% Plot the heatmap

imagesc(ts1, 0 , zall_array{1, 946});

% Apply the custom colormap
colormap(custom_colormap);
% Restrict the color axis range to [2, 0]
clim([-1 1]);
% Add a colorbar
colorbar;

hold off;

% % Create subplot for the mean and raw data
subplot(2, 1, 2);
hold on;

% Plot the mean as a thick black line
meanData = mean(zall_array{1, 946});
plot(ts1, meanData, 'r', 'LineWidth', 2, 'Color', custom_colormap(end, :));


% Plot the raw data in grey with transparency
for trial = 1:size(zall_array{1, 946}, 1)
    plot(ts1, zall_array{1, 946}(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
    hold on;
end
