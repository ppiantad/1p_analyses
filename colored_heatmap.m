% Define the custom colormap from white to orange
custom_colormap = [
    1, 1, 1; % white
    1, 0.9, 0.8;
    1, 0.8, 0.6;
    1, 0.7, 0.4;
    1, 0.6, 0.2;
    1, 0.5, 0; % orange
];

% Generate more intermediate colors for a smoother transition
n = 256; % Number of colors
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));

% Plot the heatmap
figure;
imagesc(ts1, 0 , zall_mouse{1, 4}{1, 8});

% Apply the custom colormap
colormap(custom_colormap);
% Restrict the color axis range to [2, 0]
clim([-1 2]);
% Add a colorbar
colorbar;
