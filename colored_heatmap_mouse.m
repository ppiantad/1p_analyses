
plot_num = 46; 

array_to_plot = 1; % depends on the structure of zall

select_mouse = 'BLA_Insc_25';

% for RDT D1 BLA_Insc_25:
%prechoice neuron num 46
%postchoice rew num 38
%consumption num 39
%shock num 11

select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';

second_session = 'RDT_D1';


%%
% in order to trim off excess (because calcium recording starts before
% behavior), you'll need the start time. unfortunately I don't save that
% variable in the "final" struct, but you can get it from the adjustment to
% some of the columns of BehavData. e.g., see below - but make sure to
% update the session etc as necessary! 

% BehavData = final.(select_mouse).(first_session).choiceTime.uv.BehavData;
BehavData = final.(select_mouse).(first_session).uv.BehavData;
% because the first trial possible is ALWAYS 60 seconds after ABET is
% issued, you can determine what adjustment has been made to this column
% (adding time to account for calcium recording starting first) by
% subtracting off 60 from the first element
stTime = BehavData.TrialPossible(1)-60; 



time2Collect = BehavData.collectionTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
trialStartTime = BehavData.stTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
median_trialStartTime = median(trialStartTime)
median_time2Collect = median(time2Collect)
xline(median_trialStartTime)
xline(median_time2Collect)
[numTrials, ~] = size(time2Collect);
Tris = [1:numTrials]';




%%





% Define the custom colormap from white to orange
% custom_colormap = [
%     1, 1, 1; % white
%     1, 0.9, 0.8;
%     1, 0.8, 0.6;
%     1, 0.7, 0.4;
%     1, 0.6, 0.2;
%     1, 0.5, 0; % orange
% ];




custom_colormap = [
    1, 1, 1;       % white
    0.9, 0.95, 0.9;
    0.8, 0.9, 0.8;
    0.6, 0.8, 0.6;
    0.4, 0.7, 0.4;
    0.2, 0.6, 0.2;
    0.13, 0.55, 0.13; % forest green
];

% custom_colormap = [
%     1, 1, 1;         % white
%     0.9, 0.95, 0.95;
%     0.8, 0.9, 0.9;
%     0.6, 0.85, 0.85;
%     0.4, 0.8, 0.8;
%     0.2, 0.8, 0.8;
%     0.0, 0.8, 0.8;   % robin's egg blue
% ];

% custom_colormap = [
%     1, 1, 1;         % white
%     0.9, 0.9, 0.95;
%     0.8, 0.8, 0.9;
%     0.6, 0.6, 0.8;
%     0.4, 0.4, 0.7;
%     0.2, 0.2, 0.6;
%     0.0, 0.0, 0.55;   % dark blue
% ];

% custom_colormap = [
%     1, 1, 1; % white
%     1, 0.9, 0.9;
%     1, 0.8, 0.8;
%     1, 0.7, 0.7;
%     1, 0.6, 0.6;
%     1, 0.5, 0.5;
%     1, 0.4, 0.4;
%     1, 0.3, 0.3;
%     1, 0.2, 0.2;
%     1, 0.1, 0.1;
%     1, 0, 0;   % red
% ];


% Generate more intermediate colors for a smoother transition
n = 256; % Number of colors
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 250, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
ax1 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, 1:size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1), zall_mouse{select_mouse_index, array_to_plot}{1, plot_num});


% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar(ax1, 'eastoutside');
set(c, 'YTick', clim); % 

ylim([0.5, size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1) + 0.5]);

xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', []);
set(gca, 'YTick', [1, size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1)]);
xline(0)
scatter(time2Collect, Tris               , 'Marker', 'p')
scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;


% Plot the raw data in grey with transparency
for trial = 1:size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1)
    plot(ts1, zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
end

% Plot the mean as a thick black line
meanData = mean(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num});
plot(ts1, meanData, 'r', 'LineWidth', 2, 'Color', 'k');

ylim([-4 4]);
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
xline(0)
yline(0)
fontsize(18, 'points')
hold off;
