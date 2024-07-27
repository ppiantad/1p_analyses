% Pre_RDT_RM & RDT D1 matched
% RDT D1
% collectionTime
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 1
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 2, 'BLOCK', 3

% actually, load overall workspace (e.g. workspace with 5 variables
% identified, then identify the one that you want to see change, (e.g.,
% collectionTime Blocks 2 & 3), then change the variable below to do that
% comparison with the appropriate array #s

arrays_to_examine = [3 10];

event_for_figures = 1; 

%%

prechoice_block_1 = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} ~= 1 & respClass_all_array{1, 3} ~= 1;
prechoice_blocks_2_and_3 = respClass_all_array{1, 8} == 1 & respClass_all_array{1, 9} ~= 1 & respClass_all_array{1, 10} ~= 1;

postchoice_reward_block_1 = respClass_all_array{1, 2} == 1 & respClass_all_array{1, 1} ~= 1 & respClass_all_array{1, 3} ~= 1;
postchoice_reward_blocks_2_and_3 = respClass_all_array{1, 9} == 1 & respClass_all_array{1, 8} ~= 1 & respClass_all_array{1, 10} ~= 1;

collect_block_1 = respClass_all_array{1, 3} == 1 & respClass_all_array{1, 1} ~= 1 & respClass_all_array{1, 2} ~= 1;
collect_blocks_2_and_3 = respClass_all_array{1, 10} == 1 & respClass_all_array{1, 8} ~= 1 & respClass_all_array{1, 9} ~= 1;

collect_conserved = collect_block_1 == event_for_figures & collect_blocks_2_and_3 == event_for_figures;
collect_conserved_sum = sum(collect_conserved)

collect_lost = collect_block_1 == event_for_figures & collect_blocks_2_and_3 ~= event_for_figures;
collect_lost_sum = sum(collect_lost)

collect_remapped = collect_block_1 ~= event_for_figures & collect_blocks_2_and_3 == event_for_figures;
collect_remapped_sum = sum(collect_remapped)

%%



collect_all_possible = sum([sum(respClass_all_array{1, arrays_to_examine(1)} ==1), sum(respClass_all_array{1, arrays_to_examine(2)} ==1)])
collect_conserved = respClass_all_array{1, arrays_to_examine(1)} == event_for_figures & respClass_all_array{1, arrays_to_examine(2)} == event_for_figures;
collect_conserved_sum = sum(collect_conserved)
collect_original_activated_sum = sum(respClass_all_array{1, arrays_to_examine(1)}==1)
collect_original_inhibited_sum = sum(respClass_all_array{1, arrays_to_examine(1)}==2)
collect_remapped = respClass_all_array{1, arrays_to_examine(1)} ~= event_for_figures & respClass_all_array{1, arrays_to_examine(2)} == event_for_figures;
collect_remapped_sum = sum(collect_remapped)
not_collect_original_sum = sum(respClass_all_array{1, arrays_to_examine(1)} == 3)

collect_lost = respClass_all_array{1, arrays_to_examine(1)} == event_for_figures & respClass_all_array{1, arrays_to_examine(2)} ~= event_for_figures;
collect_lost_sum = sum(collect_lost)
collect_all_possible = sum([sum(respClass_all_array{1, arrays_to_examine(1)} == 1), sum(respClass_all_array{1, arrays_to_examine(2)} ==1)])
collect_conserved_percent = collect_conserved_sum/neuron_num
collect_remapped_percent = collect_remapped_sum/neuron_num
collect_lost_percent = collect_lost_sum/neuron_num
collect_new_activated_sum = sum(respClass_all_array{1, arrays_to_examine(2)}==1)
collect_new_inhibited_sum = sum(respClass_all_array{1, arrays_to_examine(2)}==2)
not_collect_new_sum = sum(respClass_all_array{1, arrays_to_examine(2)} ==3)

figure; pie([collect_original_activated_sum collect_original_inhibited_sum not_collect_original_sum])
title({"Event: " + strrep(full_filter_string{arrays_to_examine(1)}, '_', '-'), "Num of events (across mice): " + num_trials_per_event(arrays_to_examine(1)), "Num of neurons (across mice): " + neuron_num}, 'FontSize', 9)
labels = {'activated: ' + string(collect_original_activated_sum), 'inhibited: ' + string(collect_original_inhibited_sum), 'neutral: ' + string(not_collect_original_sum)};
legend(labels)

figure; pie([collect_new_activated_sum collect_new_inhibited_sum not_collect_new_sum])
title({"Event: " + strrep(full_filter_string{arrays_to_examine(2)}, '_', '-'), "Num of events (across mice): " + num_trials_per_event(arrays_to_examine(2)), "Num of neurons (across mice): " + neuron_num}, 'FontSize', 9)
labels = {'activated: ' + string(collect_new_activated_sum), 'inhibited: ' + string(collect_new_inhibited_sum), 'neutral: ' + string(not_collect_new_sum)};
legend(labels)


figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(respClass_all_array{1, arrays_to_examine(1)} == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(1)}(respClass_all_array{1, arrays_to_examine(1)} == event_for_figures, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(respClass_all_array{1, arrays_to_examine(2)} == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(respClass_all_array{1, arrays_to_examine(2)} == event_for_figures, :)), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-8 8]);
ylim([-0.6 0.8])
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

inner_pie = [collect_conserved_sum/neuron_num,...
            
            collect_lost_sum/neuron_num,...

            collect_remapped_sum/neuron_num,...
           
            (neuron_num- [collect_conserved_sum+collect_remapped_sum+collect_lost_sum]) / neuron_num];


labels = ["conserved", "lost", "new (remapped)",  "other"];

figure; donutchart(inner_pie, labels, 'InnerRadius', 0.5)

figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
% xtickformat('%.2f');
ytickformat('%.2f');
N = 5;
axes('ColorOrder',brewermap(N,'Pastel2'),'NextPlot','replacechildren')
plot(ts1, mean(neuron_mean_array{1, 1}(collect_conserved  == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, arrays_to_examine(1)}(collect_remapped == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, arrays_to_examine(2)}(collect_remapped == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, arrays_to_examine(2)}(collect_lost == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, arrays_to_examine(1)}(collect_lost == 1, :)))


figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
% xtickformat('%.2f');
ytickformat('%.2f');
hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(collect_conserved  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(collect_conserved  ==1, :)), 'lineProps', {'color', acton(1,:)});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(collect_remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(collect_remapped  ==1, :)), 'lineProps', {'color', acton(50,:)});
h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(collect_remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(2)}(collect_remapped  ==1, :)), 'lineProps', {'color', acton(150,:)});
% h(4) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,2}(collect_lost  ==1, :)), nanmean(neuron_sem_array{1, 2}(collect_lost  ==1, :)), 'lineProps', {'color', acton(200,:)});
% h(5) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,1}(collect_lost  ==1, :)), nanmean(neuron_sem_array{1, 1}(collect_lost  ==1, :)), 'lineProps', {'color', acton(250,:)});
% legend([h(1).mainLine h(2).mainLine h(3).mainLine h(4).mainLine h(5).mainLine], 'conserved', 'pre-remapped (safe blocks)', 'post-remapped (risky blocks)')
legend([h(1).mainLine h(2).mainLine h(3).mainLine], 'conserved', 'pre-remapped (safe blocks)', 'post-remapped (risky blocks)')
xlim([-8 8]);
ylim([-0.6 0.8])

% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');

%%
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(collect_conserved  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(collect_conserved  ==1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(collect_conserved  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(2)}(collect_conserved  ==1, :)), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], 'conserved (safe block)', 'conserved (risky blocks)')
xlim([-8 8]);
ylim([-0.6 0.7])

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
ax1 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(1)}(collect_conserved  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(1)}(collect_conserved  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(1)}(collect_conserved  ==1), 2)]);
xline(0)
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(2)}(collect_conserved  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(2)}(collect_conserved  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(2)}(collect_conserved  ==1), 2)]);
xline(0)
hold off;

%%
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(collect_lost  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(collect_lost  ==1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(collect_lost  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(2)}(collect_lost  ==1, :)), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], 'lost (safe block)', 'lost (risky blocks)')
xlim([-8 8]);
ylim([-0.6 0.7])

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
ax1 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(1)}(collect_lost  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(1)}(collect_lost  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(1)}(collect_lost  ==1), 2)]);
xline(0)
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(2)}(collect_lost  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(2)}(collect_lost  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(2)}(collect_lost  ==1), 2)]);
xline(0)
hold off;


%%
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(collect_remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(collect_remapped  ==1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(collect_remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(2)}(collect_remapped  ==1, :)), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], 'new (safe block)', 'new (risky blocks)')
xlim([-8 8]);
ylim([-0.6 0.7])

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
ax1 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(1)}(collect_remapped  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(1)}(collect_remapped  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(1)}(collect_remapped  ==1), 2)]);
xline(0)
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(2)}(collect_remapped  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(2)}(collect_remapped  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(2)}(collect_remapped  ==1), 2)]);
xline(0)
hold off;