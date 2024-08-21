% Pre_RDT_RM & RDT D1 matched
% RDT D1
% collectionTime
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 1
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 2, 'BLOCK', 3

% actually, load overall workspace (e.g. workspace with 5 variables
% identified, then identify the one that you want to see change, (e.g.,
% collectionTime Blocks 2 & 3), then change the variable below to do that
% comparison with the appropriate array #s

% use this to specify which arrays should be compared below. this is
% critical - you must know which arrays correspond to variables to be
% compared, i..e, prechoice B1 vs prechoice B2/3, which may correspond to
% respClass_all_arrays 1 and 4, in the case of arrays of size 6. this code
% will break and could make invalid comparisons if you do not adjust these values appropriately 
if size(respClass_all_array, 2) == 10
    comparison_arrays = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays = [1 2 3; 4 5 6]
end



arrays_to_examine = [1 8];

inhib_or_excite = 1;

event_for_figures = 1; 

%%

% prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite;
prechoice_blocks_2_and_3 = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;

% postchoice_reward_block_1 = respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
postchoice_reward_block_1 = respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite;
postchoice_reward_blocks_2_and_3 = respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;

% collect_block_1 = respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
collect_block_1 = respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite;
collect_blocks_2_and_3 = respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite;

% shk_neurons = respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} == inhib_or_excite;
shk_neurons = 0;

rest_of_neurons = neuron_num - [sum(prechoice_block_1)+sum(postchoice_reward_block_1)+sum(collect_block_1)+sum(shk_neurons)];
figure; pie([sum(prechoice_block_1), sum(postchoice_reward_block_1), sum(collect_block_1), sum(shk_neurons), rest_of_neurons])

block_1_pre_and_post = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite;
sum(block_1_pre_and_post)
block_1_post_and_consumption = respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite;
sum(block_1_post_and_consumption)
block_1_pre_and_consumption = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite;
sum(block_1_pre_and_consumption)


block_2_and_3_pre_and_post = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;
sum(block_2_and_3_pre_and_post)
block_2_and_3_post_and_consumption = respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite;
sum(block_2_and_3_post_and_consumption)
block_2_and_3_pre_and_consumption = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite;
sum(block_2_and_3_pre_and_consumption)


%% 
% Run one cell below, then do all the figure generation in the following
% code. Then run the next cell below, and do the same figure generation

%%
arrays_to_examine = [comparison_arrays(1, 1) comparison_arrays(2, 1)];
conserved = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
conserved_sum(1) = sum(conserved)

lost = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 ~= event_for_figures;
lost_sum(1) = sum(lost)

remapped = prechoice_block_1 ~= event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
remapped_sum(1) = sum(remapped)

vars_to_use = {'prechoice_block_1', 'prechoice_blocks_2_and_3'};

%%
arrays_to_examine = [comparison_arrays(1, 2) comparison_arrays(2, 2)];
conserved = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
conserved_sum(2) = sum(conserved)

lost = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 ~= event_for_figures;
lost_sum(2) = sum(lost)

remapped = postchoice_reward_block_1 ~= event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
remapped_sum(2) = sum(remapped)

vars_to_use = {'postchoice_reward_block_1', 'postchoice_reward_blocks_2_and_3'};

%%
arrays_to_examine = [comparison_arrays(1, 3) comparison_arrays(2, 3)];
conserved = collect_block_1 == event_for_figures & collect_blocks_2_and_3 == event_for_figures;
conserved_sum(3) = sum(conserved)

lost = collect_block_1 == event_for_figures & collect_blocks_2_and_3 ~= event_for_figures;
lost_sum(3) = sum(lost)

remapped = collect_block_1 ~= event_for_figures & collect_blocks_2_and_3 == event_for_figures;
remapped_sum(3) = sum(remapped)

vars_to_use = {'collect_block_1', 'collect_blocks_2_and_3'};

%% DO NOT USE IF USING CODE ABOVE

collect_all_possible = sum([sum(respClass_all_array{1, arrays_to_examine(1)} ==1), sum(respClass_all_array{1, arrays_to_examine(2)} ==1)])
conserved = respClass_all_array{1, arrays_to_examine(1)} == event_for_figures & respClass_all_array{1, arrays_to_examine(2)} == event_for_figures;
conserved_sum = sum(conserved)
collect_original_activated_sum = sum(respClass_all_array{1, arrays_to_examine(1)}==1)
collect_original_inhibited_sum = sum(respClass_all_array{1, arrays_to_examine(1)}==2)
remapped = respClass_all_array{1, arrays_to_examine(1)} ~= event_for_figures & respClass_all_array{1, arrays_to_examine(2)} == event_for_figures;
remapped_sum = sum(remapped)
not_collect_original_sum = sum(respClass_all_array{1, arrays_to_examine(1)} == 3)

lost = respClass_all_array{1, arrays_to_examine(1)} == event_for_figures & respClass_all_array{1, arrays_to_examine(2)} ~= event_for_figures;
lost_sum = sum(lost)
collect_all_possible = sum([sum(respClass_all_array{1, arrays_to_examine(1)} == 1), sum(respClass_all_array{1, arrays_to_examine(2)} ==1)])
collect_conserved_percent = conserved_sum/neuron_num
collect_remapped_percent = remapped_sum/neuron_num
collect_lost_percent = lost_sum/neuron_num
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

%% START HERE IF USING CODE SEQUENCES ABOVE

figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(eval(vars_to_use{1, 1})  == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(1)}(eval(vars_to_use{1, 1})  == event_for_figures, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-8 8]);
ylim([-0.6 0.8])
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');


inner_pie = [conserved_sum/neuron_num,...
            
            lost_sum/neuron_num,...

            remapped_sum/neuron_num,...
           
            (neuron_num- [conserved_sum+remapped_sum+lost_sum]) / neuron_num];


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
plot(ts1, mean(neuron_mean_array{1, 1}(conserved  == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, arrays_to_examine(1)}(remapped == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, arrays_to_examine(2)}(remapped == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, arrays_to_examine(2)}(lost == 1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, arrays_to_examine(1)}(lost == 1, :)))


figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
% xtickformat('%.2f');
ytickformat('%.2f');
hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(conserved  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(conserved  ==1, :)), 'lineProps', {'color', acton(1,:)});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(remapped  ==1, :)), 'lineProps', {'color', acton(50,:)});
h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(2)}(remapped  ==1, :)), 'lineProps', {'color', acton(150,:)});
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
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(conserved  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(conserved  ==1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(conserved  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(2)}(conserved  ==1, :)), 'lineProps', {'color', 'b'});
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
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(1)}(conserved  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(1)}(conserved  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(1)}(conserved  ==1), 2)]);
xline(0)
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(2)}(conserved  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(2)}(conserved  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(2)}(conserved  ==1), 2)]);
xline(0)
hold off;

%%
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(lost  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(lost  ==1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(lost  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(2)}(lost  ==1, :)), 'lineProps', {'color', 'b'});
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
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(1)}(lost  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(1)}(lost  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(1)}(lost  ==1), 2)]);
xline(0)
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(2)}(lost  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(2)}(lost  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(2)}(lost  ==1), 2)]);
xline(0)
hold off;


%%
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(remapped  ==1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(2)}(remapped  ==1, :)), 'lineProps', {'color', 'b'});
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
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(1)}(remapped  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(1)}(remapped  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(1)}(remapped  ==1), 2)]);
xline(0)
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, [], neuron_mean_array{1,arrays_to_examine(2)}(remapped  ==1, :))

% % Apply the custom colormap
% colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(neuron_mean_array{1,arrays_to_examine(2)}(remapped  ==1), 2)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(neuron_mean_array{1,arrays_to_examine(2)}(remapped  ==1), 2)]);
xline(0)
hold off;

%% MULTIPLEX NEURONS
%%
arrays_to_examine = [1 8];
conserved = block_1_pre_and_post  == event_for_figures & block_2_and_3_pre_and_post == event_for_figures;
conserved_sum = sum(conserved)

lost = block_1_pre_and_post == event_for_figures & block_2_and_3_pre_and_post ~= event_for_figures;
lost_sum = sum(lost)

remapped = block_1_pre_and_post ~= event_for_figures & block_2_and_3_pre_and_post == event_for_figures;
remapped_sum = sum(remapped)

vars_to_use = {'block_1_pre_and_post', 'block_2_and_3_pre_and_post'};

%%
arrays_to_examine = [2 9];
conserved = block_1_pre_and_post == event_for_figures & block_2_and_3_pre_and_post == event_for_figures;
conserved_sum = sum(conserved)

lost = block_1_pre_and_post == event_for_figures & block_2_and_3_pre_and_post ~= event_for_figures;
lost_sum = sum(lost)

remapped = block_1_pre_and_post ~= event_for_figures & block_2_and_3_pre_and_post == event_for_figures;
remapped_sum = sum(remapped)

vars_to_use = {'block_1_pre_and_post', 'block_2_and_3_pre_and_post'};

%%
arrays_to_examine = [3 10];
conserved = block_1_pre_and_post == event_for_figures & block_2_and_3_pre_and_post == event_for_figures;
conserved_sum = sum(conserved)

lost = block_1_pre_and_post == event_for_figures & block_2_and_3_pre_and_post ~= event_for_figures;
lost_sum = sum(lost)

remapped = block_1_pre_and_post ~= event_for_figures & block_2_and_3_pre_and_post == event_for_figures;
remapped_sum = sum(remapped)

vars_to_use = {'block_1_pre_and_post', 'block_2_and_3_pre_and_post'};

%%
arrays_to_examine = [1 8];
conserved = block_1_post_and_consumption  == event_for_figures & block_2_and_3_post_and_consumption == event_for_figures;
conserved_sum = sum(conserved)

lost = block_1_post_and_consumption == event_for_figures & block_2_and_3_post_and_consumption ~= event_for_figures;
lost_sum = sum(lost)

remapped = block_1_post_and_consumption ~= event_for_figures & block_2_and_3_post_and_consumption == event_for_figures;
remapped_sum = sum(remapped)

vars_to_use = {'block_1_post_and_consumption', 'block_2_and_3_post_and_consumption'};

%%
arrays_to_examine = [2 9];
conserved = block_1_post_and_consumption == event_for_figures & block_2_and_3_post_and_consumption == event_for_figures;
conserved_sum = sum(conserved)

lost = block_1_post_and_consumption == event_for_figures & block_2_and_3_post_and_consumption ~= event_for_figures;
lost_sum = sum(lost)

remapped = block_1_post_and_consumption ~= event_for_figures & block_2_and_3_post_and_consumption == event_for_figures;
remapped_sum = sum(remapped)

vars_to_use = {'block_1_post_and_consumption', 'block_2_and_3_post_and_consumption'};

%%
arrays_to_examine = [3 10];
conserved = block_1_post_and_consumption == event_for_figures & block_2_and_3_post_and_consumption == event_for_figures;
conserved_sum = sum(conserved)

lost = block_1_post_and_consumption == event_for_figures & block_2_and_3_post_and_consumption ~= event_for_figures;
lost_sum = sum(lost)

remapped = block_1_post_and_consumption ~= event_for_figures & block_2_and_3_post_and_consumption == event_for_figures;
remapped_sum = sum(remapped)

vars_to_use = {'block_1_post_and_consumption', 'block_2_and_3_post_and_consumption'};

%% for plotting changes on a donut (specific) and pie (broad) charts
all_conserved_sum = sum(conserved_sum)
all_lost_sum = sum(lost_sum)
all_remapped_sum = sum(remapped_sum)
remaining_neurons = neuron_num - (all_conserved_sum + all_lost_sum +all_remapped_sum);

figure;
piechart([all_conserved_sum/neuron_num, all_lost_sum/neuron_num, all_remapped_sum/neuron_num, remaining_neurons/neuron_num])

figure;
donutchart([conserved_sum/neuron_num, lost_sum/neuron_num, remapped_sum/neuron_num, remaining_neurons/neuron_num])

all_neurons = [conserved_sum; lost_sum; remapped_sum]
%all_neurons_2 = [conserved_sum; lost_sum; remapped_sum]

%% chi-square test on all_neurons. load in data for 1 set of all_neurons above, then load in second set and create all_neurons_2, then run code below

% Given arrays
% all_neurons = [ ... ]; % Your 3x3 array
% all_neurons_2 = [ ... ]; % Your other 3x3 array

% Initialize a matrix to store p-values
p_values = zeros(size(all_neurons));

% Perform chi-square test for each corresponding cell
for i = 1:size(all_neurons, 1)
    for j = 1:size(all_neurons, 2)
        observed = [all_neurons(i,j), all_neurons_2(i,j)];
        expected = mean(observed);
        chi2_stat = sum((observed - expected).^2 ./ expected);
        p_values(i,j) = 1 - chi2cdf(chi2_stat, 1); % degrees of freedom = 1
    end
end

% Display the p-values
disp('P-values for each corresponding cell:');
disp(p_values);
