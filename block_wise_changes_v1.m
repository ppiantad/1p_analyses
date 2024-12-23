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
if size(respClass_all_array, 2) == 10 | size(respClass_all_array, 2) == 11 | size(respClass_all_array, 2) == 12
    comparison_arrays = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays = [1 2 3; 4 5 6]
elseif size(respClass_all_array, 2) == 7
    comparison_arrays = [1 2 3; 5 6 7]
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
block_1_pre_and_post_and_consumption = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite;
sum(block_1_pre_and_post_and_consumption)



block_2_and_3_pre_and_post = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;
sum(block_2_and_3_pre_and_post)
block_2_and_3_post_and_consumption = respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite;
sum(block_2_and_3_post_and_consumption)
block_2_and_3_pre_and_consumption = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite;
sum(block_2_and_3_pre_and_consumption)


lost_prechoice = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 ~= event_for_figures;
lost_postchoice = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 ~= event_for_figures;
lost_consumption = collect_block_1 == event_for_figures & collect_blocks_2_and_3 ~= event_for_figures;
lost_all = lost_prechoice == 1 | lost_postchoice == 1 | lost_consumption == 1; 


conserved_prechoice = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
conserved_postchoice = postchoice_reward_block_1 == event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
conserved_consumption = collect_block_1 == event_for_figures & collect_blocks_2_and_3 == event_for_figures;
conserved_all = conserved_prechoice == 1 | conserved_postchoice == 1 | conserved_consumption == 1;


remapped_prechoice = prechoice_block_1 ~= event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
remapped_postchoice = postchoice_reward_block_1 ~= event_for_figures & postchoice_reward_blocks_2_and_3 == event_for_figures;
remapped_consumption = collect_block_1 ~= event_for_figures & collect_blocks_2_and_3 == event_for_figures;
remapped_all = remapped_prechoice == 1 | remapped_postchoice == 1 | remapped_consumption == 1;


sum(lost_all)
intersect(lost_prechoice, lost_postchoice)
true_neutral = respClass_all_array{1, comparison_arrays(1, 1)} == 3 & respClass_all_array{1, comparison_arrays(1, 2)} ~= 3 & respClass_all_array{1, comparison_arrays(1, 3)} ~= 3;

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

%%
for kk = 1:size(animalIDs, 1)
    % prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
    prechoice_block_1_mouse{kk, :} = respClass_all_array_mouse{kk, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 3)} ~= inhib_or_excite;
    prechoice_blocks_2_and_3_mouse{kk, :} = respClass_all_array_mouse{kk, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 2)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 3)} ~= inhib_or_excite;

    % postchoice_reward_block_1 = respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
    postchoice_reward_block_1_mouse{kk, :}  = respClass_all_array_mouse{kk, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 3)} ~= inhib_or_excite;
    postchoice_reward_blocks_2_and_3_mouse{kk, :}  = respClass_all_array_mouse{kk, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 3)} ~= inhib_or_excite;

    % collect_block_1 = respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
    collect_block_1_mouse{kk, :}  = respClass_all_array_mouse{kk, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 2)} ~= inhib_or_excite;
    collect_blocks_2_and_3_mouse{kk, :}  = respClass_all_array_mouse{kk, comparison_arrays(2, 3)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 2)} ~= inhib_or_excite;

    num_cells_mouse(kk) = size(respClass_all_array_mouse{kk, 1}, 2);
    true_neutral_block_1_mouse{kk, :} = respClass_all_array_mouse{kk, comparison_arrays(1, 1)} == 3 & respClass_all_array_mouse{kk, comparison_arrays(1, 2)} == 3 & respClass_all_array_mouse{kk, comparison_arrays(1, 3)} == 3;
end




%%
for kk = 1:size(animalIDs, 1)
    % prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
    prechoice_conserved_mouse{kk, :} = prechoice_block_1_mouse{kk, :} == event_for_figures & prechoice_blocks_2_and_3_mouse{kk, :} == event_for_figures;
    prechoice_lost_mouse{kk, :} = prechoice_block_1_mouse{kk, :} == event_for_figures & prechoice_blocks_2_and_3_mouse{kk, :} ~= event_for_figures;
    prechoice_remapped_mouse{kk, :} = prechoice_block_1_mouse{kk, :} ~= event_for_figures & prechoice_blocks_2_and_3_mouse{kk, :} == event_for_figures;

    postchoice_conserved_mouse{kk, :} = postchoice_reward_block_1_mouse{kk, :} == event_for_figures & postchoice_reward_blocks_2_and_3_mouse{kk, :} == event_for_figures;
    postchoice_lost_mouse{kk, :} = postchoice_reward_block_1_mouse{kk, :} == event_for_figures & postchoice_reward_blocks_2_and_3_mouse{kk, :} ~= event_for_figures;
    postchoice_remapped_mouse{kk, :} = postchoice_reward_block_1_mouse{kk, :} ~= event_for_figures & postchoice_reward_blocks_2_and_3_mouse{kk, :} == event_for_figures;

    collection_conserved_mouse{kk, :} = collect_block_1_mouse{kk, :} == event_for_figures & collect_blocks_2_and_3_mouse{kk, :} == event_for_figures;
    collection_lost_mouse{kk, :} = collect_block_1_mouse{kk, :} == event_for_figures & collect_blocks_2_and_3_mouse{kk, :} ~= event_for_figures;
    collection_remapped_mouse{kk, :} = collect_block_1_mouse{kk, :} ~= event_for_figures & collect_blocks_2_and_3_mouse{kk, :} == event_for_figures;

end

%%
for kk = 1:size(animalIDs, 1)
    
    
    conserved_sum(kk) = sum([sum(prechoice_conserved_mouse{kk, :}), sum(postchoice_conserved_mouse{kk, :}), sum(collection_conserved_mouse{kk, :})]);
    lost_sum(kk) = sum([sum(prechoice_lost_mouse{kk, :}), sum(postchoice_lost_mouse{kk, :}), sum(collection_lost_mouse{kk, :})]);
    remapped_sum(kk) = sum([sum(prechoice_remapped_mouse{kk, :}), sum(postchoice_remapped_mouse{kk, :}), sum(collection_remapped_mouse{kk, :})]);
    
    conserved_ratio(kk) = conserved_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);
    lost_ratio(kk) = lost_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);
    remapped_ratio(kk) = remapped_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);

    conserved_prechoice_sum(kk) = sum(prechoice_conserved_mouse{kk, :});
    remapped_prechoice_sum(kk) = sum(prechoice_remapped_mouse{kk, :});
    lost_prechoice_sum(kk) = sum(prechoice_lost_mouse{kk, :});

    conserved_prechoice_ratio(kk) = conserved_prechoice_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);
    remapped_prechoice_ratio(kk) = remapped_prechoice_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);
    lost_prechoice_ratio(kk) = lost_prechoice_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);


    conserved_postchoice_sum(kk) = sum(postchoice_conserved_mouse{kk, :});
    remapped_postchoice_sum(kk) = sum(postchoice_remapped_mouse{kk, :});
    lost_postchoice_sum(kk) = sum(postchoice_lost_mouse{kk, :});

    conserved_postchoice_ratio(kk) = conserved_postchoice_sum(kk)/size(postchoice_conserved_mouse{kk, :}, 2);
    remapped_postchoice_ratio(kk) = remapped_postchoice_sum(kk)/size(postchoice_conserved_mouse{kk, :}, 2);
    lost_postchoice_ratio(kk) = lost_postchoice_sum(kk)/size(postchoice_conserved_mouse{kk, :}, 2);

    conserved_collection_sum(kk) = sum(collection_conserved_mouse{kk, :});
    remapped_collection_sum(kk) = sum(collection_remapped_mouse{kk, :});
    lost_collection_sum(kk) = sum(collection_lost_mouse{kk, :});

    conserved_collection_ratio(kk) = conserved_collection_sum(kk)/size(collection_conserved_mouse{kk, :}, 2);
    remapped_collection_ratio(kk) = remapped_collection_sum(kk)/size(collection_conserved_mouse{kk, :}, 2);
    lost_collection_ratio(kk) = lost_collection_sum(kk)/size(collection_conserved_mouse{kk, :}, 2);
end



%% change the "y" variable below to broadly check correlations w/ variables defined in the "variables" section
% Define a cell array of variable names that end with "ratio"
variables = {'conserved_ratio', 'lost_ratio', 'remapped_ratio', ...
             'conserved_prechoice_ratio', 'remapped_prechoice_ratio', 'lost_prechoice_ratio', ...
             'conserved_postchoice_ratio', 'remapped_postchoice_ratio', 'lost_postchoice_ratio', ...
             'conserved_collection_ratio', 'remapped_collection_ratio', 'lost_collection_ratio'};

% Initialize an empty table to store the results
correlation_results = table();

% Iterate over each variable
for i = 1:length(variables)
    % Extract the variable by its name
    % x = eval(variables{i})';  % Get the corresponding data for the current variable
    % y = risk_table.Mean_1_to_3;  % This is the column you want to correlate with

    % for BLA-NAcSh data since there are so few cells uncomment below
    x = eval(variables{i});
    x = x(num_cells_mouse > 30)';
    y = risk_table.Var11(num_cells_mouse > 30);
    % Compute the correlation coefficient
    [r, pval] = corrcoef(x, y);
    
    % Store the results in the table
    correlation_results = [correlation_results; table({variables{i}}, r(2), pval(2), 'VariableNames', {'Variable', 'Correlation', 'PValue'})];
end

%%
% use this to get mean activity of prechoice ensemble across mice for
% Block1 and Blocks2/3

for gg = 1:size(animalIDs, 1) %size(animalIDs, 1)
    % if gg ~= 2 && gg ~= 3 && gg ~= 6
    select_mouse = animalIDs{gg};

    select_mouse_index = find(strcmp(animalIDs, select_mouse));
    BehavData = final_behavior.(select_mouse).(first_session).uv.BehavData;


    ca = final.(select_mouse).(first_session).CNMFe_data.(ca_data_type);
    % ca = ca(prechoice_indices_for_PV{aa, select_mouse_index} == 1, :);
    % ca_zscored = zscore(ca, [], 2);
    ca_zscored = normalize(ca, 2);
    time_array = final.(select_mouse).(first_session).time;
    % time_array_all{aa, gg} = time_array;
    % Reshape the data and take the mean of each bin
    [n_neurons, n_samples] = size(ca_zscored);
    ca_binned = squeeze(mean(reshape(ca_zscored(:, 1:floor(n_samples/bin_factor)*bin_factor), n_neurons, bin_factor, []), 2));

    % Squeeze the data to remove the singleton dimension
    % ca_binned = squeeze(ca_binned);
    ca = ca_binned;
    % Assuming time_array is a 1D array of time points corresponding to each sample
    time_array_binned = mean(reshape(time_array(1:floor(length(time_array)/bin_factor)*bin_factor), bin_factor, []), 1);
    time_array = time_array_binned;





    
    % prechoice_block_1_column_1_data = {};
    % prechoice_mean_mouse = [];
    % prechoice_block_1_column_1_data = zall_mouse{select_mouse_index, 11}(prechoice_block_1_mouse{select_mouse_index, 1} == 1);
    % for cc = 1:size(prechoice_block_1_column_1_data, 2)
    %     prechoice_data = prechoice_block_1_column_1_data{1, cc};
    %     prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    % end
    % prechoice_block_1_mean_mouse_array{gg} = prechoice_mean_mouse;
    % 
    % prechoice_block_2_3_column_1_data = {};
    % prechoice_mean_mouse = [];
    % prechoice_block_2_3_column_1_data = zall_mouse{select_mouse_index, 11}(prechoice_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);
    % 
    % for cc = 1:size(prechoice_block_2_3_column_1_data, 2)
    %     prechoice_data = prechoice_block_2_3_column_1_data{1, cc};
    %     prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    % end
    % prechoice_block_2_3_mean_mouse_array{gg} = prechoice_mean_mouse;
    % 
    % prechoice_block_all_column_1_data = {};
    % prechoice_mean_mouse = [];
    % prechoice_block_all_column_1_data = zall_mouse{select_mouse_index, 11};
    % 
    % for cc = 1:size(prechoice_block_all_column_1_data, 2)
    %     prechoice_data = prechoice_block_all_column_1_data{1, cc};
    %     prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    % end
    % prechoice_all_mean_mouse_array{gg} = prechoice_mean_mouse;


    prechoice_remapped_column_1_data = {};
    prechoice_mean_mouse = [];
    prechoice_remapped_column_1_data = zall_mouse{select_mouse_index, 8}(prechoice_remapped_mouse{select_mouse_index, 1} == 1);
    for cc = 1:size(prechoice_remapped_column_1_data, 2)
        prechoice_data = prechoice_remapped_column_1_data{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    end
    prechoice_remapped_mean_mouse_block_2_3_array{gg} = prechoice_mean_mouse;

    prechoice_remapped_column_1_data = {};
    prechoice_mean_mouse = [];
    prechoice_remapped_column_1_data = zall_mouse{select_mouse_index, 1}(prechoice_remapped_mouse{select_mouse_index, 1} == 1);
    for cc = 1:size(prechoice_remapped_column_1_data, 2)
        prechoice_data = prechoice_remapped_column_1_data{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    end
    prechoice_remapped_mean_mouse_block_1_array{gg} = prechoice_mean_mouse;

    
    % postchoice_block_1_column_1_data = {};
    % postchoice_mean_mouse = [];
    % postchoice_block_1_column_1_data = zall_mouse{select_mouse_index, 11}(postchoice_reward_block_1_mouse{select_mouse_index, 1} == 1);
    % for cc = 1:size(postchoice_block_1_column_1_data, 2)
    %     postchoice_data = postchoice_block_1_column_1_data{1, cc};
    %     postchoice_mean_mouse(:, cc) = mean(postchoice_data (:, ts1 >= 0 & ts1 <= 2), 2);
    % end
    % postchoice_block_1_mean_mouse_array{gg} = postchoice_mean_mouse;
    % 
    % postchoice_block_2_3_column_1_data = {};
    % postchoice_mean_mouse = [];
    % postchoice_block_2_3_column_1_data = zall_mouse{select_mouse_index, 11}(postchoice_reward_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);
    % 
    % for cc = 1:size(postchoice_block_2_3_column_1_data, 2)
    %     postchoice_data = postchoice_block_2_3_column_1_data{1, cc};
    %     postchoice_mean_mouse(:, cc) = mean(postchoice_data (:, ts1 >= 0 & ts1 <= 2), 2);
    % end
    % postchoice_block_2_3_mean_mouse_array{gg} = postchoice_mean_mouse;
    % 
    % postchoice_block_all_column_1_data = {};
    % postchoice_mean_mouse = [];
    % postchoice_block_all_column_1_data = zall_mouse{select_mouse_index, 11};
    % 
    % for cc = 1:size(postchoice_block_all_column_1_data, 2)
    %     postchoice_data = postchoice_block_all_column_1_data{1, cc};
    %     postchoice_mean_mouse(:, cc) = mean(postchoice_data (:, ts1 >= 0 & ts1 <= 2), 2);
    % end
    % postchoice_all_mean_mouse_array{gg} = postchoice_mean_mouse;

    
    % collect_block_1_column_1_data = {};
    % collect_mean_mouse = [];
    % collect_block_1_column_1_data = zall_mouse{select_mouse_index, 12}(collect_block_1_mouse{select_mouse_index, 1} == 1);
    % for cc = 1:size(collect_block_1_column_1_data, 2)
    %     collect_data = collect_block_1_column_1_data{1, cc};
    %     collect_mean_mouse(:, cc) = mean(collect_data (:, ts1 >= 1 & ts1 <= 3), 2);
    % end
    % collect_block_1_mean_mouse_array{gg} = collect_mean_mouse;
    % 
    % collect_block_2_3_column_1_data = {};
    % collect_mean_mouse = [];
    % collect_block_2_3_column_1_data = zall_mouse{select_mouse_index, 12}(collect_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);
    % 
    % for cc = 1:size(collect_block_2_3_column_1_data, 2)
    %     collect_data = collect_block_2_3_column_1_data{1, cc};
    %     collect_mean_mouse(:, cc) = mean(collect_data (:, ts1 >= 1 & ts1 <= 3), 2);
    % end
    % collect_block_2_3_mean_mouse_array{gg} = collect_mean_mouse;
    % 
    % collecte_block_all_column_1_data = {};
    % collect_mean_mouse = [];
    % collect_block_all_column_1_data = zall_mouse{select_mouse_index, 12};
    % 
    % for cc = 1:size(collect_block_all_column_1_data, 2)
    %     collect_data = collect_block_all_column_1_data{1, cc};
    %     collect_mean_mouse(:, cc) = mean(collect_data (:, ts1 >= 1 & ts1 <= 3), 2);
    % end
    % collect_all_mean_mouse_array{gg} = collect_mean_mouse;

end

for ff = 1:size(prechoice_remapped_mean_mouse_block_1_array, 2)

    mean_prechoice_remapped_mouse_block_1(ff) = mean(mean(prechoice_remapped_mean_mouse_block_1_array{1, ff}))
    mean_prechoice_remapped_mouse_block_2_3(ff) = mean(mean(prechoice_remapped_mean_mouse_block_2_3_array{1, ff}))

end

diff_risky_vs_safe_prechoice_remapped = mean_prechoice_remapped_mouse_block_2_3-mean_prechoice_remapped_mouse_block_1;


%%
x = diff_risky_vs_safe_prechoice_remapped';
y = risk_table.Mean_1_to_3;
% x = remapped_prechoice_ratio(num_cells_mouse > 30)';
% y = risk_table.Mean_1_to_3(num_cells_mouse > 30);
% Create a new figure with specific dimensions
figure;
width = 250; % Width of the figure
height = 350; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]

hold on;
scatter(x, y, 100)
% Set the axis labels to have 2 decimal places
xtickformat('%.2f');
ytickformat('%.2f');
% Add a regression line (You can keep this part unchanged)
coefficients = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);
plot(x_fit, y_fit, 'r');

% Calculate R-squared value (You can keep this part unchanged)
y_pred = polyval(coefficients, x);
ssr = sum((y_pred - mean(y)).^2);
sst = sum((y - mean(y)).^2);
r_squared = ssr / sst;

% can also calculate the r-squared this way
% Calculate the R^2 value
[r, pval] = corrcoef(x, y); % Compute correlation coefficient matrix
rsq = r(1, 2)^2; % Extract and square the correlation coefficient
% Get axes limits
ax = gca;
xLimits = xlim(ax);
yLimits = ylim(ax);

% Set text position relative to axes limits
xPos = xLimits(2) - 0.05 * range(xLimits); % Slightly inside the top-right
yPos = yLimits(2) - 0.05 * range(yLimits); % Slightly inside the top-right

% Add R-squared value to the plot (You can keep this part unchanged)
text(xPos, yPos, ...
    {['R^2 = ' num2str(r_squared, '%.2f')], ...
     ['p = ' num2str(pval(2), '%.2f')]}, ...
    'FontSize', 12, ...
    'Color', 'blue', ...
    'HorizontalAlignment', 'right', ... % Align text to the right
    'VerticalAlignment', 'top');       % Align text to the top

hold off;


%% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav
% this outputs a ever so slightly wonky diagram. a few nodes that do not
% actually overlap minimally overlap (but intersections are 0), and 1 node
% that has 1 overlap does not overlap at all. 
conserved_prechoice_ind = find(conserved_prechoice == 1);
remapped_prechoice_ind = find(remapped_prechoice == 1);
shk_ind = find(respClass_all_array{1,4} == 1);
pre_choice_active_ind = find(respClass_all_array{1,1} == 1);
pre_choice_active_block_2_3_ind = find(respClass_all_array{1,8} == 1);
consum_active_ind = find(respClass_all_array{1,3} == 1);
consum_active_block_2_3 = find(respClass_all_array{1,10} == 1);
post_choice_active_ind = find(respClass_all_array{1,2} == 1);
aa_active_ind = find(respClass_all_array{1,11} == 1);
% consum_inhibited_ind = find(all_consum_inhibited == 1);
setListData = {shk_ind, conserved_prechoice_ind, remapped_prechoice_ind, aa_active_ind};
setLabels = ["Shk excited", "Prechoice conserved", "Prechoice new", "Approach-Abort excited"];
figure;
ve_diagram = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

ve_diagram.ShowIntersectionCounts = true;
ve_diagram.ShowIntersectionAreas = true;
% h.SetLabels = [];

%% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav
% this outputs a ever so slightly wonky diagram. a few nodes that do not
% actually overlap minimally overlap (but intersections are 0), and 1 node
% that has 1 overlap does not overlap at all. 
shk_ind_mouse = [];
aa_ind_mouse = [];
shk_ind_mouse = find(respClass_all_array_mouse{6, 4} == 1);
aa_ind_mouse = find(respClass_all_array_mouse{6, 11} == 1);
% remapped_prechoice_ind = find(remapped_prechoice == 1);
% shk_ind = find(respClass_all_array{1,4} == 1);
% pre_choice_active_ind = find(respClass_all_array{1,1} == 1);
% pre_choice_active_block_2_3_ind = find(respClass_all_array{1,8} == 1);
% consum_active_ind = find(respClass_all_array{1,3} == 1);
% consum_active_block_2_3 = find(respClass_all_array{1,10} == 1);
% post_choice_active_ind = find(respClass_all_array{1,2} == 1);

% consum_inhibited_ind = find(all_consum_inhibited == 1);
setListData = {shk_ind_mouse, aa_ind_mouse};
setLabels = ["Shk excited", "Approach-Abort excited"];
figure;
ve_diagram = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

ve_diagram.ShowIntersectionCounts = true;
ve_diagram.ShowIntersectionAreas = true;
% h.SetLabels = [];

for dd = 1:size(respClass_all_array_mouse, 1)
    shk_ind_mouse = find(respClass_all_array_mouse{dd, 4} == 1);
    aa_ind_mouse = find(respClass_all_array_mouse{dd, 11} == 1);
    overlap_between_shk_abort_neurons(dd) = size(intersect(shk_ind_mouse, aa_ind_mouse), 2)/size(respClass_all_array_mouse{dd, 4} == 1, 2);

end


