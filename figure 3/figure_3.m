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



shk_mean = mean(neuron_mean_array{1, 4}(:, ts1 > 0 & ts1 < 2),  2);

% [peak_values, time_of_peak_activity] = max(neuron_mean_array{1, 1}, [], 2);
[~, sort_indices] = sort(shk_mean);
neuron_mean_sorted = neuron_mean_array{1, 4}(sort_indices, :);


% Sort the rows of activated_neuron_mean based on peak_times.
% [~, sort_indices] = sort(time_of_peak_activity);
% activated_neuron_mean_sorted = activated_rows(sort_indices, :);

% Now, activated_neuron_mean_sorted contains the rows of neuron_mean filtered by respClass_all == 1
% and sorted by the time of peak activity.

figure;
% Generate the heatmap
imagesc(ts1, 1, neuron_mean_sorted);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

%%
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
% 

custom_colormap = [
    1, 1, 1;   % white
    1, 0.95, 0.9;
    1, 0.9, 0.8;
    1, 0.85, 0.7;
    1, 0.75, 0.55;
    1, 0.65, 0.4;
    1, 0.55, 0.3;
    1, 0.45, 0.2;
    1, 0.35, 0.1;
    1, 0.25, 0.05;
    1, 0.15, 0;  % orange
];

% Generate more intermediate colors for a smoother transition
n = 256; % Number of colors
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 350, 600]); % [left, bottom, width, height]
hold on
% Create a tiled layout with 2 rows and 1 column
% tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
% ax1 = nexttile;
% hold on;

% Plot the heatmap
imagesc(ts1, 1, neuron_mean_sorted);

% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 
ylim([1, neuron_num]);
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
set(gca, 'YTick', [1, neuron_num]);
xline(0)
% scatter(time2Collect, Tris               , 'Marker', 'p')
% scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

%% plot shk activated neurons

mean_data_array = {neuron_mean_array{1, 4}(respClass_all_array{1,4} == 1, :)};
sem_data_array = {neuron_sem_array{1, 4}(respClass_all_array{1,4} == 1, :), neuron_sem_array{1, 10}(conserved_consumption ~= 1 & shk_activated == 1, :)};

[comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1)


%%

exclusive_shk_activated = respClass_all_array{1,4} == 1 & respClass_all_array{1,1} == 3 &respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 3;
exclusive_collection_activated = respClass_all_array{1,4} == 3 & respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 1;
shk_event = respClass_all_array{1,4} == 1;
post_choice_rew_event = respClass_all_array{1,2} == 1;
post_choice_both_excited = respClass_all_array{1,2} == 1 & respClass_all_array{1,4} == 1;
% this is the start of checking if neurons are MORE active than during
% other events, i.e. if you wanted to check if REW activated neurons are
% significantly differentially activated by SHK. preliminary poking around
% seems to suggest that few large reward active neurons have their activity
% increase in response to SHK
co_activated_indices = find(post_choice_both_excited(1,:) == 1);
co_activated_indices_sum = numel(co_activated_indices);
clear h p ci stats
for qq = 1:size(co_activated_indices, 2)
    [h(qq),p(qq),ci{qq},stats{qq}] = ttest(neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx),neuron_mean_array{1, 4}(co_activated_indices(qq),sub_window_idx));
    mean_diff(qq) = mean(neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx) - mean(neuron_mean_array{1, 4}(co_activated_indices(qq),sub_window_idx)));
end

sig_increase_shk_from_large = co_activated_indices(p < 0.05 & mean_diff < 0);
sig_increase_shk_from_large_sum = numel(sig_increase_shk_from_large);
sig_increase_shk_from_large_ind = zeros(1, size(respClass_all_array{1,4}, 2));

sig_increase_shk_from_large_ind(:, sig_increase_shk_from_large) = 1;

no_sig_increase_shk_from_large = co_activated_indices(p > 0.05);
no_sig_increase_shk_from_large_sum = numel(no_sig_increase_shk_from_large);
no_sig_increase_shk_from_large_ind = zeros(1, size(respClass_all_array{1,4}, 2));
no_sig_increase_shk_from_large_ind(:, no_sig_increase_shk_from_large ) = 1;

shk_activated = respClass_all_array{1,4} == 1 & no_sig_increase_shk_from_large_ind ~= 1
shk_activated_sum = sum(shk_activated);

figure; plot(ts1, nanmean(neuron_mean_array{1, 4}(shk_activated,:))); hold on; plot(ts1,  nanmean(neuron_mean_array{1, 2}(respClass_all_array{1,2} == 1,:)));


figure; shadedErrorBar


total_modulated = [(sum(post_choice_rew_event)/neuron_num)*100 (sum(shk_event)/neuron_num)*100];

A = total_modulated;
I = (co_activated_indices_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end

total_modulated = [(sum(shk_event)/neuron_num)*100 (sum(post_choice_rew_event)/neuron_num)*100];
A = total_modulated;
I = (co_activated_indices_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end




shk_event = respClass_all_array{1,4} == 1;

consumption_event = respClass_all_array{1, 10} == 1;

total_modulated = [(sum(shk_event)/neuron_num)*100 (sum(consumption_event)/neuron_num)*100];

shk_and_consum_both_excited = respClass_all_array{1,10} == 1 & respClass_all_array{1,4} == 1;
% this is the start of checking if neurons are MORE active than during
% other events, i.e. if you wanted to check if REW activated neurons are
% significantly differentially activated by SHK. preliminary poking around
% seems to suggest that few large reward active neurons have their activity
% increase in response to SHK
co_activated_indices = find(shk_and_consum_both_excited(1,:) == 1);
co_activated_indices_sum = numel(co_activated_indices);


A = total_modulated;
I = (co_activated_indices_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end

% Calculate different groups as percentages
only_shk = (sum(shk_event)/neuron_num)*100 - (co_activated_indices_sum/neuron_num)*100; % SHK only
only_consumption = (sum(consumption_event)/neuron_num)*100 - (co_activated_indices_sum/neuron_num)*100; % Consumption only
both = (co_activated_indices_sum/neuron_num)*100; % Overlap between SHK and Consumption
not_modulated = 100 - (only_shk + only_consumption + both); % Unmodulated neurons

% Data for the stacked bar (bottom to top: Not modulated, Consumption only, Both, SHK only)
data_for_bar_plot = [not_modulated, only_consumption, both, only_shk];

% Create the stacked bar plot
figure;
bar(1, data_for_bar_plot, 'stacked'); % Single bar at x = 1
colormap([0.7 0.7 0.7; 0 0 1; 0.8 0.8 0; 1 0 0]); % Grey (Unmodulated), Blue (Consumption), Yellow (Both), Red (SHK)

% Formatting
xticks(1); % Single bar on x-axis
xticklabels({'Neuron Modulation'});
ylabel('Percentage of Neurons (%)');
title('Neuron Modulation by Events');
ylim([0 100]); % Ensure the bar always reaches 100%
legend({'Not Modulated', 'Consumption Only', 'Both', 'SHK Only'}, 'Location', 'eastoutside');

% Adding percentage labels
y_offset = 0; % Start at the base of the bar
for i = 1:length(data_for_bar_plot)
    if data_for_bar_plot(i) > 0  % Only add labels for non-zero sections
        text(1, y_offset + data_for_bar_plot(i)/2, sprintf('%.1f%%', data_for_bar_plot(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', 'white', 'FontWeight', 'bold');
        y_offset = y_offset + data_for_bar_plot(i); % Move up for next section
    end
end

for zz = 1:size(respClass_all_array_mouse, 1)
    exclusive_shk_activated_mouse{zz} = respClass_all_array_mouse{zz,4} == 1 & respClass_all_array_mouse{zz,1} == 3 & respClass_all_array_mouse{zz,2} == 3 & respClass_all_array_mouse{zz,3} == 3;
    sum_exclusive_shk_activated_mouse(zz) = sum(exclusive_shk_activated_mouse{zz});
    percent_exclusive_shk_activated_mouse(zz) = sum_exclusive_shk_activated_mouse(zz)/size(exclusive_shk_activated_mouse{zz}, 2);
    max_activity_exclusive_shk(zz) = max(mean(neuron_mean_mouse{zz, 4}(exclusive_shk_activated_mouse{1, zz}   == 1, :)));
end

% Calculate different groups as percentages
only_shk_sum = sum(shk_event) - co_activated_indices_sum; % SHK only
only_consumption_sum = sum(consumption_event) -co_activated_indices_sum; % Consumption only
percent_shk_only = (only_shk_sum/neuron_num)*100
percent_consum_only = (only_consumption_sum/neuron_num)*100

%% for heatmap, change "plot_num" (what neuron to plot) and array_to_plot (# corresponds to which dataset)

plot_num = 147; 

array_to_plot = 4; % depends on the structure of zall

select_mouse = 'BLA_Insc_40';

% for RDT D1 BLA_Insc_25:
%prechoice neuron num 46
%postchoice rew num 38
%consumption num 39
%shock num 11


% for RDT D1 BLA_Insc_40:
%prechoice neuron num 12
%postchoice rew num 70
%consumption num 10
%shock num 11


select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';

second_session = 'RDT_D1';



% in order to trim off excess (because calcium recording starts before
% behavior), you'll need the start time. unfortunately I don't save that
% variable in the "final" struct, but you can get it from the adjustment to
% some of the columns of BehavData. e.g., see below - but make sure to
% update the session etc as necessary! 

% BehavData = final.(select_mouse).(first_session).choiceTime.uv.BehavData;
% BehavData = final.(select_mouse).(first_session).uv.BehavData;
% because the first trial possible is ALWAYS 60 seconds after ABET is
% issued, you can determine what adjustment has been made to this column
% (adding time to account for calcium recording starting first) by
% subtracting off 60 from the first element
% stTime = BehavData.TrialPossible(1)-60; 



% time2Collect = BehavData.collectionTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
% trialStartTime = BehavData.stTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
% median_trialStartTime = median(trialStartTime)
% median_time2Collect = median(time2Collect)
% xline(median_trialStartTime)
% xline(median_time2Collect)
% [numTrials, ~] = size(time2Collect);
% Tris = [1:numTrials]';

% Define the custom colormap from white to orange
% custom_colormap = [
%     1, 1, 1; % white
%     1, 0.9, 0.8;
%     1, 0.8, 0.6;
%     1, 0.7, 0.4;
%     1, 0.6, 0.2;
%     1, 0.5, 0; % orange
% ];


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

% custom_colormap = [
%     1, 1, 1;         % white
%     0.9, 0.9, 0.95;
%     0.8, 0.8, 0.9;
%     0.6, 0.6, 0.8;
%     0.4, 0.4, 0.7;
%     0.2, 0.2, 0.6;
%     0.0, 0.0, 0.55;   % dark blue
% ];

custom_colormap = [
    1, 1, 1; % white
    1, 0.9, 0.9;
    1, 0.8, 0.8;
    1, 0.7, 0.7;
    1, 0.6, 0.6;
    1, 0.5, 0.5;
    1, 0.4, 0.4;
    1, 0.3, 0.3;
    1, 0.2, 0.2;
    1, 0.1, 0.1;
    1, 0, 0;   % red
];


% Generate more intermediate colors for a smoother transition
n = 256; % Number of colors
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 200, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
ax1 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, 1:size(zall_array{array_to_plot, plot_num}, 1), zall_array{array_to_plot, plot_num});


% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar(ax1, 'eastoutside');
set(c, 'YTick', clim); % 

ylim([0.5, size(zall_array{array_to_plot, plot_num}, 1) + 0.5]);

xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', []);
set(gca, 'YTick', [1, size(zall_array{array_to_plot, plot_num}, 1)]);
xline(0)
% scatter(time2Collect, Tris               , 'Marker', 'p')
% scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;


% Plot the raw data in grey with transparency
for trial = 1:size(zall_array{array_to_plot, plot_num}, 1)
    plot(ts1, zall_array{array_to_plot, plot_num}(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
end

% Plot the mean as a thick black line
meanData = mean(zall_array{array_to_plot, plot_num});
plot(ts1, meanData, 'r', 'LineWidth', 2, 'Color', 'k');

ylim([-3 4]);
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
xline(0)
yline(0)
fontsize(18, 'points')
hold off;

%% are consumption neurons inhibited by Shk? are Shk neurons inhibited during consumption?



figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
ylim([-0.5 1.2]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8], 'YTick', [-0.5 0 0.5 1]);
% shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(respClass_all_array{1, 4}==1, :)), nanmean(neuron_sem_array{1, 4}(respClass_all_array{1, 4}==1, :)), 'lineProps', {'color', 'r'});
% shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(respClass_all_array{1, 4}==1, :)), std(neuron_mean_array{1, 4}(respClass_all_array{1, 4}==1, :))/sqrt(size(neuron_mean_array{1, 4}(respClass_all_array{1, 4}==1, :), 1)), 'lineProps', {'color', 'r'});
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(respClass_all_array{1, 4}==1, :)), std(neuron_mean_array{1, 4}(respClass_all_array{1, 4}==1, :)), 'lineProps', {'color', 'r'});
% hold on;shadedErrorBar(ts1, nanmean(zall_mean_all_array{1, 12}(respClass_all_array{1, 11}==1, :)), nanmean(sem_all_array{1, 12}(respClass_all_array{1, 11}==1, :)), 'lineProps', {'color', 'k'});

xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')

hold off
%%

%These data can be used to plot the median or mean choice
% time on a PCA graph, for example

behav_tbl_iter_single = behav_tbl_iter(1);

% Initialize the concatenated table
concatenatedTable = table();

% Iterate through the 3x1 cell array
for i = 1:numel(behav_tbl_iter_single)
    % Assuming each cell contains a 12x1 cell array of tables
    twelveByOneCellArray = behav_tbl_iter_single{i};
    
    % Initialize a temporary table to store the concatenated tables for this cell
    tempTable = table();
    
    % Iterate through the 12x1 cell array
    for j = 1:numel(twelveByOneCellArray)
        % Assuming each cell in the 12x1 cell array contains a table
        currentTable = twelveByOneCellArray{j};
        
        % Concatenate the current table to the temporary table vertically
        tempTable = vertcat(tempTable, currentTable);
    end
    
    % Concatenate the temporary table to the overall concatenated table vertically
    concatenatedTable = vertcat(concatenatedTable, tempTable);
end

median_choice_time_block_1 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_choice_time_block_2 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_choice_time_block_3 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));

median_collect_time_block_1 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_collect_time_block_2 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_collect_time_block_3 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));





median_start_time_from_choice = median(concatenatedTable.stTime - concatenatedTable.choiceTime);
median_collect_time_from_choice = median(concatenatedTable.collectionTime - concatenatedTable.choiceTime);


%%


figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(prechoice_block_1==1, :)), nanmean(neuron_sem_array{1, 8}(prechoice_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 9}(postchoice_reward_block_1==1, :)), nanmean(neuron_sem_array{1, 9}(postchoice_reward_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 8}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 10}(collect_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(respClass_all_array{1,4}==1, :)), nanmean(neuron_sem_array{1, 10}(respClass_all_array{1,4}==1, :)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')
ylim([-0.6 1.0]);
hold off

%%
prechoice_to_shk = prechoice_block_1 == 1 & respClass_all_array{1,4} == 1;
sum(prechoice_to_shk)
postchoice_reward_to_shk = postchoice_reward_block_1 == 1 & respClass_all_array{1,4} == 1;
sum(postchoice_reward_to_shk)
collect_to_shk = collect_block_1 == 1 & respClass_all_array{1,4} == 1;
sum(collect_to_shk)

nonresp_all = prechoice_block_1 == 0 & postchoice_reward_block_1 == 0 & collect_block_1 == 0;
sum(nonresp_all)

sum([sum(nonresp_all) sum(collect_block_1) sum(postchoice_reward_block_1) sum(prechoice_block_1)])
nonresp_to_shk = nonresp_all == 1 & respClass_all_array{1,4} == 1;
sum(nonresp_to_shk)

%% run eventRelated w/ SHK, 1 and LOSS_PLUS_ONE 1, then run code below
for q = 1:length (behav_tbl_iter{4, 1})
    nestedCellArray_1 = behav_tbl_iter{4, 1}{q};
    if ~isempty(nestedCellArray_1)
        nestedCellArray_2 = behav_tbl_iter{4, 1}{q};
        if size(nestedCellArray_1, 1) > size(nestedCellArray_2, 1)
            delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime(1:end-1,:);
        else
            delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
        end


        for zz = 1:size(nestedCellArray_1, 1)
            valid_start_times = nestedCellArray_1.stTime(2:end);
            valid_choice_times = nestedCellArray_1.choiceTime(1:end-1);
            % delay_to_initiation = valid_start_times - valid_choice_times;
            trial_types = nestedCellArray_1.bigSmall;
            trial_types_2 = nestedCellArray_2.bigSmall;
        end

        trial_choice_times = nestedCellArray_1.choiceTime - nestedCellArray_1.stTime;
        % delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
        delay_to_collect_post_shk = nestedCellArray_1.collectionTime - nestedCellArray_1.choiceTime;
        trial_choice_times_by_mouse{q} = trial_choice_times;
        delay_to_initiation_by_mouse{q} = delay_to_initiation;
        delay_to_collect_post_shk_by_mouse{q} = delay_to_collect_post_shk;
        trial_types_by_mouse{q} = trial_types;
        trial_types_second_var_by_mouse{q} = trial_types_2;
        clear trial_choice_times delay_to_initiation delay_to_collect_post_shk trial_types trial_types_2
    end


end

trial_choice_times_concat = cat(1, trial_choice_times_by_mouse{:});
rew_collect_times_concat = cat(1, delay_to_collect_post_shk_by_mouse{:});
delay_to_initiation_concat = cat(1, delay_to_initiation_by_mouse{:});
delay_to_collect_concat = cat(1, delay_to_collect_post_shk_by_mouse{:});

bar_separation_value = 3;

figure;
width = 100; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]

swarmchart(ones(1, length(delay_to_collect_concat)), delay_to_collect_concat);
hold on;
yline(mean(delay_to_collect_concat), 'k', 'LineWidth', 2); % Black horizontal line at the mean

figure;
width = 100; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]

swarmchart(ones(1, length(delay_to_initiation_concat)) * bar_separation_value, delay_to_initiation_concat);
yticks([0 50 100 150 200])
hold on;
yline(mean(delay_to_initiation_concat), 'k', 'LineWidth', 2); % Black horizontal line at the mean

yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off;


variable_to_correlate = delay_to_collect_post_shk_by_mouse;

%%
array_for_means = 4; 

% Initialize the new cell array to store the mean values
meanZallMouse = cell(size(zall_mouse, 2), 1);

% Define the time range for 0 to 2 seconds
% timeRange = (ts1 >= -4) & (ts1 <= 0);
timeRange = (ts1 >= 0) & (ts1 <= 5); % use a longer time window here so that the entirety of the shk activity can be used
% timeRange = (ts1 >= 1) & (ts1 <= 3);


% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_mouse)
    % Get the current nested cell array
    nestedCellArray_1 = zall_mouse{i, array_for_means};
    nestedCellArray_2 = zall_mouse{i, 2};
    % Initialize the nested cell array for storing mean values
    meanNestedCellArray = cell(size(nestedCellArray_1));
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(nestedCellArray_1)
        % Get the current double array
        currentArray = nestedCellArray_1{j};
        comparisonArray_for_size = nestedCellArray_2{j};
        
        if isequal(variable_to_correlate, delay_to_collect_post_shk_by_mouse)
            currentArray = nestedCellArray_1{j};
        else
            if size(currentArray, 1) > size(comparisonArray_for_size, 1)
                currentArray = currentArray(1:end-1,:);
            else

            end
        end
        % uncomment below if you want to mean center
        % currentArray_mean = mean(currentArray, 2);
        % currentArray = currentArray-currentArray_mean;
        % Compute the mean activity for each row in the time range 0 to 2 seconds
        meanValues = mean(currentArray(:, timeRange), 2);
        % meanValues = max(currentArray(:, timeRange), [], 2);
        
        % Store the mean values in the corresponding cell of the nested cell array
        meanNestedCellArray{j} = meanValues;
    end
    
    % Store the nested cell array of mean values in the corresponding cell of the main cell array
    meanZallMouse{i} = meanNestedCellArray;
end

% Now, meanZallMouse contains the mean activity for each row in the time period 0 to 2 seconds
% Each cell in meanZallMouse contains a nested cell array with the

%%
%%

% Assuming the following variables are defined:
% meanZallMouse: 14x1 cell array where each cell contains another cell array with mean values
% trial_choice_times_by_mouse: 1x11 cell array containing values to correlate with

% Initialize the new cell array to store the correlation results
correlationResults = cell(size(meanZallMouse));
correlationResults_sig = cell(size(meanZallMouse));


% Iterate through each level of meanZallMouse
for i = 1:length(meanZallMouse)
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{i};
    
    % Initialize the nested cell array for storing correlation results
    correlationNestedArray = zeros(size(meanNestedCellArray));
    corr_sig_NestedArray = zeros(size(meanNestedCellArray));
    
    % Determine the corresponding index in trial_choice_times_by_mouse
    % Adjust this logic based on how the indices are mapped
    trialIndex = mod(i-1, length(variable_to_correlate)) + 1;
    
    % Get the corresponding trial choice times
    trialChoiceTimes = variable_to_correlate{i};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        % Get the current mean values array
        meanValues = meanNestedCellArray{j};
        
        % Check if trialChoiceTimes has the same length as meanValues
        if length(trialChoiceTimes) == length(meanValues)
            % Compute the correlation
            [correlationCoeff, corr_sig_vals] = corr(meanValues, trialChoiceTimes(:));
        elseif length(trialChoiceTimes) < length(meanValues)
            [correlationCoeff, corr_sig_vals] = corr(meanValues(1:end-1), trialChoiceTimes(:));
        else
            % If lengths do not match, handle the mismatch (e.g., set correlation to NaN)
            [correlationCoeff, corr_sig_vals] = NaN;
        end
        
        % Store the correlation coefficient in the nested cell array
        correlationNestedArray(j) = correlationCoeff;
        corr_sig_NestedArray(j) = corr_sig_vals;
    end
    clear meanValues
    % Store the nested cell array of correlation coefficients in the main cell array
    correlationResults{i} = correlationNestedArray;
    correlationResults_sig{i} = corr_sig_NestedArray;
end

%%
% Assuming correlationResults is defined and contains the correlation coefficients

% Initialize an empty array to collect all correlation coefficients
allCorrelations = [];

% Iterate through each level of correlationResults
for i = 1:length(correlationResults)
    % Get the current nested cell array of correlation coefficients
    correlationNestedArray = correlationResults{i};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(correlationNestedArray)
        % Get the current correlation coefficient
        correlationCoeff = correlationNestedArray(j);
        
        % Check if the coefficient is not NaN (if applicable)
        if ~isnan(correlationCoeff)
            % Append the coefficient to the allCorrelations array
            allCorrelations = [allCorrelations; correlationCoeff];
        end
    end
end

% Now, allCorrelations contains all the correlation coefficients
% Create a histogram of the correlation coefficients
figure;
histogram(allCorrelations);
xlabel('Correlation Coefficient');
ylabel('Frequency');
title('Histogram of Correlation Coefficients');

% Optionally, you can add a vertical line at 0 for reference
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;


%%
%% SHK responsive neurons assumed to be stored in respClass_all_array{1, 1} for this purpose - change as necessary
only_shk_responsive_corrs = allCorrelations(respClass_all_array{1, 4}  ==1);
not_shk_responsive_corrs = allCorrelations(respClass_all_array{1, 4}  ==3);
% Now, allCorrelations contains all the correlation coefficients
% Create a histogram of the correlation coefficients
figure;
histogram(only_shk_responsive_corrs);
xlabel('Correlation Coefficient');
ylabel('Frequency');
title('Histogram of Correlation Coefficients');

% Optionally, you can add a vertical line at 0 for reference
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;
%%
% Assuming the following variables are defined:
% allCorrelations: array containing correlation coefficients from correlationResults
% only_shk_responsive_corrs: array containing correlation coefficients from a different variable

% Calculate means
mean_only_shk = mean(only_shk_responsive_corrs);
mean_not_shk = mean(not_shk_responsive_corrs);

% Create a histogram for allCorrelations
figure;
width = 250; % Width of the figure
height = 250; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
histogram(not_shk_responsive_corrs , 'Normalization', 'probability', 'FaceColor', 'blue','BinWidth', 0.05,'LineStyle','none');
hold on;

% Create a histogram for only_shk_responsive_corrs on the same figure
histogram(only_shk_responsive_corrs, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05, 'LineStyle','none');
xline(mean_only_shk, 'r')
xline(mean_not_shk, 'g')
% Add labels and title
xlabel('Correlation Coefficient');
ylabel('Probability');
% title('Histograms of Correlation Coefficients');
% legend('All Correlations', 'Only SHK Responsive Correlations');

% Optionally, you can add a vertical line at 0 for reference
yLimits = ylim;
plot([0 0], yLimits, 'k', 'LineWidth', 2);
xtickformat('%.2f');
ytickformat('%.2f');
hold off;

% Perform a Kolmogorov-Smirnov test to compare the two distributions
[h, p, k] = kstest2(not_shk_responsive_corrs , only_shk_responsive_corrs)

% % Display the results of the statistical test
% fprintf('Kolmogorov-Smirnov test result:\n');
% fprintf('h = %d (0 means the null hypothesis cannot be rejected, 1 means it can be rejected)\n', h);
% fprintf('p-value = %.4f\n', p);

[h,p,ci,stats] = ttest2(not_shk_responsive_corrs , only_shk_responsive_corrs)
%%





% % Plot means as bars
% figure;
% hold on;
% bar(1, mean_only_shk, 'FaceColor', 'r'); % Red bar for 'only_shk_responsive_corrs'
% bar(2, mean_not_shk, 'FaceColor', 'b'); % Blue bar for 'not_shk_responsive_corrs'
% 
% swarmchart(1, only_shk_responsive_corrs, 5)
% 
% % Scatter individual data points
% scatter(ones(size(only_shk_responsive_corrs)), only_shk_responsive_corrs, 'r', 'filled', 'jitter', 'on', 'jitterAmount', 0.15); % Red points with jitter
% scatter(2 * ones(size(not_shk_responsive_corrs)), not_shk_responsive_corrs, 'b', 'filled', 'jitter', 'on', 'jitterAmount', 0.15); % Blue points with jitter
% 
% % Customize plot
% xlim([0.5, 2.5]);
% xticks([1 2]);
% xticklabels({'Only Shk Responsive', 'Not Shk Responsive'});
% ylabel('Correlation Values');
% title('Correlation Values and Means');
% grid on;
% 
% % Add legend
% legend({'Mean Only Shk', 'Mean Not Shk', 'Individual Only Shk', 'Individual Not Shk'}, 'Location', 'best');
% 
% hold off;

bar_separation_value = 3;

figure;
width = 250; % Width of the figure
height = 250; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
swarmchart(ones(1, length(only_shk_responsive_corrs)), only_shk_responsive_corrs)
hold on
swarmchart(ones(1, length(not_shk_responsive_corrs))*bar_separation_value, not_shk_responsive_corrs)

% yline(mean(only_shk_responsive_corrs), ones(length(only_shk_responsive_corrs)))
plot([0.5; 1.5], [mean(only_shk_responsive_corrs); mean(only_shk_responsive_corrs)], 'LineWidth',3)
plot([bar_separation_value-.5; bar_separation_value+.5], [mean(not_shk_responsive_corrs); mean(not_shk_responsive_corrs)], 'LineWidth',3)
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off


%%
high_corrs = find(allCorrelations > 0.5)



get_mouse_ids_for_high_corrs = mouse_cells(1, high_corrs);


select_mouse_index = find(strcmp(animalIDs, select_mouse));


%% plot scatters for individual neurons & behav variables


%shk representative:
find(correlationResults{5, 1} > 0.5)
start_time = 0;% sub-window start time
end_time = 5; % sub-window end time
sub_window_idx = ts1 >= start_time & ts1 <= end_time;
sub_window_activity_session_1 = zall_mouse{5, 4}{1, 104}(:, sub_window_idx);
r_val_for_representative = correlationResults{5, 1}(1, 104)
p_val_for_representative = correlationResults_sig{5, 1}(1, 104)
choice_times_mouse = variable_to_correlate{1, 5};
% trial_types = trial_types_by_mouse{1, 5};
trial_types = trial_types_second_var_by_mouse{1, 5};



%CREATE SCATTER PLOT BASED ON SPECIFIC EVENTS - ASSUMING THEY ARE IN PAIRS.
%CHECK AND UPDATE START & END TIME DEPENDING ON EVENT OF INTEREST
paired_neurons = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 1;
% start_time = 0;% sub-window start time
% end_time = 2; % sub-window end time

% Extract the corresponding columns from neuron_mean


% sub_window_activity_session_1 = zall_mouse{5, 1}{1, 104}(:, sub_window_idx);
% choice_times_mouse = trial_choice_times_by_mouse{1, 5};
% trial_types = trial_types_by_mouse{1, 5};

% % Assume A and B are your 143x21 arrays
% correlation_coefficients = arrayfun(@(i) corr(sub_window_activity_session_1 (i, :)', sub_window_activity_session_2 (i, :)'), 1:size(sub_window_activity_session_1 , 1));


mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);


x = mean_sub_window_activity_session_1;
y = choice_times_mouse;


% x = mean_sub_window_activity_session_1(trial_types == 1.2);
% y = choice_times_mouse(trial_types == 1.2);
% trial_types = trial_types(trial_types == 1.2);


% Define colors based on trial types
colors = repmat([0.5, 0.5, 0.5], length(trial_types), 1); % Default to gray
colors(trial_types == 1.2, :) = repmat([0, 0, 1], sum(trial_types == 1.2), 1); % Blue for trial_types == 1.2
colors(trial_types == 0.3, :) = repmat([1, 0, 0], sum(trial_types == 0.3), 1); % Red for trial_types == 0.3

% Create scatter plot
figure;
set(gcf, 'Position', [100, 100, 200, 200]); % Adjust figure position and size
scatter(x, y, 36, colors, 'filled', 'MarkerEdgeColor', 'k'); % Use 'colors' for MarkerFaceColor

hold on;

% Add a regression line
coefficients = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);
plot(x_fit, y_fit, 'r');

% Calculate R-squared value
y_pred = polyval(coefficients, x);
ssr = sum((y_pred - mean(y)).^2);
sst = sum((y - mean(y)).^2);
r_squared = ssr / sst;

% Add R-squared value to the plot
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared)], 'FontSize', 12);

% Add labels and a legend if needed
xlabel('Mean Sub-window Activity Session 1');
ylabel('Choice Times Mouse');
title('Scatter Plot with Regression Line and R^2 Value');
hold off;


%% are consumption neurons inhibited by Shk? are Shk neurons inhibited during consumption?



figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
ylim([-0.5 1.2]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8], 'YTick', [-0.5 0 0.5 1]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(collect_block_1==1, :)), nanmean(neuron_sem_array{1, 4}(collect_block_1==1, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 10}(respClass_all_array{1, 4}==1, :)), nanmean(neuron_sem_array{1, 10}(respClass_all_array{1, 4}==1, :)), 'lineProps', {'color', 'k'});

xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')

hold off