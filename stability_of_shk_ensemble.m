%% load matched RDT_D1 & RDT_D2 data
% filter on SHK, 1

shk_consistent = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1;
sum(shk_consistent)
figure; plot(ts1, neuron_mean_array{1,1}(shk_consistent, :)); 
figure; plot(ts1, neuron_mean_array{1,2}(shk_consistent, :)); 


figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(shk_consistent==1, :)), nanmean(neuron_sem_array{1, 1}(shk_consistent==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 2}(shk_consistent==1, :)), nanmean(neuron_sem_array{1, 2}(shk_consistent==1, :)), 'lineProps', {'color', batlowW(iter,:)});
% xtickformat('%.2f');
ytickformat('%.2f');
xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
% legend({'pre-choice active', 'post-choice reward active', 'consumption', 'not active'}, 'Location','northwest')

%%

