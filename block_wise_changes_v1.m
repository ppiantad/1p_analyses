% Pre_RDT_RM & RDT D1 matched
% RDT D1
% collectionTime
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 1
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 2, 'BLOCK', 3

collect_all_possible = sum([sum(respClass_all_array{1, 1} ==1), sum(respClass_all_array{1, 2} ==1)])
collect_conserved = respClass_all_array{1, 1}==1 & respClass_all_array{1, 2}==1;
collect_conserved_sum = sum(collect_conserved)
collect_original_sum = sum(respClass_all_array{1, 1}==1)
collect_remapped = respClass_all_array{1, 1} ~=1 & respClass_all_array{1, 2}==1;
collect_remapped_sum = sum(collect_remapped)
not_collect_sum = sum(respClass_all_array{1, 1} ~=1)

collect_lost = respClass_all_array{1, 1} ==1 & respClass_all_array{1, 2} ~=1;
collect_lost_sum = sum(collect_lost)
collect_all_possible = sum([sum(respClass_all_array{1, 1} ==1), sum(respClass_all_array{1, 2} ==1)])
collect_conserved_percent = collect_conserved_sum/neuron_num
collect_remapped_percent = collect_remapped_sum/neuron_num
collect_lost_percent = collect_lost_sum/neuron_num

inner_pie = [collect_conserved_sum/neuron_num,...
            
            collect_lost_sum/neuron_num,...

            collect_remapped_sum/neuron_num,...
           
            (neuron_num- [collect_conserved_sum+collect_remapped_sum+collect_lost_sum]) / neuron_num];


labels = ["conserved", "lost", "new (remapped)",  "other"];

figure; donutchart(inner_pie, labels, 'InnerRadius', 0.5)

figure; 
N = 5;
axes('ColorOrder',brewermap(N,'Pastel2'),'NextPlot','replacechildren')
plot(ts1, mean(neuron_mean_array{1, 1}(collect_conserved  ==1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, 1}(collect_remapped  ==1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, 2}(collect_remapped  ==1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, 2}(collect_lost  ==1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, 1}(collect_lost  ==1, :)))


figure; hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,1}(collect_conserved  ==1, :)), nanmean(neuron_sem_array{1, 1}(collect_conserved  ==1, :)), 'lineProps', {'color', acton(1,:)});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,1}(collect_remapped  ==1, :)), nanmean(neuron_sem_array{1, 1}(collect_remapped  ==1, :)), 'lineProps', {'color', acton(50,:)});
h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,2}(collect_remapped  ==1, :)), nanmean(neuron_sem_array{1, 2}(collect_remapped  ==1, :)), 'lineProps', {'color', acton(150,:)});
% h(4) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,2}(collect_lost  ==1, :)), nanmean(neuron_sem_array{1, 2}(collect_lost  ==1, :)), 'lineProps', {'color', acton(200,:)});
% h(5) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,1}(collect_lost  ==1, :)), nanmean(neuron_sem_array{1, 1}(collect_lost  ==1, :)), 'lineProps', {'color', acton(250,:)});
legend([h(1).mainLine h(2).mainLine h(3).mainLine h(4).mainLine h(5).mainLine], 'conserved', 'pre-remapped (safe blocks)', 'post-remapped (risky blocks)')
xlim([-8 8]);
ylim([-0.6 0.8])

% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');

%%
figure; hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,1}(collect_conserved  ==1, :)), nanmean(neuron_sem_array{1, 1}(collect_conserved  ==1, :)), 'lineProps', {'color', acton(1,:)});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,2}(collect_conserved  ==1, :)), nanmean(neuron_sem_array{1, 2}(collect_conserved  ==1, :)), 'lineProps', {'color', acton(200,:)});
legend([h(1).mainLine h(2).mainLine], 'conserved (safe block)', 'conserved (risky blocks)')
xlim([-8 8]);
ylim([-0.6 0.7])

%%
figure; hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,1}(collect_lost  ==1, :)), nanmean(neuron_sem_array{1, 1}(collect_lost  ==1, :)), 'lineProps', {'color', acton(1,:)});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,2}(collect_lost  ==1, :)), nanmean(neuron_sem_array{1, 2}(collect_lost  ==1, :)), 'lineProps', {'color', acton(200,:)});
legend([h(1).mainLine h(2).mainLine], 'lost (safe block)', 'lost (risky blocks)')
xlim([-8 8]);
ylim([-0.6 0.7])

%%
figure; hold on;
h(1) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,1}(collect_remapped  ==1, :)), nanmean(neuron_sem_array{1, 1}(collect_remapped  ==1, :)), 'lineProps', {'color', acton(1,:)});
h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,2}(collect_remapped  ==1, :)), nanmean(neuron_sem_array{1, 2}(collect_remapped  ==1, :)), 'lineProps', {'color', acton(200,:)});
legend([h(1).mainLine h(2).mainLine], 'new (safe block)', 'new (risky blocks)')
xlim([-8 8]);
ylim([-0.6 0.7])