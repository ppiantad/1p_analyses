

%% Fig. 4A 
% for heatmap, change "plot_num" (what neuron to plot) and array_to_plot (# corresponds to which dataset)

plot_num = 18 % 81 / 31 or 58 or 70 / 2

array_to_plot = [1 8]; % depends on the structure of zall

select_mouse = 'BLA_Insc_34';


% for RDT D1 BLA_Insc_34:
%prechoice neuron num 18
%postchoice rew num 104
%consumption num 2


select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';

second_session = 'RDT_D1';



BehavData = final_behavior.(select_mouse).(first_session).uv.BehavData;



time2Collect = BehavData.collectionTime(trials_per_mouse{select_mouse_index, array_to_plot(1)}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot(1)});
trialStartTime = BehavData.stTime(trials_per_mouse{select_mouse_index, array_to_plot(1)}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot(1)});
median_trialStartTime = median(trialStartTime)
median_time2Collect = median(time2Collect)
xline(median_trialStartTime)
xline(median_time2Collect)
[numTrials, ~] = size(time2Collect);
Tris = [1:numTrials]';

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

custom_colormap = [
    1, 1, 1;         % white
    0.9, 0.95, 0.95;
    0.8, 0.9, 0.9;
    0.6, 0.85, 0.85;
    0.4, 0.8, 0.8;
    0.2, 0.8, 0.8;
    0.0, 0.8, 0.8;   % robin's egg blue
];

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

n = 256; 
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));


figure('Position', [100, 100, 250, 600]); 
t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');


ax1 = nexttile(t);
hold on;
xlim([-8 8]);
set(gca, 'XTick', []);
shadedErrorBar(ts1, nanmean(zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num}), nanmean(neuron_sem_array{1, 1}(prechoice_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num}), nanmean(neuron_sem_array{1, 1}(postchoice_reward_block_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});

xline(0);
ylim([-1 1]);
hold off

ax2 = nexttile(t);
hold on;

imagesc(ts1, 1:size(zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num}, 1), zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num});

colormap(custom_colormap);

clim([-1 1]);



ylim([0.5, size(zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num}, 1) + 0.5]);

xlim([-8 8]);

set(gca, 'XTick', []);
set(gca, 'YTick', [1, size(zall_mouse{select_mouse_index, array_to_plot(1)}{1, plot_num}, 1)]);
xline(0)
xline(median_trialStartTime, 'g')
xline(median_time2Collect, 'r')
scatter(time2Collect, Tris               , 'Marker', 'p')
scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

ax3 = nexttile(t);
hold on;

time2Collect = BehavData.collectionTime(trials_per_mouse{select_mouse_index, array_to_plot(2)}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot(2)});
trialStartTime = BehavData.stTime(trials_per_mouse{select_mouse_index, array_to_plot(2)}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot(2)});
median_trialStartTime = median(trialStartTime)
median_time2Collect = median(time2Collect)
xline(median_trialStartTime, 'g')
xline(median_time2Collect, 'r')
[numTrials, ~] = size(time2Collect);
Tris = [1:numTrials]';



hold on;

imagesc(ts1, 1:size(zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num}, 1), zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num});

colormap(custom_colormap);

clim([-1 1]);

c = colorbar(ax2, 'eastoutside');
set(c, 'YTick', clim); % 

ylim([0.5, size(zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num}, 1) + 0.5]);

xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', []);
set(gca, 'YTick', [1, size(zall_mouse{select_mouse_index, array_to_plot(2)}{1, plot_num}, 1)]);
xline(0)
scatter(time2Collect, Tris               , 'Marker', 'p')
scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

c = colorbar(ax3, 'eastoutside');
set(c, 'YTick', clim); % 



%% Fig. 4B
 mean_data_array = {neuron_mean_array{1, 1}(prechoice_block_1==1, :), neuron_mean_array{1, 8}(prechoice_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 1}(prechoice_block_1==1, :), neuron_sem_array{1, 8}(prechoice_block_1==1, :)};

 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 1], [0 0.6], 3);



 %
 mean_data_array = {neuron_mean_array{1, 2}(postchoice_reward_block_1==1, :), neuron_mean_array{1, 9}(postchoice_reward_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 2}(postchoice_reward_block_1==1, :), neuron_sem_array{1, 9}(postchoice_reward_block_1==1, :)};


 
 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 3], [-0.1 0.8], 3)

 %
 mean_data_array = {neuron_mean_array{1, 3}(collect_block_1==1, :), neuron_mean_array{1, 10}(collect_block_1==1, :)};
 sem_data_array = {neuron_sem_array{1, 3}(collect_block_1==1, :), neuron_sem_array{1, 10}(collect_block_1==1, :)};


 [comparison] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [0 4], [-0.2 0.8], 3);





%% Fig. 4C
all_conserved_sum = sum(conserved_sum)
all_lost_sum = sum(lost_sum)
all_remapped_sum = sum(remapped_sum)
remaining_neurons = neuron_num - (all_conserved_sum + all_lost_sum +all_remapped_sum);

figure;
piechart([all_conserved_sum/neuron_num, all_lost_sum/neuron_num, all_remapped_sum/neuron_num, remaining_neurons/neuron_num])

%% Fig. 4D


conserved_prechoice = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
conserved_prechoice_sum = sum(conserved_prechoice)

lost_prechoice = prechoice_block_1 == event_for_figures & prechoice_blocks_2_and_3 ~= event_for_figures;
lost_prechoice_sum = sum(lost_prechoice)

remapped_prechoice = prechoice_block_1 ~= event_for_figures & prechoice_blocks_2_and_3 == event_for_figures;
remapped_prechoice_sum = sum(remapped_prechoice)


mean_data_array = {neuron_mean_array{1, 8}(conserved_prechoice  ==1, :), neuron_mean_array{1, 8}(remapped_prechoice  ==1, :), neuron_mean_array{1, 8}(lost_prechoice  ==1, :)}
sem_data_array = {neuron_sem_array{1, 8}(conserved_prechoice  ==1, :), neuron_sem_array{1, 8}(remapped_prechoice  ==1, :), neuron_sem_array{1, 8}(lost_prechoice  ==1, :)}


[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 5], [-0.5 0.7], 3)


%% Fig. 4E - see PV code

%% Fig. 4F

% all possible combos
variables = {'conserved_ratio', 'lost_ratio', 'remapped_ratio', ...
             'conserved_prechoice_ratio', 'remapped_prechoice_ratio', 'lost_prechoice_ratio', ...
             'conserved_postchoice_ratio', 'remapped_postchoice_ratio', 'lost_postchoice_ratio', ...
             'conserved_collection_ratio', 'remapped_collection_ratio', 'lost_collection_ratio'};

correlation_results = table();


for i = 1:length(variables)
    x = eval(variables{i})';  
    y = risk_table.Mean_1_to_3;  
    [r, pval] = corrcoef(x, y);
    correlation_results = [correlation_results; table({variables{i}}, r(2), pval(2), 'VariableNames', {'Variable', 'Correlation', 'PValue'})];
end


%% Fig. 4G

 mean_data_array = {neuron_mean_array{1, 1}(remapped_prechoice==1, :), neuron_mean_array{1, 8}(remapped_prechoice==1, :)};
 sem_data_array = {neuron_sem_array{1, 1}(remapped_prechoice==1, :), neuron_sem_array{1, 8}(remapped_prechoice==1, :)};


 [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-5 1], [-0.05 0.4], 3);



%% Fig. 4H

remapped_postchoice_was_neutral = remapped_postchoice == 1 & (collect_block_1 == 0 & prechoice_block_1 == 0);
remapped_postchoice_was_prechoice = remapped_postchoice == 1 & prechoice_block_1 == 1;
remapped_postchoice_was_collect = remapped_postchoice == 1 & collect_block_1 == 1;

remapped_consumption_was_neutral = remapped_consumption == 1 & (postchoice_reward_block_1 == 0 & prechoice_block_1 == 0);
remapped_consumption_was_prechoice = remapped_consumption == 1 & prechoice_block_1 == 1;
remapped_consumption_was_postchoice = remapped_consumption == 1 & postchoice_reward_block_1 == 1;

remapped_prechoice_was_neutral = remapped_prechoice == 1 & (collect_block_1 == 0 & postchoice_reward_block_1 == 0);
remapped_prechoice_was_postchoice = remapped_prechoice == 1 & postchoice_reward_block_1 == 1;
remapped_prechoice_was_collect = remapped_prechoice == 1 & collect_block_1 == 1;


remapped_prechoice_retrospective_stack = [sum(remapped_prechoice_was_neutral),...
                                          0,... %for prechoice
                                          sum(remapped_prechoice_was_postchoice), ...
                                          sum(remapped_prechoice_was_collect)];

remapped_postchoice_retrospective_stack = [sum(remapped_postchoice_was_neutral), ...
                                           sum(remapped_postchoice_was_prechoice), ...
                                           0,... %for postchoice
                                           sum(remapped_postchoice_was_collect)];

remapped_consumption_retrospective_stack = [sum(remapped_consumption_was_neutral), ...
                                            sum(remapped_consumption_was_prechoice), ...
                                            sum(remapped_consumption_was_postchoice),...
                                            0]; %for collect


stacked_data = [remapped_prechoice_retrospective_stack; ...
                remapped_postchoice_retrospective_stack; ...
                remapped_consumption_retrospective_stack];


categories = {'New Prechoice', 'New Postchoice', 'New Consumption'};
subcategories = {'Neutral', 'Prechoice', 'Post-choice', 'Collect'};


figure;
bar(stacked_data, 'stacked');
colormap(parula); 
legend(subcategories, 'Location', 'northeastoutside');
xticks(1:length(categories));
xticklabels(categories);
ylabel('# of neurons');


%% For Fig. 4J-M - use behav_figs.m



