

%% Fig. 6F


behav_tbl_iter_single = behav_tbl_iter(1);


concatenatedTable = table();


for i = 1:numel(behav_tbl_iter_single) 
    twelveByOneCellArray = behav_tbl_iter_single{i};  
    tempTable = table();
    for j = 1:numel(twelveByOneCellArray)
        currentTable = twelveByOneCellArray{j};
        tempTable = vertcat(tempTable, currentTable);
    end
    concatenatedTable = vertcat(concatenatedTable, tempTable);
end

median_start_time_from_choice = median(concatenatedTable.stTime - concatenatedTable.choiceTime);
median_collect_time_from_choice = median(concatenatedTable.collectionTime - concatenatedTable.choiceTime);

median_start_time_from_choice_large = median(concatenatedTable.stTime(concatenatedTable.bigSmall == 1.2) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 1.2));
median_start_time_from_choice_small = median(concatenatedTable.stTime(concatenatedTable.bigSmall == 0.3) - concatenatedTable.choiceTime(concatenatedTable.bigSmall == 0.3));

figure;
hold on


width = 350; 
height = 300; 
set(gcf, 'Position', [50, 25, width, height]); 
xlim([-8 8]);

set(gca, 'XTick', [-8, 0, 8], 'YTick', [-0.6, 0, 0.6]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(prechoice_block_1==1, :)), std(neuron_mean_array{1, 1}(prechoice_block_1==1, :)/sqrt(size(neuron_mean_array{1, 1}(prechoice_block_1==1, :), 1))), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(postchoice_reward_block_1==1, :)), std(neuron_mean_array{1, 1}(postchoice_reward_block_1==1, :)/sqrt(size(neuron_mean_array{1, 1}(postchoice_reward_block_1==1, :), 1))), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(collect_block_1==1, :)), std(neuron_mean_array{1, 1}(collect_block_1==1, :)/sqrt(size(neuron_mean_array{1, 1}(collect_block_1==1, :), 1))), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
ylim([-0.6 0.6]);

hold off

%% Fig. 6G

prechoice_b1_percent_bla = sum(prechoice_block_1_bla)/size(prechoice_block_1_bla, 2)*100;
prechoice_b1_percent_bla_nacs = sum(prechoice_block_1_bla_nacs)/size(prechoice_block_1_bla_nacs, 2)*100;

non_prechoice_b1_percent_bla = 100 - prechoice_b1_percent_bla;

non_prechoice_b1_percent_bla_nacs = 100 - prechoice_b1_percent_bla_nacs;


outer_pie = [non_prechoice_b1_percent_bla prechoice_b1_percent_bla];
inner_pie = [non_prechoice_b1_percent_bla_nacs prechoice_b1_percent_bla_nacs];

C = {outer_pie,... 
     inner_pie};


nested_pie(C, 'LegendOrder', [2 1],...
    'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);


n1 = sum(prechoice_block_1_bla); 
n2 = sum(prechoice_block_1_bla_nacs);    
x1 = size(prechoice_block_1_bla, 2); 
x2 = size(prechoice_block_1_bla_nacs, 2); 


p1 = x1 / n1;
p2 = x2 / n2;

% proportion
p_combined = (x1 + x2) / (n1 + n2);

% Z
z = (p1 - p2) / sqrt(p_combined * (1 - p_combined) * (1/n1 + 1/n2));

z_val = abs(z)


p_value = 2 * (1 - normcdf(abs(z)));

%  results
fprintf('Z-statistic: %.4f\n', z);
fprintf('P-value: %.4f\n', p_value);

%
postchoice_b1_percent_bla = sum(postchoice_reward_block_1_bla)/size(postchoice_reward_block_1_bla, 2)*100;
postchoice_b1_percent_bla_nacs = sum(postchoice_reward_block_1_bla_nacs)/size(postchoice_reward_block_1_bla_nacs, 2)*100;

non_postchoice_reward_block_1_bla = 100 - postchoice_b1_percent_bla;

non_postchoice_reward_block_1_bla_nacs = 100 - postchoice_b1_percent_bla_nacs;


outer_pie = [ non_postchoice_reward_block_1_bla postchoice_b1_percent_bla];
inner_pie = [ non_postchoice_reward_block_1_bla_nacs postchoice_b1_percent_bla_nacs];

C = {...
    outer_pie,... 
     inner_pie};


nested_pie(C, 'LegendOrder', [2, 1],...
    'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);


title('Nested Pie Chart');

% Variables
n1 = sum(postchoice_reward_block_1_bla); 
n2 = sum(postchoice_reward_block_1_bla_nacs);     
x1 = size(postchoice_reward_block_1_bla, 2); 
x2 = size(postchoice_reward_block_1_bla_nacs, 2);     


p1 = x1 / n1;
p2 = x2 / n2;

%  proportion
p_combined = (x1 + x2) / (n1 + n2);

% Z
z = (p1 - p2) / sqrt(p_combined * (1 - p_combined) * (1/n1 + 1/n2));

z_val = abs(z)


p_value = 2 * (1 - normcdf(abs(z)));

%  results
fprintf('Z-statistic: %.4f\n', z);
fprintf('P-value: %.4f\n', p_value);

%
collect_b1_percent_bla = sum(collect_block_1_bla)/size(collect_block_1_bla, 2)*100;
collect_b1_percent_bla_nacs = sum(collect_block_1_bla_nacs)/size(collect_block_1_bla_nacs, 2)*100;

non_collect_b1_percent_bla = 100 - collect_b1_percent_bla;

non_collect_b1_percent_bla_nacs = 100 - collect_b1_percent_bla_nacs;


outer_pie = [ non_collect_b1_percent_bla collect_b1_percent_bla];
inner_pie = [ non_collect_b1_percent_bla_nacs collect_b1_percent_bla_nacs];

C = {...
    outer_pie,... 
     inner_pie};

nested_pie(C, 'LegendOrder', [2, 1],...
    'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);

n1 = sum(collect_block_1_bla); 
n2 = sum(collect_block_1_bla_nacs);   
x1 = size(collect_block_1_bla, 2); 
x2 = size(collect_block_1_bla_nacs, 2);     


p1 = x1 / n1;
p2 = x2 / n2;

%  proportion
p_combined = (x1 + x2) / (n1 + n2);

% Z
z = (p1 - p2) / sqrt(p_combined * (1 - p_combined) * (1/n1 + 1/n2));

z_val = abs(z)

p_value = 2 * (1 - normcdf(abs(z)));

% results
fprintf('Z-statistic: %.4f\n', z);
fprintf('P-value: %.4f\n', p_value);


%% Fig. 6H

non_shk_only_BLA = 100 - percent_shk_only_BLA;

non_shk_only_bla_nacs = 100 - percent_shk_only_bla_nacs;



outer_pie = [percent_shk_only_BLA non_shk_only_BLA];
inner_pie = [percent_shk_only_bla_nacs non_shk_only_bla_nacs];

C = {...
    inner_pie,... % Inner to outer layer
    outer_pie};


nested_pie(C,'LegendOrder', [2 1], 'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);


total_BLA_neurons = only_shk_sum_BLA/(percent_shk_only_BLA/100);
total_bla_nacs_neurons = only_shk_sum_bla_nacs/(percent_shk_only_bla_nacs/100);


n1 = only_shk_sum_BLA; 
n2 = only_shk_sum_bla_nacs;     
x1 = total_BLA_neurons; 
x2 =total_bla_nacs_neurons;    

p1 = x1 / n1;
p2 = x2 / n2;

% proportion
p_combined = (x1 + x2) / (n1 + n2);

% Z
z = (p1 - p2) / sqrt(p_combined * (1 - p_combined) * (1/n1 + 1/n2));

z_val = abs(z)


p_value = 2 * (1 - normcdf(abs(z)));

% results
fprintf('Z-statistic: %.4f\n', z);
fprintf('P-value: %.4f\n', p_value);




%
non_consum_only_BLA = 100 - percent_consum_only_BLA;

non_consum_only_bla_nacs = 100 - percent_consum_only_bla_nacs;


outer_pie = [percent_consum_only_BLA non_consum_only_BLA];
inner_pie = [percent_consum_only_bla_nacs non_consum_only_bla_nacs];

C = {...
    inner_pie,... % Inner to outer layer
    outer_pie};


nested_pie(C,'LegendOrder', [2 1], 'PercentStatus', {'on', 'on'}, 'RhoLower', 0.4);

total_BLA_neurons = only_consumption_sum_BLA/(percent_consum_only_BLA/100);
total_bla_nacs_neurons = only_consumption_sum_bla_nacs/(percent_consum_only_bla_nacs/100);


n1 = only_consumption_sum_BLA; 
n2 = only_consumption_sum_bla_nacs;     
x1 = total_BLA_neurons; 
x2 =total_bla_nacs_neurons;     


p1 = x1 / n1;
p2 = x2 / n2;

% proportion
p_combined = (x1 + x2) / (n1 + n2);

% Z
z = (p1 - p2) / sqrt(p_combined * (1 - p_combined) * (1/n1 + 1/n2));

z_val = abs(z)


p_value = 2 * (1 - normcdf(abs(z)));

% results
fprintf('Z-statistic: %.4f\n', z);
fprintf('P-value: %.4f\n', p_value);



%% Fig. 6I
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


 %% For Fig. 6K-N - use behav_figs.m