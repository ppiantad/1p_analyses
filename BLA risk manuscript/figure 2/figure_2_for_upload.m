

%% For Figure 2C

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

custom_colormap = [
    1, 1, 1;         % white
    0.9, 0.9, 0.95;
    0.8, 0.8, 0.9;
    0.6, 0.6, 0.8;
    0.4, 0.4, 0.7;
    0.2, 0.2, 0.6;
    0.0, 0.0, 0.55;   % dark blue
];

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



% this is done for each ensemble

% for RDT D1 BLA_Insc_40:
%prechoice neuron num 12
%postchoice rew num 70
%consumption num 10
%shock num 11

pre_choice_window = [-4 0];     
post_choice_window = [0 2];   
consumption_window = [2 5];    

windows = {pre_choice_window, post_choice_window, consumption_window};
plot_num = [15, 70, 10];

array_to_plot = [1, 1, 1]; % depends on the structure of zall

select_mouse = 'BLA_Insc_40';

% for RDT D1 BLA_Insc_25:
%prechoice neuron num 46
%postchoice rew num 38
%consumption num 39
%shock num 11


% for Pre RDT RM BLA_Insc_40:
%prechoice neuron num 12
%postchoice rew num 70
%consumption num 10
%shock num 11


select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'Pre_RDT_RM';

second_session = 'RDT_D1';

num_neurons = length(plot_num);

figure;

for neuron_idx = 1:num_neurons
    current_neuron = plot_num(neuron_idx);
    current_window = windows{neuron_idx};
    array_to_plot_current = array_to_plot(neuron_idx);
    neuron_data = zall_mouse{select_mouse_index, array_to_plot_current}{1, current_neuron};
    neuron_data_s_array = caTraceTrials_spikes_mouse{select_mouse_index, array_to_plot_current}{1, current_neuron};
    time_indices = ts1 >= current_window(1) & ts1 <= current_window(2);
    restricted_data = neuron_data(:, time_indices);
    mean_activity = mean(restricted_data, 2);
    
    [~, sorted_trial_indices] = sort(mean_activity, 'descend');
    selected_trials = sorted_trial_indices(1:min(5, size(neuron_data, 1)));
    
    y_limits = [min(neuron_data(:)), max(neuron_data(:))];
    
    for subplot_idx = 1:length(selected_trials)
        subplot_idx_global = (neuron_idx - 1) * 5 + subplot_idx;
        subplot(num_neurons, 5, subplot_idx_global);
        hold on;
        
        trial = selected_trials(subplot_idx);
        
        plot(ts1, neuron_data(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
        hold on; plot(ts1, neuron_data_s_array(trial, :), 'Color', [custom_colormap(end, :), 0.5]);

        ylim(y_limits);
        xlim([-8 8]);
        set(gca, 'XTick', [-8, 0, 8]);
        

        xline(0, 'k--');
        yline(0, 'k--');
        set(gca, 'FontSize', 12);
        
        hold off;
    end
    
end


%% Fig. 2E


pre_choice_neurons = neuron_mean_array{1, 1}(prechoice_block_1, :);
post_choice_reward_neurons = neuron_mean_array{1, 1}(postchoice_reward_block_1, :);
consumption_neurons = neuron_mean_array{1, 1}(collect_block_1, :);

only_active_array_stacked = [pre_choice_neurons; post_choice_reward_neurons; consumption_neurons];

% sort by peak time
[peak_values, time_of_peak_activity] = max(pre_choice_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
pre_choice_neurons_sorted = pre_choice_neurons(sort_indices, :);


[peak_values, time_of_peak_activity] = max(post_choice_reward_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
post_choice_reward_neurons_sorted = post_choice_reward_neurons(sort_indices, :);


[peak_values, time_of_peak_activity] = max(consumption_neurons, [], 2);
[~, sort_indices] = sort(time_of_peak_activity);
consumption_neurons_sorted = consumption_neurons(sort_indices, :);

sorted_only_active_array_stacked = [pre_choice_neurons_sorted; post_choice_reward_neurons_sorted; consumption_neurons_sorted];



figure;
imagesc(ts1, 1, sorted_only_active_array_stacked);
colorbar;
xlabel('Time (s)');
ylabel('Neuron');
set(gca, 'YDir', 'reverse');
clim([-1 1])
xline(0);
colormap(gray);


caxis([-1 1]); 
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 


%% Fig. 2F - separate decoding scripts


%% For Figure 2G 
% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav

figure;
pre_choice_active_ind = find(respClass_all_array{1,1} == 1);
consum_active_ind = find(respClass_all_array{1,3} == 1);
post_choice_active_ind = find(respClass_all_array{1,2} == 1);
setListData = {pre_choice_active_ind, consum_active_ind, post_choice_active_ind};
setLabels = ["Pre-choice excited", "Consumption excited", "Post-choice excited"];

h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

h.ShowIntersectionCounts = true;
h.ShowIntersectionAreas = true;
% h.SetLabels = [];


%% Fig. 2H


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

% For Figure 2H (top)
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


% For Figure 2H (bottom)
pre_choice_window = [-4 0];  
post_choice_window = [0 2];   
consumption_window = [1 3];     


pre_choice_neuron_count = 0;

for i = 1:size(neuron_mean_array{1,1}, 1)
    if prechoice_block_1(i) == 1
        pre_choice_neuron_count = pre_choice_neuron_count+1;
        selected_data = neuron_mean_array{1, 1}(i, :);

        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        action_auc_pre_choice(pre_choice_neuron_count) = trapz(selected_data(pre_choice_indices));
        action_auc_post_choice(pre_choice_neuron_count) = trapz(selected_data(post_choice_indices));
        action_auc_consumption(pre_choice_neuron_count) = trapz(selected_data(consumption_indices));
    else

    end
end

post_choice_neuron_count = 0;
for i = 1:size(neuron_mean_array{1,1}, 1)
    if postchoice_reward_block_1(i) == 1
        post_choice_neuron_count = post_choice_neuron_count+1;
        selected_data = neuron_mean_array{1, 1}(i, :);

        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        post_choice_reward_auc_pre_choice(post_choice_neuron_count) = trapz(selected_data(pre_choice_indices));
        post_choice_reward_auc_post_choice(post_choice_neuron_count) = trapz(selected_data(post_choice_indices));
        post_choice_reward_auc_consumption(post_choice_neuron_count) = trapz(selected_data(consumption_indices));
    else

    end
end


consumption_neuron_count = 0;
for i = 1:size(neuron_mean_array{1,3}, 1)
    if collect_block_1(i) == 1
        consumption_neuron_count = consumption_neuron_count+1;
        selected_data = neuron_mean_array{1, 3}(i, :);

        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        consumption_auc_pre_choice(consumption_neuron_count) = trapz(selected_data(pre_choice_indices));
        consumption_auc_post_choice(consumption_neuron_count) = trapz(selected_data(post_choice_indices));
        consumption_auc_consumption(consumption_neuron_count) = trapz(selected_data(consumption_indices));
    else

    end
end

mean_pre_choice = [mean(action_auc_pre_choice(:)), mean(post_choice_reward_auc_pre_choice(:)), mean(consumption_auc_pre_choice(:))];
sem_pre_choice = [std(action_auc_pre_choice(:))/sqrt(numel(action_auc_pre_choice)), std(post_choice_reward_auc_pre_choice(:))/sqrt(numel(post_choice_reward_auc_pre_choice)), std(consumption_auc_pre_choice(:))/sqrt(numel(consumption_auc_pre_choice))];

mean_post_choice = [mean(action_auc_post_choice(:)), mean(post_choice_reward_auc_post_choice(:)), mean(consumption_auc_post_choice(:))];
sem_post_choice = [std(action_auc_post_choice(:))/sqrt(numel(action_auc_post_choice)), std(post_choice_reward_auc_post_choice(:))/sqrt(numel(post_choice_reward_auc_post_choice)), std(consumption_auc_post_choice(:))/sqrt(numel(consumption_auc_post_choice))];

mean_consumption = [mean(action_auc_consumption(:)), mean(post_choice_reward_auc_consumption(:)), mean(consumption_auc_consumption(:))];
sem_consumption = [std(action_auc_consumption(:))/sqrt(numel(action_auc_consumption)), std(post_choice_reward_auc_consumption(:))/sqrt(numel(post_choice_reward_auc_consumption)), std(consumption_auc_consumption(:))/sqrt(numel(consumption_auc_consumption))];


figure;
width = 350; 
height = 200; 
set(gcf, 'Position', [50, 25, width, height]); 
set(gca, 'YTick', [-10, 0, 10]);
bar_groups = 1:3; 
bar_width = 0.3; 
hold on;
bar(bar_groups - bar_width, mean_pre_choice, bar_width, 'b');
bar(bar_groups, mean_post_choice, bar_width, 'g');
bar(bar_groups + bar_width, mean_consumption, bar_width, 'r');

errorbar(bar_groups - bar_width, mean_pre_choice, sem_pre_choice, 'k', 'LineStyle', 'none');
errorbar(bar_groups, mean_post_choice, sem_post_choice, 'k', 'LineStyle', 'none');
errorbar(bar_groups + bar_width, mean_consumption, sem_consumption, 'k', 'LineStyle', 'none');

ylabel('Mean AUC');
legend('Pre-choice', 'Post-choice', 'Consumption');
xticks(bar_groups);
grid off;
hold off;

action_data = [action_auc_pre_choice(:); action_auc_post_choice(:); action_auc_consumption(:)];

groups = [ones(size(action_auc_consumption(:)));  
          2*ones(size(action_auc_post_choice(:))); 
          3*ones(size(action_auc_pre_choice(:)))]; 


postchoice_data = [post_choice_reward_auc_pre_choice(:); post_choice_reward_auc_post_choice(:); post_choice_reward_auc_consumption(:)];

groups = [ones(size(post_choice_reward_auc_consumption(:)));  
          2*ones(size(post_choice_reward_auc_post_choice(:))); 
          3*ones(size(post_choice_reward_auc_pre_choice(:)))]; 




%% For Figure 2I (left)
% create correlation matrix heatmap with exclusively active cells
test = [];
test = [neuron_mean_array{1, 1}(prechoice_block_1 == 1, :)];
test = [test; neuron_mean_array{1, 1}(postchoice_reward_block_1 == 1, :)];
test = [test; neuron_mean_array{1, 1}(collect_block_1 == 1,:)];

pre_choice_index = [1:sum(prechoice_block_1)];
post_choice_index = [pre_choice_index(end)+1:pre_choice_index(end)+sum(postchoice_reward_block_1)];
consumption_index = [post_choice_index(end)+1:post_choice_index(end)+sum(collect_block_1)];
neutral_index = [consumption_index(end)+1:consumption_index(end)+sum(respClass_all_array{1, 2}~=1 & respClass_all_array{1, 1}~=1 & respClass_all_array{1,3}~=1)];


data = test;

alpha = 0.0001;


correlation_matrix = zeros(size(data, 1));
p_value_matrix = zeros(size(data, 1));


for i = 1:size(data, 1)
    for j = 1:size(data, 1)
        [corr_coeff, p_value] = corrcoef(data(i, :)', data(j, :)');
        correlation_matrix(i, j) = corr_coeff(1, 2); % Store correlation coefficient
        p_value_matrix(i, j) = p_value(1, 2); % Store p-value
    end
end


figure;
imagesc(correlation_matrix);

axis square; 

ylim([1  size(test, 1)])

xticks([1  size(test, 1)]);
yticks([1  size(test, 1)]);

colormap(gray);


caxis([-1 1]); % assuming correlations range from -1 to 1

c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

%% For Figure 2I (right)
action_p_value_matrix = p_value_matrix(pre_choice_index, pre_choice_index);
action_correl_matrix = correlation_matrix(pre_choice_index, pre_choice_index);

n = size(pre_choice_index, 2); 
k = 2;   



action_positive_count = 0;
action_negative_count = 0;
action_no_correlation_count = 0;


matrix_size = size(action_correl_matrix, 1);
uu = 1;

for i = 1:matrix_size
    for j = i+1:matrix_size 

        action_ensemble_corr_overall(uu) = action_correl_matrix(i, j);
        uu = uu+1;
        % 
        if action_p_value_matrix(i, j) < alpha
            if action_correl_matrix(i, j) > 0
                action_positive_count = action_positive_count + 1;
            elseif action_correl_matrix(i, j) < 0
                action_negative_count = action_negative_count + 1;
            end
        else
            action_no_correlation_count = action_no_correlation_count + 1;
        end
    end
end

action_comparisons_possible = [action_positive_count + action_negative_count + action_no_correlation_count];


action_data = [(action_positive_count/action_comparisons_possible)*100, (action_negative_count/action_comparisons_possible)*100, (action_no_correlation_count/action_comparisons_possible)*100];



%
post_choice_p_value_matrix = p_value_matrix(post_choice_index, post_choice_index);
post_choice_correl_matrix = correlation_matrix(post_choice_index, post_choice_index);

n = size(post_choice_index, 2); 
k = 2;   



post_choice_positive_count = 0;
post_choice_negative_count = 0;
post_choice_no_correlation_count = 0;


matrix_size = size(post_choice_correl_matrix, 1);
uu = 1

for i = 1:matrix_size
    for j = i+1:matrix_size % Start from i+1 to exclude the diagonal
        post_choice_ensemble_corr_overall(uu) = post_choice_correl_matrix(i, j);
        uu = uu+1;
        % Check if p-value is less than 0.01
        if post_choice_p_value_matrix(i, j) < alpha
            if post_choice_correl_matrix(i, j) > 0
                post_choice_positive_count = post_choice_positive_count + 1;
            elseif post_choice_correl_matrix(i, j) < 0
                post_choice_negative_count = post_choice_negative_count + 1;
            end
        else
            post_choice_no_correlation_count = post_choice_no_correlation_count + 1;
        end
    end
end

post_choice_comparisons_possible = [post_choice_positive_count + post_choice_negative_count + post_choice_no_correlation_count];

post_choice_data = [(post_choice_positive_count/post_choice_comparisons_possible)*100, (post_choice_negative_count/post_choice_comparisons_possible)*100, (post_choice_no_correlation_count/post_choice_comparisons_possible)*100];

consumption_p_value_matrix = p_value_matrix(consumption_index, consumption_index);
consumption_correl_matrix = correlation_matrix(consumption_index, consumption_index);

n = size(consumption_index, 2); 
k = 2;  


consumption_positive_count = 0;
consumption_negative_count = 0;
consumption_no_correlation_count = 0;


matrix_size = size(consumption_correl_matrix, 1);
uu = 1

for i = 1:matrix_size
    for j = i+1:matrix_size 
        consumption_ensemble_corr_overall(uu) = consumption_correl_matrix(i, j);
        uu = uu+1;
        % 
        if consumption_p_value_matrix(i, j) < alpha
            if consumption_correl_matrix(i, j) > 0
                consumption_positive_count = consumption_positive_count + 1;
            elseif consumption_correl_matrix(i, j) < 0
                consumption_negative_count = consumption_negative_count + 1;
            end
        else
            consumption_no_correlation_count = consumption_no_correlation_count + 1;
        end
    end
end

consumption_comparisons_possible = [consumption_positive_count + consumption_negative_count + consumption_no_correlation_count];



consumption_data = [(consumption_positive_count/consumption_comparisons_possible)*100, (consumption_negative_count/consumption_comparisons_possible)*100, (consumption_no_correlation_count/consumption_comparisons_possible)*100];

action_post_choice_p_value_matrix = p_value_matrix(pre_choice_index, post_choice_index);
action_post_choice_correl_matrix = correlation_matrix(pre_choice_index, post_choice_index);

n1 = size(action_post_choice_p_value_matrix, 1); 
n2 = size(action_post_choice_p_value_matrix, 2); 
k = 2;    

action_post_choice_positive_count = 0;
action_post_choice_negative_count = 0;
action_post_choice_no_correlation_count = 0;

[num_neurons_1, num_neurons_2] = size(action_post_choice_correl_matrix);


for i = 1:num_neurons_1
    for j = 1:num_neurons_2
        correlation = action_post_choice_correl_matrix(i, j);
        p_value = action_post_choice_p_value_matrix(i, j);
        if p_value < alpha
            if correlation > 0
                action_post_choice_positive_count = action_post_choice_positive_count + 1;
            elseif correlation < 0
                action_post_choice_negative_count = action_post_choice_negative_count + 1;
            end
        else
            action_post_choice_no_correlation_count = action_post_choice_no_correlation_count + 1;
        end
    end
end

action_post_choice_comparisons_possible = [action_post_choice_positive_count+action_post_choice_negative_count+action_post_choice_no_correlation_count];



action_post_choice_data = [(action_post_choice_positive_count/action_post_choice_comparisons_possible)*100, (action_post_choice_negative_count/action_post_choice_comparisons_possible)*100, (action_post_choice_no_correlation_count/action_post_choice_comparisons_possible)*100];

action_consumption_p_value_matrix = p_value_matrix(pre_choice_index, consumption_index);
action_consumption_correl_matrix = correlation_matrix(pre_choice_index, consumption_index);

n1 = size(action_consumption_p_value_matrix, 1); 
n2 = size(action_consumption_p_value_matrix, 2);
k = 2;   


action_consumption_positive_count = 0;
action_consumption_negative_count = 0;
action_consumption_no_correlation_count = 0;


[num_neurons_1, num_neurons_2] = size(action_consumption_correl_matrix);

for i = 1:num_neurons_1
    for j = 1:num_neurons_2
        correlation = action_consumption_correl_matrix(i, j);
        p_value = action_consumption_p_value_matrix(i, j);
        if p_value < alpha
            if correlation > 0
                action_consumption_positive_count = action_consumption_positive_count + 1;
            elseif correlation < 0
                action_consumption_negative_count = action_consumption_negative_count + 1;
            end
        else
            action_consumption_no_correlation_count = action_consumption_no_correlation_count + 1;
        end
    end
end


action_consumption_comparisons_possible = [action_consumption_positive_count+action_consumption_negative_count+action_consumption_no_correlation_count];


action_consumption_data = [(action_consumption_positive_count/action_consumption_comparisons_possible)*100, (action_consumption_negative_count/action_consumption_comparisons_possible)*100, (action_consumption_no_correlation_count/action_consumption_comparisons_possible)*100];



post_choice_consumption_p_value_matrix = p_value_matrix(post_choice_index, consumption_index);
post_choice_consumption_correl_matrix = correlation_matrix(post_choice_index, consumption_index);

n1 = size(post_choice_consumption_p_value_matrix, 1);
n2 = size(post_choice_consumption_p_value_matrix, 2); 
k = 2;    



post_choice_consumption_positive_count = 0;
post_choice_consumption_negative_count = 0;
post_choice_consumption_no_correlation_count = 0;


[num_neurons_1, num_neurons_2] = size(post_choice_consumption_correl_matrix);

for i = 1:num_neurons_1
    for j = 1:num_neurons_2
        correlation = post_choice_consumption_correl_matrix(i, j);
        p_value = post_choice_consumption_p_value_matrix(i, j);
        if p_value < alpha
            if correlation > 0
                post_choice_consumption_positive_count = post_choice_consumption_positive_count + 1;
            elseif correlation < 0
                post_choice_consumption_negative_count = post_choice_consumption_negative_count + 1;
            end
        else
            post_choice_consumption_no_correlation_count = post_choice_consumption_no_correlation_count + 1;
        end
    end
end


post_choice_consumption_comparisons_possible = [post_choice_consumption_positive_count+post_choice_consumption_negative_count+post_choice_consumption_no_correlation_count];



post_choice_consumption_data = [(post_choice_consumption_positive_count/post_choice_consumption_comparisons_possible)*100, (post_choice_consumption_negative_count/post_choice_consumption_comparisons_possible)*100, (post_choice_consumption_no_correlation_count/post_choice_consumption_comparisons_possible)*100];

x = 1:6;
figure;
barh(x, [action_data; post_choice_data; consumption_data; action_post_choice_data; action_consumption_data; post_choice_consumption_data], 'stacked');
legend('Positive', 'Negative', 'No sig correl');



%% Figure 2J, Figure 2K, and Figure 2L

array_for_means = 3; 


for q = 1:length (behav_tbl_iter{array_for_means, 1})
    nestedCellArray_1 = behav_tbl_iter{array_for_means, 1}{q};
    if ~isempty(nestedCellArray_1)

        for zz = 1:size(nestedCellArray_1, 1)
            valid_start_times = nestedCellArray_1.stTime(2:end);
            valid_choice_times = nestedCellArray_1.choiceTime(1:end-1);
            delay_to_initiation = valid_start_times - valid_choice_times;
            trial_types = nestedCellArray_1.bigSmall;
        end

        trial_choice_times = nestedCellArray_1.choiceTime - nestedCellArray_1.stTime;
        % delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
        delay_to_collect_post_shk = nestedCellArray_1.collectionTime - nestedCellArray_1.choiceTime;
        
        trial_choice_times_by_mouse{q} = trial_choice_times;
        delay_to_initiation_by_mouse{q} = delay_to_initiation;
        delay_to_collect_post_shk_by_mouse{q} = delay_to_collect_post_shk;
        trial_types_by_mouse{q} = trial_types;
        consum_times = nestedCellArray_1.collectionTime_end - nestedCellArray_1.collectionTime;
        consum_times_by_mouse{q} = consum_times;
        clear trial_choice_times delay_to_initiation delay_to_collect_post_shk trial_types consum_times
    end


end

trial_types_concat = cat(1, trial_types_by_mouse{:});
trial_choice_times_concat = cat(1, trial_choice_times_by_mouse{:});
rew_collect_times_concat = cat(1, delay_to_collect_post_shk_by_mouse{:});
consum_times_concat = cat(1, consum_times_by_mouse{:});

bar_separation_value = 3;

figure;


width = 400; 
height = 500; 
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]


subplot(1, 3, 1); 
hold on;
colors = repmat([0.5, 0.5, 0.5], length(trial_choice_times_concat), 1); 
colors(trial_types_concat == 1.2, :) = repmat([0, 0, 1], sum(trial_types_concat == 1.2), 1); 
colors(trial_types_concat == 0.3, :) = repmat([1, 0, 0], sum(trial_types_concat == 0.3), 1); 


s = swarmchart(ones(1, length(trial_choice_times_concat)), trial_choice_times_concat, 36, colors, 'filled');
s.MarkerEdgeColor = 'k'; 


mean_12 = mean(trial_choice_times_concat(trial_types_concat == 1.2));
mean_03 = mean(trial_choice_times_concat(trial_types_concat == 0.3));
plot([0.8, 1.2], [mean_12, mean_12], 'b-', 'LineWidth', 2);
plot([0.8, 1.2], [mean_03, mean_03], 'r-', 'LineWidth', 2); 


[h_choice_times p_choice_times, ~, stats_choice_times] = ttest2(trial_choice_times_concat(trial_types_concat == 1.2), trial_choice_times_concat(trial_types_concat == 0.3))

% tests for normality for reviewers
[h_choice_times_large_normality, p_choicetimes_large_normality] = kstest(trial_choice_times_concat(trial_types_concat == 1.2))
[h,p,adstat,cv] = adtest(trial_choice_times_concat(trial_types_concat == 1.2))
[h_choice_times_small_normality, p_choicetimes_small_normality] = kstest(trial_choice_times_concat(trial_types_concat == 0.3))
[h,p,adstat,cv] = adtest(trial_choice_times_concat(trial_types_concat == 0.3))

normalitytest(trial_choice_times_concat(trial_types_concat == 1.2)')
normalitytest(trial_choice_times_concat(trial_types_concat == 0.3)')

[p_rank_sum_choice_times h_rank_sum_choice_times stats_rank_sum_choice_times] = ranksum(trial_choice_times_concat(trial_types_concat == 1.2), trial_choice_times_concat(trial_types_concat == 0.3))

U_choice_times = (stats_rank_sum_choice_times.ranksum  ) - [size(trial_choice_times_concat(trial_types_concat == 1.2), 1)*(size(rew_collect_times_concat(trial_types_concat == 1.2), 1)+1)/2];


normalitytest(rew_collect_times_concat(trial_types_concat == 1.2)')
normalitytest(rew_collect_times_concat(trial_types_concat == 0.3)')

[p_rank_sum_collect_times h_rank_sum_collect_times stats_rank_sum_collect_times] = ranksum(rew_collect_times_concat(trial_types_concat == 1.2), rew_collect_times_concat(trial_types_concat == 0.3))

U_collect_times = (stats_rank_sum_collect_times.ranksum  ) - [size(rew_collect_times_concat(trial_types_concat == 1.2), 1)*(size(rew_collect_times_concat(trial_types_concat == 1.2), 1)+1)/2];



hold off;
xlabel('Choice latency');
ylabel('Choice latency (s)');
yline(0); 
xtickformat('%.1f');
ytickformat('%.1f');


subplot(1, 3, 2); 
hold on;

colors = repmat([0.5, 0.5, 0.5], length(rew_collect_times_concat), 1); %
colors(trial_types_concat == 1.2, :) = repmat([0, 0, 1], sum(trial_types_concat == 1.2), 1); 
colors(trial_types_concat == 0.3, :) = repmat([1, 0, 0], sum(trial_types_concat == 0.3), 1); 


s = swarmchart(ones(1, length(rew_collect_times_concat)) * bar_separation_value, rew_collect_times_concat, 36, colors, 'filled');
s.MarkerEdgeColor = 'k'; 


mean_12 = mean(rew_collect_times_concat(trial_types_concat == 1.2));
mean_03 = mean(rew_collect_times_concat(trial_types_concat == 0.3));
plot([bar_separation_value - 0.2, bar_separation_value + 0.2], [mean_12, mean_12], 'b-', 'LineWidth', 2); 
plot([bar_separation_value - 0.2, bar_separation_value + 0.2], [mean_03, mean_03], 'r-', 'LineWidth', 2); 

[h_rew_times p_rew_times, ~, stats_rew_times] = ttest2(rew_collect_times_concat(trial_types_concat == 1.2), rew_collect_times_concat(trial_types_concat == 0.3))

hold off;
xlabel('Reward collection');
ylabel('Collection latency (s)');
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');


variable_to_correlate = consum_times_by_mouse;


%



meanZallMouse = cell(size(zall_mouse, 2), 1);

% Define the time range
% timeRange = (ts1 >= -4) & (ts1 <= 0);
% timeRange = (ts1 >= 0) & (ts1 <= 2);
timeRange = (ts1 >= 1) & (ts1 <= 3);



for i = 1:length(zall_mouse)

    nestedCellArray_1 = zall_mouse{i, array_for_means};
    nestedCellArray_2 = zall_mouse{i, 1};

    meanNestedCellArray = cell(size(nestedCellArray_1));
    
    for j = 1:length(nestedCellArray_1)
     
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
        
        meanNestedCellArray{j} = meanValues;
    end
    
    meanZallMouse{i} = meanNestedCellArray;
end





%


correlationResults = cell(size(meanZallMouse));
correlationResults_sig = cell(size(meanZallMouse));



for i = 1:length(meanZallMouse)

    meanNestedCellArray = meanZallMouse{i};
    

    correlationNestedArray = zeros(size(meanNestedCellArray));
    corr_sig_NestedArray = zeros(size(meanNestedCellArray));
    trialIndex = mod(i-1, length(variable_to_correlate)) + 1;
    

    trialChoiceTimes = variable_to_correlate{i}';
        

    for j = 1:length(meanNestedCellArray)

        meanValues = meanNestedCellArray{j};
        

        



        if length(trialChoiceTimes) == length(meanValues)

            [correlationCoeff, corr_sig_vals] = corr(meanValues, trialChoiceTimes(:));
        elseif length(trialChoiceTimes) < length(meanValues)
            [correlationCoeff, corr_sig_vals] = corr(meanValues(1:end-1), trialChoiceTimes(:));
        else
          
            [correlationCoeff, corr_sig_vals] = NaN;
        end
        
        
        correlationNestedArray(j) = correlationCoeff;
        corr_sig_NestedArray(j) = corr_sig_vals;
    end
    clear meanValues

    correlationResults{i} = correlationNestedArray;
    correlationResults_sig{i} = corr_sig_NestedArray;
end



allCorrelations = [];


for i = 1:length(correlationResults)

    correlationNestedArray = correlationResults{i};
    for j = 1:length(correlationNestedArray)
        correlationCoeff = correlationNestedArray(j);
        if ~isnan(correlationCoeff)
            allCorrelations = [allCorrelations; correlationCoeff];
        end
    end
end

figure;
histogram(allCorrelations);
xlabel('Correlation coefficient');
ylabel('Frequency');

hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;

% only_shk_responsive_corrs = allCorrelations(kmeans_idx' == 3);
% only_shk_responsive_corrs = allCorrelations(prechoice_block_1 == 1);
% only_shk_responsive_corrs = allCorrelations(postchoice_reward_block_1 == 1);
only_shk_responsive_corrs = allCorrelations(collect_block_1 == 1);
% only_shk_responsive_corrs = allCorrelations(prechoice_blocks_2_and_3 == 1);
% not_shk_responsive_corrs = allCorrelations(prechoice_block_1 ~=1);
% not_shk_responsive_corrs = allCorrelations(kmeans_idx' ~= 3);c
not_shk_responsive_corrs = allCorrelations(true_neutral ==1);

figure;
histogram(only_shk_responsive_corrs);
xlabel('Correlation coefficient');

hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;


mean_only_shk = mean(only_shk_responsive_corrs);
mean_not_shk = mean(not_shk_responsive_corrs);

figure;
width = 250;
height = 250; 
set(gcf, 'Position', [50, 25, width, height]); 
histogram(not_shk_responsive_corrs , 'Normalization', 'probability', 'FaceColor', 'blue','BinWidth', 0.05,'LineStyle','none');
hold on;


histogram(only_shk_responsive_corrs, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05, 'LineStyle','none');
xline(mean_only_shk, 'r')
xline(mean_not_shk, 'g')
xlabel('Correlation coefficient');
ylabel('Probability');


yLimits = [0 0.15];
plot([0 0], yLimits, 'k', 'LineWidth', 2);
xtickformat('%.2f');
ytickformat('%.2f');
hold off;
[h, p, k] = kstest2(not_shk_responsive_corrs , only_shk_responsive_corrs)

% stats
fprintf('Kolmogorov-Smirnov test result:\n');
fprintf('h = %d (0 means the null hypothesis cannot be rejected, 1 means it can be rejected)\n', h);
fprintf('p-value = %.4f\n', p);

[h,p,ci,stats] = ttest2(not_shk_responsive_corrs , only_shk_responsive_corrs)

bar_separation_value = 3;

figure;
width = 250; 
height = 250; 
set(gcf, 'Position', [50, 25, width, height]); 
swarmchart(ones(1, length(only_shk_responsive_corrs)), only_shk_responsive_corrs)
hold on
swarmchart(ones(1, length(not_shk_responsive_corrs))*bar_separation_value, not_shk_responsive_corrs)

plot([0.5; 1.5], [mean(only_shk_responsive_corrs); mean(only_shk_responsive_corrs)], 'LineWidth',3)
plot([bar_separation_value-.5; bar_separation_value+.5], [mean(not_shk_responsive_corrs); mean(not_shk_responsive_corrs)], 'LineWidth',3)
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off

%% Fig. 2K and Fig. 2L
% plot scatters for individual neurons & behav variables



% % prechoice representative: 
% find(correlationResults{5, 1} < -0.3)
% start_time = -4;% sub-window start time
% end_time = 0; % sub-window end time
% sub_window_idx = ts1 >= start_time & ts1 <= end_time;
% sub_window_activity_session_1 = zall_mouse{5, 1}{1, 44}(:, sub_window_idx);
% r_val_for_representative = correlationResults{5, 1}(1, 44)
% p_val_for_representative = correlationResults_sig{5, 1}(1, 44)
% choice_times_mouse = trial_choice_times_by_mouse{1, 5};
% trial_types = trial_types_by_mouse{1, 5};

% postchoice representative:
find(correlationResults{5, 1} < -0.3)
start_time = 0;% sub-window start time
end_time = 2; % sub-window end time
sub_window_idx = ts1 >= start_time & ts1 <= end_time;
sub_window_activity_session_1 = zall_mouse{5, 1}{1, 65}(:, sub_window_idx);
r_val_for_representative = correlationResults{5, 1}(1, 65)
p_val_for_representative = correlationResults_sig{5, 1}(1, 65)
choice_times_mouse = trial_choice_times_by_mouse{1, 5};
trial_types = trial_types_by_mouse{1, 5};

%consumption representative:
% find(correlationResults{5, 1} > 0.3)
% start_time = 1;% sub-window start time
% end_time = 3; % sub-window end time
% sub_window_idx = ts1 >= start_time & ts1 <= end_time;
% sub_window_activity_session_1 = zall_mouse{5, 3}{1, 1}(:, sub_window_idx);
% choice_times_mouse = consum_times_by_mouse{1, 5};
% trial_types = trial_types_by_mouse{1, 5};


mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);


x = mean_sub_window_activity_session_1;
y = choice_times_mouse;
size(y)

colors = repmat([0.5, 0.5, 0.5], length(trial_types), 1); % Default to gray
colors(trial_types == 1.2, :) = repmat([0, 0, 1], sum(trial_types == 1.2), 1); % Blue for trial_types == 1.2
colors(trial_types == 0.3, :) = repmat([1, 0, 0], sum(trial_types == 0.3), 1); % Red for trial_types == 0.3

figure;
set(gcf, 'Position', [100, 100, 200, 200]); % Adjust figure position and size
scatter(x, y, 36, colors, 'filled', 'MarkerEdgeColor', 'k'); % Use 'colors' for MarkerFaceColor

hold on;

x_large = x(trial_types == 1.2);
y_large = y(trial_types == 1.2);
coefficients_large = polyfit(x_large, y_large, 1);
x_fit_large = linspace(min(x_large), max(x_large), 100);
y_fit_large = polyval(coefficients_large, x_fit_large);
plot(x_fit_large, y_fit_large, 'b', 'LineWidth', 2);

y_pred_large = polyval(coefficients_large, x_large);
ssr_large = sum((y_pred_large - mean(y_large)).^2);
sst_large = sum((y_large - mean(y_large)).^2);
r_squared_large = ssr_large / sst_large;

x_small = x(trial_types == 0.3);
y_small = y(trial_types == 0.3);
coefficients_small = polyfit(x_small, y_small, 1);
x_fit_small = linspace(min(x_small), max(x_small), 100);
y_fit_small = polyval(coefficients_small, x_fit_small);
plot(x_fit_small, y_fit_small, 'r', 'LineWidth', 2);

y_pred_small = polyval(coefficients_small, x_small);
ssr_small = sum((y_pred_small - mean(y_small)).^2);
sst_small = sum((y_small - mean(y_small)).^2);
r_squared_small = ssr_small / sst_small;

text(min(x) + 0.1, max(y) - 0.1, ['Large R^2 = ' num2str(r_squared_large, '%.3f')], 'FontSize', 10, 'Color', 'b');
text(min(x) + 0.1, max(y) - 0.3, ['Small R^2 = ' num2str(r_squared_small, '%.3f')], 'FontSize', 10, 'Color', 'r');

xlabel('Mean dF/F activity');
ylabel('Latency');
hold off;


