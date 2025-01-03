
load('BLA_C_raw_no_additional_filtering_Pre_RDT_RM_only_completed_sessions_zall_window_base_workspace_3_categories.mat')



%%
% arrays_to_examine = [1 8];

inhib_or_excite = 1;

event_for_figures = 1; 

%%


for kk = 1:size(respClass_all_array_mouse, 1)
    respClass_data = respClass_all_array_mouse(kk, :);
    
    prechoice{kk, 1} = respClass_data{1, 1} == inhib_or_excite & respClass_data{1, 2} ~= inhib_or_excite & respClass_data{1, 3} ~= inhib_or_excite;
    prechoice_inds{kk, 1} = find(respClass_data{1, 1} == inhib_or_excite & respClass_data{1, 2} ~= inhib_or_excite & respClass_data{1, 3} ~= inhib_or_excite);

    % prechoice_no_consumption{kk, 1} = respClass_data{1, 1} == inhib_or_excite & respClass_data{1, 2} ~= inhib_or_excite & respClass_data{1, 3} ~= inhib_or_excite;
    % prechoice_inds{kk, 1} = find(respClass_data{1, 1} == inhib_or_excite & respClass_data{1, 2} ~= inhib_or_excite & respClass_data{1, 3} ~= inhib_or_excite);


    % 
    postchoice{kk, 1} = respClass_data{1, 2} == inhib_or_excite & respClass_data{1, 1} ~= inhib_or_excite & respClass_data{1, 3} ~= inhib_or_excite;
    postchoice_inds{kk, 1} = find(respClass_data{1, 2} == inhib_or_excite & respClass_data{1, 1} ~= inhib_or_excite & respClass_data{1, 3} ~= inhib_or_excite);
    % 
    collect{kk, 1} = respClass_data{1, 3} == inhib_or_excite & respClass_data{1, 1} ~= inhib_or_excite & respClass_data{1, 2} ~= inhib_or_excite;
    collect_inds{kk, 1} = find(respClass_data{1, 3} == inhib_or_excite & respClass_data{1, 1} ~= inhib_or_excite & respClass_data{1, 2} ~= inhib_or_excite);

    prechoice_no_consumption_mouse{kk, 1} = respClass_data{1, 1} == inhib_or_excite & respClass_data{1, 2} ~= inhib_or_excite & respClass_data{1, 3} == 3;
    prechoice_no_consumption_mouse_inds{kk, 1} = find(respClass_data{1, 1} == inhib_or_excite & respClass_data{1, 2} ~= inhib_or_excite & respClass_data{1, 3} == 3);
end



%%

relevant_behav_array = behav_tbl_iter{1, 1};


for zz = 1:size(relevant_behav_array, 1)
    BehavData = relevant_behav_array{zz};
    
    % Create a logi0.10cal index array based on your conditions
    logical_index_long_ITI{zz} = BehavData.stTime - BehavData.TrialPossible >= 10 & BehavData.stTime - BehavData.TrialPossible <= 50;
    logical_index_short_ITI{zz} = BehavData.stTime - BehavData.TrialPossible <= 10 & BehavData.stTime - BehavData.TrialPossible <= 50;
    % Use the logical index array to subset BehavData
    % BehavData = BehavData(logical_index,: );
    % trials = trials(logical_index);

end

%%
% only_prechoice_neurons = zall_array(1, prechoice == 1);

only_choice_aligned_array =  zall_mouse(:, 1);

mean_data_from_cells = []; 

count = 1; 
count_prechoice = 1;
count_postchoice = 1;
count_collect = 1; 
count_prechoice_no_collect = 1; 
for pp = 1:size(only_choice_aligned_array, 1)
    mouse_data = only_choice_aligned_array{pp, 1};
    mouse_data_prechoice = mouse_data(prechoice{pp, 1} == 1);
    mouse_data_postchoice = mouse_data(postchoice{pp, 1} == 1);
    mouse_data_collect = mouse_data(collect{pp, 1} == 1);
    mouse_data_prechoice_no_consumption = mouse_data(prechoice_no_consumption_mouse{pp, 1} == 1);
    for gg = 1:size(mouse_data_prechoice, 2)
        % cell_data = mouse_data{gg};
        % cell_data = cell_data(logical_index{1, pp}, :);
        % mean_data_from_cells(count, :) = mean(cell_data);


        % prechoice_cell_data = mouse_data_prechoice{gg};
        prechoice_cell_data = mouse_data_prechoice{gg};
        prechoice_cell_data_long = prechoice_cell_data(logical_index_long_ITI{1, pp}, :);
        prechoice_cell_data_short = prechoice_cell_data(logical_index_short_ITI{1, pp}, :);
        prechoice_mean_data_from_cells_long(count_prechoice, :) = mean(prechoice_cell_data_long);
        prechoice_sem_data_from_cells_long(count_prechoice, :)= nanstd(prechoice_cell_data_long,1)/(sqrt(size(prechoice_cell_data_long, 1)));
        prechoice_mean_data_from_cells_short(count_prechoice, :) = mean(prechoice_cell_data_short);
        prechoice_sem_data_from_cells_short(count_prechoice, :)= nanstd(prechoice_cell_data_short,1)/(sqrt(size(prechoice_cell_data_short, 1)));
        count_prechoice = count_prechoice+1;
    end

    for gg = 1:size(mouse_data_postchoice, 2)
        postchoice_cell_data = mouse_data_postchoice{gg};
        postchoice_cell_data_long = postchoice_cell_data(logical_index_long_ITI{1, pp}, :);
        postchoice_cell_data_short = postchoice_cell_data(logical_index_short_ITI{1, pp}, :);
        postchoice_mean_data_from_cells_long(count_postchoice, :) = mean(postchoice_cell_data_long);
        postchoice_sem_data_from_cells_long(count_postchoice, :)= nanstd(postchoice_cell_data_long,1)/(sqrt(size(postchoice_cell_data_long, 1)));
        postchoice_mean_data_from_cells_short(count_postchoice, :) = mean(postchoice_cell_data_short);
        postchoice_sem_data_from_cells_short(count_postchoice, :)= nanstd(postchoice_cell_data_short,1)/(sqrt(size(postchoice_cell_data_short, 1)));
        count_postchoice = count_postchoice+1;
    end


    for gg = 1:size(mouse_data_collect, 2)
        collect_cell_data = mouse_data_collect{gg};
        collect_cell_data_long = collect_cell_data(logical_index_long_ITI{1, pp}, :);
        collect_cell_data_short = collect_cell_data(logical_index_short_ITI{1, pp}, :);
        collect_mean_data_from_cells_long(count_collect, :) = mean(collect_cell_data_long);
        collect_sem_data_from_cells_long(count_collect, :)= nanstd(collect_cell_data_long,1)/(sqrt(size(collect_cell_data_long, 1)));
        collect_mean_data_from_cells_short(count_collect, :) = mean(collect_cell_data_short);
        collect_sem_data_from_cells_short(count_collect, :)= nanstd(collect_cell_data_short,1)/(sqrt(size(collect_cell_data_short, 1)));
        count_collect = count_collect+1;

    end

    for gg = 1:size(mouse_data_prechoice_no_consumption, 2)
        prechoice_no_consumption_data = mouse_data_prechoice_no_consumption{gg};
        prechoice_no_consumption_data_long = prechoice_no_consumption_data(logical_index_long_ITI{1, pp}, :);
        prechoice_no_consumption_data_short = prechoice_no_consumption_data(logical_index_short_ITI{1, pp}, :);
        prechoice_no_consumption_mean_data_from_cells_long(count_prechoice_no_collect, :) = mean(prechoice_no_consumption_data_long);
        prechoice_no_consumption_sem_data_from_cells_long(count_prechoice_no_collect, :)= nanstd(prechoice_no_consumption_data_long,1)/(sqrt(size(prechoice_no_consumption_data_long, 1)));
        prechoice_no_consumption_mean_data_from_cells_short(count_prechoice_no_collect, :) = mean(prechoice_no_consumption_data_short);
        prechoice_no_consumption_sem_data_from_cells_short(count_prechoice_no_collect, :)= nanstd(prechoice_no_consumption_data_short,1)/(sqrt(size(prechoice_no_consumption_data_short, 1)));
        count_prechoice_no_collect = count_prechoice_no_collect+1;

    end
    
end


%%


figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(prechoice_mean_data_from_cells_long), nanmean(prechoice_sem_data_from_cells_long), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(prechoice_mean_data_from_cells_short), nanmean(prechoice_sem_data_from_cells_short), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], 'long ITI', 'short ITI interval')
xlim([-8 8]);
ylim([-0.6 0.8])
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');
%%


figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(prechoice_no_consumption_mean_data_from_cells_long), nanmean(prechoice_no_consumption_sem_data_from_cells_long), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(prechoice_no_consumption_mean_data_from_cells_short), nanmean(prechoice_no_consumption_sem_data_from_cells_short), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], 'long ITI', 'short ITI interval')
xlim([-8 8]);
ylim([-0.6 0.8])
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');


%%


figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(postchoice_mean_data_from_cells_long), nanmean(postchoice_sem_data_from_cells_long), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(postchoice_mean_data_from_cells_short), nanmean(postchoice_sem_data_from_cells_short), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-8 8]);
ylim([-0.6 0.8])
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

%%


figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(collect_mean_data_from_cells_long), nanmean(collect_sem_data_from_cells_long), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(collect_mean_data_from_cells_short), nanmean(collect_sem_data_from_cells_short), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-8 8]);
ylim([-0.6 0.8])
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

%%
prechoice_no_consumption = respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 3; 
prechoice_no_consumption_inds = find(respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 3); 



%%


figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(collect_mean_data_from_cells_long), nanmean(collect_sem_data_from_cells_long), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(collect_mean_data_from_cells_short), nanmean(collect_sem_data_from_cells_short), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-8 8]);
ylim([-0.6 0.8])
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');