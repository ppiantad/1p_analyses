

mouse_to_check = 9; 

prechoice_block_2_3_column_1_data = zall_mouse{mouse_to_check, 1}(prechoice_remapped_mouse{mouse_to_check, 1} == 1);
prechoice_block_2_3_column_2_data = zall_mouse{mouse_to_check, 5} (prechoice_remapped_mouse{mouse_to_check, 1} == 1);

behav_part_1 = behav_tbl_iter{1, 1}{mouse_to_check, 1}  
behav_part_2 = behav_tbl_iter{5, 1}{mouse_to_check, 1} 
concat_behav = vertcat(behav_part_1, behav_part_2 )

% Initialize the concatenated cell array
concatenated_columns = cell(1, size(prechoice_block_2_3_column_1_data, 2));


% Iterate through each column and concatenate the data
for i = 1:size(prechoice_block_2_3_column_1_data, 2)
    concatenated_columns{i} = vertcat(prechoice_block_2_3_column_1_data{i}, ...
                                      prechoice_block_2_3_column_2_data{i});
end

% Initialize an array to store the resulting mean values
[numRows, numCols] = size(concatenated_columns{1}); % Assume all cells have the same dimensions
mean_array = zeros(numRows, numCols); % Resulting double array

% Loop through each cell array
for i = 1:length(concatenated_columns)
    % Add the values from the current cell array to the mean calculation
    mean_array = mean_array + concatenated_columns{i};
end

% Divide by the number of cell arrays to calculate the mean
mean_array = mean_array / length(concatenated_columns);

figure; imagesc(ts1, [], mean_array)

% figure; plot(concat_behav.bigSmall)



%% TO RUN CODE BELOW, load 10x variables then run eventRelatedActivityAndClassification_PTP_v4:
% 1. choiceTime 'OMITALL', 0, 'BLANK_TOUCH', 0
% 2. collectionTime'OMITALL', 0, 'BLANK_TOUCH', 0
%%
if size(respClass_all_array, 2) == 10 | size(respClass_all_array, 2) == 11 | size(respClass_all_array, 2) == 12
    comparison_arrays_full = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays_full = [1 2 3; 4 5 6]
elseif size(respClass_all_array, 2) == 7
    comparison_arrays_full = [1 2 3; 5 6 7]
end

% ca_data_type = uv.ca_data_type

% Assume 'ca' is your calcium imaging data, sampled every 100 ms
% bin_size is the desired bin size in ms
bin_size = 100;  % Set this as needed
bin_factor = bin_size / (uv.dt*1000);  % Determine how many samples to bin together

first_session = session_to_analyze;

ts1 = (uv.evtWin(1):bin_factor/10:uv.evtWin(2)-0.1);
%%
% for aa = 1:size(comparison_arrays_full, 1) %size(comparison_arrays_full, 1)
% 
%     comparison_arrays = comparison_arrays_full(aa, :)
%     first_session = session_to_analyze;

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





    %%
    prechoice_block_1_column_1_data = {};
    prechoice_mean_mouse = [];
    prechoice_sd_mouse = [];
    prechoice_sem_mouse = [];
    prechoice_iqr = [];
    prechoice_block_1_column_1_data = zall_mouse{select_mouse_index, 11}(prechoice_block_1_mouse{select_mouse_index, 1} == 1);
    for cc = 1:size(prechoice_block_1_column_1_data, 2)
        prechoice_data = prechoice_block_1_column_1_data{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
        prechoice_sd_mouse(:, cc) = std(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 0, 2);
        prechoice_sem_mouse(:, cc) = std(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 0, 2)/sqrt(size(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2));
        prechoice_iqr(:, cc) = iqr(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
        
    end
    prechoice_block_1_mean_mouse_array{gg} = prechoice_mean_mouse;
    prechoice_block_1_sd_mouse_array{gg} = prechoice_sd_mouse;

    prechoice_block_2_3_column_1_data = {};
    prechoice_mean_mouse = [];
    prechoice_block_2_3_column_1_data = zall_mouse{select_mouse_index, 11}(prechoice_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);

    for cc = 1:size(prechoice_block_2_3_column_1_data, 2)
        prechoice_data = prechoice_block_2_3_column_1_data{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    end
    prechoice_block_2_3_mean_mouse_array{gg} = prechoice_mean_mouse;

    prechoice_block_all_column_1_data = {};
    prechoice_mean_mouse = [];
    prechoice_block_all_column_1_data = zall_mouse{select_mouse_index, 11};

    for cc = 1:size(prechoice_block_all_column_1_data, 2)
        prechoice_data = prechoice_block_all_column_1_data{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    end
    prechoice_all_mean_mouse_array{gg} = prechoice_mean_mouse;


    prechoice_remapped_column_1_data = {};
    prechoice_mean_mouse = [];
    prechoice_remapped_column_1_data = zall_mouse{select_mouse_index, 11}(prechoice_remapped_mouse{select_mouse_index, 1} == 1);
    for cc = 1:size(prechoice_remapped_column_1_data, 2)
        prechoice_data = prechoice_remapped_column_1_data{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    end
    prechoice_remapped_mean_mouse_array{gg} = prechoice_mean_mouse;



    %%
    postchoice_block_1_column_1_data = {};
    postchoice_mean_mouse = [];
    postchoice_block_1_column_1_data = zall_mouse{select_mouse_index, 11}(postchoice_reward_block_1_mouse{select_mouse_index, 1} == 1);
    for cc = 1:size(postchoice_block_1_column_1_data, 2)
        postchoice_data = postchoice_block_1_column_1_data{1, cc};
        postchoice_mean_mouse(:, cc) = mean(postchoice_data (:, ts1 >= 0 & ts1 <= 2), 2);
    end
    postchoice_block_1_mean_mouse_array{gg} = postchoice_mean_mouse;

    postchoice_block_2_3_column_1_data = {};
    postchoice_mean_mouse = [];
    postchoice_block_2_3_column_1_data = zall_mouse{select_mouse_index, 11}(postchoice_reward_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);

    for cc = 1:size(postchoice_block_2_3_column_1_data, 2)
        postchoice_data = postchoice_block_2_3_column_1_data{1, cc};
        postchoice_mean_mouse(:, cc) = mean(postchoice_data (:, ts1 >= 0 & ts1 <= 2), 2);
    end
    postchoice_block_2_3_mean_mouse_array{gg} = postchoice_mean_mouse;

    postchoice_block_all_column_1_data = {};
    postchoice_mean_mouse = [];
    postchoice_block_all_column_1_data = zall_mouse{select_mouse_index, 11};

    for cc = 1:size(postchoice_block_all_column_1_data, 2)
        postchoice_data = postchoice_block_all_column_1_data{1, cc};
        postchoice_mean_mouse(:, cc) = mean(postchoice_data (:, ts1 >= 0 & ts1 <= 2), 2);
    end
    postchoice_all_mean_mouse_array{gg} = postchoice_mean_mouse;

    %%
    collect_block_1_column_1_data = {};
    collect_mean_mouse = [];
    collect_block_1_column_1_data = zall_mouse{select_mouse_index, 12}(collect_block_1_mouse{select_mouse_index, 1} == 1);
    for cc = 1:size(collect_block_1_column_1_data, 2)
        collect_data = collect_block_1_column_1_data{1, cc};
        collect_mean_mouse(:, cc) = mean(collect_data (:, ts1 >= 1 & ts1 <= 3), 2);
    end
    collect_block_1_mean_mouse_array{gg} = collect_mean_mouse;

    collect_block_2_3_column_1_data = {};
    collect_mean_mouse = [];
    collect_block_2_3_column_1_data = zall_mouse{select_mouse_index, 12}(collect_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);

    for cc = 1:size(collect_block_2_3_column_1_data, 2)
        collect_data = collect_block_2_3_column_1_data{1, cc};
        collect_mean_mouse(:, cc) = mean(collect_data (:, ts1 >= 1 & ts1 <= 3), 2);
    end
    collect_block_2_3_mean_mouse_array{gg} = collect_mean_mouse;

    collecte_block_all_column_1_data = {};
    collect_mean_mouse = [];
    collect_block_all_column_1_data = zall_mouse{select_mouse_index, 12};

    for cc = 1:size(collect_block_all_column_1_data, 2)
        collect_data = collect_block_all_column_1_data{1, cc};
        collect_mean_mouse(:, cc) = mean(collect_data (:, ts1 >= 1 & ts1 <= 3), 2);
    end
    collect_all_mean_mouse_array{gg} = collect_mean_mouse;

end
% end

%%
for hh = 1:size(prechoice_all_mean_mouse_array, 2)

    extracted = prechoice_all_mean_mouse_array{hh};
    session_long_prechoice(:, hh) = mean(extracted, 2);

end

mean_session_long_prechoice = mean(session_long_prechoice, 2)



for hh = 1:size(prechoice_remapped_mean_mouse_array, 2)

    extracted = prechoice_remapped_mean_mouse_array{hh};
    session_long_prechoice_block_1(:, hh) = mean(extracted, 2);

end

mean_session_long_prechoice = mean(session_long_prechoice_block_1, 2)

for hh = 1:size(prechoice_block_2_3_mean_mouse_array, 2)

    extracted = prechoice_block_2_3_mean_mouse_array{hh};
    session_long_prechoice_block_2_3(:, hh) = mean(extracted, 2);

end

mean_session_long_prechoice = mean(session_long_prechoice_block_2_3, 2)

%%
% Create the figure
figure;

% Plot the first dataset on the left Y-axis
yyaxis left;
plot(mean_session_long_prechoice, '-b', 'LineWidth', 1.5); % Blue line for the first dataset
ylabel('Pre-choice ensemble PV'); % Label for the left Y-axis
xlabel('Index'); % Common X-axis label

% Plot the second dataset on the right Y-axis
yyaxis right;
plot(meanData_large, '-r', 'LineWidth', 1.5); % Red line for the second dataset
ylabel('% large rew choice'); % Label for the right Y-axis

%%

first_data_to_plot = prechoice_block_1_mean_mouse_array;
second_data_to_plot = prechoice_remapped_mean_mouse_array;
behav_sequence_data_to_plot = large_sequences_mouse;



% Calculate the mean across rows for each cell in prechoice_remapped_mean_mouse_array
mean_first_data = cellfun(@(x) mean(x, 2), first_data_to_plot, 'UniformOutput', false);

mean_second_data = cellfun(@(x) mean(x, 2), second_data_to_plot, 'UniformOutput', false);

% Calculate the overall mean across all cells for prechoice_block_1_mean_mouse_array
overall_mean_first_data = mean(cell2mat(mean_first_data), 2); % Combine all cells and take the mean

overall_mean_second_data = mean(cell2mat(mean_second_data), 2); % Combine all cells and take the mean

% Calculate the overall mean across all rows for large_sequences_mouse
overall_mean_sequence = mean(behav_sequence_data_to_plot, 1)'; % Mean across rows


% Number of subplots (equal to the number of cells in the array)
num_subplots = numel(first_data_to_plot);

% Create a figure
figure;

% Loop through each cell and create a subplot
for i = 1:num_subplots
    % Extract the mean data for the current cell
    mean_data = mean_first_data{i}; % 90x1 double
    
    % Extract the corresponding row from large_sequences_mouse
    sequence_data = behav_sequence_data_to_plot(i, :)'; % 90x1 double
    
    % Create a subplot
    subplot(ceil(sqrt(num_subplots)), ceil(sqrt(num_subplots)), i);
    
    % Create a left Y-axis plot for the large_sequences_mouse data
    yyaxis right;
    plot(sequence_data, 'r', 'DisplayName', 'Sequence Data');
    ylabel('Large Rew choice');
    
    % Create a right Y-axis plot for the prechoice_remapped_mean_mouse_array data
    yyaxis left;
    plot(mean_data, 'b', 'DisplayName', 'Mean Data');
    ylabel('Mean PV');
    ylim([-0.1 0.5]);

    % Add title and labels
    title(['Mouse ' num2str(i)]);
    xlabel('Trial');
end



% Create a new figure
figure;

% Create a left Y-axis plot for the large_sequences_mouse mean
yyaxis right;
h1 = plot(overall_mean_sequence, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Overall Large Rew choice'); % Red solid line
ylabel('Overall Large Rew choice');

% Create a right Y-axis plot for the prechoice_block_1_mean_mouse_array mean
yyaxis left;
h2 = plot(overall_mean_first_data, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Overall Mean PV Block 1'); % Blue solid line
hold on;
h3 = plot(overall_mean_second_data, 'b', 'LineWidth', 1.5, 'DisplayName', 'Overall Mean PV Blocks 2/3'); % Green solid line
ylabel('Overall Mean PV');
ylim([-0.1 0.5]);

xlabel('Trial');

% Create the legend using the plot handles
legend([h1, h2, h3], {'Overall Large Rew choice', 'Overall Mean PV Block 1', 'Overall Mean PV Blocks 2/3'}, 'Location', 'best');

%%
concatenatedData = [];
concatenatedData = horzcat(prechoice_block_1_mean_mouse_array{:});
concatenatedData_2 = horzcat(prechoice_block_2_3_mean_mouse_array{:});
concatenatedData_3 = horzcat(prechoice_remapped_mean_mouse_array{:});

% Assuming concatenatedData is a 90x163 array
rowsPerPart = size(concatenatedData, 1) / 3; % Determine number of rows per part
part1 = concatenatedData(1:rowsPerPart, :); % First third
part2 = concatenatedData(rowsPerPart+1:2*rowsPerPart, :); % Second third
part3 = concatenatedData(2*rowsPerPart+1:end, :); % Third third

% Calculate the standard deviation for each column in each part
std_part1 = std(part1);
std_part2 = std(part2);
std_part3 = std(part3);

% Combine the results into a single array
std_per_column = [std_part1; std_part2; std_part3];

std_trials = std(concatenatedData, 0, 2);

figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
% xtickformat('%.2f');
ytickformat('%.2f');
hold on;
h(1) = shadedErrorBar(1:90, nanmean(concatenatedData, 2), iqr(concatenatedData, 2), 'lineProps', {'color', acton(1,:)});
ylim([-0.6 1])
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(1)}(remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(1)}(remapped  ==1, :)), 'lineProps', {'color', acton(50,:)});

% Transpose std_per_column for row-wise boxplot
data_for_boxplot = std_per_column';

% Create a figure
figure;
hold on;

% Boxplot for each row (across all columns)
boxplot(data_for_boxplot, 'Labels', 1:size(std_per_column, 1));

% % Overlay individual data points
% for i = 1:size(data_for_boxplot, 1) % Loop through each column of std_per_column
%     plot(1:size(std_per_column, 1), data_for_boxplot(i, :), '-o', 'LineWidth', 1.5);
% end

% Customize the plot
xlabel('Row Index (Part)');
ylabel('Standard Deviation');
title('Row-wise Boxplots with Individual Data Points');
legend('show', 'Location', 'bestoutside');
grid on;
hold off;
%%
% Perform Kruskal-Wallis test
[p, tbl, stats] = kruskalwallis(data_for_boxplot, [], 'off');

% Display the p-value
disp(['Kruskal-Wallis Test p-value: ', num2str(p)]);

% Perform pairwise comparisons
pairwise_results = multcompare(stats);

% Display pairwise results
disp('Pairwise Comparisons Results:');
disp(array2table(pairwise_results, ...
    'VariableNames', {'Group1', 'Group2', 'LowerLimit', 'Estimate', 'UpperLimit', 'pValue'}));

%%
% trying to find the best way to track, on a trial x trial basis,
% "ensemble" membership
%% TO RUN CODE BELOW, load 10x variables then run eventRelatedActivityAndClassification_PTP_v4:
% 1. choiceTime 'OMITALL', 0, 'BLANK_TOUCH', 0
% 2. collectionTime'OMITALL', 0, 'BLANK_TOUCH', 0


pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [1 3];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

% Find indices corresponding to each time window
pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

indices_to_use = pre_choice_indices;

% zall(11, :) should contain choice aligned data, all 90 trials
session_long_data = zall_array(1, :);
% session_long_data = zall_array(12, :);


data_for_shuffling = zall_array(1, :);
% data_for_shuffling = zall_array(10, :);
% data_for_shuffling = zall_array(12, :);

% get the mean of the first 30 trials for each neuron - this is how Block 1
% neurons are determined

for gg = 1:size(session_long_data, 2)
    current_session_long_data = session_long_data{1, gg};
    current_data_for_shuffling = data_for_shuffling{1, gg};
    % if doing ALL neurons, subselect 30 trials from the larger array
    % current_data_for_shuffling = current_data_for_shuffling(randperm(size(current_data_for_shuffling, 1), 30), :);

    % 9/13/2024
    % alternate way of shuffling, maybe worth trying
    [trial_num, sample_num] = size(current_data_for_shuffling);
    % [trial_num, sample_num] = size(caTraceTrials);
    % shift_val = randi(sample_num)
    for g = 1:uv.resamples                                              %for each resampling of the data
        %sort the data index to create a new list of indices
        % here I am using the data that were sub-selected to be the same
        % size as Block 1 (30 trials). You could also replace
        % current_data_for_shuffling with current_session_long_data and use
        % 1:30 for Block 1, and 31:90 for Blocks 2/3
        for t = 1:size(current_data_for_shuffling, 1) %31:90 1:30

            shift_val = randi(sample_num); %for each trial
            shuffledTrace(t,:) = circshift(current_data_for_shuffling(t,:), shift_val,2);     %shuffle the calcium trace
            %         shuffledEvtRate(t,:) = caEvtRateTrials(t,shuffledIDX(t,:)); %shuffle the event rate
        end
        nullDistTrace(g,:) = nanmean(shuffledTrace);                    %calculate the NaN mean of the shuffled traces
        %     nullDistEvtRate(g,:) = nanmean(shuffledEvtRate);                %calculate the NaN mean of the shuffled event rates
    end
    clear shuffled* g t trialCt
    nullDist = nullDistTrace;                                       %direct transfer
    % empiricalTrialWin = nanmean...
    %     (current_session_long_data(:,ts1 >= -4 & ts1 <= 0));                   %NaN mean of the fluorescent response across trials, within the time window. this gets the within trial mean.
    % empiricalSEM = nansem...
    %     (current_session_long_data(:,:));

    % empiricalWinAvg = nanmean(empiricalTrialWin);                   %across trial mean
    % empiricalSEMAvg = nanmean(empiricalSEM);

    sdNull = mean(nanstd(nullDist(:, indices_to_use)));   
    
    upperSD = nanmean(nullDist(:)) + (uv.sigma*sdNull);                    %calculate upper limit of interval around the mean
    % lowerSD = nanmean(nullDist(:)) - (uv.sigma*sdNull);

    % prechoice_block_1_mean_mouse(:, gg) = mean(current_session_long_data (:, ts1 >= -4 & ts1 <= 0), 2);
    prechoice_all_blocks_mean_mouse(:, gg) = mean(current_session_long_data (:, indices_to_use), 2);
    
    prechoice_logical_ind(:, gg) = prechoice_all_blocks_mean_mouse(:, gg) > upperSD ;

    % if prechoice_all_blocks_mean_mouse(:, gg) > upperSD  %empiricalWinAvg > upperSD & empiricalWinAvg > otherEventWinAvg1 & empiricalWinAvg > otherEventWinAvg2
    %     prechoice_logical_ind(:, gg) = 1;
    %     % respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args)(u,1) = 1;
    %     % respClass_all_array_mouse{ii, iter}(u) = 1;
    % % elseif empiricalWinAvg < lowerSD  % empiricalWinAvg < lowerSD & empiricalWinAvg < otherPeriodWinAvg & empiricalWinAvg < otherEventWinAvg2
    % %     respClass_all(neuron_num) = 2;
    % %     respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args)(u,1) = 2;
    % %     respClass_all_array_mouse{ii, iter}(u) = 2;
    % else
    %     prechoice_logical_ind(:, gg) = 0;
    %     % respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args)(u,1) = 3;
    %     % respClass_all_array_mouse{ii, iter}(u) = 3;
    % end
    % 
    % % prechoice_logical_ind(:, gg) = prechoice_all_blocks_mean_mouse(:, gg) > prechoice_block_1_mean_mouse(:, gg)
    % 


end

only_prechoice_block_1 = prechoice_logical_ind(:, prechoice_block_1 == 1);
% only_prechoice_block_1 = prechoice_logical_ind;

hold on
figure; plot(1:90, mean(only_prechoice_block_1, 2))
ylim([0.3 .8])
hold off