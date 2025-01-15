
load('acton.mat')
load('batlowW.mat')

% match the number of neurons in each group (if necessary - data will stay
% the same if sizes are equal
clear idx temp
for i = 1:length(zall_mean_all_array)

    temp(i) = (size(zall_mean_all_array{1, i},1));
end

[minSize,minIdx] = min(temp);

for i = 1:length(zall_mean_all_array)
    if i ~= minIdx
        idx_for_subsample = randperm(temp(i),minSize);
        idx_for_subsample = sort(idx_for_subsample);
        zdataTemp{i} = zall_mean_all_array{1, i}(idx_for_subsample,:);
    else
        zdataTemp{i} = zall_mean_all_array{1, i};

    end
end



%DO NOT CHANGE FROM THESE SETTINGS UNLESS WILLING TO REVERT BACK

% MY APPROACH
% % use zall array if you want to check how trials compare across block
neuron_mean_concat = horzcat(zdataTemp{:});
% neuron_mean_concat = horzcat(zall_mean_all_array{1, 11}, zall_mean_all_array{1, 12});
%12/14/2024 MAYBE I DONT NEED THIS LINE BELOW? MAYBE IT ARTIFICALLY DEFLATES
%DIFFERENCES?

neuron_mean_concat = zscore(neuron_mean_concat, 0 , 2);


% neuron_mean_concat = horzcat(reformat_zall_mean_array{:});
% neuron_mean_concat = horzcat(neuron_mean_array{:});

% RUAIRI APPROACH
% % use neuron_mean_all_unnorm if you want to check how things differ across
% % time
% neuron_mean_concat = horzcat(neuron_mean_all_unnormalized{:});
% 
% neuron_mean_concat = mean_center_columnwise(neuron_mean_concat);



% % use below to normalize THEN concatenate, which might be a better
% approach?
% Iterate through each cell and apply zscore
% for i = 1:size(neuron_mean_all_unnormalized, 2)
%     % Apply zscore to the data in each cell
%     neuron_mean_all_normalized{i} = zscore(neuron_mean_all_unnormalized{i}, 0, 2);
%     neuron_mean_reformat{i} = neuron_mean_all_unnormalized{i}';
% end
% % % 
% neuron_mean_concat = horzcat(neuron_mean_all_normalized{:});
% neuron_mean_concat = zscore(neuron_mean_concat, 0 , 2);

% neuron_mean_reformat_concat = vertcat(neuron_mean_reformat{:});
% neuron_mean_reformat_concat_normalized = zscore(neuron_mean_reformat_concat, 0 , 1);


%%
% Make behav_tbl_iters the same size (RDT contains one extra column)

% Loop through the first level of behav_tbl_iter
for i = 1:length(behav_tbl_iter)
    % Get the 12x1 cell array from the current cell
    subArray = behav_tbl_iter{i};
    
    % Loop through the second level of the 12x1 cell array
    for j = 1:length(subArray)
        % Check if the current cell contains a table
        if istable(subArray{j})
            % Get the table from the current cell
            tbl = subArray{j};
            
            % Check if the 'type_binary' column exists
            if ismember('type_binary', tbl.Properties.VariableNames)
                % Delete the 'type_binary' column
                tbl.type_binary = [];
            end
            
            % Update the subArray with the modified table
            subArray{j} = tbl;
        end
    end
    
    % Update the behav_tbl_iter with the modified subArray
    behav_tbl_iter{i} = subArray;
end

%%

[concatenatedTable_all, concatenate_all_tables] = get_median_choice_and_collect_fn(behav_tbl_iter);

[full_table_all] = get_median_choice_and_collect_fn(behav_tbl_iter);
full_table = vertcat(full_table_all{:});

median_start_time_large = median(full_table.stTime(full_table.bigSmall == 1.2) - full_table.choiceTime(full_table.bigSmall == 1.2));
median_start_time_small = median(full_table.stTime(full_table.bigSmall == 0.3) - full_table.choiceTime(full_table.bigSmall == 0.3));

median_start_time_block_1 = median(full_table.stTime(full_table.Block == 1) - full_table.choiceTime(full_table.Block == 1));
median_start_time_block_2 = median(full_table.stTime(full_table.Block == 2) - full_table.choiceTime(full_table.Block == 2));
median_start_time_block_3 = median(full_table.stTime(full_table.Block == 3) - full_table.choiceTime(full_table.Block == 3));

median_choice_time_large = median(full_table.choiceTime(full_table.bigSmall == 1.2) - full_table.stTime(full_table.bigSmall == 1.2));
median_choice_time_small = median(full_table.choiceTime(full_table.bigSmall == 0.3) - full_table.stTime(full_table.bigSmall == 0.3));

median_choice_time_block_1 = median(full_table.choiceTime(full_table.Block == 1) - full_table.stTime(full_table.Block == 1));
median_choice_time_block_2 = median(full_table.choiceTime(full_table.Block == 2) - full_table.stTime(full_table.Block == 2));
median_choice_time_block_3 = median(full_table.choiceTime(full_table.Block == 3) - full_table.stTime(full_table.Block == 3));

median_collect_time_large = median(full_table.collectionTime(full_table.bigSmall == 1.2) - full_table.stTime(full_table.bigSmall == 1.2));
median_collect_time_small = median(full_table.collectionTime(full_table.bigSmall == 0.3) - full_table.stTime(full_table.bigSmall == 0.3));

median_collect_time_block_1 = median(full_table.collectionTime(full_table.Block == 1) - full_table.stTime(full_table.Block == 1));
median_collect_time_block_2 = median(full_table.collectionTime(full_table.Block == 2) - full_table.stTime(full_table.Block == 2));
median_collect_time_block_3 = median(full_table.collectionTime(full_table.Block == 3) - full_table.stTime(full_table.Block == 3));

[~, closest_index_start_time_large] = min(abs(ts1 - median_start_time_large));
[~, closest_index_start_time_small] = min(abs(ts1 - median_start_time_small));

[~, closest_index_start_time_block_1] = min(abs(ts1 - median_start_time_block_1));
[~, closest_index_start_time_block_2] = min(abs(ts1 - median_start_time_block_2));
[~, closest_index_start_time_block_3] = min(abs(ts1 - median_start_time_block_3));

[~, closest_index_collect_time_large] = min(abs(ts1 - median_collect_time_large));
[~, closest_index_collect_time_small] = min(abs(ts1 - median_collect_time_small));

[~, closest_index_collect_time_block_1] = min(abs(ts1 - median_collect_time_block_1));
[~, closest_index_collect_time_block_2] = min(abs(ts1 - median_collect_time_block_2));
[~, closest_index_collect_time_block_3] = min(abs(ts1 - median_collect_time_block_3));

for zz = 1:size(concatenatedTable_all, 2)
    median_start_time(1, zz) = median(concatenatedTable_all{1, zz}.stTime - concatenatedTable_all{1, zz}.choiceTime)
    median_collect_time(1, zz) = median(concatenatedTable_all{1, zz}.collectionTime - concatenatedTable_all{1, zz}.choiceTime)
    [~, closest_index_start(zz)] = min(abs(ts1 - median_start_time(1, zz)));
    [~, closest_index_collect(zz)] = min(abs(ts1 - median_collect_time(1, zz)));
    [~, closest_index_zero(zz)] = min(abs(ts1 - 0));
end



%% PCA


numNeuronsPerCondition = neuron_num;


numConditions = size(neuron_mean_concat, 2)/numMeasurements;


eventIdx = 1:numConditions;

NumPC = 4; %2

array_size = size(neuron_mean_concat, 2);

result = cell(size(varargin_list));

% Loop through the input cell array
for i = 1:numel(varargin_list)
    % Initialize an empty string to store the concatenated values
    concat_str = '';
    % Loop through the elements in the 1x4 cell array
    for j = 1:numel(varargin_list{i})
        % Convert doubles to strings and concatenate with underscores
        if isnumeric(varargin_list{i}{j})
            concat_str = [concat_str, num2str(varargin_list{i}{j}), '_'];
        else
            concat_str = [concat_str, varargin_list{i}{j}, '_'];
        end
    end
    % Remove the trailing underscore
    concat_str = concat_str(1:end-1);
    % Store the concatenated string in the result cell array
    result{i} = concat_str;
end

% Convert the result cell array into a 3x1 string array
eventNames = string(result);

disp(result)


for i = 1:numConditions
    start_index = (i - 1) * numMeasurements + 1;
    end_index = i * numMeasurements;
%     start_index = (i - 1) * numNeuronsPerCondition + 1;
%     end_index = i * numNeuronsPerCondition;    
    % Handle the last condition which may have a different end index
    if i == numConditions
        end_index = array_size;
    end
    condition_ranges{i} = {start_index,end_index};
end


for qq = 1:numConditions

    condition_data{qq} = neuron_mean_concat(:, condition_ranges{1, qq}{1}:condition_ranges{1, qq}{2});
    % condition_data{qq} = neuron_mean_concat(condition_ranges{1, qq}{1}:condition_ranges{1, qq}{2}, :);
    
end


[coef,score, latent, ~, explained, ~] = pca(neuron_mean_concat');

for i = eventIdx
    temp = condition_data{i};

    PCScore{i} = coef(1:size(temp,1), 1: NumPC)'*temp;

end


figure;
bar(explained(1:40, :), 'FaceColor', [0.2 0.6 0.8]); % Creates a bar plot with custom color
title('Scree Plot');
xlabel('Principal Component');
ylabel('Variance Explained (%)');
grid on; % Optional: adds a grid to the plot

%% if applying PCA to individual trials (will only work if each mouse has the same # of trials. might be best to investigate how to do this within a mouse

% Example input: zall_array is a 1xN cell array where each cell contains a double array.
num_cells = numel(zall_array);
num_rows = size(zall_array{1}, 1); % Number of rows in the double arrays
num_cols = size(zall_array{1}, 2); % Number of columns in the double arrays

% Initialize the output cell array
zall_array_trial_concat = cell(num_rows, 1);

% Loop through each row of the double arrays
for row_idx = 1:num_rows
    % Initialize an array to concatenate rows from all cells
    concatenated_rows = [];
    
    % Loop through each cell in zall_array
    for cell_idx = 1:num_cells
        % Extract the current row and concatenate it
        concatenated_rows = [concatenated_rows; zall_array{cell_idx}(row_idx, :)];
    end
    
    % Store the concatenated rows in the output cell array
    zall_array_trial_concat{row_idx} = concatenated_rows;
end

%ALSO TRY SPLITTING SESSION INTO 9 PARTS (10 TRIALS EACH) AND RUN PCA 10
%TIMES

%run PCA separately on risky vs. non-risky mice

for hh = 1:size(zall_array_trial_concat, 1)
    current_zall_array_trial_concat = zall_array_trial_concat{hh};

    [coef,score, latent, ~, explained, ~] = pca(current_zall_array_trial_concat');
    coef_by_trial{hh} = coef;
    score_by_trial{hh} = score;
    latent_by_trial{hh} = latent;
    explained_by_trial{hh} = explained;
    clear coef score latent explained


end



for kk = 1:size(explained_by_trial, 2)
    current_explained_by_trial = explained_by_trial{kk};

    variance_explained_by_first_10_PCs(kk) = sum(current_explained_by_trial(1:10));

end
figure; plot(variance_explained_by_first_10_PCs)

% concat_all_trials_all_neurons = horzcat(zall_array_trial_concat{:});
% 
% [coef,score, latent, ~, explained, ~] = pca(concat_all_trials_all_neurons');


% Initialize array to store the number of PCs for each trial
num_PCs_for_90_variance = zeros(1, size(explained_by_trial, 2));

for kk = 1:size(explained_by_trial, 2)
    current_explained_by_trial = explained_by_trial{kk};
    
    % Compute cumulative variance
    cumulative_variance = cumsum(current_explained_by_trial);
    
    % Find the number of PCs needed to reach or exceed 90% variance
    num_PCs_for_90_variance(kk) = find(cumulative_variance >= 90, 1, 'first');
    num_PCs_for_80_variance(kk) = find(cumulative_variance >= 80, 1, 'first');
    num_PCs_for_70_variance(kk) = find(cumulative_variance >= 70, 1, 'first');
    num_PCs_for_60_variance(kk) = find(cumulative_variance >= 60, 1, 'first');
    num_PCs_for_50_variance(kk) = find(cumulative_variance >= 50, 1, 'first');
    num_PCs_for_40_variance(kk) = find(cumulative_variance >= 40, 1, 'first');
end

% Plot the result
figure;
plot(num_PCs_for_90_variance);
hold on; plot(num_PCs_for_80_variance);
hold on; plot(num_PCs_for_70_variance);
hold on; plot(num_PCs_for_60_variance);
hold on; plot(num_PCs_for_50_variance);
hold on; plot(num_PCs_for_40_variance);
xlabel('Trial #');
ylabel('Number of PCs required for variance threshold');
ylim([10 20]);
legend({'90%', '80%', '70%', '60%', '50%', '40%'})





%%
% apply PCA to each trial
for q = 1:size(zall_array_trial_concat, 1)
    temp_trial = zall_array_trial_concat{q, 1};
    temp_trial = zscore(temp_trial, 0, 2);
    PCScore_trials{q, 1} = coef(1:size(temp,1), 1: NumPC)'*temp_trial;


end

for zz = 1:size(PCScore_trials, 1)
    PC_1_all_trials(zz, :) = PCScore_trials{zz, 1}(1, :);
    PC_2_all_trials(zz, :) = PCScore_trials{zz, 1}(2, :);
    PC_3_all_trials(zz, :) = PCScore_trials{zz, 1}(3, :);
    PC_4_all_trials(zz, :) = PCScore_trials{zz, 1}(4, :);
end


% Define the bin size
num_bins_for_pcs = 30;

% Initialize arrays to store the binned data
PC_1_binned = [];
PC_2_binned = [];
PC_3_binned = [];
PC_4_binned = [];

% Determine the number of bins
num_rows = size(PC_1_all_trials, 1);
num_bins = floor(num_rows / num_bins_for_pcs);

% Bin the data for each array
for i = 1:num_bins
    % Compute the row indices for the current bin
    row_start = (i-1) * num_bins_for_pcs + 1;
    row_end = row_start + num_bins_for_pcs - 1;
    
    % Calculate the mean across the rows for each principal component array
    PC_1_binned(i, :) = mean(PC_1_all_trials(row_start:row_end, :), 1);
    PC_2_binned(i, :) = mean(PC_2_all_trials(row_start:row_end, :), 1);
    PC_3_binned(i, :) = mean(PC_3_all_trials(row_start:row_end, :), 1);
    PC_4_binned(i, :) = mean(PC_4_all_trials(row_start:row_end, :), 1);
end

% Check that ts1 matches the number of columns in PC_1_binned
if length(ts1) ~= size(PC_1_binned, 2)
    error('The length of ts1 must match the number of columns in the binned data.');
end

% Define an array of PC data for easier looping
PC_binned_data = {PC_1_binned, PC_2_binned, PC_3_binned, PC_4_binned};
titles = {'PC 1', 'PC 2', 'PC 3', 'PC 4'};

% Loop through each PC binned data array
for pc_idx = 1:length(PC_binned_data)
    % Get the current PC binned data
    PC_binned = PC_binned_data{pc_idx};
    
    % Create a new figure
    figure;
    hold on;
    
    % Plot each bin
    for bin_idx = 1:size(PC_binned, 1)
        plot(ts1, PC_binned(bin_idx, :), 'DisplayName', sprintf('Bin #%d', bin_idx));
    end
    
    % Add legend, title, and labels
    legend('show', 'Location', 'best');
    title(titles{pc_idx});
    xlabel('Time (s)');
    ylabel('PC data');
    

    hold off;
end


%% plot for publication 3D trajectories


% % Create a figure and plot the initial state of the lines
% figure;
% % d1 = smoothforward(PCScore{1,1}(:,1:5:end), [1,size(PCScore{1,1},2);], 5, 15, 'mono_dir');
% % d2 = smoothforward(PCScore{1,2}(:,1:5:end), [1,size(PCScore{1,2},2);], 5, 15, 'mono_dir');
% % d3 = smoothforward(PCScore{1,3}(:,1:5:end), [1,size(PCScore{1,3},2);], 5, 15, 'mono_dir');
% 
% % d1 = PCScore{1,1};
% % d2 = PCScore{1,2};
% % d3 = PCScore{1,3};
% d1 = smoothforward(PCScore{1,1}, [1,size(PCScore{1,1},2);], 5, 15, 'mono_dir');
% d2 = smoothforward(PCScore{1,2}, [1,size(PCScore{1,2},2);], 5, 15, 'mono_dir');
% % d3 = smoothforward(PCScore{1,3}, [1,size(PCScore{1,3},2);], 5, 15, 'mono_dir');
% 
% downsampled_array_factor = numMeasurements/size(d1, 2);
% downsampled_start_time = size(d1, 2)/2;
% downsampled_median_choice_time_block_1 = downsampled_start_time+ (median_choice_time_block_1/downsampled_array_factor);
% downsampled_median_choice_time_block_2 = downsampled_start_time+ (median_choice_time_block_2/downsampled_array_factor);
% downsampled_median_choice_time_block_3 = downsampled_start_time+ (median_choice_time_block_3/downsampled_array_factor);
% 
% downsampled_median_collect_time_block_1 = downsampled_start_time+ (median_collect_time_block_1/downsampled_array_factor);
% downsampled_median_collect_time_block_2 = downsampled_start_time+ (median_collect_time_block_2/downsampled_array_factor);
% downsampled_median_collect_time_block_3 = downsampled_start_time+ (median_collect_time_block_3/downsampled_array_factor);
% 
% find(ts1 == ceil(median_start_time_block_1))
% 
% 
% d_marker_loc = ceil([downsampled_median_choice_time_block_1, downsampled_median_choice_time_block_2, downsampled_median_choice_time_block_3]);
% 
% 
% 
%     hold on;
%     p1 = plot(d1(1, :), d1(2, :), 'DisplayName', d_legend{1});
%     p1.Color(1: 3) = l_color{1}; p1.Color(4) = l_opacity; p1.LineWidth = l_width;
%     p1.Marker = '.'; p1.MarkerFaceColor = p_color{1}; p1.MarkerEdgeColor = p_color{1}; p1.MarkerIndices = [1: p_freq: size(d1, 2)]; p1.MarkerSize = p_size;
%     e11 = scatter(d1(1, closest_index_start_time_large), d1(2, closest_index_start_time_large), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
%     e12 = scatter(d1(1, ts1 == 0), d1(2, ts1 == 0), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
%     e13 = scatter(d1(1, closest_index_collect_time_large), d1(2, closest_index_collect_time_large), d_marker_size, [0, 0, 0]/255, '^', 'filled', 'HandleVisibility', 'off');
% %     e13 = scatter(d1(1, d_marker_loc(3)), d1(2, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
% 
%     p2 = plot(d2(1, :), d2(2, :), 'DisplayName', d_legend{2});
%     p2.Color(1: 3) = l_color{2}; p2.Color(4) = l_opacity; p2.LineWidth = l_width;
%     p2.Marker = '.'; p2.MarkerFaceColor = p_color{2}; p2.MarkerEdgeColor = p_color{2}; p2.MarkerIndices = [1: p_freq: size(d2, 2)]; p2.MarkerSize = p_size;
%     e14 = scatter(d2(1, closest_index_start_time_small), d2(2, closest_index_start_time_small), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
%     e15 = scatter(d2(1, ts1 == 0), d2(2, ts1 == 0), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
%     e16 = scatter(d2(1, closest_index_collect_time_small), d2(2, closest_index_collect_time_small), d_marker_size, [0, 0, 0]/255, '^', 'filled', 'HandleVisibility', 'off');
% %    
%     p3 = plot(d3(1, :), d3(2, :), 'DisplayName', d_legend{3});
%     p3.Color(1: 3) = l_color{3}; p3.Color(4) = l_opacity; p3.LineWidth = l_width;
%     p3.Marker = '.'; p3.MarkerFaceColor = p_color{3}; p3.MarkerEdgeColor = p_color{3}; p3.MarkerIndices = [1: p_freq: size(d3, 2)]; p3.MarkerSize = p_size;
%     e17 = scatter(d3(1, closest_index_start_time_block_3), d3(2, closest_index_start_time_block_3), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
%     e18 = scatter(d3(1, ts1 == 0), d3(2, ts1 == 0), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
%     e19 = scatter(d3(1, closest_index_collect_time_block_3), d3(2, closest_index_collect_time_block_3), d_marker_size, [0, 0, 0]/255, '^', 'filled', 'HandleVisibility', 'off');
% %    
%     hold off;
% 
% 
% 
% lines = {d1, d2, d3};
% 
% line_1 = d1; 
% line_2 = d2;
% line_3 = d3;
% 
% 
% % Parameters for animation
% numLines = size(lines, 2); % Total number of lines
% numFrames = size(d1, 2); % Total number of frames (assuming 40 data points)
% pauseTime = 0.1; % Pause time between frames (adjust as needed)
% 

d_legend = eventNames;
d_marker_loc = [1, 11, 21];
d_marker_size = 60;
d_marker_legend = {'Start', 'CS Onset', 'End'};
l_color = {[120, 114, 176]/255, [1, 1, 1]/255, [227, 124, 39]/255};
l_opacity = 0.6;
l_width = 5;
p_color = ["black",  "black",  "black"];
p_size = 5;
p_freq = 1;

line_color_space = round(linspace(1, size(acton, 1), numConditions));


figure;
for ff = 1:numConditions
    % PCA_traj{ff} = smoothforward(PCScore{1,ff}, [1,size(PCScore{1,ff},2);], 5, 15, 'mono_dir');
    PCA_traj{ff} = PCScore{1, ff};
    hold on;
    p1 = plot(PCA_traj{ff}(1, :), PCA_traj{ff}(2, :), 'DisplayName', d_legend{ff});
    p1.Color(1: 3) = acton(line_color_space(ff), :); 
    p1.Color(4) = l_opacity; 
    p1.LineWidth = l_width;
    p1.Marker = '.'; p1.MarkerFaceColor = p_color{1}; p1.MarkerEdgeColor = p_color{1}; p1.MarkerIndices = [1: p_freq: size(PCA_traj{ff}, 2)]; p1.MarkerSize = p_size;
    e11 = scatter(PCA_traj{ff}(1, closest_index_start(ff)), PCA_traj{ff}(2, closest_index_start(ff)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e12 = scatter(PCA_traj{ff}(1, closest_index_zero(ff)), PCA_traj{ff}(2, closest_index_zero(ff)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    e13 = scatter(PCA_traj{ff}(1, closest_index_collect(ff)), PCA_traj{ff}(2, closest_index_collect(ff)), d_marker_size, [0, 0, 0]/255, '^', 'filled', 'HandleVisibility', 'off');
%     e13 = scatter(d1(1, d_marker_loc(3)), d1(2, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    legend

end

%%

% Initialize a cell array to store the path lengths
path_lengths = cell(1, 2);

% Loop through each cell in PCA_traj
for i = 1:length(PCA_traj)
    % Get the current 4x80 double array
    data = PCA_traj{i};
    
    % Extract the X and Y coordinates (first 2 rows)
    X = data(1, :);
    Y = data(2, :);
    
    % Calculate the differences between consecutive points
    deltaX = diff(X);
    deltaY = diff(Y);
    
    % Calculate the Euclidean distances between consecutive points
    distances = sqrt(deltaX.^2 + deltaY.^2);
    
    % Sum the distances to get the total path length
    path_lengths{i} = sum(distances);
end

% Display the results
disp('Path lengths for each trajectory:');
disp(path_lengths);



%% based on ChatGPT to find the angle over the first 10 samples of the PCScore arrays

% Initialize an array to store the angles
angles = zeros(1, 2); % Since PCScore is a 1x2 cell array

for i = 1:2
    % Extract the first two rows (X and Y coordinates) for each cell
    X = PCScore{i}(1, :); % First row contains X coordinates
    Y = PCScore{i}(2, :); % Second row contains Y coordinates
    
    % Compute the vector between the first index and the 10th index
    deltaX = X(10) - X(1);
    deltaY = Y(10) - Y(1);
    
    % Calculate the angle in radians using the arctangent
    angle = atan2(deltaY, deltaX); % atan2 gives the angle relative to X-axis
    
    % Convert the angle to degrees for easier interpretation
    angles(i) = rad2deg(angle);
end

% Display the calculated angles for each trajectory
disp(angles);

%%
% Extract the first two PCs for each dataset
PC1_set1 = PCScore{1,1}(1:2, :); % 2x160 for dataset 1
PC1_set2 = PCScore{1,2}(1:2, :); % 2x160 for dataset 2



% Calculate differences between consecutive samples
diff_set1 = diff(PC1_set1, 1, 2); % Take differences along columns (time)

% Compute Euclidean distances
distances_set1 = sqrt(sum(diff_set1.^2, 1)); % Sum squares along rows (PCs) and take square root

% Sum the distances to get the total length
trajectory_length_set1 = sum(distances_set1);


% Same calculation for dataset 2
diff_set2 = diff(PC1_set2, 1, 2); % Differences along columns
distances_set2 = sqrt(sum(diff_set2.^2, 1)); % Euclidean distances
trajectory_length_set2 = sum(distances_set2);


% Calculate the difference between corresponding points
diff_trajectories = PC1_set1 - PC1_set2; % Element-wise difference between trajectories

% Compute the Euclidean distance for each bin
euclidean_distances = sqrt(sum(diff_trajectories.^2, 1)); % Sum squares across PCs, then take square root

figure;
plot(ts1, euclidean_distances, 'LineWidth', 2);
xlabel('Bin (Time Point)');
ylabel('Euclidean Distance');
title('Distance Between PCA Trajectories (Bin-by-Bin)');
grid on;

diff_trajectories_PC1_only = PC1_set1(1, :) - PC1_set2(1, :); % Element-wise difference between trajectories
% Compute the Euclidean distance for each bin
euclidean_distances = sqrt(sum(diff_trajectories_PC1_only.^2, 1)); % Sum squares across PCs, then take square root

figure;
plot(ts1, euclidean_distances, 'LineWidth', 2);
xlabel('Bin (Time Point)');
ylabel('Euclidean Distance');
title('Distance Between PCA Trajectories (Bin-by-Bin)');
grid on;

figure; plot(ts1, PC1_set1(1, :))
hold on; plot(ts1, PC1_set2(1, :))


diff_trajectories_PC2_only = PC1_set1(2, :) - PC1_set2(2, :); % Element-wise difference between trajectories
% Compute the Euclidean distance for each bin
euclidean_distances = sqrt(sum(diff_trajectories_PC2_only.^2, 1)); % Sum squares across PCs, then take square root

figure;
plot(ts1, euclidean_distances, 'LineWidth', 2);
xlabel('Bin (Time Point)');
ylabel('Euclidean Distance');
title('Distance Between PCA Trajectories (Bin-by-Bin)');
grid on;

figure; plot(ts1, PC1_set1(2, :))
hold on; plot(ts1, PC1_set2(2, :))

%% Hao code data setup

[uniqueStrings, ~, indices] = unique(mouse_cells(1, :));

disp(indices);

% indices = each cell labeled with the mouse index
% NumPC = number of PCs from above
% neuron_num = total # of neurons in the sample
% numMeasurements = # of samples in window

data_LOO_first_variable = pca2LOO_new(coef, neuron_mean_concat(:, condition_ranges{1, 1}{1}:condition_ranges{1, 1}{2}), indices, NumPC, neuron_num, numMeasurements);

data_LOO_second_variable = pca2LOO_new(coef, neuron_mean_concat(:, condition_ranges{1, 2}{1}:condition_ranges{1, 2}{2}), indices, NumPC, neuron_num, numMeasurements);


%% leave one out for length and distance
% match the number of neurons
e = ephys.zdata(ephys.cell_labels(:,1)==1,:);
c = ephys.zdata(ephys.cell_labels(:,1)==0,:);
idx = randperm(size(c,1),size(e,1));
idx = sort(idx);
c = c(idx,:);
%find neutral tirals that are 0   103:180 just arbitrary
nc = sum(c(:,105:180),2)==0;
ne = sum(e(:,105:180),2)==0;

% finding out which neurons coming from which mice
einfo = ephys.info(ephys.cell_labels(:,1)==1,1);
cinfo = ephys.info(ephys.cell_labels(:,1)==0,:);
cinfo = cinfo(idx,:);
file = [cinfo(~nc,1) ; einfo(~ne,1)];

% data = [c;e];
% file = [cinfo(:,1) ; einfo(:,1)];

file = string(file);
f = 0;
for i = 2: size(file,1)
    f(i) = strcmp(file(i),file(i-1));
end
f = 1-f';
f_c = find(f(1:col)==1);
f_e = find(f(col+1:end)==1);

if length(f_c)<length(f_e)
    K = length(f_c);
else
    K = length(f_e);
end
f_final = find(f==1);

f_c(length(f_c)+1) = size(c(~nc,:),1);
f_e(length(f_e)+1) = size(e(~ne,:),1);

cellID_c = [];
for i = 1:length(f_c)-1
    cellID_c(f_c(i):f_c(i+1),1) = i;
end


cellID_e = [];
for i = 1:length(f_e)-1
    cellID_e(f_e(i):f_e(i+1),1) = i;
end

% leave one out
num_pcs = find(cumsum(explained)>=90, 1);
data_LOO_c = pca2LOO_new(coef(1:size(c(~nc,:),1),:), c(~nc,1:binSize*3), cellID_c, num_pcs, min(diff(f_c)), floor(median(diff(f_c))));
data_LOO_e = pca2LOO_new(coef(size(c(~nc,:),1)+1:end,:), e(~ne,1:binSize*3), cellID_e, num_pcs, min(diff(f_e)), floor(median(diff(f_e))));

%calculate length and distance using vecnorm
l_c_reward = zeros(1, size(data_LOO_first_variable, 3));
l_c_neutral = zeros(1, size(data_LOO_c, 3));
l_c_shock = zeros(1, size(data_LOO_c, 3));
l_e_reward = zeros(1, size(data_LOO_second_variable, 3));
l_e_neutral = zeros(1, size(data_LOO_e, 3));
l_e_shock = zeros(1, size(data_LOO_e, 3));


for i = 1: size(data_LOO_first_variable, 3)
    for j = 1: size(data_LOO_first_variable, 2)
        l_c_reward(1, i) = l_c_reward(1, i) + vecnorm(data_LOO_first_variable(:,j,i) - data_LOO_first_variable(:, j, i), 2, 1);

    end
end

for i = 1: size(data_LOO_second_variable, 3)
    for j = 1: size(data_LOO_second_variable, 2)
        l_e_reward(1, i) = l_e_reward(1, i) + vecnorm(data_LOO_second_variable(:,j,i) - data_LOO_second_variable(:, j, i), 2, 1);

    end
end

% distance
d_c_reward = zeros(101, nchoosek(size(data_LOO_c, 3), 2));
counter = 1;
for i = 1: size(data_LOO_c, 3)
    for j = i+1: size(data_LOO_c, 3)
        for k = 1: 101
            d_c_reward(k, counter) = vecnorm(data_LOO_c(:, k, i) - data_LOO_c(:, k, j));
        end
        counter = counter + 1;
    end
end

%% BRITT EucDistance calculation

eucD1to7= [];
%here, p = s1 and q = s7
for i = 1:size(PCScore{1, 1}  ,2)
    p1 = PCScore{1, 1}(1,i);
    p2 = PCScore{1, 1}(2,i);
    p3 = PCScore{1, 1}(3,i);

    q1 = PCScore{1, 2}(1,i);
    q2 = PCScore{1, 2}(2,i);
    q3 = PCScore{1, 2}(3,i); 

    eucDtemp = sqrt((q1-p1).^2 + (q2-p2).^2 + (q3-p3).^2);
    % eucDtemp = sqrt((q1-p1).^2 + (q2-p2).^2);
eucD1to7(1,i) = eucDtemp;
end 

figure; plot(ts1, eucD1to7)

%%
% Data
s1_data = PCScore{1,1}; % Data for s1
s2_data = PCScore{1,2}; % Data for s2

% Time points (adjust based on your sample rate)
time_points = linspace(-4, 2, size(s1_data, 2));

% Choice and collect times
choice_time = 0; % Choice happens at time 0
collect_time_s1 = 0.8790; % Collection time for s1
collect_time_s2 = 1.1940; % Collection time for s2

% Create the figure
figure;
hold on;

% Plot s1 trajectory
s1_pre_choice = time_points < choice_time;
plot3(s1_data(1, s1_pre_choice), s1_data(2, s1_pre_choice), s1_data(3, s1_pre_choice), 'b--'); % Dotted line pre-choice
plot3(s1_data(1, ~s1_pre_choice), s1_data(2, ~s1_pre_choice), s1_data(3, ~s1_pre_choice), 'b-', 'LineWidth', 1.5); % Solid line post-choice
scatter3(s1_data(1, time_points == collect_time_s1), ...
         s1_data(2, time_points == collect_time_s1), ...
         s1_data(3, time_points == collect_time_s1), ...
         100, 's', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); % Collect time marker

% Plot s2 trajectory
s2_pre_choice = time_points < choice_time;
plot3(s2_data(1, s2_pre_choice), s2_data(2, s2_pre_choice), s2_data(3, s2_pre_choice), 'r--'); % Dotted line pre-choice
plot3(s2_data(1, ~s2_pre_choice), s2_data(2, ~s2_pre_choice), s2_data(3, ~s2_pre_choice), 'r-', 'LineWidth', 1.5); % Solid line post-choice
scatter3(s2_data(1, time_points == collect_time_s2), ...
         s2_data(2, time_points == collect_time_s2), ...
         s2_data(3, time_points == collect_time_s2), ...
         100, 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); % Collect time marker

% Configure the 3D plot
grid on;
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('3D Trajectories of PC Scores');
legend({'s1 Pre-Choice', 's1 Post-Choice', 's1 Collect', ...
        's2 Pre-Choice', 's2 Post-Choice', 's2 Collect'}, ...
        'Location', 'best');
view(3); % Adjust view angle for better visualization

%%
% Data
s1_data = PCScore{1,1}; % Data for s1
s2_data = PCScore{1,2}; % Data for s2

% Time points (adjust based on your sample rate)
time_points = linspace(-4, 2, size(s1_data, 2));

% Choice and collect times
choice_time = 0; % Choice happens at time 0
collect_time_s1 = 0.8790; % Collection time for s1
collect_time_s2 = 1.1940; % Collection time for s2

collect_idx_s1 = find(abs(time_points - collect_time_s1) == min(abs(time_points - collect_time_s1)), 1);
collect_idx_s2 = find(abs(time_points - collect_time_s2) == min(abs(time_points - collect_time_s2)), 1);


% Create the figure
figure;
hold on;
grid on;
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('Animated 3D Trajectories of PC Scores');
view(3); % 3D perspective

% Initialize plot handles
s1_line_pre = plot3(nan, nan, nan, 'b--', 'LineWidth', 1.5); % Dotted pre-choice line for s1
s1_line_post = plot3(nan, nan, nan, 'b-', 'LineWidth', 1.5); % Solid post-choice line for s1
% s1_collect_marker = scatter3(nan, nan, nan, 100, 's', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); % s1 collect marker

s2_line_pre = plot3(nan, nan, nan, 'r--', 'LineWidth', 1.5); % Dotted pre-choice line for s2
s2_line_post = plot3(nan, nan, nan, 'r-', 'LineWidth', 1.5); % Solid post-choice line for s2
% s2_collect_marker = scatter3(nan, nan, nan, 100, 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); % s2 collect marker

for frame = 1:size(s1_data, 2)
    % Determine current frame
    current_frame = min(frame, size(s1_data, 2));
    
    % Update s1 trajectory
    pre_choice_frame = time_points < choice_time;

    % Update s1 trajectory
    valid_indices_pre_s1 = find(pre_choice_frame & (1:length(time_points) <= current_frame));
    set(s1_line_pre, 'XData', s1_data(1, valid_indices_pre_s1), ...
        'YData', s1_data(2, valid_indices_pre_s1), ...
        'ZData', s1_data(3, valid_indices_pre_s1));

    valid_indices_post_s1 = find(~pre_choice_frame & (1:length(time_points) <= current_frame));
    set(s1_line_post, 'XData', s1_data(1, valid_indices_post_s1), ...
        'YData', s1_data(2, valid_indices_post_s1), ...
        'ZData', s1_data(3, valid_indices_post_s1));

    set(s1_collect_marker, 'XData', s1_data(1, collect_idx_s1), ...
        'YData', s1_data(2, collect_idx_s1), ...
        'ZData', s1_data(3, collect_idx_s1));

    % Update s2 trajectory
    valid_indices_pre_s2 = find(pre_choice_frame & (1:length(time_points) <= current_frame));
    set(s2_line_pre, 'XData', s2_data(1, valid_indices_pre_s2), ...
        'YData', s2_data(2, valid_indices_pre_s2), ...
        'ZData', s2_data(3, valid_indices_pre_s2));

    valid_indices_post_s2 = find(~pre_choice_frame & (1:length(time_points) <= current_frame));
    set(s2_line_post, 'XData', s2_data(1, valid_indices_post_s2), ...
        'YData', s2_data(2, valid_indices_post_s2), ...
        'ZData', s2_data(3, valid_indices_post_s2));

    set(s2_collect_marker, 'XData', s2_data(1, collect_idx_s2), ...
        'YData', s2_data(2, collect_idx_s2), ...
        'ZData', s2_data(3, collect_idx_s2));

    % Update s1 collect marker
    if current_frame >= collect_idx_s1
        scatter3(s1_data(1, collect_idx_s1), s1_data(2, collect_idx_s1), s1_data(3, collect_idx_s1), ...
            100, 'filled', 'MarkerFaceColor', 'b');
    else
        set(s1_collect_marker, ...
            'XData', nan, 'YData', nan, 'ZData', nan, 'Visible', 'off'); % Hide marker
    end

    % Update s2 collect marker
    if current_frame >= collect_idx_s2
        scatter3(s2_data(1, collect_idx_s2), s2_data(2, collect_idx_s2), s2_data(3, collect_idx_s2), ...
            100, 'filled', 'MarkerFaceColor', 'r');
    else
        set(s2_collect_marker, ...
            'XData', nan, 'YData', nan, 'ZData', nan, 'Visible', 'off'); % Hide marker
    end

    % Refresh the plot
   

    % Update axis limits dynamically
    xlim([min([s1_data(1,:), s2_data(1,:)]) max([s1_data(1,:), s2_data(1,:)])]);
    ylim([min([s1_data(2,:), s2_data(2,:)]) max([s1_data(2,:), s2_data(2,:)])]);
    zlim([min([s1_data(3,:), s2_data(3,:)]) max([s1_data(3,:), s2_data(3,:)])]);
    
    % Pause for animation effect
    pause(0.1);
end


% Optional: Save animation to video
% video = VideoWriter('3D_Trajectory_Animation.mp4', 'MPEG-4');
% open(video);
% writeVideo(video, frames);
% close(video);
