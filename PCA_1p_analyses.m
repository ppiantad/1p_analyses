

% % use zall array if you want to check how trials compare across block
% neuron_mean_concat = horzcat(zall_mean_all_array{:});

% neuron_mean_concat = horzcat(neuron_mean_array{:});

% % use neuron_mean_all_unnorm if you want to check how things differ across
% % time
% neuron_mean_concat = horzcat(neuron_mean_all_unnormalized{:});

% neuron_mean_concat = zscore(neuron_mean_concat, 0 , 2);

% % use below to normalize THEN concatenate, which might be a better
% approach?
% Iterate through each cell and apply zscore
for i = 1:size(neuron_mean_all_unnormalized, 2)
    % Apply zscore to the data in each cell
    neuron_mean_all_normalized{i} = zscore(neuron_mean_all_unnormalized{i}, 0, 2);
    neuron_mean_reformat{i} = neuron_mean_all_unnormalized{i}';
end
% % 
neuron_mean_concat = horzcat(neuron_mean_all_normalized{:});
neuron_mean_concat = zscore(neuron_mean_concat, 0 , 2);

neuron_mean_reformat_concat = vertcat(neuron_mean_reformat{:});
neuron_mean_reformat_concat_normalized = zscore(neuron_mean_reformat_concat, 0 , 1);


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
    median_collect_time(1, zz) = median(concatenatedTable_all{1, zz}.collectionTime - concatenatedTable_all{1, zz}.stTime)
end

%% PCA
% Load your data if not already loaded
% load('neuron_mean_concat.mat');

numNeuronsPerCondition = neuron_num;
% change depending on the number of behaviors to decode!
% numConditions = size(neuron_mean_concat, 1)/numNeuronsPerCondition;

% numConditions = size(neuron_mean_concat, 2)/numNeuronsPerCondition;

numConditions = size(neuron_mean_concat, 2)/numMeasurements;

% numConditions = numNeuronsPerCondition/numMeasurements;

eventIdx = 1:numConditions;

NumPC = 4; %2

array_size = size(neuron_mean_concat, 2);



% Initialize an empty cell array for the result
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

% Display the result
disp(result)

% Generate the ranges & then use this to loop through and conduct PCA STILL
% NEED TO FINISH UPDATING CODE BELOW LINE 32
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


[coef, ~, ~, ~, explained, ~] = pca(neuron_mean_concat');

for i = eventIdx
    temp = condition_data{i};

    PCScore{i} = coef(1:size(temp,1), 1: NumPC)'*temp;

end

% Assuming explained is a double array with variance explained for each principal component
% explained = [var1, var2, var3, ...]; % Replace this with your actual explained array

% Create the scree plot
figure;
bar(explained, 'FaceColor', [0.2 0.6 0.8]); % Creates a bar plot with custom color
title('Scree Plot');
xlabel('Principal Component');
ylabel('Variance Explained (%)');
grid on; % Optional: adds a grid to the plot


%% plot for publication 3D trajectories
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

figure ();
plot3traj_3Dgaussian_hao(PCScore{1,1},PCScore{1,2},PCScore{1,3}, d_legend, d_marker_loc, d_marker_size, d_marker_legend, l_color, l_opacity, l_width, p_color, p_size, p_freq)
% xlim([-5 10])
% ylim([-1 2])
hold on
%d_marker_size = 20;
%p_size = 1;
% plot3traj_3Dgaussian_hao(PCScore{2,1},PCScore{2,2},PCScore{2,3}, d_legend, d_marker_loc, d_marker_size, d_marker_legend, l_color, l_opacity, l_width, p_color, p_size, p_freq)
%  xlim([0 1.5])
%  ylim([0 2])


%%
% Create a figure and plot the initial state of the lines
figure;
% d1 = smoothforward(PCScore{1,1}(:,1:5:end), [1,size(PCScore{1,1},2);], 5, 15, 'mono_dir');
% d2 = smoothforward(PCScore{1,2}(:,1:5:end), [1,size(PCScore{1,2},2);], 5, 15, 'mono_dir');
% d3 = smoothforward(PCScore{1,3}(:,1:5:end), [1,size(PCScore{1,3},2);], 5, 15, 'mono_dir');

% d1 = PCScore{1,1};
% d2 = PCScore{1,2};
% d3 = PCScore{1,3};
d1 = smoothforward(PCScore{1,1}, [1,size(PCScore{1,1},2);], 5, 15, 'mono_dir');
d2 = smoothforward(PCScore{1,2}, [1,size(PCScore{1,2},2);], 5, 15, 'mono_dir');
% d3 = smoothforward(PCScore{1,3}, [1,size(PCScore{1,3},2);], 5, 15, 'mono_dir');

downsampled_array_factor = numMeasurements/size(d1, 2);
downsampled_start_time = size(d1, 2)/2;
downsampled_median_choice_time_block_1 = downsampled_start_time+ (median_choice_time_block_1/downsampled_array_factor);
downsampled_median_choice_time_block_2 = downsampled_start_time+ (median_choice_time_block_2/downsampled_array_factor);
downsampled_median_choice_time_block_3 = downsampled_start_time+ (median_choice_time_block_3/downsampled_array_factor);

downsampled_median_collect_time_block_1 = downsampled_start_time+ (median_collect_time_block_1/downsampled_array_factor);
downsampled_median_collect_time_block_2 = downsampled_start_time+ (median_collect_time_block_2/downsampled_array_factor);
downsampled_median_collect_time_block_3 = downsampled_start_time+ (median_collect_time_block_3/downsampled_array_factor);

find(ts1 == ceil(median_start_time_block_1))


d_marker_loc = ceil([downsampled_median_choice_time_block_1, downsampled_median_choice_time_block_2, downsampled_median_choice_time_block_3]);



    hold on;
    p1 = plot(d1(1, :), d1(2, :), 'DisplayName', d_legend{1});
    p1.Color(1: 3) = l_color{1}; p1.Color(4) = l_opacity; p1.LineWidth = l_width;
    p1.Marker = '.'; p1.MarkerFaceColor = p_color{1}; p1.MarkerEdgeColor = p_color{1}; p1.MarkerIndices = [1: p_freq: size(d1, 2)]; p1.MarkerSize = p_size;
    e11 = scatter(d1(1, closest_index_start_time_large), d1(2, closest_index_start_time_large), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e12 = scatter(d1(1, ts1 == 0), d1(2, ts1 == 0), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    e13 = scatter(d1(1, closest_index_collect_time_large), d1(2, closest_index_collect_time_large), d_marker_size, [0, 0, 0]/255, '^', 'filled', 'HandleVisibility', 'off');
%     e13 = scatter(d1(1, d_marker_loc(3)), d1(2, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    
    p2 = plot(d2(1, :), d2(2, :), 'DisplayName', d_legend{2});
    p2.Color(1: 3) = l_color{2}; p2.Color(4) = l_opacity; p2.LineWidth = l_width;
    p2.Marker = '.'; p2.MarkerFaceColor = p_color{2}; p2.MarkerEdgeColor = p_color{2}; p2.MarkerIndices = [1: p_freq: size(d2, 2)]; p2.MarkerSize = p_size;
    e14 = scatter(d2(1, closest_index_start_time_small), d2(2, closest_index_start_time_small), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e15 = scatter(d2(1, ts1 == 0), d2(2, ts1 == 0), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    e16 = scatter(d2(1, closest_index_collect_time_small), d2(2, closest_index_collect_time_small), d_marker_size, [0, 0, 0]/255, '^', 'filled', 'HandleVisibility', 'off');
%    
    p3 = plot(d3(1, :), d3(2, :), 'DisplayName', d_legend{3});
    p3.Color(1: 3) = l_color{3}; p3.Color(4) = l_opacity; p3.LineWidth = l_width;
    p3.Marker = '.'; p3.MarkerFaceColor = p_color{3}; p3.MarkerEdgeColor = p_color{3}; p3.MarkerIndices = [1: p_freq: size(d3, 2)]; p3.MarkerSize = p_size;
    e17 = scatter(d3(1, closest_index_start_time_block_3), d3(2, closest_index_start_time_block_3), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e18 = scatter(d3(1, ts1 == 0), d3(2, ts1 == 0), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    e19 = scatter(d3(1, closest_index_collect_time_block_3), d3(2, closest_index_collect_time_block_3), d_marker_size, [0, 0, 0]/255, '^', 'filled', 'HandleVisibility', 'off');
%    
    hold off;



lines = {d1, d2, d3};

line_1 = d1; 
line_2 = d2;
line_3 = d3;


% Parameters for animation
numLines = size(lines, 2); % Total number of lines
numFrames = size(d1, 2); % Total number of frames (assuming 40 data points)
pauseTime = 0.1; % Pause time between frames (adjust as needed)





%% animate PCA lines

% Calculate the minimum and maximum x and y values from all lines
minX = min(cellfun(@(line) min(line(1, :)), lines));
maxX = max(cellfun(@(line) max(line(1, :)), lines));
minY = min(cellfun(@(line) min(line(2, :)), lines));
maxY = max(cellfun(@(line) max(line(2, :)), lines));


% Create a figure and axes
figure;
ax = gca;
xlabel('X');
ylabel('Y');
grid on;

% Set fixed dimensions for the plot based on the minimum and maximum values
axis([floor(minX) ceil(maxX) floor(minY) ceil(maxY)]);

% Parameters for animation
pauseTime = 0.1; % Pause time between animations (adjust as needed)


% Initialize empty lines
lineObj = cell(1, numLines);
for i = 1:numLines
    lineObj{i} = line('XData', [], 'YData', [], 'LineWidth', 2);
end

% Iterate over lines
for i = 1:numLines
    % Iterate over frames for the current line
    for frame = 1:numFrames
        % Update the current line's data
        set(lineObj{i}, 'XData', lines{1, i}(1, 1:frame), 'YData', lines{1, i}(2, 1:frame));
        
        % Update the figure
        drawnow;
        
        % Pause to control the animation speed
        pause(pauseTime);
    end
end



