animalIDs = fieldnames(final_DLC);

session_to_analyze = 'D3';

recorded_fps = 30;



% Parameters
stimulus_duration = 2 * 60; % 2 minutes in seconds
num_repeats = 3;
total_stimuli = 2;

% Initialize variables
stimulus_times = cell(total_stimuli, 1);
current_time = 0;

% Loop through each stimulus alternately and calculate start and end times
for j = 1:num_repeats
    for i = 1:total_stimuli
        start_time = current_time;
        end_time = start_time + stimulus_duration;
        
        if isempty(stimulus_times{i})
            stimulus_times{i} = [start_time, end_time];
        else
            stimulus_times{i} = [stimulus_times{i}; start_time, end_time];
        end
        
        current_time = end_time;
    end
end


% Multiply every value in stimulus_times by the FPS
for i = 1:total_stimuli
    stimulus_frames{i} = stimulus_times{i} * recorded_fps;
end



%%
clear aversive_context safe_context mean_aversive mean_safe mean_aversive_context mean_safe_context


mouse_count = 0;
for gg = 1:size(animalIDs, 1)
    current_mouse = animalIDs{gg};
    
    if strcmp(final_DLC.(current_mouse).experimental_grp, 'No Shock')
        mouse_count = mouse_count+1;
        mouse_in_cond(mouse_count, :) = current_mouse;
        mouse_data = final_DLC.(current_mouse).(session_to_analyze);
        fields = fieldnames(mouse_data); % Get the field names
        matchIdx = find(contains(fields, 'DLC')); % Find indices of matches
        


        DLC_data_mouse = final_DLC.(current_mouse).(session_to_analyze).movement_data;
        if any(strcmp('freeze', DLC_data_mouse.Properties.VariableNames))
            freeze_data = DLC_data_mouse.freeze; 
        elseif any(strcmp('freeze_status', DLC_data_mouse.Properties.VariableNames))
            freeze_data = DLC_data_mouse.freeze_status;
        elseif any(strcmp('was_freezing', DLC_data_mouse.Properties.VariableNames))
            freeze_data = DLC_data_mouse.was_freezing;
        end
        % freeze_data = final_DLC.(current_mouse).(session_to_analyze).DLC_data_raw.freeze; 
        frames_data = DLC_data_mouse.frame;
        
        for qq = 1:num_repeats
            safe_context{mouse_count, qq} = freeze_data(stimulus_frames{1, 1}(qq,1)+1:stimulus_frames{1, 1}(qq,2));
        end
        
        
        for pp = 1:size(safe_context, 2)
            % mean_safe_context(mouse_count, pp) = sum(freeze_data(stimulus_frames{1, 1}(pp,1)+1:stimulus_frames{1, 1}(pp,2)))/size(safe_context, 1);
            % Calculate the proportion of 1s in the column
            b = sum(safe_context{mouse_count,pp}) / size(safe_context{mouse_count,pp}, 1);

            % Store the proportion in mean_aversive_context
            mean_safe_context(mouse_count, pp) = b;

            % Calculate the number of rows in the column
            n = size(safe_context, 1);

            % Calculate the standard error of the proportion
            standard_error_safe(mouse_count, pp) = sqrt(b * (1 - b) / n);
        end


        for qq = 1:num_repeats
            if qq <= 2
                aversive_context{mouse_count, qq} = freeze_data(stimulus_frames{1, 2}(qq,1)+1:stimulus_frames{1, 2}(qq,2));
            elseif qq > 2
                aversive_context{mouse_count, qq} = freeze_data(stimulus_frames{1, 2}(qq,1)+1:end);
            end
        end


        for pp = 1:size(aversive_context, 2)
            % Calculate the proportion of 1s in the column
            p = sum(aversive_context{mouse_count,pp}) / size(aversive_context{mouse_count,pp}, 1);

            % Store the proportion in mean_aversive_context
            mean_aversive_context(mouse_count, pp) = p;

            % Calculate the number of rows in the column
            n = size(mean_aversive_context, 1);

            % Calculate the standard error of the proportion
            standard_error_aversive(mouse_count, pp) = sqrt(p * (1 - p) / n);
        end

    end
    



end


% Calculate the mean of each column
mean_safe = mean(mean_safe_context);
mean_aversive = mean(mean_aversive_context);

mean_sem_safe = mean(standard_error_safe);
mean_sem_aversive = mean(standard_error_aversive);

% Initialize the interleaved array
interleaved_means = zeros(1, 6);
interleaved_sems = zeros(1,6);
interleaved_raw = zeros(size(mean_aversive_context, 1),6);

% Interleave the means
interleaved_means(1:2:end) = mean_safe;
interleaved_means(2:2:end) = mean_aversive;

% Interleave the SEMs
interleaved_sems(1:2:end) = mean_sem_safe;
interleaved_sems(2:2:end) = mean_sem_aversive;

interleaved_raw(1:2:end) = mean_safe_context;
interleaved_raw(2:2:end) = mean_aversive_context;

% Display the interleaved means
disp('Interleaved Means:');
disp(interleaved_means);

% figure; plot(interleaved_means);

test_interleave_mean = zeros(size(mean_safe_context, 1), 6);
test_interleave_mean(:, [1 3 5]) = mean_safe_context;
test_interleave_mean(:, [2 4 6]) = mean_aversive_context;

test_interleave_sem = zeros(size(mean_safe_context, 1), 6);
test_interleave_sem(:, [1 3 5]) = standard_error_safe;
test_interleave_sem(:, [2 4 6]) = standard_error_aversive;

figure; shadedErrorBar(1:size(mean(test_interleave_mean), 2), mean(test_interleave_mean), mean(test_interleave_sem));
% hold on; plot(1:size(interleaved_means, 2), interleaved_raw)
% hold on;shadedErrorBar(1:size(interleaved_means, 2), interleaved_means, interleaved_sems);
hold off; 


% Initialize a new cell array to store the concatenated results
combined_context = cell(size(aversive_context, 1), 1);

for i = 1:size(aversive_context, 1)
    concatenated_array = [];
    
    for j = 1:size(aversive_context, 2)
        % Concatenate the data from the corresponding cells of aversive_context and safe_context
        concatenated_array = [concatenated_array; aversive_context{i, j}; safe_context{i, j}];
    end
    
    % Store the concatenated array in the new cell array
    combined_context{i} = concatenated_array;
    max_size(i) = size(combined_context{i}, 1);
end

trimmed_combined_context = []; 

for zz = 1:size(combined_context, 1)
    trimmed_combined_context(zz, :) = combined_context{zz, 1}(1:min(max_size), 1); 



end

session_long_mean = mean(trimmed_combined_context); 



%%
% Plotting the data
figure;
plot(final_DLC.C68604.D3.movement_data.frame, final_DLC.C68604.D3.movement_data.was_freezing);
hold on;

% Define colors and transparency
color_1 = [0, 0, 1]; % Blue for stimulus 1
color_2 = [1, 0, 0]; % Red for stimulus 2
transparency = 0.3; % Transparency level

% Plot rectangles for each stimulus presentation
for i = 1:total_stimuli
    for j = 1:num_repeats
        % Get start and end frames
        start_frame = stimulus_frames{i}(j, 1);
        end_frame = stimulus_frames{i}(j, 2);
        
        % Determine the Y-axis range
        y_limits = ylim;
        
        % Plot the rectangle
        rectangle('Position', [start_frame, y_limits(1), end_frame-start_frame, y_limits(2)-y_limits(1)], ...
                  'FaceColor', [color_1 transparency] * (i == 1) + [color_2 transparency] * (i == 2), ...
                  'EdgeColor', 'none');
    end
end

% Finalize the plot

hold off;

%% conditioning
% experimental_grps = readtable('E:\MATLAB\my_repo\context fear\organize_SLEAP_data\full_pilot_mice.xlsx');
experimental_grps = readtable('I:\MATLAB\my_repo\context fear\organize_DLC_data\pilot groups.xlsx');

% Define parameters
threshold = 1; % Velocity threshold
sample_duration = 0.03; % Duration of each sample in seconds
min_duration = 2; % Minimum duration to trigger labeling in seconds

% Calculate the minimum number of consecutive rows needed
min_samples = min_duration / sample_duration;

animalIDs = fieldnames(final_DLC);

session_to_analyze = 'D2_Afternoon';

group_to_analyze = 'No Shock';

mouse_count = 0;
for gg = 1:size(animalIDs, 1)
    current_mouse = animalIDs{gg};
    
    if isfield(final_DLC.(current_mouse), session_to_analyze)
        mouse_count = mouse_count+1;
        DLC_data_mouse = final_DLC.(current_mouse).(session_to_analyze).movement_data;
        
        body_velocity = [];
        labels = [];
        % Get the body_velocity column
        body_velocity = final_DLC.(current_mouse).(session_to_analyze).movement_data.body_velocity;

        % Initialize the new column
        labels = zeros(size(body_velocity));

        % Find consecutive segments where body_velocity < threshold
        below_threshold = body_velocity < threshold;
        start_idx = find(diff([0; below_threshold]) == 1); % Start indices
        end_idx = find(diff([below_threshold; 0]) == -1); % End indices

        % Iterate through each segment and label
        for i = 1:length(start_idx)
            segment_length = end_idx(i) - start_idx(i) + 1;
            if segment_length >= min_samples
                labels(start_idx(i):end_idx(i)) = 1;
            end
        end

        % Add the labels as a new column to the table
        % final_DLC.B46837.D1_Afternoon.movement_data.freeze_label = labels;
        freeze_data(mouse_count, :) = labels(1:21590)';
        experimental_grps_updated(mouse_count, :) = experimental_grps(gg, :);
    end
end

% samples = final_DLC.(current_mouse).(session_to_analyze).movement_data.frame(1:21590)';
% mean_freeze_experimental = mean(freeze_data(strcmp(experimental_grps.group, group_to_analyze), :));
% std_freeze_experimental = std(freeze_data(strcmp(experimental_grps.group, group_to_analyze), :));
% sem_freeze_experimental = std_freeze_experimental/sqrt(size(freeze_data(strcmp(experimental_grps.group, group_to_analyze)), 1));
% 
% 
% mean_freeze_experimental_percent = mean_freeze_experimental*100;
% figure; plot(mean_freeze_experimental_percent);

%%
% Number of bins
num_bins = 48;

% Original number of columns
num_columns = size(freeze_data, 2);

% Bin size (number of columns per bin)
bin_size = floor(num_columns / num_bins);

% Preallocate binned data array
binned_data = zeros(size(freeze_data, 1), num_bins);

% Loop through each bin and calculate the mean for each mouse
for bin_idx = 1:num_bins
    % Determine the start and end columns for the current bin
    start_col = (bin_idx - 1) * bin_size + 1;
    if bin_idx == num_bins
        % Ensure the last bin includes any remaining columns
        end_col = num_columns;
    else
        end_col = bin_idx * bin_size;
    end
    
    % Average data within the current bin
    binned_data(:, bin_idx) = mean(freeze_data(:, start_col:end_col), 2);
end

binned_data_std = std(binned_data(strcmp(experimental_grps_updated.group, group_to_analyze), :));
binned_data_sem = binned_data_std/sqrt(size(binned_data(strcmp(experimental_grps_updated.group, group_to_analyze)), 1));

% %%
% figure; plot(mean(binned_data(strcmp(experimental_grps.group, 'Experimental'), :))); 
% hold on; plot(mean(binned_data(strcmp(experimental_grps.group, 'One Context'), :)));
% hold on; plot(mean(binned_data(strcmp(experimental_grps.group, 'No Shock'), :)));

%%
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(1:48, mean(binned_data(strcmp(experimental_grps_updated.group, 'Experimental'), :)), binned_data_std/sqrt(size(binned_data(strcmp(experimental_grps_updated.group, 'Experimental')), 1)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(1:48, mean(binned_data(strcmp(experimental_grps_updated.group, 'One Context'), :)), binned_data_std/sqrt(size(binned_data(strcmp(experimental_grps_updated.group, 'One Context')), 1)), 'lineProps', {'color', 'k'});
h(2) = shadedErrorBar(1:48, mean(binned_data(strcmp(experimental_grps_updated.group, 'No Shock'), :)), binned_data_std/sqrt(size(binned_data(strcmp(experimental_grps_updated.group, 'No Shock')), 1)), 'lineProps', {'color', 'b'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(remapped  ==1, :)), nanmean(neuron_sem_array{1, arrays_to_examine(2)}(remapped  ==1, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], 'new (safe block)', 'new (risky blocks)')


binned_data_mean_experimental = mean(binned_data(strcmp(experimental_grps_updated.group, 'Experimental'), :));
binned_data_mean_one_context = mean(binned_data(strcmp(experimental_grps_updated.group, 'One Context'), :));
binned_data_mean_no_shock = mean(binned_data(strcmp(experimental_grps_updated.group, 'No Shock'), :));

binned_data_sem_experimental = binned_data_std/sqrt(size(binned_data(strcmp(experimental_grps_updated.group, 'Experimental')), 1));
binned_data_sem_one_context = binned_data_std/sqrt(size(binned_data(strcmp(experimental_grps_updated.group, 'One Context')), 1));
binned_data_sem_no_shock = binned_data_std/sqrt(size(binned_data(strcmp(experimental_grps_updated.group, 'No Shock')), 1));