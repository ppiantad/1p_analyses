animalIDs = fieldnames(final_DLC);

session_to_analyze = 'D4';

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
% experimental_grps = readtable('E:\MATLAB\my_repo\context fear\organize_SLEAP_data\full_pilot_mice.xlsx');
experimental_grps = readtable('I:\MATLAB\my_repo\context fear\organize_DLC_data\pilot groups.xlsx');


%%
clear aversive_context safe_context mean_aversive mean_safe mean_aversive_context mean_safe_context


mouse_count = 0;
for gg = 1:size(animalIDs, 1)
    current_mouse = animalIDs{gg};
    
    if strcmp(final_DLC.(current_mouse).experimental_grp, 'Experimental')
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

experimental_grps = readtable('I:\MATLAB\my_repo\context fear\organize_DLC_data\pilot groups.xlsx');

% experimental_grps = readtable('I:\MATLAB\my_repo\context fear\organize_DLC_data\PFC mice.xlsx');

% experimental_grps = readtable('I:\MATLAB\my_repo\context fear\organize_SLEAP_data\full_pilot_mice.xlsx');

% Define parameters
threshold = 1; % Velocity threshold
sample_duration = 0.03; % Duration of each sample in seconds
min_duration = 2; % Minimum duration to trigger labeling in seconds
frame_rate = 30; % Frames per second
% Calculate the minimum number of consecutive rows needed
min_samples = min_duration / sample_duration;

animalIDs = fieldnames(final_DLC);

session_to_analyze = 'D1_Afternoon';

mouse_count = 0;
for gg = 1:size(animalIDs, 1)
    current_mouse = animalIDs{gg};
    if strcmp(session_to_analyze, 'D1_Morning') & strcmp(current_mouse, 'B57417')
        continue
    else
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

            % Find the row index where the 'mouse' column matches 'current_mouse'
            row_idx = strcmp(experimental_grps.mouse, current_mouse);

            freeze_data(mouse_count, :) = labels(1:21590)';
            % freeze_data(mouse_count, :) = DLC_data_mouse.was_freezing(1:21590);
            experimental_grps_updated(mouse_count, :) = experimental_grps(row_idx, :);
        end
    end
end

% Calculate the total number of samples from body_velocity
num_samples = size(freeze_data, 2); % Assuming body_velocity is a table column

% Generate time array in seconds
time_array = (0:num_samples-1) / frame_rate;

% Convert time array to minutes if needed
time_array_minutes = time_array / 60;

no_shk_period_ind = time_array_minutes < 4;

avg_freeze_no_shk_period = sum(freeze_data(:, no_shk_period_ind), 2)/sum(no_shk_period_ind == 1);

avg_freeze_shk_period = sum(freeze_data(:, ~no_shk_period_ind), 2)/sum(~no_shk_period_ind == 1);

experimental_data_no_shk_period = avg_freeze_no_shk_period(strcmp(experimental_grps_updated.group, 'Experimental'), :);
experimental_sem_no_shk_period = std(experimental_data_no_shk_period)/sqrt(size(experimental_data_no_shk_period, 1));
experimental_mice_no_shk_period = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental'), :);

one_context_data_no_shk_period = avg_freeze_no_shk_period(strcmp(experimental_grps_updated.group, 'One Context'), :);
one_context_sem_no_shk_period = std(one_context_data_no_shk_period)/sqrt(size(one_context_data_no_shk_period, 1));
one_context_mice_no_shk_period = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'One Context'), :);

no_shock_data_no_shk_period = avg_freeze_no_shk_period(strcmp(experimental_grps_updated.group, 'No Shock'), :);
no_shock_sem_no_shk_period = std(no_shock_data_no_shk_period)/sqrt(size(no_shock_data_no_shk_period, 1));
no_shock_mice_no_shk_period = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'No Shock'), :);

experimental_data_shk_period = avg_freeze_shk_period(strcmp(experimental_grps_updated.group, 'Experimental'), :);
experimental_sem_shk_period = std(experimental_data_shk_period)/sqrt(size(experimental_data_shk_period, 1));
experimental_mice_shk_period = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental'), :);

one_context_data_shk_period = avg_freeze_shk_period(strcmp(experimental_grps_updated.group, 'One Context'), :);
one_context_sem_shk_period = std(one_context_data_shk_period)/sqrt(size(one_context_data_shk_period, 1));
one_context_mice_shk_period = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'One Context'), :);

no_shock_data_shk_period = avg_freeze_shk_period(strcmp(experimental_grps_updated.group, 'No Shock'), :);
no_shock_sem_shk_period = std(no_shock_data_shk_period)/sqrt(size(no_shock_data_shk_period, 1));
no_shock_mice_shk_period = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'No Shock'), :);

if any("sex" == string(experimental_grps.Properties.VariableNames))
    experimental_data_no_shk_period_males = avg_freeze_no_shk_period(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_sem_no_shk_period_males = std(experimental_data_no_shk_period_males)/sqrt(size(experimental_data_no_shk_period_males, 1));
    experimental_mice_no_shk_period_males = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);

    experimental_data_no_shk_period_females = avg_freeze_no_shk_period(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_sem_no_shk_period_females = std(experimental_data_no_shk_period_females)/sqrt(size(experimental_data_no_shk_period_females, 1));
    experimental_mice_no_shk_period_females = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);

    experimental_data_shk_period_males = avg_freeze_shk_period(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_sem_shk_period_males = std(experimental_data_shk_period_males)/sqrt(size(experimental_data_shk_period_males, 1));
    experimental_mice_shk_period_males = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);

    experimental_data_shk_period_females = avg_freeze_shk_period(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_sem_shk_period_females = std(experimental_data_shk_period_females)/sqrt(size(experimental_data_shk_period_females, 1));
    experimental_mice_shk_period_females = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
end

%%
% Compute mean data for each group for no_shk_period and shk_period

if any("sex" == string(experimental_grps.Properties.VariableNames))

    no_shk_data = [experimental_data_no_shk_period_males, experimental_data_no_shk_period_females];
    no_shk_means = [mean(experimental_data_no_shk_period_males); ...
        mean(experimental_data_no_shk_period_females)];

    no_shk_sems = [experimental_sem_no_shk_period_males; ...
        experimental_sem_no_shk_period_females];

    shk_data = [experimental_data_shk_period_males, experimental_data_shk_period_females];
    shk_means = [mean(experimental_data_shk_period_males); ...
        mean(experimental_data_shk_period_females)];

    shk_sems = [experimental_sem_shk_period_males; ...
        experimental_sem_shk_period_females];

    % Group labels for the x-axis
    group_labels = {'male experimental', 'female experimental'};

    % Define custom bar positions
    no_shk_x = [1, 2]; % Closer together for no_shk
    shk_x = [5, 6];    % Closer together for shk

    % Define colors for groups
    group_colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]}; % Blue, Orange, Yellow

else

    no_shk_data = {experimental_data_no_shk_period, one_context_data_no_shk_period, no_shock_data_no_shk_period};
    no_shk_means = [mean(experimental_data_no_shk_period); ...
        mean(one_context_data_no_shk_period); ...
        mean(no_shock_data_no_shk_period)];

    no_shk_sems = [experimental_sem_no_shk_period; ...
        one_context_sem_no_shk_period; ...
        no_shock_sem_no_shk_period];

    shk_data = {experimental_data_shk_period, one_context_data_shk_period, no_shock_data_shk_period};
    shk_means = [mean(experimental_data_shk_period); ...
        mean(one_context_data_shk_period); ...
        mean(no_shock_data_shk_period)];

    shk_sems = [experimental_sem_shk_period; ...
        one_context_sem_shk_period; ...
        no_shock_sem_shk_period];

    % Group labels for the x-axis
    group_labels = {'experimental', 'one_context', 'no_shock'};

    % Define custom bar positions
    no_shk_x = [1, 2, 3]; % Closer together for no_shk
    shk_x = [5, 6, 7];    % Closer together for shk

    % Define colors for groups
    group_colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]}; % Blue, Orange, Yellow

end


% Create a figure
figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
hold on;

% Plot bars
bar_width = 0.6; % Adjust bar width
for i = 1:length(group_labels)
    % Plot "No Shock" bars
    bar(no_shk_x(i), no_shk_means(i), bar_width, 'FaceColor', group_colors{i}, 'EdgeColor', 'none');
    % Plot "Shock" bars
    bar(shk_x(i), shk_means(i), bar_width, 'FaceColor', group_colors{i}, 'EdgeColor', 'none');
end

% Add error bars
errorbar(no_shk_x, no_shk_means, no_shk_sems, 'k', 'LineStyle', 'none', 'LineWidth', 1);
errorbar(shk_x, shk_means, shk_sems, 'k', 'LineStyle', 'none', 'LineWidth', 1);

jitter_amount = 0.2;
for i = 1:size(no_shk_data, 2) % Loop over groups (columns)
    no_shk_data_temp = [];
    shk_data_temp = [];

    no_shk_data_temp = no_shk_data{1, i};
    shk_data_temp = shk_data{1, i};

    % Scatter points for no_shk_data
    scatter(no_shk_x(i) + jitter_amount * (rand(size(no_shk_data_temp, 1), 1) - 0.5), ...
            no_shk_data_temp, 30, group_colors{i}, 'filled');
    % Scatter points for shk_data
    scatter(shk_x(i) + jitter_amount * (rand(size(shk_data_temp, 1), 1) - 0.5), ...
            shk_data_temp, 30, group_colors{i}, 'filled');
end

% Customization
xticks([mean(no_shk_x), mean(shk_x)]); % Set ticks at the center of each group

xlim([0, 8]); % Adjust x-axis limits to provide spacing

legend(group_labels, 'Location', 'NorthWest');
ylim([0 .7])

hold off;

%%
% Number of bins
num_bins = 24;

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




if any("sex" == string(experimental_grps.Properties.VariableNames))
    experimental_data_male = binned_data(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_sem_male = std(experimental_data_male)/sqrt(size(experimental_data_male, 1));
    experimental_mice_male = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);

    experimental_data_female = binned_data(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_sem_female = std(experimental_data_female)/sqrt(size(experimental_data_female, 1));
    experimental_mice_female = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);

    figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
    hold on;

    h(1) = shadedErrorBar(1:num_bins, mean(experimental_data_male), experimental_sem_male, 'lineProps', {'color', 'r'});
    h(2) = shadedErrorBar(1:num_bins, mean(experimental_data_female), experimental_sem_female, 'lineProps', {'color', 'k'});

    % legend([h(1).mainLine h(2).mainLine], 'new (safe block)', 'new (risky blocks)')
    % Adjust x-axis ticks and labels
    % Adjust x-axis ticks and labels
    xlim([1 num_bins]); % Set x-axis limits to match the data range
    xticks([1:4:num_bins, num_bins]); % Add the last tick explicitly
    xticklabels([0:2:12]); % Label ticks with corresponding time in minutes

    ylim([0 1]); % Set y-axis limits

else

    experimental_data = binned_data(strcmp(experimental_grps_updated.group, 'Experimental'), :);
    experimental_sem = std(experimental_data)/sqrt(size(experimental_data, 1));
    experimental_mice = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental'), :);

    one_context_data = binned_data(strcmp(experimental_grps_updated.group, 'One Context'), :);
    one_context_sem = std(one_context_data)/sqrt(size(one_context_data, 1));
    one_context_mice = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'One Context'), :);

    no_shock_data = binned_data(strcmp(experimental_grps_updated.group, 'No Shock'), :);
    no_shock_sem = std(no_shock_data)/sqrt(size(no_shock_data, 1));
    no_shock_mice = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'No Shock'), :);

    figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
    hold on;

    h(1) = shadedErrorBar(1:num_bins, mean(experimental_data), experimental_sem, 'lineProps', {'color', 'r'});
    h(2) = shadedErrorBar(1:num_bins, mean(one_context_data), one_context_sem, 'lineProps', {'color', 'k'});
    h(2) = shadedErrorBar(1:num_bins, mean(no_shock_data), no_shock_sem, 'lineProps', {'color', 'b'});

    % legend([h(1).mainLine h(2).mainLine], 'new (safe block)', 'new (risky blocks)')
    % Adjust x-axis ticks and labels
    % Adjust x-axis ticks and labels
    xlim([1 num_bins]); % Set x-axis limits to match the data range
    xticks([1:4:num_bins, num_bins]); % Add the last tick explicitly
    xticklabels([0:2:12]); % Label ticks with corresponding time in minutes

    ylim([0 1]); % Set y-axis limits


end




%% plot individual data from a given session - make sure to update variables and indices if using!
% Load the data
body_velocity = final_DLC.B46851        .D1_Afternoon.movement_data.body_velocity; % Assuming this is a table column
freeze_data_extracted = freeze_data(6,:); % Get the first row of freeze_data

% Trim the data
start_index = 101; % Omit first 100 indices
end_index = 21590; % Limit to index 21590
% Create time vector in minutes

% Create time vector in minutes
time_vector = (start_index:end_index) / (frame_rate * 60);

body_velocity_trimmed = body_velocity(start_index:end_index);
freeze_data_trimmed = freeze_data_extracted(start_index:end_index);

% Find the maximum velocity for plotting rectangles
max_velocity = max(body_velocity_trimmed);

% Plot the body velocity
figure;
plot(time_vector, body_velocity_trimmed, 'b', 'LineWidth', 1.5);
hold on;

% Find runs of 1s in freeze_data_trimmed
freeze_start_indices = find(diff([0, freeze_data_trimmed]) == 1);
freeze_end_indices = find(diff([freeze_data_trimmed, 0]) == -1);

% Add rectangles for freeze periods
for i = 1:length(freeze_start_indices)
    % No additional adjustment needed for freeze indices relative to time_vector
    x_start = time_vector(freeze_start_indices(i)); 
    x_end = time_vector(freeze_end_indices(i)); 
    rectangle('Position', [x_start, 0, x_end - x_start, max_velocity], ...
              'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', 'none');
end

% Format the plot
hold off;
xlabel('Time (min)');
ylabel('Velocity (cm/s)');

ylim([0, 300]);
x_ticks = 0:1:max(time_vector); % 1-minute increments
xticks(x_ticks);
xlim([min(time_vector), max(time_vector)]);
% legend('Body Velocity', 'Freeze Periods');

%%
uv.evtWin = [-2 4]; %what time do y
frame_rate = 30; % Frames per second
% Parameters
shock_start_time = 4 * 60; % First shock in seconds
shock_interval = 60; % Interval between shocks in seconds
shock_duration = 2; % Duration of each shock in seconds
num_shocks = 6;

% Initialize variables
shk_on = zeros(1, num_shocks);
shk_off = zeros(1, num_shocks);

% Calculate on and off times for each shock
for i = 0:(num_shocks-1)
    shk_on(i+1) = shock_start_time + i * shock_interval; % Shock start time in seconds
    shk_off(i+1) = shk_on(i+1) + shock_duration; % Shock end time in seconds
end

ts1 = (uv.evtWin(1):1/frame_rate:uv.evtWin(2)-1/frame_rate);
animalIDs = (fieldnames(final_DLC));
mouse_count = 0;
for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    if isfield(final_DLC.(currentanimal), session_to_analyze)
        mouse_count = mouse_count+1;


        body_velocity = [];
        labels = [];
        % Get the body_velocity column
        body_velocity = final_DLC.(currentanimal).(session_to_analyze).movement_data.body_velocity;

        eTS = shk_on'; %get time stamps

        %calculate time windows for each event
        evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
        % Define parameters
        

        % Calculate the total number of samples from body_velocity
        num_samples = height(body_velocity); % Assuming body_velocity is a table column

        % Generate time array in seconds
        time_array = (0:num_samples-1) / frame_rate;

        % Convert time array to minutes if needed
        time_array_minutes = time_array / 60;
        for t = 1:size(eTS,1)
            % set each trial's temporal boundaries
            timeWin = [eTS(t)+uv.evtWin(1,1):1/frame_rate:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
            % BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
            % unitTrace_zscored = zscore(unitTrace);


            if min(timeWin) > min(time_array) && max(timeWin) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                % get unit event counts in trials
                % get unit ca traces in trials
                idx = time_array >= min(timeWin) & time_array < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);

              
                velocity_for_shocks(t,:) = body_velocity(idx);
                 

            end

        end
        velocity_for_shocks_mouse{ii} = velocity_for_shocks;
        mean_velocity_for_shocks(ii, :) = mean(velocity_for_shocks);
        sem_velocity_for_shocks(ii, :) = std(velocity_for_shocks)/sqrt(size(velocity_for_shocks, 1));
        clear velocity_for_shocks

    end

end

if any("sex" == string(experimental_grps.Properties.VariableNames))
    experimental_data_males = mean_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_sems_mice_males = sem_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_sem_males = std(experimental_data_males)/sqrt(size(experimental_data_males, 1));
    experimental_mice_males = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);

    experimental_data_females = mean_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_sems_mice_females = sem_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_sem_females = std(experimental_data_females)/sqrt(size(experimental_data_females, 1));
    experimental_mice_females = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    
    figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
    hold on;

    h(1) = shadedErrorBar(ts1, mean(experimental_data_males), experimental_sem_males, 'lineProps', {'color', 'r'});
    h(2) = shadedErrorBar(ts1, mean(experimental_data_females), experimental_sem_females, 'lineProps', {'color', 'k'});

    % legend([h(1).mainLine h(2).mainLine], 'new (safe block)', 'new (risky blocks)')
    ylim([0 90]);

    mean_data_array = {experimental_data_males, experimental_data_females};
    sem_data_array = {experimental_sems_mice_males, experimental_sems_mice_females};

    % need to make sure consec_thresh in perm_and_bCI is set to 10!
    [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1);
else

    experimental_data = mean_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'Experimental'), :);
    experimental_sems_mice = sem_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'Experimental'), :);
    experimental_sem = std(experimental_data)/sqrt(size(experimental_data, 1));
    experimental_mice = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental'), :);

    one_context_data = mean_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'One Context'), :);
    one_context_sems_mice = sem_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'One Context'), :);
    one_context_sem = std(one_context_data)/sqrt(size(one_context_data, 1));
    one_context_mice = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'One Context'), :);

    no_shock_data = mean_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'No Shock'), :);
    no_shock_sems_mice = sem_velocity_for_shocks(strcmp(experimental_grps_updated.group, 'No Shock'), :);
    no_shock_sem = std(no_shock_data)/sqrt(size(no_shock_data, 1));
    no_shock_mice = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'No Shock'), :);



    figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
    hold on;

    h(1) = shadedErrorBar(ts1, mean(experimental_data), experimental_sem, 'lineProps', {'color', 'r'});
    h(2) = shadedErrorBar(ts1, mean(one_context_data), one_context_sem, 'lineProps', {'color', 'k'});
    h(2) = shadedErrorBar(ts1, mean(no_shock_data), no_shock_sem, 'lineProps', {'color', 'b'});

    ylim([0 90]);

    mean_data_array = {experimental_data, one_context_data, no_shock_data};
    sem_data_array = {experimental_sems_mice, one_context_sems_mice, no_shock_sems_mice};

    % need to make sure consec_thresh in perm_and_bCI is set to 10!
    [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1);

end
%% discrimination sessions
clear aversive_context safe_context mean_aversive mean_safe mean_aversive_context mean_safe_context

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


% experimental_grps = readtable('I:\MATLAB\my_repo\context fear\organize_DLC_data\PFC mice.xlsx');
% experimental_grps = readtable('E:\MATLAB\my_repo\context fear\organize_SLEAP_data\full_pilot_mice.xlsx');
experimental_grps = readtable('I:\MATLAB\my_repo\context fear\organize_DLC_data\pilot groups.xlsx');

% Define parameters
threshold = 1; % Velocity threshold
sample_duration = 0.03; % Duration of each sample in seconds
min_duration = 2; % Minimum duration to trigger labeling in seconds
frame_rate = 30; % Frames per second
% Calculate the minimum number of consecutive rows needed
min_samples = min_duration / sample_duration;

animalIDs = fieldnames(final_DLC);

session_to_analyze = 'D3';

mouse_count = 0;
for gg = 1:size(animalIDs, 1)
    current_mouse = animalIDs{gg};
    if strcmp(session_to_analyze, 'D1_Morning') & strcmp(current_mouse, 'B57417')
        continue
    elseif strcmp(session_to_analyze, 'D3') & strcmp(current_mouse, 'B51618')
        continue
    else
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

            % Find the row index where the 'mouse' column matches 'current_mouse'
            row_idx = strcmp(experimental_grps.mouse, current_mouse);

            freeze_data(mouse_count, :) = labels(1:21590)';
            % freeze_data(mouse_count, :) = DLC_data_mouse.was_freezing(1:21590)';
            experimental_grps_updated(mouse_count, :) = experimental_grps(row_idx, :);
        end
    end
end



if strcmp(session_to_analyze, 'D3')
    for gg = 1:size(freeze_data, 1)
        freeze_data_session = freeze_data(gg, :);

        for qq = 1:num_repeats
            safe_context{gg, qq} = freeze_data_session(stimulus_frames{1, 1}(qq,1)+1:stimulus_frames{1, 1}(qq,2));
        end

        for qq = 1:num_repeats
            if qq <= 2
                aversive_context{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:stimulus_frames{1, 2}(qq,2));
            elseif qq > 2
                aversive_context{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:end);
            end
        end
    end

elseif strcmp(session_to_analyze, 'D4')

    for gg = 1:size(freeze_data, 1)
        freeze_data_session = freeze_data(gg, :);

        for qq = 1:num_repeats
            aversive_context{gg, qq} = freeze_data_session(stimulus_frames{1, 1}(qq,1)+1:stimulus_frames{1, 1}(qq,2));
        end

        for qq = 1:num_repeats
            if qq <= 2
                safe_context{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:stimulus_frames{1, 2}(qq,2));
            elseif qq > 2
                safe_context{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:end);
            end
        end
    end
end




% for gg = 1:size(freeze_data, 1)
%     freeze_data_session = freeze_data(gg, :);
% 
%     for qq = 1:num_repeats
%         safe_context{gg, qq} = freeze_data_session(stimulus_frames{1, 1}(qq,1)+1:stimulus_frames{1, 1}(qq,2));
%     end
% 
%     for qq = 1:num_repeats
%         if qq <= 2
%             aversive_context{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:stimulus_frames{1, 2}(qq,2));
%         elseif qq > 2
%             aversive_context{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:end);
%         end
%     end
% end

for pp = 1:size(safe_context, 1)
    for ff = 1:size(safe_context, 2)
        % mean_safe_context(mouse_count, pp) = sum(freeze_data(stimulus_frames{1, 1}(pp,1)+1:stimulus_frames{1, 1}(pp,2)))/size(safe_context, 1);
        % Calculate the proportion of 1s in the column
        b = sum(safe_context{pp,ff}) / size(safe_context{pp,ff}, 2);

        % Store the proportion in mean_aversive_context
        mean_safe_context(pp, ff) = b;

        sem_safe_context(pp, ff) = b/sqrt(size(safe_context{pp,ff}, 2));

        % % Calculate the number of rows in the column
        % n = size(safe_context{pp,ff}, 2);
        % 
        % % Calculate the standard error of the proportion
        % standard_error_safe(pp, ff) = sqrt(b * (1 - b) / n);
    end
end


for pp = 1:size(aversive_context, 1)
    for ff = 1:size(aversive_context, 2)
        % mean_safe_context(mouse_count, pp) = sum(freeze_data(stimulus_frames{1, 1}(pp,1)+1:stimulus_frames{1, 1}(pp,2)))/size(safe_context, 1);
        % Calculate the proportion of 1s in the column
        b = sum(aversive_context{pp,ff}) / size(aversive_context{pp,ff}, 2);

        % Store the proportion in mean_aversive_context
        mean_aversive_context(pp, ff) = b;

        sem_aversive_context(pp, ff) = b/sqrt(size(aversive_context{pp,ff}, 2));

        % % Calculate the number of rows in the column
        % n = size(aversive_context{pp,ff}, 2);
        % 
        % % Calculate the standard error of the proportion
        % standard_error_safe(pp, ff) = sqrt(b * (1 - b) / n);
    end
end

if any("sex" == string(experimental_grps.Properties.VariableNames))


    experimental_data_safe_male = mean_safe_context(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_mean_safe_male = mean(experimental_data_safe_male);
    experimental_mean_safe_collapsed_male = mean(experimental_data_safe_male, 2);
    experimental_sem_safe_male = std(experimental_data_safe_male)/sqrt(size(experimental_data_safe_male, 1));
    experimental_mice_safe_male = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);

    experimental_data_aversive_male = mean_aversive_context(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_mean_aversive_male = mean(experimental_data_aversive_male);
    experimental_mean_aversive_collapsed_male = mean(experimental_data_aversive_male, 2);
    experimental_sem_aversive_male = std(experimental_data_aversive_male)/sqrt(size(experimental_data_aversive_male, 1));
    experimental_mice_aversive_male = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);

    experimental_data_safe_female = mean_safe_context(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_mean_safe_female = mean(experimental_data_safe_female);
    experimental_mean_safe_collapsed_female = mean(experimental_data_safe_female, 2);
    experimental_sem_safe_female = std(experimental_data_safe_female)/sqrt(size(experimental_data_safe_female, 1));
    experimental_mice_safe_female = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);

    experimental_data_aversive_female = mean_aversive_context(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_mean_aversive_female = mean(experimental_data_aversive_female);
    experimental_mean_aversive_collapsed_female = mean(experimental_data_aversive_female, 2);
    experimental_sem_aversive_female = std(experimental_data_aversive_female)/sqrt(size(experimental_data_aversive_female, 1));
    experimental_mice_aversive_female = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);


else

    experimental_data_safe = mean_safe_context(strcmp(experimental_grps_updated.group, 'Experimental'), :);
    experimental_mean_safe = mean(experimental_data_safe);
    experimental_mean_safe_collapsed = mean(experimental_data_safe, 2);
    experimental_sem_safe = std(experimental_data_safe)/sqrt(size(experimental_data_safe, 1));
    experimental_mice_safe = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental'), :);

    one_context_data_safe = mean_safe_context(strcmp(experimental_grps_updated.group, 'One Context'), :);
    one_context_mean_safe = mean(one_context_data_safe);
    one_context_mean_safe_collapsed = mean(one_context_data_safe, 2);
    one_context_sem_safe = std(one_context_data_safe)/sqrt(size(one_context_data_safe, 1));
    one_context_mice_safe = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'One Context'), :);

    no_shock_data_safe = mean_safe_context(strcmp(experimental_grps_updated.group, 'No Shock'), :);
    no_shock_mean_safe = mean(no_shock_data_safe);
    no_shock_mean_safe_collapsed = mean(no_shock_data_safe, 2);
    no_shock_sem_safe = std(no_shock_data_safe)/sqrt(size(no_shock_data_safe, 1));
    no_shock_mice_safe = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'No Shock'), :);

    experimental_data_aversive = mean_aversive_context(strcmp(experimental_grps_updated.group, 'Experimental'), :);
    experimental_mean_aversive = mean(experimental_data_aversive);
    experimental_mean_aversive_collapsed = mean(experimental_data_aversive, 2);
    experimental_sem_aversive = std(experimental_data_aversive)/sqrt(size(experimental_data_aversive, 1));
    experimental_mice_aversive = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental'), :);

    one_context_data_aversive = mean_aversive_context(strcmp(experimental_grps_updated.group, 'One Context'), :);
    one_context_mean_aversive = mean(one_context_data_aversive);
    one_context_mean_aversive_collapsed = mean(one_context_data_aversive, 2);
    one_context_sem_aversive = std(one_context_data_aversive)/sqrt(size(one_context_data_aversive, 1));
    one_context_mice_aversive = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'One Context'), :);

    no_shock_data_aversive = mean_aversive_context(strcmp(experimental_grps_updated.group, 'No Shock'), :);
    no_shock_mean_aversive = mean(no_shock_data_aversive);
    no_shock_mean_aversive_collapsed = mean(no_shock_data_aversive, 2);
    no_shock_sem_aversive = std(no_shock_data_aversive)/sqrt(size(no_shock_data_aversive, 1));
    no_shock_mice_aversive = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'No Shock'), :);
end

%%

if strcmp(session_to_analyze, 'D3')
    
    if any("sex" == string(experimental_grps.Properties.VariableNames))

        interleave_mean_experimental_male = zeros(size(experimental_data_safe_male, 1), 6);
        interleave_mean_experimental_male(:, [1 3 5]) = experimental_data_safe_male;
        interleave_mean_experimental_male(:, [2 4 6]) = experimental_data_aversive_male;

        interleave_mean_experimental_female = zeros(size(experimental_data_safe_female, 1), 6);
        interleave_mean_experimental_female(:, [1 3 5]) = experimental_data_safe_female;
        interleave_mean_experimental_female(:, [2 4 6]) = experimental_data_aversive_female;

    else

        interleave_mean_experimental = zeros(size(experimental_data_safe, 1), 6);
        interleave_mean_experimental(:, [1 3 5]) = experimental_data_safe;
        interleave_mean_experimental(:, [2 4 6]) = experimental_data_aversive;

        interleave_mean_one_context = zeros(size(one_context_data_safe, 1), 6);
        interleave_mean_one_context(:, [1 3 5]) = one_context_data_safe;
        interleave_mean_one_context(:, [2 4 6]) = one_context_data_aversive;

        interleave_mean_no_shock = zeros(size(no_shock_data_safe, 1), 6);
        interleave_mean_no_shock(:, [1 3 5]) = no_shock_data_safe;
        interleave_mean_no_shock(:, [2 4 6]) = no_shock_data_aversive;

        % interleave_sem_experimental = zeros(size(mean_safe_context, 1), 6);
        % interleave_sem_experimental(:, [1 3 5]) = standard_error_safe;
        % interleave_sem_experimental(:, [2 4 6]) = standard_error_aversive;

    end

elseif strcmp(session_to_analyze, 'D4')
    
    if any("sex" == string(experimental_grps.Properties.VariableNames))


        interleave_mean_experimental_male = zeros(size(experimental_data_safe_male, 1), 6);
        interleave_mean_experimental_male(:, [1 3 5]) = experimental_data_aversive_male;
        interleave_mean_experimental_male(:, [2 4 6]) = experimental_data_safe_male;

        interleave_mean_experimental_female = zeros(size(experimental_data_safe_female, 1), 6);
        interleave_mean_experimental_female(:, [1 3 5]) = experimental_data_aversive_female;
        interleave_mean_experimental_female(:, [2 4 6]) = experimental_data_safe_female;
    else

        interleave_mean_experimental = zeros(size(experimental_data_safe, 1), 6);
        interleave_mean_experimental(:, [1 3 5]) = experimental_data_aversive;
        interleave_mean_experimental(:, [2 4 6]) = experimental_data_safe;

        interleave_mean_one_context = zeros(size(one_context_data_safe, 1), 6);
        interleave_mean_one_context(:, [1 3 5]) = one_context_data_aversive;
        interleave_mean_one_context(:, [2 4 6]) = one_context_data_safe;

        interleave_mean_no_shock = zeros(size(no_shock_data_safe, 1), 6);
        interleave_mean_no_shock(:, [1 3 5]) = no_shock_data_aversive;
        interleave_mean_no_shock(:, [2 4 6]) = no_shock_data_safe;

    end
end

if any("sex" == string(experimental_grps.Properties.VariableNames))
    figure; shadedErrorBar(1:size(mean(interleave_mean_experimental_male), 2), mean(interleave_mean_experimental_male), std(interleave_mean_experimental_male));
    hold on; shadedErrorBar(1:size(mean(interleave_mean_experimental_female), 2), mean(interleave_mean_experimental_female), std(interleave_mean_experimental_female));
    hold off;
    % xticks([1:4:num_bins, num_bins]); % Add the last tick explicitly
    % xticklabels([0:2:12]); % Label ticks with corresponding time in minutes




else

    figure; shadedErrorBar(1:size(mean(interleave_mean_experimental), 2), mean(interleave_mean_experimental), std(interleave_mean_experimental));
    hold on; shadedErrorBar(1:size(mean(interleave_mean_one_context), 2), mean(interleave_mean_one_context), std(interleave_mean_one_context));
    hold on; shadedErrorBar(1:size(mean(interleave_mean_no_shock), 2), mean(interleave_mean_no_shock), std(interleave_mean_no_shock));
    hold off;
    % xticks([1:4:num_bins, num_bins]); % Add the last tick explicitly
    % xticklabels([0:2:12]); % Label ticks with corresponding time in minutes

end



%%
% Number of bins
num_bins = 24;

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


if any("sex" == string(experimental_grps.Properties.VariableNames))

    experimental_data_male = binned_data(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_sem_male = std(experimental_data_male)/sqrt(size(experimental_data_male, 1));
    experimental_mice_male = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);

    experimental_data_female = binned_data(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_sem_female = std(experimental_data_female)/sqrt(size(experimental_data_female, 1));
    experimental_mice_female = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
 
    figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
    hold on;
    h(1) = shadedErrorBar(1:num_bins, mean(experimental_data_male), experimental_sem_male, 'lineProps', {'color', 'r'});
    h(2) = shadedErrorBar(1:num_bins, mean(experimental_data_female), experimental_sem_female, 'lineProps', {'color', 'k'});

    % legend([h(1).mainLine h(2).mainLine], 'new (safe block)', 'new (risky blocks)')

    % Adjust x-axis ticks and labels
    xlim([1 num_bins]); % Set x-axis limits to match the data range
    xticks([1:4:num_bins, num_bins]); % Add the last tick explicitly
    xticklabels([0:2:12]); % Label ticks with corresponding time in minutes

    ylim([0 0.9]); % Set y-axis limits
    % Add rectangles for freeze periods
    % x_tick_indices = [1:4:num_bins]
    % for i = 1:length(x_tick_indices)
    %     % No additional adjustment needed for freeze indices relative to time_vector
    %     x_start = x_tick_indices(i);
    %     x_end = x_tick_indices(i+1);
    %     rectangle('Position', [x_start, 0, x_end - x_start, 1], ...
    %               'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', 'none');
    % end
    %
    %


else


    experimental_data = binned_data(strcmp(experimental_grps_updated.group, 'Experimental'), :);
    experimental_sem = std(experimental_data)/sqrt(size(experimental_data, 1));
    experimental_mice = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental'), :);

    one_context_data = binned_data(strcmp(experimental_grps_updated.group, 'One Context'), :);
    one_context_sem = std(one_context_data)/sqrt(size(one_context_data, 1));
    one_context_mice = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'One Context'), :);

    no_shock_data = binned_data(strcmp(experimental_grps_updated.group, 'No Shock'), :);
    no_shock_sem = std(no_shock_data)/sqrt(size(no_shock_data, 1));
    no_shock_mice = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'No Shock'), :);

    figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
    hold on;

    h(1) = shadedErrorBar(1:num_bins, mean(experimental_data), experimental_sem, 'lineProps', {'color', 'r'});
    h(2) = shadedErrorBar(1:num_bins, mean(one_context_data), one_context_sem, 'lineProps', {'color', 'k'});
    h(2) = shadedErrorBar(1:num_bins, mean(no_shock_data), no_shock_sem, 'lineProps', {'color', 'b'});


    % legend([h(1).mainLine h(2).mainLine], 'new (safe block)', 'new (risky blocks)')

    % Adjust x-axis ticks and labels
    xlim([1 num_bins]); % Set x-axis limits to match the data range
    xticks([1:4:num_bins, num_bins]); % Add the last tick explicitly
    xticklabels([0:2:12]); % Label ticks with corresponding time in minutes

    ylim([0 0.9]); % Set y-axis limits
    % Add rectangles for freeze periods
    % x_tick_indices = [1:4:num_bins]
    % for i = 1:length(x_tick_indices)
    %     % No additional adjustment needed for freeze indices relative to time_vector
    %     x_start = x_tick_indices(i);
    %     x_end = x_tick_indices(i+1);
    %     rectangle('Position', [x_start, 0, x_end - x_start, 1], ...
    %               'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', 'none');
    % end
    %
    %


end




