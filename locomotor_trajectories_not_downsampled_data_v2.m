
animalIDs = (fieldnames(final_SLEAP));

select_mouse = 'BLA_Insc_26';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

session_to_analyze = 'RDT_D1';

shapeData = final_SLEAP.(select_mouse).(session_to_analyze).shapeData;

onset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.stTime';
choice_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.choiceTime';
offset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.collectionTime';
fs_cam = 30; %set sampling rate according to camera, this is hard coded for now
time_ranges_trials = [onset_trials; choice_trials; offset_trials];



SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;

X_data = SLEAP_data.corrected_x_pix;
Y_data = SLEAP_data.corrected_y_pix;

% velocity_data = final_SLEAP.(select_mouse).(session_to_analyze).zscored_SLEAP_data_velocity';

velocity_data = zscore(SLEAP_data.vel_cm_s)';

% SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data;


BehavData = final_SLEAP.(select_mouse).(session_to_analyze).BehavData;
adjusted_start_time = BehavData.TrialPossible(1)-60;
SLEAP_data.idx_time = SLEAP_data.idx_time+adjusted_start_time;
% gcamp_normalized = ((Y_dF_all_session)-mean(Y_dF_all_session))/std(Y_dF_all_session);
% SLEAP_data_vel_filtered_session_normalized = ((SLEAP_data_vel_filtered_session)-mean(SLEAP_data_vel_filtered_session))/std(SLEAP_data_vel_filtered_session);


% trial_starts_array = BehavData.stTime-BehavData.choiceTime;
% trial_ends_array = BehavData.collectionTime - BehavData.choiceTime;





 %%
        % FILTER ALL EXISTING DATA ON THESE TIME RANGES
        % filter streams
        if ~isempty(SLEAP_data)
            filtered_motion = [];
            max_ind = SLEAP_data.idx_frame(end);
            idx_frame_redo = 1:1:size(SLEAP_data, 1);
            good_index = 1;
            for j = 1:size(time_ranges_trials,2)
                onset = round(time_ranges_trials(1,j)*fs_cam)+1;
                choice = round(time_ranges_trials(2,j)*fs_cam)+1;
                offset = round(time_ranges_trials(3,j)*fs_cam)+1;

                % throw it away if onset or offset extends beyond recording window
                if isinf(offset)
                    if onset <= max_ind && onset > 0
                        filtered_motion{good_index} = SLEAP_data(onset:end);
                        break %return
                    end
                else
                    if offset <= max_ind && offset > 0 && onset <= max_ind && onset > 0
                        % buffering this by adding +1 to the end time for now, for
                        % some reason the array seems too short without?
                        % after some extensive checking, it seems like the
                        % strangeness where the body @ start and @ end does not
                        % overlap often comes from the fact that the mouse's tail
                        % can trigger the IR beam in the food cup on a non-trivial
                        % # of trials
                        filtered_motion{j}= [X_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j))'; Y_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j))']; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
                        choice_times{j} = [X_data(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'; Y_data(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'];
                        filtered_velocity{j}= velocity_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j));
                        % filtered_gcamp{j}= Y_dF_all_session(gcamp_samples(onset:offset));
                        % filtered{j} = Y_data_filtered(SLEAP_data.idx_frame(onset:offset))'; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
                        good_index = good_index + 1;
                    end
                end
            end
            % if KEEPDATA
            %     data.streams.Motion.filtered = filtered;
            % else
            %     data.streams.Motion.data = filtered;
            %     data.streams.Motion.filtered = [];
            % end
        end

%%
figure;
hold on
% Adding labels and title (customize as needed)
xlabel('X-axis Label');
ylabel('Y-axis Label');
title('Plotting the Line');

% Display grid if desired
grid on;

for ii = 1:size(filtered_motion, 2)
    % Assuming filtered{1, 1} contains your X and Y coordinates
    X = filtered_motion{1, ii}(1, :);
    Y = filtered_motion{1, ii}(2, :);

    % Plotting the line
    plot(X, Y, '-x');  % You can use different markers or line styles if needed


end
hold off
%%
% Combine X and Y coordinates into a 2-column matrix
% coordinates = [X', Y'];

% Calculate pairwise distances between consecutive points
for qq = 1:size(filtered_motion, 2)
    coordinates = filtered_motion{1, qq}';
    distances = pdist2(coordinates, coordinates);

    % Sum up the distances to get the total path length
    path_length = sum(diag(distances, 1));

    disp(['Path Length: ', num2str(path_length)]);
    
    distances_matrix{qq} = distances;
    path_length_array(qq) = path_length;
    clear coordinates distances path_length
end
%%
mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 1))
mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 2))
mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 3))

mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 1))
mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 2))
mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 3))

%% this codeblock attempts to plot a bunch of trajectories, with the color of the lines representing the calcium activity, and the size of the lines representing the velocity


for i = 1:size(filtered_velocity,2)
    original_array = filtered_velocity{1,i};
    if ~isempty(original_array)
        % Assuming filtered_gcamp{1,1} is your array

        % Normalize the array to be between 0 and 1
        min_value = min(original_array(:));
        max_value = max(original_array(:));
        normalized_velocity_bounded{i} = (original_array - min_value) / (max_value - min_value);
        clear original_array min_value max_value
    elseif isempty(original_array)
        normalized_velocity_bounded{i} = [];
        disp('The array is empty. Skipping normalization.');
    end
end

% large_block_1_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 1 & BehavData.ForceFree == 0;
% large_block_2_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 2 & BehavData.ForceFree == 0;
% large_block_3_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 3 & BehavData.ForceFree == 0;

large_block_1_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 1;
large_block_2_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 2;
large_block_3_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 3;


small_block_1_ind = BehavData.bigSmall == 0.3 & BehavData.Block == 1;
small_block_2_ind = BehavData.bigSmall == 0.3 & BehavData.Block == 2;
small_block_3_ind = BehavData.bigSmall == 0.3 & BehavData.Block == 3;


% Find the indices where large_block_1_ind is true
large_block_1_true_indices = find(large_block_1_ind);
large_block_2_true_indices = find(large_block_2_ind);
large_block_3_true_indices = find(large_block_3_ind);

small_block_1_true_indices = find(small_block_1_ind);
small_block_2_true_indices = find(small_block_2_ind);
small_block_3_true_indices = find(small_block_3_ind);


%%
% Loop through each line and plot
figure;
for j = 1:size(large_block_1_true_indices , 1) %num_lines
    plot_index = large_block_1_true_indices (j);
    x = [];
    y = [];
    % figure; %comment me out if you want everything plotted on the same fig! 
    
    % Extract X-Y coordinates
    x = filtered_motion{plot_index}(1, :);
    y = filtered_motion{plot_index}(2, :);
    
    % Extract velocity data
    velocity_values = normalized_velocity_bounded{plot_index};

    % Normalize velocity_values to the range [0, 1] for colormap mapping
    normalized_velocity = (velocity_values - min(velocity_values)) / (max(velocity_values) - min(velocity_values));

    % Set marker colors based on normalized_velocity
    marker_colors = colormap('jet');
    marker_colors_mapped = interp1(linspace(0, 1, size(marker_colors, 1)), marker_colors, normalized_velocity);

    % Plot the line using scatter with varying marker colors
    scatter(x, y, 50, marker_colors_mapped, 'filled');
    hold on;
    plot(x, y, 'Color', 'k');
    % Add a black square at the start of the line
    scatter(x(1), y(1), 200, 'k', 's', 'filled');

    % Add a black circke at the choice time
    scatter(choice_times{1, large_block_1_true_indices(j)}(1), choice_times{1, large_block_1_true_indices(j)}(2), 200, 'k', 'o', 'filled');


    % Add a black triangle at the end of the line
    scatter(x(end), y(end), 200, 'k', '^', 'filled');
end



% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end

hold on;
for j = 1:size(small_block_1_true_indices , 1) %num_lines
    plot_index = small_block_1_true_indices (j);
    x = [];
    y = [];
    % figure; %comment me out if you want everything plotted on the same fig! 
    
    % Extract X-Y coordinates
    x = filtered_motion{plot_index}(1, :);
    y = filtered_motion{plot_index}(2, :);
    
    % Extract velocity data
    velocity_values = normalized_velocity_bounded{plot_index};

    % Normalize velocity_values to the range [0, 1] for colormap mapping
    normalized_velocity = (velocity_values - min(velocity_values)) / (max(velocity_values) - min(velocity_values));

    % Set marker colors based on normalized_velocity
    marker_colors = colormap('jet');
    marker_colors_mapped = interp1(linspace(0, 1, size(marker_colors, 1)), marker_colors, normalized_velocity);

    % Plot the line using scatter with varying marker colors
    scatter(x, y, 50, marker_colors_mapped, 'filled');
    hold on;
    plot(x, y, 'Color', 'k');
    % Add a black square at the start of the line
    scatter(x(1), y(1), 200, 'k', 's', 'filled');

    % Add a black circke at the choice time
    scatter(choice_times{1, small_block_1_true_indices(j)}(1), choice_times{1, small_block_1_true_indices(j)}(2), 200, 'k', 'o', 'filled');


    % Add a black triangle at the end of the line
    scatter(x(end), y(end), 200, 'k', '^', 'filled');
end


hold off

%% Plot random large reward trial(s)

lines_to_plot = 1;
colors = distinguishable_colors(lines_to_plot);


% Randomly select indices
large_block_1_selected_indices = large_block_1_true_indices(randperm(length(large_block_1_true_indices), lines_to_plot ));
large_block_2_selected_indices = large_block_2_true_indices(randperm(length(large_block_2_true_indices), lines_to_plot ));
large_block_3_selected_indices = large_block_3_true_indices(randperm(length(large_block_3_true_indices), lines_to_plot ));

% Initialize variables to store min and max values across all sets of data
x_min_all = inf;
x_max_all = -inf;
y_min_all = inf;
y_max_all = -inf;

% Loop through large_block_1_selected_indices to find min and max values
for j = 1:length(large_block_1_selected_indices)
    x = filtered_motion{large_block_1_selected_indices(j)}(1, :);
    y = filtered_motion{large_block_1_selected_indices(j)}(2, :);

    x_min_all = min(x_min_all, min(x));
    x_max_all = max(x_max_all, max(x));
    y_min_all = min(y_min_all, min(y));
    y_max_all = max(y_max_all, max(y));
end

% Loop through large_block_2_selected_indices to find min and max values
for j = 1:length(large_block_2_selected_indices)
    x = filtered_motion{large_block_2_selected_indices(j)}(1, :);
    y = filtered_motion{large_block_2_selected_indices(j)}(2, :);

    x_min_all = min(x_min_all, min(x));
    x_max_all = max(x_max_all, max(x));
    y_min_all = min(y_min_all, min(y));
    y_max_all = max(y_max_all, max(y));
end

% Loop through large_block_3_selected_indices to find min and max values
for j = 1:length(large_block_3_selected_indices)
    x = filtered_motion{large_block_3_selected_indices(j)}(1, :);
    y = filtered_motion{large_block_3_selected_indices(j)}(2, :);

    x_min_all = min(x_min_all, min(x));
    x_max_all = max(x_max_all, max(x));
    y_min_all = min(y_min_all, min(y));
    y_max_all = max(y_max_all, max(y));
end

% Create one figure with three separate subplots
figure;

% Plot for large_block_1_selected_indices
subplot(3, 1, 1);
hold on;
loop_num = 1;
for j = 1:length(large_block_1_selected_indices)
    plot_line(loop_num,  colors, choice_times{large_block_1_selected_indices(j)}, filtered_motion, normalized_velocity_bounded, large_block_1_selected_indices(j), 'Large Block 1');
    loop_num = loop_num+1;
end
% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end
hold off;
title('Large Block 1');
xlim([x_min_all, x_max_all]);
ylim([y_min_all, y_max_all]);

% Plot for large_block_2_selected_indices
subplot(3, 1, 2);
hold on;
loop_num = 1;
for j = 1:length(large_block_2_selected_indices)
    plot_line(loop_num,  colors, choice_times{large_block_2_selected_indices(j)}, filtered_motion, normalized_velocity_bounded, large_block_2_selected_indices(j), 'Large Block 2');
    loop_num = loop_num+1;
end
% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end
hold off;
title('Large Block 2');
xlim([x_min_all, x_max_all]);
ylim([y_min_all, y_max_all]);

% Plot for large_block_3_selected_indices
subplot(3, 1, 3);
hold on;
loop_num = 1;
for j = 1:length(large_block_3_selected_indices)
    plot_line(loop_num,  colors, choice_times{large_block_3_selected_indices(j)}, filtered_motion, normalized_velocity_bounded, large_block_3_selected_indices(j), 'Large Block 3');
    loop_num = loop_num+1;
end
% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end
hold off;
title('Large Block 3');
xlim([x_min_all, x_max_all]);
ylim([y_min_all, y_max_all]);






%% Plot random small reward trial(s)

lines_to_plot = 1;
colors = distinguishable_colors(lines_to_plot);


% Randomly select indices
small_block_1_selected_indices = small_block_1_true_indices(randperm(length(small_block_1_true_indices), lines_to_plot ));
small_block_2_selected_indices = small_block_2_true_indices(randperm(length(small_block_2_true_indices), lines_to_plot ));
small_block_3_selected_indices = small_block_3_true_indices(randperm(length(small_block_3_true_indices), lines_to_plot ));




% Initialize variables to store min and max values across all sets of data
x_min_all = inf;
x_max_all = -inf;
y_min_all = inf;
y_max_all = -inf;

% Loop through large_block_1_selected_indices to find min and max values
for j = 1:length(small_block_1_selected_indices)
    x = filtered_motion{small_block_1_selected_indices(j)}(1, :);
    y = filtered_motion{small_block_1_selected_indices(j)}(2, :);

    x_min_all = min(x_min_all, min(x));
    x_max_all = max(x_max_all, max(x));
    y_min_all = min(y_min_all, min(y));
    y_max_all = max(y_max_all, max(y));
end

% Loop through large_block_2_selected_indices to find min and max values
for j = 1:length(small_block_2_selected_indices)
    x = filtered_motion{small_block_2_selected_indices(j)}(1, :);
    y = filtered_motion{small_block_2_selected_indices(j)}(2, :);

    x_min_all = min(x_min_all, min(x));
    x_max_all = max(x_max_all, max(x));
    y_min_all = min(y_min_all, min(y));
    y_max_all = max(y_max_all, max(y));
end

% Loop through large_block_3_selected_indices to find min and max values
for j = 1:length(small_block_3_selected_indices)
    x = filtered_motion{small_block_3_selected_indices(j)}(1, :);
    y = filtered_motion{small_block_3_selected_indices(j)}(2, :);

    x_min_all = min(x_min_all, min(x));
    x_max_all = max(x_max_all, max(x));
    y_min_all = min(y_min_all, min(y));
    y_max_all = max(y_max_all, max(y));
end

% Create one figure with three separate subplots
figure;

% Plot for large_block_1_selected_indices
subplot(3, 1, 1);
hold on;
loop_num = 1;
for j = 1:length(small_block_1_selected_indices)
    plot_line(loop_num,  colors, choice_times{small_block_1_selected_indices(j)}, filtered_motion, normalized_velocity_bounded, small_block_1_selected_indices(j), 'Small Block 1');
    loop_num = loop_num+1;
end
% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end
hold off;
title('Small Block 1');
xlim([x_min_all, x_max_all]);
ylim([y_min_all, y_max_all]);

% Plot for large_block_2_selected_indices
subplot(3, 1, 2);
hold on;
loop_num = 1;
for j = 1:length(small_block_2_selected_indices)
    plot_line(loop_num,  colors, choice_times{small_block_2_selected_indices(j)}, filtered_motion, normalized_velocity_bounded, small_block_2_selected_indices(j), 'Small Block 2');
    loop_num = loop_num+1;
end
% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end
hold off;
title('Small Block 2');
xlim([x_min_all, x_max_all]);
ylim([y_min_all, y_max_all]);

% Plot for large_block_3_selected_indices
subplot(3, 1, 3);
hold on;
loop_num = 1;
for j = 1:length(small_block_3_selected_indices)
    plot_line(loop_num,  colors, choice_times{small_block_3_selected_indices(j)}, filtered_motion, normalized_velocity_bounded, small_block_3_selected_indices(j), 'Small Block 3');
    loop_num = loop_num+1;
end
% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end
hold off;
title('Small Block 3');
xlim([x_min_all, x_max_all]);
ylim([y_min_all, y_max_all]);


%%
% Loop through each line and plot
figure;
for j = 1:size(large_block_1_true_indices , 1) %num_lines
    plot_index = large_block_1_true_indices (j);
    x = [];
    y = [];
    % figure; %comment me out if you want everything plotted on the same fig! 
    
    % Extract X-Y coordinates
    x = filtered_motion{plot_index}(1, :);
    y = filtered_motion{plot_index}(2, :);
    

    % Plot the line using scatter with varying marker colors
    % scatter(x, y, 50, ">", 'filled');
    hold on;
    plot(x, y, 'Color', 'b');
    % % Add a black square at the start of the line
    scatter(x(1), y(1), 200, 'b', 's');

    % Add a black circke at the choice time
    scatter(choice_times{1, large_block_1_true_indices(j)}(1), choice_times{1, large_block_1_true_indices(j)}(2), 200, 'b', 'o');


    % Add a black triangle at the end of the line
    scatter(x(end), y(end), 200, 'b', '^');
end



% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end

hold off;

figure;
for j = 1:size(small_block_1_true_indices , 1) %num_lines
    plot_index = small_block_1_true_indices (j);
    x = [];
    y = [];
    % figure; %comment me out if you want everything plotted on the same fig! 
    
    % Extract X-Y coordinates
    x = filtered_motion{plot_index}(1, :);
    y = filtered_motion{plot_index}(2, :);
    


    % % Plot the line using scatter with varying marker colors
    % scatter(x, y, 50, marker_colors_mapped, 'filled');
    hold on;
    plot(x, y, 'Color', 'r');
    % Add a black square at the start of the line
    scatter(x(1), y(1), 200, 'r', 's');

    % Add a black circke at the choice time
    scatter(choice_times{1, small_block_1_true_indices(j)}(1), choice_times{1, small_block_1_true_indices(j)}(2), 200, 'r', 'o');


    % Add a black triangle at the end of the line
    scatter(x(end), y(end), 200, 'r', '^');
end
% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end

hold off



%%
% Loop through each line and plot
figure;
for j = 1:size(large_block_3_true_indices , 1) %num_lines
    plot_index = large_block_3_true_indices (j);
    x = [];
    y = [];
    % figure; %comment me out if you want everything plotted on the same fig! 
    
    % Extract X-Y coordinates
    x = filtered_motion{plot_index}(1, :);
    y = filtered_motion{plot_index}(2, :);
    

    % Plot the line using scatter with varying marker colors
    % scatter(x, y, 50, ">", 'filled');
    hold on;
    plot(x, y, 'Color', 'b');
    % % Add a black square at the start of the line
    scatter(x(1), y(1), 200, 'b', 's');

    % Add a black circke at the choice time
    scatter(choice_times{1, large_block_3_true_indices(j)}(1), choice_times{1, large_block_3_true_indices(j)}(2), 200, 'b', 'o');


    % Add a black triangle at the end of the line
    scatter(x(end), y(end), 200, 'b', '^');
end



% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end

hold off;


figure

hold on;
for j = 1:size(small_block_3_true_indices , 1) %num_lines
    plot_index = small_block_3_true_indices (j);
    x = [];
    y = [];
    % figure; %comment me out if you want everything plotted on the same fig! 
    
    % Extract X-Y coordinates
    x = filtered_motion{plot_index}(1, :);
    y = filtered_motion{plot_index}(2, :);
    


    % % Plot the line using scatter with varying marker colors
    % scatter(x, y, 50, marker_colors_mapped, 'filled');
    hold on;
    plot(x, y, 'Color', 'r');
    % Add a black square at the start of the line
    scatter(x(1), y(1), 200, 'r', 's');

    % Add a black circke at the choice time
    scatter(choice_times{1, small_block_3_true_indices(j)}(1), choice_times{1, small_block_3_true_indices(j)}(2), 200, 'r', 'o');


    % Add a black triangle at the end of the line
    scatter(x(end), y(end), 200, 'r', '^');
end


% Plot circles and squares from shapeData
for k = 1:numel(shapeData)
    if strcmp(shapeData{k}.Type, 'Circle')
        viscircles(shapeData{k}.Center, shapeData{k}.Radius);
    elseif strcmp(shapeData{k}.Type, 'Square')
        if strcmp(shapeData{k}.Location, 'left screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'b');
        elseif strcmp(shapeData{k}.Location, 'right screen')
            square_center = shapeData{k}.Center;
            square_size = shapeData{k}.Size;
            % square_rotation = shapeData{k}.Rotation;
            rectangle('Position', [square_center - square_size / 2, square_size], 'EdgeColor', 'r');
        end
    end
end

xlim([-20 120])

hold off

