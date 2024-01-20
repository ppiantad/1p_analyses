%% 
animalIDs = (fieldnames(final_SLEAP));

select_mouse = 'BLA_Insc_29';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

session_to_analyze = 'RDT_D1';

onset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.stTime'; 
offset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.collectionTime';
fs_cam = 10; %set sampling rate according to camera, this is hard coded for now
time_ranges_trials = [onset_trials; offset_trials];

% gcamp_samples = 1:1:size(Y_dF_all_session, 2);

% gcamp_time = (0:length(F405_downsampled_data)-1)/fs_cam;


velocity_data = final_SLEAP.(select_mouse).(session_to_analyze).zscored_SLEAP_data_velocity';

SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data;


BehavData = final_SLEAP.(select_mouse).(session_to_analyze).BehavData;
adjusted_start_time = BehavData.TrialPossible(1)-60;
% SLEAP_data.idx_time = SLEAP_data.idx_time-adjusted_start_time;
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
        offset = round(time_ranges_trials(2,j)*fs_cam)+1;

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
                filtered_motion{j}= [SLEAP_data.x_pix(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(2,j))'; SLEAP_data.y_pix(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(2,j))']; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
                filtered_velocity{j}= velocity_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(2,j));
                % filtered_gcamp{j}= Y_dF_all_session(gcamp_samples(onset:offset));
                % filtered{j} = SLEAP_data.y_pix_filtered(SLEAP_data.idx_frame(onset:offset))'; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
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

% for i = 1:size(filtered_gcamp,2)
% %     BL_shifted(pp,:)=[BASELINE_PER(1)+time2EPOC(i) BASELINE_PER(2)+time2EPOC(i)]; %BL_shifted(pp,:)=[BASELINE_PER(1)+(-1*time2Collect(i)) BASELINE_PER(2)+time2Collect(i)];
% %     ind = ts2(1,:) < BL_shifted(pp,2) & ts2(1,:) > BL_shifted(pp,1);
%     ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
% 
%     %use if you want to take the Z-score using the entire window mean
%     zb_gcamp = mean(filtered_gcamp{i}); % baseline period mean
%     zsd_gcamp = std(filtered_gcamp{i}); % baseline period stdev
%     zb_motion = mean(filtered_velocity{i}); % baseline period mean
%     zsd_motion = std(filtered_velocity{i}); % baseline period stdev    
%     %use if you want to calculate the Z-score using your specified baseline
% %     zb = mean(Y_dF_all(i,ind)); % baseline period mean
% %     zbmedian = median(Y_dF_all(i,length(ts1)));
% %     zsd = std(Y_dF_all(i,ind)); % baseline period stdev
%     array_sz = size(filtered_gcamp{i}, 2);
% 
% 
%     tmp = 0;
%     pp=pp+1;
%     for j = 1:array_sz % Z score per bin
%         tmp = tmp + 1;
%         zall_gcamp{i}(1, tmp)=(filtered_gcamp{1,i}(j) - zb_gcamp)/zsd_gcamp;
%         % zall_velocity{i}(1,tmp)=(filtered_velocity{1,i}(j) - zb_motion)/zsd_motion;
%     end
%     tmp=0;
% 
%     original_array = filtered_velocity{1,i};
%     if ~isempty(original_array)
%         % Assuming filtered_gcamp{1,1} is your array
% 
%         % Normalize the array to be between 0 and 1
%         min_value = min(original_array(:));
%         max_value = max(original_array(:));
%         normalized_velocity_bounded{i} = (original_array - min_value) / (max_value - min_value);
%         clear original_array min_value max_value
%     elseif isempty(original_array)
%         normalized_velocity_bounded{i} = [];
%         disp('The array is empty. Skipping normalization.');
%     end
% end

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


% Find the indices where large_block_1_ind is true
large_block_1_true_indices = find(large_block_1_ind);
large_block_2_true_indices = find(large_block_2_ind);
large_block_3_true_indices = find(large_block_3_ind);





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

    % Add a black triangle at the end of the line
    scatter(x(end), y(end), 200, 'k', '^', 'filled');
end



%%


% Randomly select 5 indices
large_block_1_selected_indices = large_block_1_true_indices(randperm(length(large_block_1_true_indices), 2));
large_block_2_selected_indices = large_block_2_true_indices(randperm(length(large_block_2_true_indices), 2));
large_block_3_selected_indices = large_block_3_true_indices(randperm(length(large_block_3_true_indices), 2));
% Example data
num_lines = size(filtered_motion, 2); % Number of lines



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
for j = 1:length(large_block_1_selected_indices)
    plot_line(filtered_motion, normalized_velocity_bounded, large_block_1_selected_indices(j), 'Large Block 1');
end
hold off;
title('Large Block 1');
xlim([x_min_all, x_max_all]);
ylim([y_min_all, y_max_all]);

% Plot for large_block_2_selected_indices
subplot(3, 1, 2);
hold on;
for j = 1:length(large_block_2_selected_indices)
    plot_line(filtered_motion, normalized_velocity_bounded, large_block_2_selected_indices(j), 'Large Block 2');
end
hold off;
title('Large Block 2');
xlim([x_min_all, x_max_all]);
ylim([y_min_all, y_max_all]);

% Plot for large_block_3_selected_indices
subplot(3, 1, 3);
hold on;
for j = 1:length(large_block_3_selected_indices)
    plot_line(filtered_motion, normalized_velocity_bounded, large_block_3_selected_indices(j), 'Large Block 3');
end
hold off;
title('Large Block 3');
xlim([x_min_all, x_max_all]);
ylim([y_min_all, y_max_all]);


%%
% Specify the path to your video file
videoPath = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-27\RDT D1\RDT D1\BLA-Insc-27_RDT_D1_2023-01-02T12_28_19.avi';

% Create a VideoReader object
videoObj = VideoReader(videoPath);

% Read a random frame from the video
randomFrameIndex = randi([1, videoObj.NumFrames]);
randomFrame = read(videoObj, randomFrameIndex);

% Create a figure and display the random frame
figure;
imshow(randomFrame);
title('Mouse over the frame to display coordinates');

% Set up a callback function for the mouse-over event
set(gcf, 'WindowButtonMotionFcn', @mouseOverCallback);



%%

% Create a figure for the average plot
figure;

% Compute the mean X-Y coordinates
avg_x = mean(cell2mat(cellfun(@(m) m(1, :), filtered_motion, 'UniformOutput', false)), 2);
avg_y = mean(cell2mat(cellfun(@(m) m(2, :), filtered_motion, 'UniformOutput', false)), 2);


% Smoothed curve for the average gcamp_values
avg_gcamp = zeros(1, max(cellfun(@length, zall_gcamp)));
avg_velocity = zeros(1, max(cellfun(@length, normalized_velocity_bounded)));
for j = 1:num_lines
    gcamp_values = zall_gcamp{j};
    avg_gcamp(1:length(gcamp_values)) = avg_gcamp(1:length(gcamp_values)) + gcamp_values;
    
    velocity_values = normalized_velocity_bounded{j};
    avg_velocity(1:length(velocity_values)) = avg_velocity(1:length(velocity_values)) + velocity_values;
end
avg_gcamp = avg_gcamp / num_lines;
avg_velocity = avg_velocity / num_lines;

% Plot the scatter plot for the average data
scatter(avg_x, avg_y, avg_velocity * 100, avg_gcamp, 'filled');

% Customize the average plot
colormap('jet');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('Average Scatter Plot with Color based on avg\_gcamp, Size based on avg\_velocity');
grid on;