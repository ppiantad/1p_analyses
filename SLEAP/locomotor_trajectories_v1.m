%% Need to run PhotometryAnalysis_single_mouse first!! 

onset_trials = BehavData.stTime'; 
offset_trials = BehavData.collectionTime';
fs_cam = 30; %set sampling rate according to camera, this is hard coded for now
time_ranges_trials = [onset_trials; offset_trials];

gcamp_samples = 1:1:size(Y_dF_all_session, 2);

gcamp_time = (0:length(F405_downsampled_data)-1)/fs_cam;


% idx_start = gcamp_time >= timeStart;
% trim_470 = F465_downsampled_data(idx_start);
% F465_downsampled_data = [];
% F465_downsampled_data = trim_470;
% trim_405 = F405_downsampled_data(idx_start);
% F405_downsampled_data = [];
% F405_downsampled_data = trim_405;
% trim_SLEAP_velocity = SLEAP_data_vel_filtered_session(idx_start);
% SLEAP_data_vel_filtered_session = [];
% SLEAP_data_vel_filtered_session = trim_SLEAP_velocity;
% 
% 
% time2 = gcamp_time(idx_start)-timeStart;
% % find last valid timestamp (sometimes mice omit final bunch of trials
% % leaving no valid collectionTime
% if BehavData.collectionTime(end,:) ~= 0
%     idx_end = time2 <= (BehavData.collectionTime(end,:)+10);
% elseif BehavData.collectionTime(end,:) == 0
%     idx_end = time2 <= (BehavData.choiceTime(end,:)+10);
% end
% 
% 
% trim_470 = F465_downsampled_data(idx_end);
% F465_downsampled_data = [];
% F465_downsampled_data = trim_470;
% trim_405 = F405_downsampled_data(idx_end);
% F405_downsampled_data = [];
% F405_downsampled_data = trim_405;
% trim_SLEAP_velocity = SLEAP_data_vel_filtered_session(idx_end);
% SLEAP_data_vel_filtered_session = [];
% SLEAP_data_vel_filtered_session = trim_SLEAP_velocity;
% 
% 
% % test3 = test1(idx_end);
% time3 = time2(idx_end);
% 
% 
% 
% % fit 405 to 465 to detrend signal (Y_dF_all)
% bls = polyfit(F405_downsampled_data(1:end), F465_downsampled_data(1:end), 1); %polyfit(F465(1:end), F405(1:end), 1);
% Y_fit_all = bls(1) .* F405_downsampled_data + bls(2);
% Y_dF_all = F465_downsampled_data - Y_fit_all;
% Y_dF_all_delta_F_over_F = Y_dF_all./Y_fit_all;



% gcamp_y=(neuron.C(neuron_num_to_model,:)-mean(neuron.C(neuron_num_to_model,:)))./std(neuron.C(neuron_num_to_model,:));
gcamp_normalized = ((Y_dF_all_session)-mean(Y_dF_all_session))/std(Y_dF_all_session);
SLEAP_data_vel_filtered_session_normalized = ((SLEAP_data_vel_filtered_session)-mean(SLEAP_data_vel_filtered_session))/std(SLEAP_data_vel_filtered_session);


trial_starts_array = BehavData.stTime-BehavData.choiceTime;
trial_ends_array = BehavData.collectionTime - BehavData.choiceTime;

%%
% FILTER ALL EXISTING DATA ON THESE TIME RANGES
% filter streams
if ~isempty(SLEAP_data)
    filtered_motion = [];
    max_ind = max(size(SLEAP_data));
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
                filtered_motion{j}= [SLEAP_data.x_pix_filtered(SLEAP_data.idx_frame(onset:offset))'; SLEAP_data.y_pix_filtered(SLEAP_data.idx_frame(onset:offset))']; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
                filtered_velocity{j}= SLEAP_data_vel_filtered_session(gcamp_samples(onset:offset));
                filtered_gcamp{j}= Y_dF_all_session(gcamp_samples(onset:offset));
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

for i = 1:size(filtered_gcamp,2)
%     BL_shifted(pp,:)=[BASELINE_PER(1)+time2EPOC(i) BASELINE_PER(2)+time2EPOC(i)]; %BL_shifted(pp,:)=[BASELINE_PER(1)+(-1*time2Collect(i)) BASELINE_PER(2)+time2Collect(i)];
%     ind = ts2(1,:) < BL_shifted(pp,2) & ts2(1,:) > BL_shifted(pp,1);
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
    
    %use if you want to take the Z-score using the entire window mean
    zb_gcamp = mean(filtered_gcamp{i}); % baseline period mean
    zsd_gcamp = std(filtered_gcamp{i}); % baseline period stdev
    zb_motion = mean(filtered_velocity{i}); % baseline period mean
    zsd_motion = std(filtered_velocity{i}); % baseline period stdev    
    %use if you want to calculate the Z-score using your specified baseline
%     zb = mean(Y_dF_all(i,ind)); % baseline period mean
%     zbmedian = median(Y_dF_all(i,length(ts1)));
%     zsd = std(Y_dF_all(i,ind)); % baseline period stdev
    array_sz = size(filtered_gcamp{i}, 2);


    tmp = 0;
    pp=pp+1;
    for j = 1:array_sz % Z score per bin
        tmp = tmp + 1;
        zall_gcamp{i}(1, tmp)=(filtered_gcamp{1,i}(j) - zb_gcamp)/zsd_gcamp;
        % zall_velocity{i}(1,tmp)=(filtered_velocity{1,i}(j) - zb_motion)/zsd_motion;
    end
    tmp=0;
    
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

large_block_1_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 1;
% Example data
num_lines = size(zall_gcamp, 2); % Number of lines

% Create a figure
figure;

% Set the colormap for the heatmap
colormap('jet');





% Loop through each line and plot
for j = 1:num_lines
    figure; %comment me out if you want everything plotted on the same fig! 
    
    % Extract X-Y coordinates
    x = filtered_motion{j}(1, :);
    y = filtered_motion{j}(2, :);
    
    % Extract calcium imaging data
    gcamp_values = zall_gcamp{j};

    % Extract velocity data
    velocity_values = normalized_velocity_bounded{j};

    % % Plot the line using scatter with colors based on gcamp_values
    % scatter(x, y, 50, gcamp_values, 'filled');
    hold on;
    
    % Loop through each point and add scatter points with size based on velocity_values
    for k = 1:length(x)
        if velocity_values(k) == 0 
            velocity_values(k) = 0.0001;
        end
        scatter(x(k), y(k), velocity_values(k) * 100, gcamp_values(k), 'filled');
    end
    % Add a black square at the start of the line
    scatter(x(1), y(1), 200, 'k', 's', 'filled');

    % Add a black triangle at the end of the line
    scatter(x(end), y(end), 200, 'k', '^', 'filled');
end

% Add colorbar for the heatmap
colorbar;

% Customize the plot as needed
xlabel('X-axis');
ylabel('Y-axis');
title('Lines with Heatmap Colors');

% Filter out empty cells before computing color limits
non_empty_cells = ~cellfun('isempty', zall_gcamp);
min_value = min(cellfun(@min, zall_gcamp(non_empty_cells)));
max_value = max(cellfun(@max, zall_gcamp(non_empty_cells)));

% Adjust color limits based on non-empty cells
caxis([min_value, max_value]);


% Show grid if needed
grid on;
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