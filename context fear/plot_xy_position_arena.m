valid_session_names = {'D1_Morning', 'D1_Afternoon', 'D2_Morning', 'D2_Afternoon', 'D3', 'D4'};

% Example inputs
struct_data = final_DLC.B46852  ; % Replace with your actual struct
% final_DLC.B57422: experimental, ok to plot
% final_DLC.B46852: No Shock, good to plot
% final_DLC.B46857: One Context, good to plot

% Get field names of the top-level struct
fields = fieldnames(struct_data);

% Filter the fields based on valid_session_names
matching_fields = fields(ismember(fields, valid_session_names));

% Initialize variables to store global min and max
global_x_min = inf;
global_x_max = -inf;
global_y_min = inf;
global_y_max = -inf;

% First pass: Calculate global min and max
for i = 1:numel(valid_session_names)
    session_name = valid_session_names{i};
    
    if isfield(struct_data, session_name) && isfield(struct_data.(session_name), 'movement_data')
        movement_data = struct_data.(session_name).movement_data;
        
        % Extract X and Y positions, omitting the first 100 rows
        x_positions = movement_data.mean_x_pix(101:end);
        y_positions = movement_data.mean_y_pix(101:end);
        
        % Update global min and max
        global_x_min = min(global_x_min, min(x_positions));
        global_x_max = max(global_x_max, max(x_positions));
        global_y_min = min(global_y_min, min(y_positions));
        global_y_max = max(global_y_max, max(y_positions));
    end
end

% Prepare figure
num_sessions = numel(matching_fields);
figure;
tiledlayout(1, num_sessions); % Use tiledlayout for 1-row subplots

% Second pass: Plot data with consistent axes
for i = 1:num_sessions
    session_name = valid_session_names{i};
    
    % Check if the session contains the expected field
    if isfield(struct_data, session_name) && isfield(struct_data.(session_name), 'movement_data')
        movement_data = struct_data.(session_name).movement_data;
        
        % Extract X and Y positions, omitting the first 100 rows
        x_positions = movement_data.body_x_pix(101:end);
        y_positions = movement_data.body_y_pix(101:end);
        
        % Create a subplot
        nexttile;
        plot(x_positions, y_positions, 'LineWidth', 1.5);
        title(sprintf('Session: %s', session_name), 'Interpreter', 'none');
        xlabel('X Position (pixels)');
        ylabel('Y Position (pixels)');
        % grid on;
        
        % Set consistent axis limits
        xlim([global_x_min, global_x_max]);
        ylim([global_y_min, global_y_max]);
    else
        warning('No "movement_data" field in session: %s', session_name);
    end
end

% Adjust the layout
sgtitle('Mouse Movement Across Sessions');

%% plot XY position during shocks

session_name = valid_session_names{2};


movement_data = struct_data.(session_name).movement_data;
velocity_data = movement_data.body_velocity;
x_positions = movement_data.body_x_pix;
y_positions = movement_data.body_y_pix;

uv.evtWin = [-2 25]; %what time do y
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



eTS = shk_on'; %get time stamps

%calculate time windows for each event
evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
% Define parameters


% Calculate the total number of samples from body_velocity
num_samples = height(velocity_data); % Assuming body_velocity is a table column

% Generate time array in seconds
time_array = (0:num_samples-1) / frame_rate;

% Convert time array to minutes if needed
time_array_minutes = time_array / 60;
velocity_by_shock = [];
xposition_by_shock = [];
yposition_by_shock = [];




for t = 1:size(eTS,1)
    % set each trial's temporal boundaries
    timeWin = [eTS(t)+uv.evtWin(1,1):1/frame_rate:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
    % BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
    % unitTrace_zscored = zscore(unitTrace);


    if min(timeWin) > min(time_array) && max(timeWin) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
        % get unit event counts in trials
        % get unit ca traces in trials
        idx = time_array >= min(timeWin) & time_array < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);


        velocity_by_shock(t,:) = velocity_data(idx);
        xposition_by_shock(t,:) = x_positions(idx);
        yposition_by_shock(t,:) = y_positions(idx);

    end

end

% smoothing the velocity data a bit to make it less noisy 
velocity_by_shock_filtered = sgolayfilt(velocity_by_shock, 5, 21, [], 2);


time_vector = [];
% Determine time vector for each trajectory
num_points = size(xposition_by_shock, 2);
time_vector = linspace(uv.evtWin(1), uv.evtWin(2), num_points);
%%
% Prepare figure
figure;
hold on;

% Define colormap
cmap = jet(256); % Jet colormap with 256 levels
num_shocks = size(xposition_by_shock, 1);

% Define velocity range for colormap scaling
velocity_range = [0, 50]; % Adjust this range as necessary

% Loop through each row of data
for i = 1:num_shocks
    x_positions = xposition_by_shock(i, :);
    y_positions = yposition_by_shock(i, :);
    velocities = velocity_by_shock_filtered(i, :);

    % Normalize velocities for colormap indexing
    norm_velocities = (velocities - min(velocities)) / (max(velocities) - min(velocities));
    color_indices = round(norm_velocities * (size(cmap, 1) - 1)) + 1; % Map to colormap indices

    % Plot each segment of the line with corresponding color
    for j = 1:(length(x_positions) - 1)
        plot(x_positions(j:j+1), y_positions(j:j+1), 'Color', cmap(color_indices(j), :), 'LineWidth', 1.5);
    end
    % Find indices for shock start and end
    [~, shock_start_idx] = min(abs(time_vector - 0));  % Time = 0
    [~, shock_end_idx] = min(abs(time_vector - 1));   % Time = 1 second

    % Add scatter symbols for shock start and end
    scatter(x_positions(shock_start_idx), y_positions(shock_start_idx), 100, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 1.5); % Shock start (green circle)
    scatter(x_positions(shock_end_idx), y_positions(shock_end_idx), 100, 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5); % Shock end (red square)
end

% Add a colorbar to indicate velocity scale
colormap(cmap);
c = colorbar;
c.Label.String = 'Velocity';

% Set consistent axis limits
xlim([global_x_min, global_x_max]); % Use the global values calculated earlier
ylim([global_y_min, global_y_max]); % Use the global values calculated earlier

% Add labels and title
xlabel('X Position (pixels)');
ylabel('Y Position (pixels)');
title('Mouse Position Colored by Velocity');
grid on;
hold off;


%% Prepare figure
figure;

% Define colormap
cmap = jet(256); % Jet colormap with 256 levels
num_shocks = size(xposition_by_shock, 1);

% Define velocity range for colormap scaling
velocity_range = [0, 50]; % Adjust this range as necessary

% Loop through each row of data
for i = 1:num_shocks
    % Create a subplot for each shock
    subplot(num_shocks, 1, i);
    hold on;
    
    x_positions = xposition_by_shock(i, :);
    y_positions = yposition_by_shock(i, :);
    velocities = velocity_by_shock_filtered(i, :);

    % Normalize velocities for colormap indexing
    norm_velocities = (velocities - min(velocities)) / (max(velocities) - min(velocities));
    color_indices = round(norm_velocities * (size(cmap, 1) - 1)) + 1; % Map to colormap indices

    % Plot each segment of the line with corresponding color
    for j = 1:(length(x_positions) - 1)
        plot(x_positions(j:j+1), y_positions(j:j+1), 'Color', cmap(color_indices(j), :), 'LineWidth', 1.5);
    end

    % Find indices for shock start and end
    [~, shock_start_idx] = min(abs(time_vector - 0));  % Time = 0
    [~, shock_end_idx] = min(abs(time_vector - 1));   % Time = 1 second

    % Add scatter symbols for shock start and end
    scatter(x_positions(shock_start_idx), y_positions(shock_start_idx), 100, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 1.5); % Shock start (green circle)
    scatter(x_positions(shock_end_idx), y_positions(shock_end_idx), 100, 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5); % Shock end (red square)

    % Set consistent axis limits
    xlim([global_x_min, global_x_max]); % Use the global values calculated earlier
    ylim([global_y_min, global_y_max]); % Use the global values calculated earlier

    % Add labels and title
    xlabel('X Position (pixels)');
    ylabel('Y Position (pixels)');
    title(['Shock ', num2str(i)]);
    grid on;
    hold off;
end

% Add a colorbar to indicate velocity scale (global, outside of subplots)
colormap(cmap);
c = colorbar('Position', [0.93 0.1 0.02 0.8]); % Adjust position as necessary
c.Label.String = 'Velocity';

%%
%% Prepare figure
figure;

% Define colormap
cmap = jet(256); % Jet colormap with 256 levels
num_shocks = size(xposition_by_shock, 1);

% Define velocity range for colormap scaling
velocity_range = [0, 200]; % Adjust this range as necessary

% Conversion factor from pixels to centimeters
pixels_per_cm = 7.01;

% Circle properties
radius_cm = 45; % Radius of the circle in cm
circle_diameter = 2 * radius_cm;

% Loop through each row of data
for i = 2:2:num_shocks
    % Create a subplot for each shock
    subplot(num_shocks, 1, i);
    hold on;
    
    % Convert positions to centimeters
    x_positions_cm = xposition_by_shock(i, :) / pixels_per_cm;
    y_positions_cm = yposition_by_shock(i, :) / pixels_per_cm;

    % Restrict velocities to the defined range
    velocities_clipped = max(min(velocities, velocity_range(2)), velocity_range(1));

    % Normalize velocities for colormap indexing
    norm_velocities = (velocities_clipped - velocity_range(1)) / (velocity_range(2) - velocity_range(1));
    color_indices = round(norm_velocities * (size(cmap, 1) - 1)) + 1; % Map to colormap indices


    % Plot circle with a diameter of 90 cm (radius 45 cm)
    center_x = mean(x_positions_cm); % Use mean X position as center
    center_y = mean(y_positions_cm); % Use mean Y position as center
    rectangle('Position', [center_x - radius_cm, center_y - radius_cm, circle_diameter, circle_diameter], ...
              'Curvature', [1, 1], 'EdgeColor', [0.7, 0.7, 0.7], 'LineWidth', 1.5, 'LineStyle', '-');

    % Plot each segment of the line with corresponding color
    for j = 1:(length(x_positions_cm) - 1)
        plot(x_positions_cm(j:j+1), y_positions_cm(j:j+1), 'Color', cmap(color_indices(j), :), 'LineWidth', 1.5);
    end

    % Find indices for shock start and end
    [~, shock_start_idx] = min(abs(time_vector - 0));  % Time = 0
    [~, shock_end_idx] = min(abs(time_vector - 1));   % Time = 1 second

    % Add scatter symbols for shock start and end
    scatter(x_positions_cm(shock_start_idx), y_positions_cm(shock_start_idx), 100, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 1.5); % Shock start (green circle)
    scatter(x_positions_cm(shock_end_idx), y_positions_cm(shock_end_idx), 100, 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5); % Shock end (red square)

    % Adjust axis limits to match the circle's border
    xlim([center_x - radius_cm, center_x + radius_cm]);
    ylim([center_y - radius_cm, center_y + radius_cm]);

    % Ensure the subplot is square and remove axes
    % axis equal;
    axis off;  % Removes the axis ticks and labels
    hold off;
end

% % Add a fixed colorbar to indicate velocity scale
% colormap(cmap);
% c = colorbar('Position', [0.93 0.1 0.02 0.8]); % Set position to prevent resizing
% c.Label.String = 'Velocity';

%%
% Create a separate figure for the colorbar
figure;

% Define colormap
cmap = jet(256); % Jet colormap with 256 levels

% Get the full velocity data (concatenate all velocities from the shocks)
all_velocities = velocity_by_shock_filtered(:); % Flatten the velocity data across all shocks

% Normalize velocities for colormap indexing (same as before)
norm_velocities = (all_velocities - min(all_velocities)) / (max(all_velocities) - min(all_velocities));

% Calculate the min and max of the normalized velocities for colorbar scaling
norm_velocity_min = min(norm_velocities);
norm_velocity_max = max(norm_velocities);

% Create a dummy plot to set the colorbar
imagesc([norm_velocity_min, norm_velocity_max], [0, 1], [0 0 0]); % Dummy data to generate a colorbar
colormap(cmap);

% Add colorbar and set limits based on normalized velocity data
c = colorbar;
c.Label.String = 'Velocity';
caxis([norm_velocity_min, norm_velocity_max]); % Set color axis to the normalized velocity range

% Optionally, save the figure
% saveas(gcf, 'colorbar_figure.png'); % Save as image (if needed)

% You can now manually insert this colorbar figure into your document.

