
% Given data
mean_x_pix = final_DLC.C68604.D1Afternoon.movement_data.body_x_pix  ; % x-coordinates in pixels
mean_y_pix = final_DLC.C68604.D1Afternoon.movement_data.body_y_pix; % y-coordinates in pixels
pixels_per_cm = 7.01; % Conversion factor from pixels to cm

% Frame rate of the video (e.g., 30 frames per second)
frame_rate = 30; % Adjust this to your video's frame rate
time_interval = 1 / frame_rate; % Time interval between frames in seconds

% Convert coordinates from pixels to centimeters
mean_x_cm = mean_x_pix / pixels_per_cm;
mean_y_cm = mean_y_pix / pixels_per_cm;

% Calculate the distance traveled between consecutive frames
delta_x = diff(mean_x_cm); % Change in x-coordinates (cm)
delta_y = diff(mean_y_cm); % Change in y-coordinates (cm)
distance_traveled = sqrt(delta_x.^2 + delta_y.^2); % Euclidean distance (cm)

% Calculate velocity (distance / time interval)
velocity = distance_traveled / time_interval; % Velocity in cm/s

% Pad with NaN or 0 to match the original data length (optional)
velocity = [NaN; velocity]; % NaN for the first frame where velocity can't be calculated

% Display or save the result
disp('Velocity (cm/s):');
disp(velocity);

% Plot the velocity over time (optional)
figure;
plot(velocity, 'LineWidth', 1.5);
xlabel('Frame Number');
ylabel('Velocity (cm/s)');
title('Mouse Velocity Over Time');
grid on;
