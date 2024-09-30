%% --- Preloaded data ---
clearvars -except mouseData final
%%

BehavData = final.BLA_Insc_24.RDT_D1.uv.BehavData;
velocity_data = final_SLEAP.BLA_Insc_24.RDT_D1.zscored_SLEAP_data_velocity;

ca_data = final.BLA_Insc_24.RDT_D1.CNMFe_data.C_raw;
ca_data_denoised = final.BLA_Insc_24.RDT_D1.CNMFe_data.C;
ca_data_spike_inf = final.BLA_Insc_24.RDT_D1.CNMFe_data.spike_prob;
ca_data_norm = normalize(ca_data, 2);
ca_data_denoised_norm = normalize(ca_data_denoised, 2);
ca_data_spike_inf_norm = normalize(ca_data_spike_inf, 2);

stTime = BehavData.TrialPossible(1)-60; 
stTime_to_frames = round(stTime*10);
m = 1;
behavior = 'open';


% --- Speed factor for playback ---
speed_factor = 1000;  % Adjust this factor to speed up playback

% --- Load the calcium imaging video from the .tiff file ---
tiff_video_file = 'BLA_INSC_24_RDT_D1__2022-08-10-14-07-58_video_green_motion_corrected.tiff';  % Replace with your .tiff file
info = imfinfo(tiff_video_file);
num_frames = numel(info);  % Number of frames in .tiff file
video_height = info(1).Height;
video_width = info(1).Width;

% Load .tiff video data
tiff_video_data = zeros(video_height, video_width, num_frames);
for i = 1:num_frames
    tiff_video_data(:, :, i) = imread(tiff_video_file, i);
end

%%

selected_neurons = [21, 1, 42];  % Modify this to include more neurons 10, 128, 55

% --- Load the .MPG video ---
mpg_video_file = 'downsampled_video_insc_24_RDT_D1.avi';  % Replace with your .MPG file
mpg_video = VideoReader(mpg_video_file);
mpg_num_frames = mpg_video.NumFrames;
mpg_frame_rate = mpg_video.FrameRate;

% Set frame rate for calcium imaging (adjusted)
calcium_frame_rate = (final.BLA_Insc_24.RDT_D1.uv.dt)*100;  % Frame rate for the calcium imaging video


% Get neural data from mouseData

contours = final.BLA_Insc_24.RDT_D1.CNMFe_data.Coor;
    

fluorescence_data = final.BLA_Insc_24.RDT_D1.CNMFe_data.C_raw;
fluorescence_data = normalize(fluorescence_data, 2);

choices_only = BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3);

% Initialize the continuous variable as a zeros array with length of calcium data
choiceVar = zeros(1, num_frames);

% Convert choice times (in seconds) to frame indices (sampling at 10 Hz)
choiceTimes_in_frames = round(choices_only * 10); % Multiply by 10 to get frame indices

% Ensure that frame indices are within the range of calcium data frames
choiceTimes_in_frames = choiceTimes_in_frames(choiceTimes_in_frames <= num_frames);

% Mark the frames where choices were made
choiceVar(choiceTimes_in_frames) = 1;

% choiceVar now contains 1s at the frames where choices were made and 0s elsewhere




% Get behavior data from mouseData
behavTrace = choiceVar;

% Convert NaNs in behavTrace to zeros
behavTrace(isnan(behavTrace)) = 0;

% Number of neurons and timepoints
num_neurons = length(selected_neurons);

% --- Scale behavTrace ---
% Remove NaNs and zeros before scaling
valid_behavTrace = behavTrace(~isnan(behavTrace) & behavTrace > 0);
if isempty(valid_behavTrace)
    error('Behavior trace contains no valid values.');
end

% Get maximum fluorescence value to scale the traces
max_fluorescence = max(fluorescence_data(selected_neurons));

% % Scale behavTrace to match the max fluorescence range
% scaled_behavTrace = behavTrace;
% scaled_behavTrace = scaled_behavTrace - nanmin(scaled_behavTrace);  % Shift to start from 0
% % Scale based on the max value found in valid behavTrace
% max_valid_behav = nanmax(scaled_behavTrace);
% scaled_behavTrace = (scaled_behavTrace / max_valid_behav) * max_fluorescence;  % Scale to max fluorescence value
% 

%%
hRedLine = [];  % Initialize an empty array to store line handles
num_frames = numel(info);  % Number of frames in .tiff file
video_height = info(1).Height;
video_width = info(1).Width;

% Moving window settings (300 frames before and after current frame)
window_size = 400;

% Prompt user for frame range for the calcium video
prompt = {'Enter the start frame:', 'Enter the end frame:'};
dlg_title = 'Frame Range';
num_lines = 1;
defaultans = {num2str(window_size), num2str(num_frames)};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
start_frame = str2double(answer{1});
end_frame = str2double(answer{2});
% scaled_behavTrace_2 = scaled_behavTrace(1, 400:end);  % Scale to max fluorescence value

% Validate frame range
if isnan(start_frame) || isnan(end_frame) || start_frame < 1 || end_frame > num_frames || start_frame > end_frame
    error('Invalid frame range. Please enter valid start and end frames.');
end

% Adjust frame range
frames = start_frame:end_frame;
num_frames = length(frames);

% Limit playback length based on behavTrace
trace_length = length(behavTrace);
max_frames = min([end_frame, mpg_num_frames, trace_length]);

% --- Calculate the MPG video start frame based on calcium imaging start frame ---
mpg_start_frame = round((start_frame / calcium_frame_rate) * mpg_frame_rate)-stTime_to_frames;  % Calculate the equivalent MPG frame
mpg_frame_count = max(1, mpg_start_frame);  % Ensure we start from a valid frame (at least frame 1)

velocity_data = velocity_data(:, mpg_start_frame:end);


% Create color map for neurons
cmap = lines(num_neurons);

% Initialize larger figure
figure('Position', [100, 100, 1200, 1000], 'Color', 'black');  % Adjust figure dimensions for more space

% --- Calcium imaging subplot ---
subplot('Position', [0.05, 0.55, 0.4, 0.4]);  % Left side, aligned to the top with 40% height
hTiffVideo = imshow(tiff_video_data(:, :, start_frame), [0, 500]);  % Set dynamic scaling to display correctly
set(gca, 'Color', 'black', 'XColor', 'w', 'YColor', 'w');  % Black background, white axes
hold on;
hContour = gobjects(length(selected_neurons), 1);
hFrameNumber = text(0.80, 0.95, '', 'Units', 'normalized', 'FontSize', 12, 'Color', 'white', 'FontWeight', 'bold');  % Frame number display
for i = 1:length(selected_neurons)
    neuron_idx = selected_neurons(i);
    hContour(i) = plot(contours{neuron_idx}(1,:), contours{neuron_idx}(end,:), 'LineWidth', 1, 'Color', cmap(i, :));
end
caxis([500 700]);
colormap(gray);  % Set colormap to gray

% --- MPG video subplot ---
subplot('Position', [0.55, 0.55, 0.4, 0.4]);  % Right side, aligned to the top with 40% height
hMpgVideo = imshow(read(mpg_video, 1), []);  % Display the first frame
set(gca, 'Color', 'black', 'XColor', 'w', 'YColor', 'w');  % Black background, white axes
hold on;

% --- Fluorescence and behavior trace plot ---
subplot('Position', [0.05, 0.05, 0.9, 0.45]);  % Bottom side, increased height to 45% for more trace space
hold on;

% Initialize behavior trace with fill() (make sure it covers the whole X range)
hBehavTrace = fill(nan, nan, 'w', 'FaceAlpha', 0.2);  % Initialize behavior trace
hFluorescence = gobjects(length(selected_neurons), 1);  % Initialize fluorescence traces

hVelocity = plot(nan, nan, 'LineWidth', 1, 'Color', cmap(i, :));  % Initialize traces

for i = 1:length(selected_neurons)
    neuron_idx = selected_neurons(i);
    hFluorescence(i) = plot(nan, nan, 'LineWidth', 1, 'Color', cmap(i, :));  % Initialize traces
end
set(gca, 'Color', 'black', 'XColor', 'w', 'YColor', 'w');  % Set black background and white axes
xlabel('Time (s)', 'Color', 'w');
set(gca, 'YTick', []);
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis

% Add a fixed vertical dotted white line in the middle
hMiddleLine = plot([0, 0], ylim, '--', 'Color', 'w');

% Create a new dotted line for the current frame
hCurrentFrameLine = plot([0, 0], [0, 0], '--', 'Color', 'w');  % Red dotted line

% --- Synchronization and loop ---
calcium_time_per_frame = 1 / calcium_frame_rate;  % Time per frame for calcium imaging
mpg_time_per_frame = 1 / mpg_frame_rate;  % Time per frame for MPG video
calcium_frame_count = start_frame;  % Start from the user-defined frame

offset = 10;

% Prepare the video writer
video_file = 'output_video';  % Output video file name
video_writer = VideoWriter(video_file, 'MPEG-4');
video_writer.Quality = 100;
open(video_writer);

% Error handling for video frame limits
try
    % Convert choices_only timestamps (seconds) to frame indices
    choices_frames = round(choices_only * calcium_frame_rate);  % Convert seconds to frame indices

    % Define the delay in terms of frames (144 frame delay)
    delay_frames = stTime_to_frames;  % Assuming stTime_to_frames is 144
    
    % Initializing the starting frame for calcium and mpg videos
    calcium_frame_count = start_frame;  % Assume start_frame is given, e.g., 400
    

    % Loop through frames for playback
    while calcium_frame_count <= max_frames && mpg_frame_count <= max_frames
        % Update calcium imaging video
        set(hTiffVideo, 'CData', tiff_video_data(:, :, calcium_frame_count));


        mpg_frame = read(mpg_video, mpg_frame_count);  % Read current frame from the video
        set(hMpgVideo, 'CData', mpg_frame);  % Update MPG video
        mpg_frame_count = mpg_frame_count + 1;  % Increment mpg video frame count

        % Update frame number display
        set(hFrameNumber, 'String', sprintf('Frame: %d', calcium_frame_count));

        % Calculate window range for the current frame
        window_start = max(1, calcium_frame_count - window_size);
        window_end = min(trace_length, calcium_frame_count + window_size);
        time_range = (window_start:window_end) - calcium_frame_count;  % Time relative to current frame for scrolling

        window_range = window_start:window_end;

        % Update fluorescence traces in the moving window
        max_fluorescence_current = 0;  % Reset for current frame
        for i = 1:length(selected_neurons)
            neuron_idx = selected_neurons(i);
            set(hFluorescence(i), 'XData', time_range, 'YData', fluorescence_data(neuron_idx, window_range) + i * offset);
            max_fluorescence_current = max(max_fluorescence_current, max(fluorescence_data(neuron_idx, window_range) + i * offset));
        end
         % Update velocity plot near the X-axis (remove any offset for velocity)
         set(hVelocity, 'XData', time_range, 'YData', velocity_data(1, window_range));  % No offset for velocity


        
        % Clear any previous red lines (if any exist)
        if exist('hRedLine', 'var')
            delete(hRedLine(:));  % Clear all previous lines
        end

        % Find all matching frames within the current window range
        matching_frames = choices_frames(choices_frames >= window_start & choices_frames <= window_end);

        % Plot vertical red lines at all matching frames within the window
        hRedLine = [];  % Initialize an empty array to store line handles
        if ~isempty(matching_frames)
            hold on;
            for k = 1:length(matching_frames)
                frame_idx = matching_frames(k);  % Current matching frame index
                x_pos = frame_idx - calcium_frame_count;  % Position relative to the moving window
                hRedLine(k) = xline(x_pos, 'r', 'LineWidth', 2);  % Plot a red vertical line for each matching frame
            end
            hold off;
        end

        % Update the current frame line position and Y-limits
        max_fluorescence_current = max(max_fluorescence_current, 0);  % Ensure it's at least 0
        set(hCurrentFrameLine, 'YData', [0, max_fluorescence_current]);  % Extend to current max fluorescence height

        % Update the middle line height
        set(hMiddleLine, 'YData', [0, max_fluorescence_current]);  % Match middle line height to current max

        % Set Y-axis limits to ensure minimum is 0
        ylim([0, max_fluorescence_current]);  % Set Y-limits for the plot

        % Pause to maintain playback speed
        pause(1 / (calcium_frame_rate * speed_factor));  % Speed up the playback

        % Increment calcium frame counter
        calcium_frame_count = calcium_frame_count + 1;

        % Write the current frame to the video
        frame_image = getframe(gcf);
        % writeVideo(video_writer, frame_image.cdata);
    end

    % Close the video writer
    close(video_writer);
    hold off;

catch ME
    disp('An error occurred during playback:');
    disp(ME.message);
end
