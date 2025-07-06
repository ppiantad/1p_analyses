%%
% Simplified MATLAB code for video playback with behavioral variables
% FIXED VERSION - Proper video recording speed control

%for RM video, use BLA_Insc_25 frames 15000 to 16000
%for AA video, use BLA_Insc_27 frames 24000 to 25000

select_mouse = 'BLA_Insc_27'; 
first_session = 'RDT_D1';

% Load behavioral data
BehavData = final_behavior.(select_mouse).(first_session).uv.BehavData;

% Set frame rate for calcium imaging (adjusted)
calcium_frame_rate = (final.(select_mouse).(first_session).uv.dt)*100;

% Load behavior timing data
stTime = BehavData.TrialPossible(1)-60; 
stTime_to_frames = round(stTime*10);

% Playback speed multiplier (1x = normal speed, 2x = double speed, etc.)
playback_speed = 1;  % Change this value: 1x, 2x, 3x, etc.

% Load the MPG video
mpg_video_file = 'downsampled_video_BLA-Insc-27_RDT_D1.avi';
mpg_video = VideoReader(mpg_video_file);
mpg_num_frames = mpg_video.NumFrames;
mpg_frame_rate = mpg_video.FrameRate;

% Create behavioral variables
choices_only = BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3);
start_only = BehavData.stTime;
consum_only = BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3);

% Check if type_binary column exists and create abort variables
has_type_binary = any(strcmp('type_binary', BehavData.Properties.VariableNames));
if has_type_binary
    % Get abort times for large aborts (type_binary == 1) and small aborts (type_binary == 2)
    large_abort_mask = BehavData.type_binary == 1;
    small_abort_mask = BehavData.type_binary == 2;
    
    % Extract abort times (assuming they correspond to choiceTime)
    large_abort_times = BehavData.choiceTime(large_abort_mask);
    small_abort_times = BehavData.choiceTime(small_abort_mask);
    
    % Remove NaN values
    large_abort_times = large_abort_times(~isnan(large_abort_times));
    small_abort_times = small_abort_times(~isnan(small_abort_times));
end

% Create binary arrays for each behavioral variable
choiceVar = zeros(1, mpg_num_frames);
startVar = zeros(1, mpg_num_frames);
consumVar = zeros(1, mpg_num_frames);

% Convert times to frame indices
choiceTimes_in_frames = round(choices_only * 10);
startTimes_in_frames = round(start_only * 10);
consumTimes_in_frames = round(consum_only * 10);

% Filter valid frame indices
choiceTimes_in_frames = choiceTimes_in_frames(choiceTimes_in_frames <= mpg_num_frames & choiceTimes_in_frames > 0);
startTimes_in_frames = startTimes_in_frames(startTimes_in_frames <= mpg_num_frames & startTimes_in_frames > 0);
consumTimes_in_frames = consumTimes_in_frames(consumTimes_in_frames <= mpg_num_frames & consumTimes_in_frames > 0);

% Mark frames for each behavioral event
choiceVar(choiceTimes_in_frames) = 1;
startVar(startTimes_in_frames) = 1;
consumVar(consumTimes_in_frames) = 1;

% Convert NaNs to zeros
choiceVar(isnan(choiceVar)) = 0;
startVar(isnan(startVar)) = 0;
consumVar(isnan(consumVar)) = 0;

% Process abort variables if type_binary exists
if has_type_binary
    % Create binary arrays for abort variables
    largeAbortVar = zeros(1, mpg_num_frames);
    smallAbortVar = zeros(1, mpg_num_frames);
    
    % Convert abort times to frame indices
    largeAbortTimes_in_frames = round(large_abort_times * 10);
    smallAbortTimes_in_frames = round(small_abort_times * 10);
    
    % Filter valid frame indices for aborts
    largeAbortTimes_in_frames = largeAbortTimes_in_frames(largeAbortTimes_in_frames <= mpg_num_frames & largeAbortTimes_in_frames > 0);
    smallAbortTimes_in_frames = smallAbortTimes_in_frames(smallAbortTimes_in_frames <= mpg_num_frames & smallAbortTimes_in_frames > 0);
    
    % Mark frames for each abort event
    if ~isempty(largeAbortTimes_in_frames)
        largeAbortVar(largeAbortTimes_in_frames) = 1;
    end
    if ~isempty(smallAbortTimes_in_frames)
        smallAbortVar(smallAbortTimes_in_frames) = 1;
    end
    
    % Convert NaNs to zeros for abort variables
    largeAbortVar(isnan(largeAbortVar)) = 0;
    smallAbortVar(isnan(smallAbortVar)) = 0;
end

% Moving window settings
window_size = 400;

% Prompt user for frame range
prompt = {'Enter the start frame:', 'Enter the end frame:'};
dlg_title = 'Frame Range';
num_lines = 1;
defaultans = {num2str(window_size), num2str(mpg_num_frames)};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
start_frame = str2double(answer{1});
end_frame = str2double(answer{2});

% Validate frame range
if isnan(start_frame) || isnan(end_frame) || start_frame < 1 || end_frame > mpg_num_frames || start_frame > end_frame
    error('Invalid frame range. Please enter valid start and end frames.');
end

% Calculate MPG video start frame
mpg_start_frame = max(1, start_frame - stTime_to_frames);

% Limit data based on available frames
max_frames = min([end_frame, mpg_num_frames, length(choiceVar)]);

% Initialize figure
figure('Position', [100, 100, 1200, 800], 'Color', 'black');

% Video subplot (larger, taking up more space)
subplot('Position', [0.1, 0.4, 0.8, 0.55]);
hMpgVideo = imshow(read(mpg_video, 1), []);
set(gca, 'Color', 'black', 'XColor', 'w', 'YColor', 'w');
title('RDT', 'Color', 'w');

% Behavioral traces subplot (smaller, at bottom, same width as video)
subplot('Position', [0.1, 0.05, 0.8, 0.3]);
hold on;

set(gca, 'Color', 'black', 'XColor', 'w', 'YColor', 'w');
xlabel('Time (frames relative to current)', 'Color', 'w');
set(gca, 'YTick', []);
ax1 = gca;
ax1.YAxis.Visible = 'off';   % Remove y-axis

% Create static legend with invisible dummy plots
h1 = plot(nan, nan, 'Color', [0.5 0 0.5], 'LineWidth', 2);  % Purple for trial initiation
h2 = plot(nan, nan, 'Color', [1 0.1 0.6], 'LineWidth', 2);  % Hot pink for choice
h3 = plot(nan, nan, 'Color', [0 0 0.5], 'LineWidth', 2);    % Dark blue for reward collection

% Add abort event legend entries if type_binary exists
if has_type_binary
    h4 = plot(nan, nan, 'Color', [1 1 1], 'LineWidth', 2);      % White for large abort
    h5 = plot(nan, nan, 'Color', [0.5 0.5 0.5], 'LineWidth', 2); % Grey for small abort
    legend([h1, h2, h3, h4, h5], {'Trial initiation', 'Choice', 'Rew. collection', 'LargeAbort', 'SmallAbort'}, ...
           'TextColor', 'w', 'Location', 'northeast', 'AutoUpdate', 'off');
else
    legend([h1, h2, h3], {'Trial initiation', 'Choice', 'Rew. collection'}, ...
           'TextColor', 'w', 'Location', 'northeast', 'AutoUpdate', 'off');
end

% Add middle reference line
hMiddleLine = plot([0, 0], [0, 1], '--', 'Color', 'w', 'LineWidth', 1);

% Current frame line
hCurrentFrameLine = plot([0, 0], [0, 1], '--', 'Color', 'y', 'LineWidth', 2);

% Convert behavioral event times to frame indices for vertical line plotting
choices_frames = round(choices_only * calcium_frame_rate);
start_frames = round(start_only * calcium_frame_rate);
consum_frames = round(consum_only * calcium_frame_rate);

% Convert abort times to frame indices if type_binary exists
if has_type_binary
    large_abort_frames = round(large_abort_times * calcium_frame_rate);
    small_abort_frames = round(small_abort_times * calcium_frame_rate);
end

delay_frames = stTime_to_frames;

% Prepare the video writer with controlled frame rate
video_file = 'behavioral_analysis_output';  % Output video file name
video_writer = VideoWriter(video_file, 'MPEG-4');
video_writer.Quality = 100;

% FIX: Set the output video frame rate based on playback speed
% This ensures the recorded video plays at the correct speed
desired_output_fps = calcium_frame_rate * playback_speed;
video_writer.FrameRate = desired_output_fps;

open(video_writer);

% Initialize counters
mpg_frame_count = mpg_start_frame;
current_frame = start_frame;

% FIX: Remove pause-based timing, use consistent frame writing
% The video speed is now controlled by the VideoWriter's FrameRate property

% Playback loop
try
    while current_frame <= max_frames && mpg_frame_count <= mpg_num_frames
        % Update video
        if mpg_frame_count <= mpg_video.NumFrames
            mpg_frame = read(mpg_video, mpg_frame_count);
            set(hMpgVideo, 'CData', mpg_frame);
        end
        
        % Calculate window for behavioral traces
        window_start = max(1, current_frame - window_size);
        window_end = min(length(choiceVar), current_frame + window_size);
        time_range = (window_start:window_end) - current_frame;  % Relative to current frame
        window_range = window_start:window_end;
        
        % Plot vertical lines at event times within current window
        delete(findobj(gca, 'Type', 'ConstantLine'));  % Clear all previous vertical lines
        
        % Purple lines for trial initiation events
        matching_start_frames = start_frames(start_frames >= window_start & start_frames <= window_end);
        for k = 1:length(matching_start_frames)
            frame_idx = matching_start_frames(k);
            x_pos = frame_idx - current_frame;
            xline(x_pos, 'Color', [0.5 0 0.5], 'LineWidth', 2);  % Purple
        end
        
        % Hot pink lines for choice events
        matching_choice_frames = choices_frames(choices_frames >= window_start & choices_frames <= window_end);
        for k = 1:length(matching_choice_frames)
            frame_idx = matching_choice_frames(k);
            x_pos = frame_idx - current_frame;
            xline(x_pos, 'Color', [1 0.1 0.6], 'LineWidth', 2);  % Hot pink
        end
        
        % Dark blue lines for reward collection events
        matching_consum_frames = consum_frames(consum_frames >= window_start & consum_frames <= window_end);
        for k = 1:length(matching_consum_frames)
            frame_idx = matching_consum_frames(k);
            x_pos = frame_idx - current_frame;
            xline(x_pos, 'Color', [0 0 0.5], 'LineWidth', 2);    % Dark blue
        end
        
        % Add abort event lines if type_binary exists
        if has_type_binary
            % White lines for large abort events
            matching_large_abort_frames = large_abort_frames(large_abort_frames >= window_start & large_abort_frames <= window_end);
            for k = 1:length(matching_large_abort_frames)
                frame_idx = matching_large_abort_frames(k);
                x_pos = frame_idx - current_frame;
                xline(x_pos, 'Color', [1 1 1], 'LineWidth', 2);    % White
            end
            
            % Grey lines for small abort events
            matching_small_abort_frames = small_abort_frames(small_abort_frames >= window_start & small_abort_frames <= window_end);
            for k = 1:length(matching_small_abort_frames)
                frame_idx = matching_small_abort_frames(k);
                x_pos = frame_idx - current_frame;
                xline(x_pos, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);  % Grey
            end
        end
        
        % Update reference lines
        set(hMiddleLine, 'YData', [0, 1]);
        set(hCurrentFrameLine, 'YData', [0, 1]);
        
        % Set axis limits
        xlim([-window_size, window_size]);
        ylim([0, 1]);
        
        % FIX: Keep pause for MATLAB display, but it won't affect video speed
        pause((1 / calcium_frame_rate) / playback_speed);
        
        % Write the current frame to the video
        % The playback speed is now controlled by video_writer.FrameRate
        frame_image = getframe(gcf);
        writeVideo(video_writer, frame_image.cdata);
        
        % Increment counters
        current_frame = current_frame + 1;
        mpg_frame_count = mpg_frame_count + 1;
    end
    
    % Close the video writer
    close(video_writer);
    fprintf('Video saved as: %s.mp4\n', video_file);
    fprintf('Output video frame rate: %.2f fps\n', desired_output_fps);
    
catch ME
    % Make sure to close video writer even if error occurs
    if exist('video_writer', 'var') && isvalid(video_writer)
        close(video_writer);
    end
    disp('An error occurred during playback:');
    disp(ME.message);
end

hold off;