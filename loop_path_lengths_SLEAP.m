animalIDs = (fieldnames(final_SLEAP));

% Step 3: Construct the base directory path
base_dir = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA';
dir_list = dir(base_dir);

for zz = 1:size(fieldnames(final_SLEAP), 1)


    select_mouse = animalIDs{zz};

    % select_mouse_index = find(strcmp(animalIDs, select_mouse));

    % Step 1: Extract substrings from select_mouse_index
    mouse_parts = strsplit(select_mouse, '_');


    session_to_analyze = 'Pre_RDT_RM';

    % Step 2: Extract substrings from session_to_analyze
    session_parts = strsplit(session_to_analyze, '_');
    for i = 1:numel(dir_list)
        if dir_list(i).isdir && ~strcmp(dir_list(i).name, '.') && ~strcmp(dir_list(i).name, '..')
            folder_name = dir_list(i).name;

            % Step 5: Check if folder_name contains all parts of select_mouse_index
            if all(ismember(mouse_parts, strsplit(folder_name, {'_', '-', ' '})))
                % Step 6: Inside this folder, search for subdirectories containing session_parts
                sub_dir = fullfile(base_dir, folder_name);
                sub_dir_list = dir(sub_dir);
                for j = 1:numel(sub_dir_list)
                    if sub_dir_list(j).isdir && ~strcmp(sub_dir_list(j).name, '.') && ~strcmp(sub_dir_list(j).name, '..')
                        sub_folder_name = sub_dir_list(j).name;

                        % Step 7: Check if sub_folder_name contains all parts of session_parts
                        if all(ismember(session_parts, strsplit(sub_folder_name, {'_', '-', ' '})))
                            % Step 8: Get the full path of the subfolder
                            sub_sub_dir = fullfile(sub_dir, sub_folder_name);

                            % Step 9: Loop through subfolders inside subfolder
                            sub_sub_dir_list = dir(sub_sub_dir);
                            for k = 1:numel(sub_sub_dir_list)
                                if sub_sub_dir_list(k).isdir && ~strcmp(sub_sub_dir_list(k).name, '.') && ~strcmp(sub_sub_dir_list(k).name, '..')
                                    sub_sub_folder_name = sub_sub_dir_list(k).name;

                                    % Step 10: Check if sub_sub_folder_name contains all parts of session_parts
                                    if all(ismember(lower(session_parts), lower(strsplit(sub_sub_folder_name, {'_', '-', ' '}))))
                                        % Step 11: Check if the sub_sub_folder contains .avi and .slp files
                                        avi_file = fullfile(sub_sub_dir, sub_sub_folder_name, '*.avi');
                                        slp_file = fullfile(sub_sub_dir, sub_sub_folder_name, '*.slp');
                                        avi_files = dir(avi_file);
                                        slp_files = dir(slp_file);
                                        if ~isempty(avi_files) && ~isempty(slp_files)
                                            disp(['Found matching folder: ', fullfile(sub_sub_dir, sub_sub_folder_name)]);

                                            % Step 12: Loop through .avi files and store filenames and full file directories
                                            avi_filenames = cell(1, numel(avi_files));
                                            avi_full_paths = cell(1, numel(avi_files));
                                            for l = 1:numel(avi_files)
                                                avi_filenames{l} = avi_files(l).name;
                                                avi_full_paths{l} = fullfile(sub_sub_dir, sub_sub_folder_name, avi_filenames{l});
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    
    % Specify the path to your video file
    videoPath = string(avi_full_paths);

    onset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.stTime';
    choice_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.choiceTime';
    offset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.collectionTime';
    fs_cam = 30; %set sampling rate according to camera, this is hard coded for now
    time_ranges_trials = [onset_trials; choice_trials; offset_trials];


    % gcamp_samples = 1:1:size(Y_dF_all_session, 2);

    % gcamp_time = (0:length(F405_downsampled_data)-1)/fs_cam;

    SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;

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

    %% Get random frame from video (if it is dark - re-run this line), draw a line at food cup which will be 4.6 cm
    % Create a VideoReader object
    videoObj = VideoReader(videoPath);

    % Read a random frame from the video
    randomFrameIndex = randi([1, videoObj.NumFrames]);
    randomFrame = read(videoObj, randomFrameIndex);

    % Display the image
    figure;
    imshow(randomFrame);
    title('Draw a line on a known distance in the image');

    % Allow the user to draw a line
    lineObj = imline;
    position = wait(lineObj);

    % Prompt the user to input the real-world length of the drawn line
    prompt = 'Enter the real-world length of the drawn line (e.g., in cm): '; %4.6 cm at feeder
    realWorldLength = input(prompt);

    % Get the pixel length of the drawn line
    pixelLength(zz) = norm(position(2, :) - position(1, :));

    % Close the figure
    close;


    %% Correct for slight angle of camera

    % User-defined variables
    cameraAngleDegrees = 30;  % Replace with your measured camera angle in degrees
    focalLength = 8;         % Replace with your camera's focal length in millimeters 50
    % arenaWidth = 100;         % Replace with your arena width in the same units as X coordinates
    % arenaHeight = 75;         % Replace with your arena height in the same units as Y coordinates

    % Load your X and Y data from the table (replace this with your actual table)
    X = SLEAP_data.x_pix;
    Y = SLEAP_data.y_pix;

    % Convert camera angle to radians
    cameraAngleRadians = deg2rad(cameraAngleDegrees);

    % Define transformation matrix for rotation
    rotationMatrix = [cos(cameraAngleRadians), -sin(cameraAngleRadians);
        sin(cameraAngleRadians), cos(cameraAngleRadians)];

    % Apply rotation to X and Y coordinates
    rotatedCoordinates = rotationMatrix * [X'; Y'];


    % Apply correction for focal length and pixel-to-real-world conversion
    scaleFactor = realWorldLength / pixelLength;
    correctedX = (rotatedCoordinates(1, :) * (focalLength * scaleFactor))';
    correctedY = (rotatedCoordinates(2, :) * (focalLength * scaleFactor))';

    SLEAP_data.x_pix = correctedX;
    SLEAP_data.y_pix = correctedY;

    % % Display the corrected X and Y coordinates
    % disp('Corrected X Coordinates:');
    % disp(correctedX');
    %
    % disp('Corrected Y Coordinates:');
    % disp(correctedY');


    %% Draw relevant touchscreen elements on video frame that was generated above
    % Display the image
    figure;
    imshow(randomFrame);
    title('Drag and drop circles for reward_receptacle, left_screen, and right_screen');

    % Initialize arrays to store circle data
    circleData = cell(3, 1);

    % Loop to create and extract data for each circle
    for i = 1:3
        % Prompt user to drag and drop a circle
        h = drawcircle;

        % Prompt user to type "yes" when done with the current circle
        userResponse = input(['Type "yes" when done with ' num2str(i) '-th circle: '], 's');

        % Extract data for the current circle
        circleData{i}.Center = h.Center;
        circleData{i}.Radius = h.Radius;

        % Delete the current circle to allow drawing the next one
        delete(h);


        % Check user response to decide whether to proceed
        if ~strcmpi(userResponse, 'yes')
            disp('User did not type "yes". Exiting.');
            break;
        end
    end


    % Iterate through each circle in circleData
    for i = 1:numel(circleData)
        % Scale the Center coordinates separately for X and Y
        scaledCenterX = circleData{i}.Center(1) * (scaleFactor * focalLength);
        scaledCenterY = circleData{i}.Center(2) * (scaleFactor * focalLength);

        % Apply rotation to the scaled Center coordinates
        rotatedCenter = rotationMatrix * [scaledCenterX; scaledCenterY];

        % Update the Center coordinates in circleData
        circleData{i}.Center = rotatedCenter';

        % Scale the Radius and correct for focal length
        circleData{i}.Radius = circleData{i}.Radius * (scaleFactor * focalLength);
    end

    % Display the extracted data
    reward_receptacle = circleData{1};
    left_screen = circleData{2};
    right_screen = circleData{3};

    disp('Data for reward_receptacle:');
    disp(reward_receptacle);

    disp('Data for left_screen:');
    disp(left_screen);

    disp('Data for right_screen:');
    disp(right_screen);
    disp('X Y data for right_screen:');
    disp(right_screen);
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
                    filtered_motion{zz, j}= [SLEAP_data.x_pix(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j))'; SLEAP_data.y_pix(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j))']; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
                    choice_times{zz, j} = [SLEAP_data.x_pix(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'; SLEAP_data.y_pix(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'];
                    filtered_velocity{zz, j}= velocity_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j));
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

    %%
    % Combine X and Y coordinates into a 2-column matrix
    % coordinates = [X', Y'];

    % Calculate pairwise distances between consecutive points
    for qq = 1:size(filtered_motion, 2)
        coordinates = filtered_motion{zz, qq}';
        distances = pdist2(coordinates, coordinates);

        % Sum up the distances to get the total path length
        path_length = sum(diag(distances, 1));

        disp(['Path Length: ', num2str(path_length)]);

        distances_matrix{zz, qq} = distances;
        path_length_array{zz}(qq) = path_length;
        clear coordinates distances path_length
    end

end