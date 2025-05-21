% Define the top-level directory
top_level_directory = 'E:\DREADD'; % Update this to your top-level directory

% Initialize a list to store matched folders
matched_folders = {};

% Get the full list of subfolders recursively
all_subfolders = genpath(top_level_directory);

% Split the list into individual folder paths
folders = strsplit(all_subfolders, pathsep);

% Loop through each folder and check for the target subfolders
for i = 1:length(folders)
    current_folder = folders{i};
    if isempty(current_folder)
        continue; % Skip empty entries
    end

    % Check if both "SLEAP_data" and "freeze_vid" exist in the current folder
    if isfolder(fullfile(current_folder, 'SLEAP_data')) && ...
       isfolder(fullfile(current_folder, 'freeze_vid')) && ...
       ~contains(current_folder, 'EPM')

        matched_folders{end+1} = current_folder; 
    end
end


% Display the matched folders
disp('Folders containing both "SLEAP_data" and "freeze_vid":');
disp(matched_folders);

%%
% Loop through each primary subfolder
for p = 1:length(matched_folders)
    % Define the path to the current primary subfolder
    primary_subfolder_path = matched_folders{p};
    file_name_strings = strsplit(primary_subfolder_path, '\');
    % update this depending on folder structure
    % session_to_append = file_name_strings{7}; %
    session_to_append = file_name_strings{4}; %

    SLEAP_folder = fullfile(primary_subfolder_path, 'SLEAP_data');
    freeze_folder = fullfile(primary_subfolder_path, 'freeze_vid');

    % Initialize an empty table to store all data for this subfolder
    SLEAP_raw_data = table();

    % Process SLEAP_data folder
    subfolders = {'body', 'L_ear', 'nose', 'R_ear', 'tail_base'};
    for i = 1:length(subfolders)
        % Define the current subfolder path within SLEAP_data
        sub_subfolder_path = fullfile(SLEAP_folder, subfolders{i});

        % Locate the .csv file in the current subfolder
        csv_file = dir(fullfile(sub_subfolder_path, '*.csv'));

        % Check if a CSV file exists in the subfolder
        if ~isempty(csv_file)
            % Read the CSV file
            data = readtable(fullfile(sub_subfolder_path, csv_file.name));

            % Extract 'x_pix' and 'y_pix' columns
            if all(ismember({'x_pix', 'y_pix'}, data.Properties.VariableNames))
                x_data = data.x_pix;
                y_data = data.y_pix;

                % Rename columns to include the subfolder name
                x_column_name = strcat(subfolders{i}, '_x_pix');
                y_column_name = strcat(subfolders{i}, '_y_pix');

                % Create a table with renamed columns
                temp_table = table(x_data, y_data, 'VariableNames', {x_column_name, y_column_name});

                % Append the new columns to SLEAP_raw_data table
                SLEAP_raw_data = [SLEAP_raw_data, temp_table];
            else
                warning('x_pix or y_pix columns not found in %s', csv_file.name);
            end
        else
            warning('No CSV file found in folder %s', subfolders{i});
        end
    end

    % Process freeze_vid folder
    freeze_csv_file = dir(fullfile(freeze_folder, '*.csv'));

    % Check if a CSV file exists in the freeze_vid folder
    if ~isempty(freeze_csv_file)
        % Read the CSV file from freeze_vid folder
        freeze_data = readtable(fullfile(freeze_folder, freeze_csv_file.name), 'Delimiter', ',');

        % Check if required columns 'frame' and 'was_freezing' are present
        if all(ismember({'frame', 'was_freezing'}, freeze_data.Properties.VariableNames))
            % Extract data from row 2 to the end for 'frame' and 'was_freezing'
            frame_data = freeze_data.frame(2:end);
            was_freezing_data = freeze_data.was_freezing(2:end);

            % Create a table with the extracted data and rename the columns
            freeze_table = table(frame_data, was_freezing_data, 'VariableNames', {'frame', 'was_freezing'});

            % Insert freeze data at the start of SLEAP_raw_data
            SLEAP_raw_data = [freeze_table, SLEAP_raw_data(1:size(freeze_table), :)];
        else
            warning('frame or was_freezing columns not found in %s', freeze_csv_file.name);
        end
    else
        warning('No CSV file found in freeze_vid folder');
    end

    % Calculate mean of "x_pix" and "y_pix" columns and add as new columns
    x_pix_columns = contains(SLEAP_raw_data.Properties.VariableNames, '_x_pix');
    y_pix_columns = contains(SLEAP_raw_data.Properties.VariableNames, '_y_pix');
    mean_x_pix = mean(SLEAP_raw_data{:, x_pix_columns}, 2, 'omitnan'); % mean across rows
    mean_y_pix = mean(SLEAP_raw_data{:, y_pix_columns}, 2, 'omitnan'); % mean across rows
    SLEAP_raw_data.mean_x_pix = mean_x_pix;
    SLEAP_raw_data.mean_y_pix = mean_y_pix;

    % use body point to calculate velocity
    body_x_pix = SLEAP_raw_data.body_x_pix; % x-coordinates in pixels
    body_y_pix = SLEAP_raw_data.body_y_pix; % y-coordinates in pixels
    % this value is deteremined based on the arena - see
    % get_distance_from_pixels
    pixels_per_cm = 7.01; % Conversion factor from pixels to cm

    % Frame rate of the video (e.g., 30 frames per second)
    frame_rate = 30; % Adjust this to your video's frame rate
    time_interval = 1 / frame_rate; % Time interval between frames in seconds

    % Convert coordinates from pixels to centimeters
    mean_x_cm = body_x_pix / pixels_per_cm;
    mean_y_cm = body_y_pix / pixels_per_cm;

    % Calculate the distance traveled between consecutive frames
    delta_x = diff(mean_x_cm); % Change in x-coordinates (cm)
    delta_y = diff(mean_y_cm); % Change in y-coordinates (cm)
    distance_traveled = sqrt(delta_x.^2 + delta_y.^2); % Euclidean distance (cm)

    % Calculate velocity (distance / time interval)
    body_velocity= distance_traveled / time_interval; % Velocity in cm/s

    % Pad with NaN or 0 to match the original data length (optional)
    body_velocity_all = [NaN; body_velocity]; % NaN for the first frame where velocity can't be calculated

    SLEAP_raw_data.body_velocity = body_velocity_all;

    % Define the filename and path for the output .csv file
    output_filename = strcat(session_to_append, '_SLEAP_and_freezing_combined.csv');
    output_file = fullfile(primary_subfolder_path, output_filename);

    % Save the SLEAP_raw_data table as a .csv file
    writetable(SLEAP_raw_data, output_file);

    % Display a message confirming the save
    disp(['SLEAP_raw_data table saved to ', output_file]);
end

