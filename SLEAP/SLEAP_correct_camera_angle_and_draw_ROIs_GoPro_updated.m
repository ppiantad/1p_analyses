
% LOAD final_SLEAP first. this code will not overwrite data if there
% already exists shapeData etc for a given session


%% Edit these uservariables with what you want to look at
% Define the directory path you want to start with
% startDirectory = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-24';

metaDirectory = 'F:\risk videos\BLA PdCO vs mCherry';
metaDirectory_subfolders = dir(metaDirectory );
metafolder_list = {};
p = get(0, "MonitorPositions");

% Loop through the list of subfolders
for i = 1:length(metaDirectory_subfolders)
    % Check if the item in subfolders is a directory (not "." or "..") or
    % one of the sets of files that I haven't analyzed yet (PR currently)
    if metaDirectory_subfolders(i).isdir && ~strcmp(metaDirectory_subfolders(i).name, '.') && ~strcmp(metaDirectory_subfolders(i).name, '..') && ~contains(metaDirectory_subfolders(i).name, 'PR') && ~contains(metaDirectory_subfolders(i).name, 'not in final dataset')
        % Get the full path of the subfolder
        metasubfolderPath = fullfile(metaDirectory, metaDirectory_subfolders(i).name);
        % Create a cell array for the subfolder path and append it
        % vertically to folder_list
        metafolder_list = vertcat(metafolder_list, {metasubfolderPath});
    end
end

for zz = 1:size(metafolder_list, 1)
    % Use the dir function to get a list of subfolders
    startDirectory = metafolder_list{zz};
    subfolders = dir(startDirectory);
    % Initialize folder_list as an empty cell array
    folder_list = {};
    % Loop through the list of subfolders
    for i = 1:length(subfolders)
        % Check if the item in subfolders is a directory (not "." or "..") or
        % one of the sets of files that I haven't analyzed yet (PR currently)
        if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..') && ~contains(subfolders(i).name, 'PR')
            % Get the full path of the subfolder
            subfolderPath = fullfile(startDirectory, subfolders(i).name);
            % Create a cell array for the subfolder path and append it
            % vertically to folder_list
            folder_list = vertcat(folder_list, {subfolderPath});
        end
    end



    for ii = 1:size(folder_list, 1)
        folder_list_string = strsplit(folder_list{ii}, '\');
        current_animal = folder_list_string{4}; % Would have to change this depending on your folder structure, but there should be an animal name folder given our current workflow.
        current_animal = matlab.lang.makeValidName(current_animal);
        current_session = folder_list_string{5};
        current_session = regexprep(current_session,{' ', '-'}, '_');
        modifiedString = lower(strrep(strrep(folder_list_string{end}, ' ', ''), '-', ''));


        list = dir(folder_list{ii});%grab a directory of the foldercontents
        list_folder_names = {list.name}.';

        folderMask = ~[list.isdir]; %find all of the folders in the directory and remove them from the list
        files = list(folderMask);  %now we have only files to work with
        clear folderMask



        % added these if/elseif statements to try to not have to redo all
        % the labeling if I add mice
        if isfield(final_SLEAP.(current_animal), current_session)
            if isfield(final_SLEAP.(current_animal).(current_session), 'shapeData')
                continue
            elseif ~isfield(final_SLEAP.(current_animal).(current_session), 'shapeData')


                folder_to_analyze = find(strcmpi(strrep(strrep(list_folder_names, ' ', ''), '-', ''), modifiedString));
                disp(['Analyzing subfolder: ' list_folder_names{folder_to_analyze,1}]);
                folder_to_analyze_Path = folder_list{ii};
                list = dir(folder_to_analyze_Path);%grab a directory of the foldercontents
                folderMask = ~[list.isdir]; %find all of the folders in the directory and remove them from the list
                files = list(folderMask);  %now we have only files to work with
                clear folderMask list

                idx = endsWith({files.name},'merged_resized_grayscaled.MP4'); %find the instances of .xlsx in the file list.
                %This command converts the name field into a cell array and searches
                %the cell array with strfind
                movie = files(idx); %build a mat file index
                clear idx files

                if isempty(movie)
                    disp('Missing behavior video file, skipping folder');
                    continue
                else
                    % Specify the path to your video file
                    videoPath = fullfile(movie.folder, movie.name);

                    %% Get random frame from video (if it is dark - re-run this line), draw a line at food cup which will be 4.6 cm
                    % Create a VideoReader object
                    videoObj = VideoReader(videoPath);

                    % Define brightness threshold
                    brightnessThreshold = 30; % Adjust this threshold as needed


                    randomFrameIndex = randi([1800, round(videoObj.NumFrames, -1)/2]);
                    randomFrame = read(videoObj, randomFrameIndex);

                    % while true
                    %     % Read a random frame from the video. get frame from only
                    %     % the 1st half of the vid to avoid getting a frame from
                    %     % after when the mouse was taken out or something odd
                    %     randomFrameIndex = randi([1800, round(videoObj.NumFrames, -1)/2]);
                    %     randomFrame = read(videoObj, randomFrameIndex);
                    %
                    %     % Calculate average brightness of the frame
                    %     avgBrightness = mean(randomFrame(:));
                    %
                    %     % Check if frame is bright enough
                    %     if avgBrightness >= brightnessThreshold
                    %         break; % Exit the loop if frame is bright enough
                    %     else
                    %         disp('Selected frame is too dark. Selecting a different frame...');
                    %     end
                    % end

                    % Display the image
                    f = figure;
                    f.Position = p(2, :); % second display
                    imshow(randomFrame);

                    title('Draw a line on a known distance in the image');
                    % Prompt the user to input the real-world length of the drawn line
                    prompt = 'Enter the real-world length of the drawn line (e.g., in cm): '; %4.6 cm at feeder

                    % Allow the user to draw a line
                    lineObj = imline;
                    position = wait(lineObj);


                    realWorldLength = input(prompt);

                    % Get the pixel length of the drawn line
                    pixelLength = norm(position(2, :) - position(1, :));

                    % Close the figure
                    close;
                    % Display the image
                    f = figure;
                    f.Position = p(2, :); % second display
                    imshow(randomFrame);
                    title('Drag and drop circles for reward_receptacle, left_screen, and right_screen');

                    % Initialize arrays to store shape data
                    shapeData = cell(3, 1);

                    % Loop to create and extract data for each shape
                    for i = 1:3
                        while true
                            if i == 1
                                % Prompt user to drag and drop a circle
                                h = drawcircle;
                                % Prompt user to type "yes" when done with the current shape
                                userResponse = input(['Type "yes" when done with circle at food cup: '], 's');
                                if strcmpi(userResponse, 'yes')
                                    % Extract data for the circle
                                    shapeData{i}.Type = 'Circle';
                                    shapeData{i}.Location = 'reward';
                                    shapeData{i}.Center = h.Center;
                                    shapeData{i}.Radius = h.Radius;
                                    shapeData{i}.BoundingBox = [h.Center - h.Radius, 2 * h.Radius, 2 * h.Radius];
                                    % Delete the current shape to allow drawing the next one
                                    % delete(h);
                                    break; % Exit the loop if 'yes' is entered
                                else
                                    % Delete the current shape to allow redrawing
                                    delete(h);
                                end
                            elseif i == 2
                                % Draw a square
                                h_square_first = drawrectangle('Rotatable', true);
                                % Prompt user to type "yes" when done with the current square
                                userResponse = input(['Type "yes" when done with square at left screen: '], 's');
                                if strcmpi(userResponse, 'yes')
                                    square_center = h_square_first.Position(1:2) + h_square_first.Position(3:4) / 2;
                                    square_size = h_square_first.Position(3:4);
                                    % Store square data
                                    shapeData{i}.Type = 'Square';
                                    shapeData{i}.Location = 'left screen';
                                    shapeData{i}.Center = square_center;
                                    shapeData{i}.Size = square_size;
                                    shapeData{i}.BoundingBox = [h_square_first.Position(1:2), h_square_first.Position(3:4)];
                                    % Delete the square after drawing
                                    % delete(h_square);
                                    break; % Exit the loop if 'yes' is entered
                                else
                                    % Delete the square to allow redrawing
                                    delete(h_square_first);
                                end
                            elseif i == 3
                                % Draw a square
                                h_square_second = drawrectangle('Position', h_square_first.Position, 'RotationAngle',h_square_first.RotationAngle);
                                % h_square = drawrectangle('Rotatable', true);
                                % Prompt user to type "yes" when done with the current square
                                userResponse = input(['Type "yes" when done with square at right screen: '], 's');
                                if strcmpi(userResponse, 'yes')
                                    square_center = h_square_second.Position(1:2) + h_square_second.Position(3:4) / 2;
                                    square_size = h_square_second.Position(3:4);
                                    % Store square data
                                    shapeData{i}.Type = 'Square';
                                    shapeData{i}.Location = 'right screen';
                                    shapeData{i}.Center = square_center;
                                    shapeData{i}.Size = square_size;
                                    shapeData{i}.BoundingBox = [h_square_second.Position(1:2),h_square_second.Position(3:4)];
                                    % Delete the square after drawing
                                    % delete(h_square);
                                    break; % Exit the loop if 'yes' is entered
                                else
                                    % Delete the square to allow redrawing
                                    delete(h_square_second);
                                end
                            end

                            % Check user response to decide whether to proceed
                            if ~strcmpi(userResponse, 'yes')
                                disp('User did not type "yes". Redrawing the shape.');
                            end
                        end
                    end

                    close
                end
                % Correct for slight angle of camera

                % User-defined variables
                cameraAngleDegrees = 30;  % Replace with your measured camera angle in degrees
                focalLength = 8;         % Replace with your camera's focal length in millimeters 50
                % arenaWidth = 100;         % Replace with your arena width in the same units as X coordinates
                % arenaHeight = 75;         % Replace with your arena height in the same units as Y coordinates
                % Convert camera angle to radians
                cameraAngleRadians = deg2rad(cameraAngleDegrees);

                % Define transformation matrix for rotation
                rotationMatrix = [cos(cameraAngleRadians), -sin(cameraAngleRadians);
                    sin(cameraAngleRadians), cos(cameraAngleRadians)];

                % Apply correction for focal length and pixel-to-real-world conversion
                scaleFactor = realWorldLength / pixelLength;


                SLEAP_data_raw = final_SLEAP.(current_animal).(current_session).SLEAP_data_raw;
                SLEAP_data_downsampled = final_SLEAP.(current_animal).(current_session).SLEAP_data;





                %
                %
                %
                %
                %
                %
                %
                %
                %
                % % Load your X and Y data from the table (replace this with your actual table)
                % X_raw = SLEAP_data_raw.x_pix;
                % Y_raw = SLEAP_data_raw.y_pix;
                %
                %
                % % Apply rotation to X and Y coordinates
                % rotatedCoordinates = rotationMatrix * [X_raw'; Y_raw'];
                %
                %
                %
                % correctedX = (rotatedCoordinates(1, :) * (focalLength * scaleFactor))';
                % correctedY = (rotatedCoordinates(2, :) * (focalLength * scaleFactor))';
                %
                % final_SLEAP.(current_animal).(current_session).SLEAP_data_raw.corrected_x_pix = correctedX;
                % final_SLEAP.(current_animal).(current_session).SLEAP_data_raw.corrected_y_pix = correctedY;
                %
                % % Load your X and Y data from the table (replace this with your actual table)
                % X_downsampled = SLEAP_data_downsampled.x_pix;
                % Y_downsampled = SLEAP_data_downsampled.y_pix;
                %
                %
                % % Apply rotation to X and Y coordinates
                % rotatedCoordinates = rotationMatrix * [X_downsampled'; Y_downsampled'];
                %
                %
                %
                % correctedX_downsampled = (rotatedCoordinates(1, :) * (focalLength * scaleFactor))';
                % correctedY_downsampled = (rotatedCoordinates(2, :) * (focalLength * scaleFactor))';
                %
                % final_SLEAP.(current_animal).(current_session).SLEAP_data.corrected_x_pix = correctedX_downsampled;
                % final_SLEAP.(current_animal).(current_session).SLEAP_data.corrected_y_pix = correctedY_downsampled;
                %
                %
                % % Iterate through each shape in shapeData
                % for i = 1:numel(shapeData)
                %     if strcmp(shapeData{i}.Type, 'Circle')
                %         % Scale the Center coordinates separately for X and Y
                %         scaledCenterX = shapeData{i}.Center(1) * (scaleFactor * focalLength);
                %         scaledCenterY = shapeData{i}.Center(2) * (scaleFactor * focalLength);
                %         % Apply rotation to the scaled Center coordinates (for circle, no rotation needed)
                %         rotatedCenter = rotationMatrix * [scaledCenterX; scaledCenterY];
                %         % Update the Center coordinates in shapeData
                %         shapeData{i}.Center = rotatedCenter';
                %         % Scale the Radius and correct for focal length
                %         shapeData{i}.Radius = shapeData{i}.Radius * (scaleFactor * focalLength);
                %         shapeData{i}.BoundingBox = shapeData{i}.BoundingBox * (scaleFactor * focalLength);
                %     else % For squares
                %         % Apply rotation to the square centers
                %         % Scale the center coordinates separately for X and Y
                %         scaledCenterX = shapeData{i}.Center(1) * (scaleFactor * focalLength);
                %         scaledCenterY = shapeData{i}.Center(2) * (scaleFactor * focalLength);
                %         % Apply rotation to the scaled center coordinates
                %         rotatedCenter = rotationMatrix * [scaledCenterX; scaledCenterY];
                %         % Update the center coordinates in shapeData
                %         shapeData{i}.Center = rotatedCenter';
                %         % Scale the size of squares and correct for focal length
                %         shapeData{i}.Size = shapeData{i}.Size * (scaleFactor * focalLength);
                %         shapeData{i}.BoundingBox = shapeData{i}.BoundingBox * (scaleFactor * focalLength);
                %     end
                % end
                final_SLEAP.(current_animal).(current_session).shapeData = shapeData;



                close all
            end
        end
    end
end