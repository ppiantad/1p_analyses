
animalIDs = (fieldnames(final_SLEAP));

select_mouse = 'BLA_Insc_26';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

session_to_analyze = 'RDT_D1';

% Specify the path to your video file
videoPath = 'i:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-26\RDT D1\RDT D1\BLA-Insc-26_RDT_D1_2023-01-16T12_20_15.avi';


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
pixelLength = norm(position(2, :) - position(1, :));

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

% Initialize arrays to store shape data
shapeData = cell(3, 1);

% Loop to create and extract data for each shape
for i = 1:3
    if i == 1
        % Prompt user to drag and drop a circle
        h = drawcircle;
        % Prompt user to type "yes" when done with the current shape
        userResponse = input(['Type "yes" when done with circle at food cup: '], 's');
        % Extract data for the circle
        shapeData{i}.Type = 'Circle';
        shapeData{i}.Location = 'reward';
        shapeData{i}.Center = h.Center;
        shapeData{i}.Radius = h.Radius;
        % Delete the current shape to allow drawing the next one
        % delete(h);
    elseif i == 2
        % Draw a square
        h_square = drawrectangle('Rotatable', true);
        % Prompt user to type "yes" when done with the current square
        userResponse = input(['Type "yes" when done with square at left screen: '], 's');
        square_center = h_square.Position(1:2) + h_square.Position(3:4) / 2;
        square_size = h_square.Position(3:4);
        % Store square data
        shapeData{i}.Type = 'Square';
        shapeData{i}.Location = 'left screen';
        shapeData{i}.Center = square_center;
        shapeData{i}.Size = square_size;
        % Delete the square after drawing
        % delete(h_square);
    elseif i == 3
        % Draw a square
        h_square = drawrectangle('Rotatable', true);
        % Prompt user to type "yes" when done with the current square
        userResponse = input(['Type "yes" when done with square at right screen: '], 's');
        square_center = h_square.Position(1:2) + h_square.Position(3:4) / 2;
        square_size = h_square.Position(3:4);
        % Store square data
        shapeData{i}.Type = 'Square';
        shapeData{i}.Location = 'right screen';
        shapeData{i}.Center = square_center;
        shapeData{i}.Size = square_size;
        % Delete the square after drawing
        % delete(h_square);
    end
    
    % Check user response to decide whether to proceed
    if ~strcmpi(userResponse, 'yes')
        disp('User did not type "yes". Exiting.');
        return;
    end
end

% Iterate through each shape in shapeData
for i = 1:numel(shapeData)
    if strcmp(shapeData{i}.Type, 'Circle')
        % Scale the Center coordinates separately for X and Y
        scaledCenterX = shapeData{i}.Center(1) * (scaleFactor * focalLength);
        scaledCenterY = shapeData{i}.Center(2) * (scaleFactor * focalLength);
        % Apply rotation to the scaled Center coordinates (for circle, no rotation needed)
        rotatedCenter = rotationMatrix * [scaledCenterX; scaledCenterY];
        % Update the Center coordinates in shapeData
        shapeData{i}.Center = rotatedCenter';
        % Scale the Radius and correct for focal length
        shapeData{i}.Radius = shapeData{i}.Radius * (scaleFactor * focalLength);
        
    else % For squares
        % Apply rotation to the square centers
        % Scale the center coordinates separately for X and Y
        scaledCenterX = shapeData{i}.Center(1) * (scaleFactor * focalLength);
        scaledCenterY = shapeData{i}.Center(2) * (scaleFactor * focalLength);
        % Apply rotation to the scaled center coordinates
        rotatedCenter = rotationMatrix * [scaledCenterX; scaledCenterY];
        % Update the center coordinates in shapeData
        shapeData{i}.Center = rotatedCenter';
        % Scale the size of squares and correct for focal length
        shapeData{i}.Size = shapeData{i}.Size * (scaleFactor * focalLength);
        
    end
end

% Display the extracted data
disp('Data for shapes:');
for i = 1:numel(shapeData)
    disp(['Data for ' shapeData{i}.Type ':']);
    disp(shapeData{i});
end
