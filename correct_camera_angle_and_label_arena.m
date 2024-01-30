% Specify the path to your video file
videoPath = 'd:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-27\RDT D1\RDT D1\BLA-Insc-27_RDT_D1_2023-01-02T12_28_19.avi';

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


%%



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

% Display the corrected X and Y coordinates
disp('Corrected X Coordinates:');
disp(correctedX');

disp('Corrected Y Coordinates:');
disp(correctedY');


%%
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