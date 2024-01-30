% Specify the path to your video file
videoPath = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-29\RDT D1\RDT D1\BLA-Insc-29_RDT_D1_2023-01-18T14_17_37.avi';

% Create a VideoReader object
videoObj = VideoReader(videoPath);

% Read a random frame from the video
randomFrameIndex = randi([1, videoObj.NumFrames]);
randomFrame = read(videoObj, randomFrameIndex);


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