%%
% Specify the path to your video file
videoPath = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-27\RDT D1\RDT D1\BLA-Insc-27_RDT_D1_2023-01-02T12_28_19.avi';

% Create a VideoReader object
videoObj = VideoReader(videoPath);

% Read a random frame from the video
randomFrameIndex = randi([1, videoObj.NumFrames]);
randomFrame = read(videoObj, randomFrameIndex);

% Create a figure and display the random frame
figure;
imshow(randomFrame);
title('Mouse over the frame to display coordinates');

% Set up a callback function for the mouse-over event
set(gcf, 'WindowButtonMotionFcn', @mouseOverCallback);