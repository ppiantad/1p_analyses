% File name
videoFile = 'C68604_CFD3_D2Afternoon_11132024_condition_two_video.avi';

% Load the video
video = VideoReader(videoFile);

% Extract a single frame
frameNumber = 700; % Modify if you want a different frame
frame = read(video, frameNumber);

% Display the frame
figure;
imshow(frame);
title('Draw a line to measure pixels, then close the line tool');

% Allow user to draw a line and measure it
h = imdistline(gca);
api = iptgetapi(h);
disp('Draw the line, then double-click to finalize the measurement.');

% Wait for user to finalize the line
uiwait(msgbox('Click OK after drawing the line and noting the pixel length.'));

% Get the measured pixel length
measuredPixels = api.getDistance();

% Maddy measured the opening of the arena, which is 17.5 inches (44.45 cm)
% diameter

% Prompt user to input the known real-world length
realWorldLength = inputdlg('Enter the real-world length of the line (in cm):', ...
    'Real-World Length', 1, {'1'}); % Default is 1 cm
realWorldLength = str2double(realWorldLength{1});

% Calculate the pixels-per-cm conversion factor
pixelsPerCm = measuredPixels / realWorldLength;

% Display the conversion factor
fprintf('Measured pixels: %.2f\n', measuredPixels);
fprintf('Real-world length: %.2f cm\n', realWorldLength);
fprintf('Pixels per cm: %.2f\n', pixelsPerCm);

% Store the conversion factor for future use
conversionFactor = pixelsPerCm; % Pixels per cm

% Example usage to convert a measured pixel value to cm
measuredPixelValue = 250; % Example pixel measurement
realLengthInCm = measuredPixelValue / conversionFactor;
fprintf('A measurement of %.2f pixels equals %.2f cm.\n', ...
    measuredPixelValue, realLengthInCm);
