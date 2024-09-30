% Load the video file
mpg_video = VideoReader('BLA-Insc-24_RDT_D12022-08-10T15_37_21.avi');

% Set the output video file properties
outputVideo = VideoWriter('downsampled_video_insc_24_RDT_D1.avi');
outputVideo.FrameRate = 10; % Set the frame rate to 10 Hz
open(outputVideo);

% Initialize a counter for frame extraction
frameCount = 0;

% Loop through the 30 Hz video and extract every 3rd frame
while hasFrame(mpg_video)
    frame = readFrame(mpg_video);
    frameCount = frameCount + 1;
    
    % Extract every third frame (downsample)
    if mod(frameCount, 3) == 1
        writeVideo(outputVideo, frame);
    end
end

% Close the output video file
close(outputVideo);

disp('Downsampling complete. The new video is saved as downsampled_video.mp4');