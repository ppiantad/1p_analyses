classdef DownsampledVideoReader
    properties
        OriginalVideo % Store the original VideoReader object
        DownsampleFactor % Factor by which to downsample frames
        NumFrames % Number of downsampled frames
        FrameRate % Frame rate (same as original)
        CurrentFrameIndex % Current frame index in the downsampled video
    end
    
    methods
        % Constructor method
        function obj = DownsampledVideoReader(videoFile, downsampleFactor)
            obj.OriginalVideo = VideoReader(videoFile);
            obj.DownsampleFactor = downsampleFactor;
            obj.NumFrames = floor(obj.OriginalVideo.NumFrames / downsampleFactor);
            obj.FrameRate = obj.OriginalVideo.FrameRate;
            obj.CurrentFrameIndex = 1;
        end
        
        % Method to read the next downsampled frame
        function frame = readFrame(obj)
            % Calculate the actual frame to read based on the downsample factor
            actualFrameIndex = (obj.CurrentFrameIndex - 1) * obj.DownsampleFactor + 1;
            
            % Seek to the actual frame in the original video
            obj.OriginalVideo.CurrentTime = (actualFrameIndex - 1) / obj.OriginalVideo.FrameRate;
            
            % Read the frame
            frame = readFrame(obj.OriginalVideo);
            
            % Increment the downsampled frame index
            obj.CurrentFrameIndex = obj.CurrentFrameIndex + 1;
        end
        
        % Check if there's another frame to read
        function tf = hasFrame(obj)
            actualFrameIndex = (obj.CurrentFrameIndex - 1) * obj.DownsampleFactor + 1;
            tf = actualFrameIndex <= obj.OriginalVideo.NumFrames;
        end
    end
end
