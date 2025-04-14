% Inputs
current_animal = 'RRD453';
current_session = 'RDT_OPTO_CHOICE';

% Load shapeData
shapeData = final_SLEAP.(current_animal).(current_session).shapeData;

% Load video (adjust path or use info from your structure if stored there)
video_path = "D:\risk videos\BLA PdCO vs mCherry\RRD453\RDT OPTO CHOICE\RRD453_RDT_OPTO_CHOICE_10022024_merged_resized_grayscaled.MP4";  % or manually: 'path/to/video.avi'
v = VideoReader(video_path);

% Select frame
frame_num = 1800;
v.CurrentTime = (frame_num - 1) / v.FrameRate;
video_frame = readFrame(v);

% Extract shape info for the selected frame
shapes = shapeData;  % 3x1 cell array

% Plot the video frame
figure;
imshow(video_frame); hold on;

% Loop over each shape
for i = 1:numel(shapes)
    this_shape = shapes{i};
    shape_type = this_shape.Type;

    switch lower(shape_type)
        case 'circle'
            center = this_shape.Center;
            radius = this_shape.Radius;
            viscircles(center, radius, 'Color', 'r', 'LineWidth', 1);

        case 'square'
            bbox = this_shape.("BoundingBox");
            rectangle('Position', bbox, 'EdgeColor', 'g', 'LineWidth', 1);

        otherwise
            warning('Unknown shape type: %s', shape_type);
    end
end

title(sprintf('%s - %s - Frame %d', current_animal, current_session, frame_num));