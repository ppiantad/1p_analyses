% Define the folder to search recursively
folder = 'D:\risk videos\BLA stGtACR vs EGFP';
% Initialize variables to store data for the table
fileNames = {};
changeTimes = [];

% Recursively traverse the folder structure
fileList = dir(fullfile(folder, '**', '*.MP4'));
for i = 1:length(fileList)
    % Check if the file name matches the pattern
    if endsWith(fileList(i).name, 'merged_resized_grayscaled.MP4', 'IgnoreCase', true)
        % Load the video
        videoPath = fullfile(fileList(i).folder, fileList(i).name);
        videoObj = VideoReader(videoPath);
        
        % Define the duration of the initial portion (in seconds)
        initialDuration = 50; % Adjust as needed
        
        % Read frames from the initial portion
        numFrames = round(videoObj.FrameRate * initialDuration);
        frames = cell(1, numFrames);
        for j = 1:numFrames
            frames{j} = readFrame(videoObj);
        end
        
        % Compute brightness for each frame
        brightness = zeros(1, numFrames);
        for j = 1:numFrames
            frame = frames{j};
            brightness(j) = mean(frame(:));
        end
        

        % Find the minimum brightness value and its index
        [minBrightness, minIndex] = min(brightness);
        % Look for significant brightness jumps after the minimum frame
        % Start from the frame after the minimum brightness frame
        startFrame = minIndex + 1;
        
        % Compute brightness changes from startFrame to the end
        brightnessChanges = diff(brightness(startFrame:end));
        
        % Find significant changes using a threshold
        threshold = 1; % Adjust as needed
        significantChangesIndices = find(abs(brightnessChanges) > threshold);
        
        % If there are significant changes, get the index of the first one
        if ~isempty(significantChangesIndices)
            % Adjust the index to account for starting from startFrame
            firstSignificantChangeIndex = significantChangesIndices(1) + startFrame - 1;

            % Display the results
            fprintf('Video: %s\n', fileList(i).name);
            fprintf('Minimum brightness value: %.2f\n', minBrightness);
            fprintf('Frame index with minimum brightness: %d\n', minIndex);
            fprintf('Frame index with first significant brightness change: %d\n', firstSignificantChangeIndex);
            % Display the frames in a subplot with 5 panels

            % Calculate the time of the first significant change
            changeTime = firstSignificantChangeIndex / videoObj.FrameRate;
            
            % Store the filename and change time
            fileNames{end+1} = fileList(i).name;
            changeTimes(end+1) = changeTime;



            % figure;
            % subplot(1, 5, 1);
            % imshow(frames{max(1, firstSignificantChangeIndex - 2)});
            % title('Frame -2');
            % 
            % subplot(1, 5, 2);
            % imshow(frames{max(1, firstSignificantChangeIndex - 1)});
            % title('Frame -1');
            % 
            % subplot(1, 5, 3);
            % imshow(frames{firstSignificantChangeIndex});
            % title('Frame 0 (Significant Change)');
            % 
            % subplot(1, 5, 4);
            % imshow(frames{min(numFrames, firstSignificantChangeIndex + 1)});
            % title('Frame +1');
            % 
            % subplot(1, 5, 5);
            % imshow(frames{min(numFrames, firstSignificantChangeIndex + 2)});
            % title('Frame +2');
            % 
            % sgtitle('Frames around the first significant brightness change');
            % You can also display the frame with minimum brightness and the frame with the first significant change using imshow if needed
            % imshow(frames{minIndex});
            % imshow(frames{firstSignificantChangeIndex});
        else
            fprintf('No significant brightness changes found after the frame with minimum brightness.\n');
        end
    end
end


% Create a table from the collected data
dataTable = table(fileNames', changeTimes', 'VariableNames', {'FileName', 'ChangeTimeInSeconds'});

% Save the table to a CSV file
% writetable(dataTable, 'brightness_change_data.csv');

%%


% Loop through each file
for i = 1:size(dataTable, 1)
    % Get the filename
    filename = dataTable.FileName{i};
    adjusted_start_time = dataTable.ChangeTimeInSeconds(i);
    % Extract the relevant information from the filename
    parts = strsplit(filename, '_merged_resized_grayscaled.MP4');
    pattern = '_\d+$';  % Match underscore followed by digits at the end of the string
    replacement = '';
    filtered_str = regexprep(parts(1), pattern, replacement);


    % Initialize a flag to indicate if the combination exists
    exists = false;

    % Loop through the first level of the structure
    fields_level1 = fieldnames(final_SLEAP);
    for qq = 1:numel(fields_level1)
        % Combine the first level field name with each second level field name
        fields_level2 = fieldnames(final_SLEAP.(fields_level1{qq}));
        for j = 1:numel(fields_level2)
            combination = [fields_level1{qq}, '_', fields_level2{j}];
            % Check if the combination matches the filtered_str
            if strcmp(combination, filtered_str)
                exists = true;
                
                BehavData = final_SLEAP.(fields_level1{qq}).(fields_level2{j}).BehavData;
                BehavData.TrialPossible(:)=BehavData.TrialPossible(:)+adjusted_start_time(1);
                BehavData.choiceTime(:)=BehavData.choiceTime(:)+adjusted_start_time(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+7.39500000000000;
                BehavData.collectionTime(:)=BehavData.collectionTime(:)+adjusted_start_time(1);
                BehavData.stTime(:)=BehavData.stTime(:)+adjusted_start_time(1);
                final_SLEAP.(fields_level1{qq}).(fields_level2{j}).BehavData = BehavData; 
                final_SLEAP.(fields_level1{qq}).(fields_level2{j}).first_frame = adjusted_start_time; 
                break; % Break the inner loop if match found
            end
        end
        if exists
            break; % Break the outer loop if match found
        end
    end

    % Output the result
    if exists
        disp('The combination exists.');
    else
        disp('The combination does not exist.');
    end
end

% Save the updated 'final' data structure
% save('final_updated.mat', 'final');