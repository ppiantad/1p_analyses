%% be very careful when using! set mainDir correctly

% Specify the main directory to start searching from
mainDir = 'D:\Context Data\Pilot\data'; % Opens a di Opens a dialog to select the folder, or replace with a string path

% Get all files in the directory tree
allFiles = dir(fullfile(mainDir, '**', '*.slp'));

% Loop through each file
for i = 1:length(allFiles)
    % Get the full path of the current file
    filePath = fullfile(allFiles(i).folder, allFiles(i).name);
    
    % Delete the file
    if exist(filePath, 'file') % Check if the file exists
        delete(filePath);
        fprintf('Deleted file: %s\n', filePath);
    end
end
