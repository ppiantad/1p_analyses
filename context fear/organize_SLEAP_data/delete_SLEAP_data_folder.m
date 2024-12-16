%% USE WITH CAUTION! MAKE SURE YOU SET mainDir correctly! 

% Specify the main directory to start searching from
mainDir = 'D:\Context Data\Pilot\data'; % Opens a dialog to select the folder, or replace with a string path

% Get all subfolders in the directory tree
subDirs = dir(fullfile(mainDir, '**', '*'));

% Filter for directories only
subDirs = subDirs([subDirs.isdir]);

% Loop through each subdirectory
for i = 1:length(subDirs)
    % Get the name of the current subdirectory
    folderName = subDirs(i).name;

    % Skip '.' and '..' directories
    if strcmp(folderName, '.') || strcmp(folderName, '..')
        continue;
    end

    % Check if the folder is named "SLEAP_data"
    if strcmp(folderName, 'SLEAP_data')
        % Get the full path of the folder
        folderPath = fullfile(subDirs(i).folder, folderName);
        
        % Delete the folder and its contents
        disp(['Deleting folder: ', folderPath]);
        rmdir(folderPath, 's');
    end
end
