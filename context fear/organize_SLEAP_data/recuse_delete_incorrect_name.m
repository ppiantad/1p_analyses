% Define the top-level directory containing all primary subfolders
top_level_directory = 'F:\Maddy Pilot\Abbreviated_Pilot';  % Replace with your actual directory path

% Search for all files named "SLEAP_raw_data.csv" in the directory and its subdirectories
files_to_delete = dir(fullfile(top_level_directory, '**', '*SLEAP_raw_data*.csv'));

% Loop through each file and delete it
for k = 1:length(files_to_delete)
    file_path = fullfile(files_to_delete(k).folder, files_to_delete(k).name);
    delete(file_path);
    disp(['Deleted file: ', file_path]);
end
