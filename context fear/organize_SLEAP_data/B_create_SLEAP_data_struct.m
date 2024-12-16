% Define the top-level directory containing all primary subfolders
top_level_directory = 'D:\Context Data\PFC Last\Raw Data\PFC alone\Raw Data';  % Replace with your actual directory path

% Get a list of primary subfolders in the top-level directory
primary_subfolders = dir(top_level_directory);
primary_subfolders = primary_subfolders([primary_subfolders.isdir] & ~ismember({primary_subfolders.name}, {'.', '..'}));

% Loop through each primary subfolder
for p = 1:length(primary_subfolders)
    % Define the path to the current primary subfolder
    primary_subfolder_path = fullfile(top_level_directory, primary_subfolders(p).name);
    
    % Get a list of sub-subfolders within the current primary subfolder
    secondary_subfolders = dir(primary_subfolder_path);
    secondary_subfolders = secondary_subfolders([secondary_subfolders.isdir] & ~ismember({secondary_subfolders.name}, {'.', '..'}));
    folder_list_string = strsplit(primary_subfolder_path, '\');

    current_animal = folder_list_string{4}; % Would have to change this depending on your folder structure, but there should be an animal name folder given our current workflow.



    % Loop through each subfolder in the current primary subfolder
    for s = 1:length(secondary_subfolders)
        % Define the path to the current sub-subfolder
        subfolder_path = fullfile(primary_subfolder_path, secondary_subfolders(s).name);
        folder_list_string = strsplit(subfolder_path, '\');

        current_animal = folder_list_string{4}; % Would have to change this depending on your folder structure, but there should be an animal name folder given our current workflow.


        current_session_raw = folder_list_string{5};

        % Locate the .csv file in the current subfolder
        csv_file = dir(fullfile(subfolder_path, 'SLEAP_and_freezing_combined*.csv'));

        % Check if a CSV file exists in the subfolder
        if ~isempty(csv_file)
            % Read the CSV file
            data = readtable(fullfile(subfolder_path, csv_file.name));
            final_DLC.(current_animal).(current_session_raw).movement_data = data;
            final_DLC.(current_animal).(current_session_raw).filename = csv_file.name;

        end
    end
end


%%
experimental_grps = readtable('E:\MATLAB\my_repo\context fear\organize_SLEAP_data\full_pilot_mice.xlsx');
animalIDs = fieldnames(final_DLC);

for dd = 1:size(experimental_grps, 1)
    current_mouse = experimental_grps.mouse{dd};
    current_mouse_condition = experimental_grps.group{dd};
    for hh = 1:size(animalIDs, 1)
        if strcmp(current_mouse, animalIDs(hh))
            final_DLC.(current_mouse).experimental_grp = current_mouse_condition;
            if any("sex" == string(experimental_grps.Properties.VariableNames))
                final_DLC.(current_mouse).sex = experimental_grps.sex{hh};
            end
        end

    end

end
