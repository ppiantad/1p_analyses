%% Use this code to create a meta-structure 'final_SLEAP' that stores movement-related data

%% Edit these uservariables with what you want to look at
uv.evtWin = [-10 10]; %what time do you want to look at around each event
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at
ts1 = (-10:(uv.dt):10-0.1);
% Define the directory path you want to start with
% startDirectory = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-24';

metaDirectory = 'D:\Context Data\PFC Last\Raw Data\PFC alone\Raw Data';
metaDirectory_subfolders = dir(metaDirectory );
metafolder_list = {};


% Loop through the list of subfolders
for i = 1:length(metaDirectory_subfolders)
    % Check if the item in subfolders is a directory (not "." or "..") or
    % one of the sets of files that I haven't analyzed yet (PR currently)
    if metaDirectory_subfolders(i).isdir && ~strcmp(metaDirectory_subfolders(i).name, '.') && ~strcmp(metaDirectory_subfolders(i).name, '..') && ~contains(metaDirectory_subfolders(i).name, 'PR') && ~contains(metaDirectory_subfolders(i).name, 'not in final dataset')
        % Get the full path of the subfolder
        metasubfolderPath = fullfile(metaDirectory, metaDirectory_subfolders(i).name);
        % Create a cell array for the subfolder path and append it
        % vertically to folder_list
        metafolder_list = vertcat(metafolder_list, {metasubfolderPath});
    end
end



for zz = 1:size(metafolder_list, 1)
    % Use the dir function to get a list of subfolders
    startDirectory = metafolder_list{zz};
    subfolders = dir(startDirectory);
    % Initialize folder_list as an empty cell array
    folder_list = {};
    % Loop through the list of subfolders
    for i = 1:length(subfolders)
        % Check if the item in subfolders is a directory (not "." or "..") or
        % one of the sets of files that I haven't analyzed yet (PR currently)
        if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..') && ~contains(subfolders(i).name, 'PR')
            % Get the full path of the subfolder
            subfolderPath = fullfile(startDirectory, subfolders(i).name);
            % Create a cell array for the subfolder path and append it
            % vertically to folder_list
            folder_list = vertcat(folder_list, {subfolderPath});
        end
    end



    for ii = 1:size(folder_list, 1)
        folder_list_string = strsplit(folder_list{ii}, '\');
        current_animal = folder_list_string{5}; % Would have to change this depending on your folder structure, but there should be an animal name folder given our current workflow.
        current_animal = matlab.lang.makeValidName(current_animal);
        current_session = folder_list_string{6};
        if contains(lower(current_session), 'aft') || contains(lower(current_session), 'noon')
            if contains(lower(current_session), '1')
                current_session_renamed = 'D1_Afternoon';

            elseif contains(lower(current_session), '2')
                current_session_renamed = 'D2_Afternoon';
            end
        elseif contains(lower(current_session), 'morn')
            if contains(lower(current_session), '1')
                current_session_renamed = 'D1_Morning';

            elseif contains(lower(current_session), '2')
                current_session_renamed = 'D2_Morning';
            end
        elseif contains(lower(current_session), '3')
            current_session_renamed = 'D3';

        elseif contains(lower(current_session), '4')
            current_session_renamed = 'D4';
        end


        list = dir(folder_list{ii});%grab a directory of the foldercontents
        list_folder_names = {list.name}.';
        sub_folder_list = {};

        % Loop through the list of subfolders
        for i = 1:length(list)
            % Check if the item in subfolders is a directory (not "." or "..") or
            % one of the sets of files that I haven't analyzed yet (PR currently)
            if list(i).isdir && ~strcmp(list(i).name, '.') && ~strcmp(list(i).name, '..') && contains(list(i).name, 'cylander_behaviour')
                % Get the full path of the subfolder
                subsubfolderPath = fullfile(startDirectory, current_session, list(i).name);
                % Create a cell array for the subfolder path and append it
                % vertically to folder_list
                % sub_folder_list = vertcat(sub_folder_list,{subsubfolderPath});
                sub_folder_list_dir = dir(subsubfolderPath);%grab a directory of the foldercontents

                for qq = 1:size(sub_folder_list_dir, 1)
                    if sub_folder_list_dir(qq).isdir && ~strcmp(sub_folder_list_dir(qq).name, '.') && ~strcmp(sub_folder_list_dir(qq).name, '..') && contains(sub_folder_list_dir(qq).name, 'exports')
                        % Get the full path of the subfolder
                        subsubsubfolderPath = fullfile(subsubfolderPath, sub_folder_list_dir(qq).name);
                        % Create a cell array for the subfolder path and append it
                        % vertically to folder_list
                        % sub_sub_folder_list = vertcat(sub_folder_list,{subsubfolderPath});
                        sub_sub_folder_list_dir = dir(subsubsubfolderPath);%grab a directory of the foldercontents

                    end
                end
            end
        end



        folderMask = ~[sub_sub_folder_list_dir.isdir]; %find all of the folders in the directory and remove them from the list
        files = sub_sub_folder_list_dir(folderMask);  %now we have only files to work with
        clear folderMask


        idx = ~cellfun('isempty',strfind({files.name},'.csv')); %find the instances of .xlsx in the file list.
        %This command converts the name field into a cell array and searches
        %the cell array with strfind
        csvFiles = files(idx); %build a mat file index
        clear idx
        csv_names = {csvFiles.name};

        for mm = 1:length(csv_names)
            % Check if the current name contains three distinct substrings
            if contains(lower(csv_names{mm}), 'dlc')
                disp(['DLC = ', csv_names{mm}])
                DLC_file = strcat(subsubsubfolderPath, '\', csv_names{mm});
            end
            if contains(csv_names{mm}, 'ABET')
                disp(['ABET File = ', csv_names{mm}])
                ABET_file = strcat(folder_list{ii}, '\', csv_names{mm});
            end

            % if ~contains(csv_names{mm}, 'ABET') | ~contains(lower(csv_names{mm}), '_gpio')
            %     disp('BLANK folder was not analyzed due to missing files! Check contents and try again');
            %
            %     % Skip the rest of the loop for this folder
            %     continue;
            % end
        end

        if exist('DLC_file','var')
            disp('Required files found, analyzing...');


            DLC_data = readtable(DLC_file);
            DLC_data.mean_x_pix = mean([DLC_data.left_ear_x, DLC_data.right_ear_x, DLC_data.tail_x], 2);
            DLC_data.mean_y_pix = mean([DLC_data.left_ear_y, DLC_data.right_ear_y, DLC_data.tail_y], 2);
            % save non-downsampled table as well - for various
            % analyses where we don't need to get on the same
            % timescale
            final_DLC.(current_animal).(current_session_renamed).DLC_data_raw = DLC_data;
            final_DLC.(current_animal).(current_session_renamed).filename = DLC_file;


            clear ABET_file SLEAP_csv DLC_file
            % clearvars -except final

        end
    end
end
%%
animalIDs = fieldnames(final_DLC);
experimental_grp_label = 'Experimental';
experimental_grps = repmat(experimental_grp_label, size(animalIDs, 1), 1);


for dd = 1:size(experimental_grps, 1)
    current_mouse = animalIDs{dd};
    label = experimental_grps(dd, :); 
    final_DLC.(current_mouse).experimental_grp = label;


end


