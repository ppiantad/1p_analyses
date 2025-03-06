%% Use this code to create a meta-structure 'final_SLEAP' that stores movement-related data

%% Edit these uservariables with what you want to look at
uv.evtWin = [-10 10]; %what time do you want to look at around each event
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at
ts1 = (-10:(uv.dt):10-0.1);
% Define the directory path you want to start with
% startDirectory = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-24';

metaDirectory = 'F:\Behavior Videos\BLA hM4Di vs mCherry';
metaDirectory_subfolders = dir(metaDirectory );
metafolder_list = {};

if exist('final_SLEAP', 'var') ~= 1
    final_SLEAP = struct;
end

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
        current_animal = folder_list_string{4}; % Would have to change this depending on your folder structure, but there should be an animal name folder given our current workflow.
        current_animal = matlab.lang.makeValidName(current_animal);
        current_session = folder_list_string{5};
        current_session = regexprep(current_session,{' ', '-'}, '_');
        modifiedString = lower(strrep(strrep(folder_list_string{end}, ' ', ''), '-', ''));
        
        
        list = dir(folder_list{ii});%grab a directory of the foldercontents
        list_folder_names = {list.name}.';

        folderMask = ~[list.isdir]; %find all of the folders in the directory and remove them from the list
        files = list(folderMask);  %now we have only files to work with
        clear folderMask

        if ~isfield(final_SLEAP, current_animal)
            final_SLEAP.(current_animal) = struct;

            if ~isfield(final_SLEAP.(current_animal), current_session)
                final_SLEAP.(current_animal).(current_session) = struct;
            end
        end

        if ~isfield(final_SLEAP.(current_animal), current_session)
            final_SLEAP.(current_animal).(current_session) = struct;

        end



        % check folder inside to get SLEAP file
        folder_to_analyze = find(strcmpi(strrep(strrep(list_folder_names, ' ', ''), '-', ''), modifiedString));
        disp(['Analyzing subfolder: ' list_folder_names{folder_to_analyze,1}]);
        folder_to_analyze_Path = fullfile(folder_list{ii}, list(folder_to_analyze).name);

        sleap_folders = dir(folder_to_analyze_Path);

        sleap_folder_list = {};
            % Loop through the list of subfolders
            for i = 1:length(sleap_folders)
                % Check if the item in subfolders is a directory (not "." or "..") or
                % one of the sets of files that I haven't analyzed yet (PR currently)
                if sleap_folders(i).isdir && ~strcmp(sleap_folders(i).name, '.') && ~strcmp(sleap_folders(i).name, '..')

                    % Get the full path of the subfolder
                    sleap_subfolderPath = fullfile(folder_to_analyze_Path, sleap_folders(i).name);

                    % Create a cell array for the subfolder path and append it
                    % vertically to folder_list
                    sleap_folder_list = vertcat(sleap_folder_list, {sleap_subfolderPath});

                end
            end
        

            for sleap_folder = 1:size(sleap_folder_list, 1)
                list = dir(sleap_folder_list{sleap_folder});%grab a directory of the foldercontents
                disp(['Analyzing subfolder: ' sleap_folder_list{sleap_folder}]);

                folderMask = ~[list.isdir]; %find all of the folders in the directory and remove them from the list
                files = list(folderMask);  %now we have only files to work with
                clear folderMask list

                idx = ~cellfun('isempty',strfind({files.name},'body_sleap_data.csv')); %find the instances of .xlsx in the file list.
                %This command converts the name field into a cell array and searches
                %the cell array with strfind
                csvFiles = files(idx); %build a mat file index
                clear idx files


                if isempty(csvFiles)
                    disp('Missing body_sleap_data.csv file, skipping folder');
                    continue
                else
                    SLEAP_csv = strcat(sleap_folder_list{sleap_folder}, '\', csvFiles.name);

                end
                % Check if all required files exist, if not, skip
                if exist('SLEAP_csv','var')
                    disp('Required files found, analyzing...');                           


                        SLEAP_data = readtable(SLEAP_csv);
                        
                        % save non-downsampled table as well - for various
                        % analyses where we don't need to get on the same
                        % timescale
                        final_SLEAP.(current_animal).(current_session).SLEAP_data_raw = SLEAP_data;

                        % Assuming 'Timestamp' is the timestamp variable in your table
                        timestamps = SLEAP_data.idx_time;

                        % Specify the original and target sampling rates
                        originalSamplingRate = 30; % Hz
                        targetSamplingRate = 10;    % Hz

                        % Calculate the downsampled indices
                        downsampledIndices = round(linspace(1, height(SLEAP_data), height(SLEAP_data) / (originalSamplingRate / targetSamplingRate)));

                        % Downsample the table
                        SLEAP_downsampled_data = SLEAP_data(downsampledIndices, :);

                        SLEAP_data = SLEAP_downsampled_data;



                        % SLEAP_data.vel_filtered = sgolayfilt(SLEAP_data.vel_cm_s, 2, 33);
                        SLEAP_data.vel_filtered_2 = sgolayfilt(SLEAP_data.vel_cm_s, 3, 25);
                        % SLEAP_data.x_pix_filtered = sgolayfilt(SLEAP_data.x_pix, 2, 33);
                        % SLEAP_data.y_pix_filtered = sgolayfilt(SLEAP_data.y_pix, 2, 33);
                        % SLEAP_data.pix_calc_2 = SLEAP_data.x_pix*(2.54/96);

                        % SLEAP_data.pix_calc_3= SLEAP_data.pix_calc_2 * (2.54/96) * (30/1);

                        % SLEAP_time = uv.dt:uv.dt:height(SLEAP_data)*uv.dt; %generate time trace

                        %adjust  time to account for the fact that Inscopix recording
                        %starts first (stTime(1);
                        % SLEAP_time = SLEAP_time + stTime(1);

                        % SLEAP_data.idx_time = SLEAP_time';

                        % if ~isempty(SLEAP_time_range_adjustment)
                        %     time_ranges = time_ranges-SLEAP_time_range_adjustment;
                        % end

                        SLEAP_data_vel_filtered_session = SLEAP_data.vel_filtered_2';

                        zscored_SLEAP_data_vel_filtered_session =(SLEAP_data_vel_filtered_session-mean(SLEAP_data_vel_filtered_session)./std(SLEAP_data_vel_filtered_session));

                        % final_SLEAP.(current_animal).(current_session).time = SLEAP_time; %final(i).time = caTime;
                        final_SLEAP.(current_animal).(current_session).uv = uv;
                        final_SLEAP.(current_animal).(current_session).SLEAP_data = SLEAP_data;
                        %Because the "neuron" data type is a pain in the ass to
                        %work with, we will instead save some of the variables
                        final_SLEAP.(current_animal).(current_session).SLEAP_data_velocity = SLEAP_data_vel_filtered_session;
                        final_SLEAP.(current_animal).(current_session).zscored_SLEAP_data_velocity = zscored_SLEAP_data_vel_filtered_session;


                        clear unitTS unitTrace unitXTrials unitAVG unitSEM i zall zb zsd zb_window zsd_window zall_window zb_session zsd_session zall_session neuron

                   
                else
                    disp('Remember: folder was not analyzed due to missing files! Check contents and try again');

                    % Skip the rest of the loop for this folder
                    continue;
                end
                
            end
            clear SLEAP_csv
            % clearvars -except final
    end
end

