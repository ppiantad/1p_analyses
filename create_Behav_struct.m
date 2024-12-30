%% Edit these uservariables with what you want to look at
% uv.evtWin = [-10 10]; %what time do you want to look at around each event
% uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at
ts1 = (-10:.1:10-0.1);
% Define the directory path you want to start with
% startDirectory = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-24';

metaDirectory = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA';
metaDirectory_subfolders = dir(metaDirectory );
metafolder_list = {};


% Loop through the list of subfolders
for i = 1:length(metaDirectory_subfolders)
    % Check if the item in subfolders is a directory (not "." or "..") or
    % one of the sets of files that I haven't analyzed yet (PR currently)
    if metaDirectory_subfolders(i).isdir && ~strcmp(metaDirectory_subfolders(i).name, '.') && ~strcmp(metaDirectory_subfolders(i).name, '..') &&  ~contains(metaDirectory_subfolders(i).name, 'not in final dataset')
        % if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..') && ~contains(lower(subfolders(i).name), 'shock')
        % if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..')
        % Get the full path of the subfolder
        metasubfolderPath = fullfile(metaDirectory, metaDirectory_subfolders(i).name);

        % Create a cell array for the subfolder path and append it
        % vertically to folder_list
        metafolder_list = vertcat(metafolder_list, {metasubfolderPath});



        % Add your analysis code here

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
        if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..')
            % if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..') && ~contains(lower(subfolders(i).name), 'shock')
            % if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..')
            % Get the full path of the subfolder
            subfolderPath = fullfile(startDirectory, subfolders(i).name);

            % Create a cell array for the subfolder path and append it
            % vertically to folder_list
            folder_list = vertcat(folder_list, {subfolderPath});



            % Add your analysis code here

        end
    end

    for ii = 1:size(folder_list, 1)
        list = dir(folder_list{ii});%grab a directory of the foldercontents
        disp(['Analyzing subfolder: ' folder_list{ii,1}]);

        % Initialize a flag to check if files were found in this folder
        filesFound = false;


        folderMask = ~[list.isdir]; %find all of the folders in the directory and remove them from the list
        files = list(folderMask);  %now we have only files to work with
        clear folderMask list
        %now we can sort the files into subtype
        idx = ~cellfun('isempty',strfind({files.name},'_final.mat')); %find the instances of .mat in the file list.
        %This command converts the name field into a cell array and searches
        %the cell array with strfind
        matFiles = files(idx); %build a mat file index

        currentMatFile = strcat(folder_list{ii}, '\', matFiles.name);
        clear idx

        idx = ~cellfun('isempty',strfind({files.name},'.csv')); %find the instances of .xlsx in the file list.
        %This command converts the name field into a cell array and searches
        %the cell array with strfind
        csvFiles = files(idx); %build a mat file index
        clear idx

        idx = ~cellfun('isempty',strfind({files.name},'.xlsx')); %find the instances of .xlsx in the file list.
        %This command converts the name field into a cell array and searches
        %the cell array with strfind
        xlsxFiles = files(idx); %build a mat file index
        clear idx files

        folder_strings = strsplit(folder_list{ii}, '\');
        %     session_strings = strsplit(folder_strings{end}, {'-', ' '});
        %     mat_strings = strsplit(char(matFiles.name),'_');
        %     date_strings = strsplit(mat_strings{4}, '-');
        csv_names = {csvFiles.name};
        current_animal = folder_strings{5}; % Would have to change this depending on your folder structure, but there should be an animal name folder given our current workflow.
        % current_session = folder_strings{6};
        current_animal = matlab.lang.makeValidName(current_animal);
        current_session = char(folder_strings(end));

        % Loop over each substring in the substrings array
        for mm = 1:length(csv_names)
            % Check if the current name contains three distinct substrings
            if contains(lower(csv_names{mm}), '_gpio')
                disp(['GPIO File = ', csv_names{mm}])
                GPIO_file = strcat(folder_list{ii}, '\', csv_names{mm});
            end
            if contains(csv_names{mm}, 'ABET')
                disp(['ABET File = ', csv_names{mm}])
                ABET_file = strcat(folder_list{ii}, '\', csv_names{mm});
            end
            if contains(csv_names{mm}, 'BORIS')
                disp(['boris File = ', csv_names{mm}])
                boris_file = strcat(folder_list{ii}, '\', csv_names{mm});
            end
            
        end

        if isempty(matFiles) || ~exist('GPIO_file', 'var') == 1 || ~exist('ABET_file', 'var')  == 1
            disp('Missing .mat or .csv files, skipping folder');
            clear ABET_file GPIO_file matFiles
        else
            currentMatFile = strcat(folder_list{ii}, '\', matFiles.name);
            filesFound = true; % Set the flag to true since .mat files were found
        end


        % Check the filesFound flag and print the final message
        if filesFound
            disp('Folder analyzed successfully');

            if strcmp(current_session, 'SHOCK TEST')
                current_session = regexprep(current_session,{' ', '-'}, '_');
                [BehavData,ABETfile]=ABET2TableFn_ShockTest(ABET_file);
                ABET_removeheader = ABETfile(2:end,:);
                tbl_ABET = cell2table(ABET_removeheader);
                tbl_ABET.Properties.VariableNames = ABETfile(1,:);
                gpio_tbl = readtable(GPIO_file);

                % Inscopix stupidly recently updated the way the GPIO pins
                % are written, so the default values that I used before are
                % no longer valid. updated this to reflect the new way
                % things are written - hopefully this works! 
                % stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Time_s_ > 0);
                stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Value > 5000);



                BehavData.TrialPossible(:)=BehavData.TrialPossible(:)+stTime(1);
                BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+7.39500000000000;
                % BehavData.collectionTime(:)=BehavData.collectionTime(:)+stTime(1);
                % BehavData.stTime(:)=BehavData.stTime(:)+stTime(1);
                % shk_times(:)=shk_times(:)+stTime(1);

                BehavData.choTime2 = BehavData.choiceTime-BehavData.TrialPossible;
                % BehavData.choTime3 = BehavData.Insc_TTL+BehavData.choTime2;

                %filter based on TrialFilter inputs (see TrialFilter.m for full list of
                %possibilities)
                BehavData=TrialFilter(BehavData,'SHK',1);


                final_behavior.(current_animal).(current_session).uv.BehavData = BehavData;
            elseif contains(current_session, 'PR')
                current_session = regexprep(current_session,{' ', '-'}, '_');
                [BehavData,ABETfile]=ABET2TableFn_PR(ABET_file);
                ABET_removeheader = ABETfile(2:end,:);
                tbl_ABET = cell2table(ABET_removeheader);
                tbl_ABET.Properties.VariableNames = ABETfile(1,:);
                gpio_tbl = readtable(GPIO_file);
                % Inscopix stupidly recently updated the way the GPIO pins
                % are written, so the default values that I used before are
                % no longer valid. updated this to reflect the new way
                % things are written - hopefully this works!
                % stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Time_s_ > 0);
                stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Value > 5000);
                frames = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName,'BNC Sync Output') & gpio_tbl.Value == 1);
                %check GPIO file to extract each TTL, since the TTL is 1000ms and is
                %sampled repeatedly. This will only extract events that are separated by >
                %8sec, so be sure to change this if the TTL or task structure changes
                %dramatically!
                pp = 2;
                ttl_filtered = stTime(1);
                for kk = 1:size(stTime,1)-1
                    if abs(stTime(kk)-stTime(kk+1)) > 8
                        ttl_filtered(pp) = stTime(kk+1);
                        pp=pp+1;
                    end
                end
                ttl_filtered = ttl_filtered';
                %Add TTL times received by Inscopix to data table, skipping omitted trials
                %which do not have a corresponding TTL due to a quirk in the behavioral
                %program
                % BehavData.Insc_TTL = zeros(length(BehavData.TrialPossible),1);
                % dd = 2;
                % for cc = 1:size(BehavData, 1)
                %     if BehavData.TrialPossible(cc) > stTime(1)
                %         BehavData.Insc_TTL(cc) = ttl_filtered(dd);
                %         dd = dd+1;
                %     elseif BehavData.TrialPossible(cc) <= stTime(1)
                %         BehavData.Insc_TTL(cc) = 0;
                %     end
                % end

                BehavData.PressTime(:)=BehavData.PressTime(:)+stTime(1);
                BehavData.rewDeliveryTime(BehavData.rewDeliveryTime > 0) = BehavData.rewDeliveryTime(BehavData.rewDeliveryTime > 0) + stTime(1);
                BehavData.collectionTime(BehavData.collectionTime > 0) = BehavData.collectionTime(BehavData.collectionTime > 0) + stTime(1);

                final_behavior.(current_animal).(current_session).uv.BehavData = BehavData;



            else


                current_session = regexprep(current_session,{' ', '-'}, '_');
                [BehavData,ABETfile,Descriptives, block_end, largeRewSide, smallRewSide, forced_trial_start, free_trial_start]=ABET2TableFn_Chamber_A_v6(ABET_file,[]);
                if exist('boris_file', 'var')
                    SLEAP_time_range_adjustment = []; %16.2733; %15.3983; %[]; %-16.5448; %[]; %[]16.2733; -1.23;
                    [BehavData, boris_Extract_tbl] = boris_to_table(boris_file, BehavData, block_end, largeRewSide, smallRewSide, SLEAP_time_range_adjustment, forced_trial_start, free_trial_start);
                end
                ABET_removeheader = ABETfile(2:end,:);
                tbl_ABET = cell2table(ABET_removeheader);
                tbl_ABET.Properties.VariableNames = ABETfile(1,:);
                gpio_tbl = readtable(GPIO_file);

                % Inscopix stupidly recently updated the way the GPIO pins
                % are written, so the default values that I used before are
                % no longer valid. updated this to reflect the new way
                % things are written - hopefully this works! 
                % stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Time_s_ > 0);
                stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Value > 5000);



                BehavData.TrialPossible(:)=BehavData.TrialPossible(:)+stTime(1);
                BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+7.39500000000000;
                BehavData.collectionTime(:)=BehavData.collectionTime(:)+stTime(1);
                BehavData.stTime(:)=BehavData.stTime(:)+stTime(1);
                BehavData.collectionTime_end(:)=BehavData.collectionTime_end(:)+stTime(1);



                % shk_times(:)=shk_times(:)+stTime(1);

                BehavData.choTime2 = BehavData.choiceTime-BehavData.TrialPossible;


                final_behavior.(current_animal).(current_session).uv.BehavData = BehavData;

                clear unitTS unitTrace unitXTrials unitAVG unitSEM i zall zb zsd zb_window zsd_window zall_window zb_session zsd_session zall_session neuron boris_file ABET_file GPIO_file
            end

        else
            disp('Remember: BLANK folder was not analyzed due to missing files! Check contents and try again');

            % Skip the rest of the loop for this folder
            continue;
        end
        clear boris_file ABET_file GPIO_file
    end
    % clearvars -except final
end
%%

