%% Use this code to create a meta-structure 'final_SLEAP' that stores movement-related data

%% Edit these uservariables with what you want to look at
uv.evtWin = [-10 10]; %what time do you want to look at around each event
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at
ts1 = (-10:(uv.dt):10-0.1);
% Define the directory path you want to start with
% startDirectory = 'I:\MATLAB\Sean CNMFe\pan-neuronal BLA\BLA-Insc-24';

metaDirectory = 'D:\risk videos\BLA hM4Di vs mCherry\RRD424';
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
        current_animal = folder_list_string{5}; % Would have to change this depending on your folder structure, but there should be an animal name folder given our current workflow.
        current_animal = matlab.lang.makeValidName(current_animal);
        current_session = folder_list_string{6};
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
            else
                continue
            end
        end

        if ~isfield(final_SLEAP.(current_animal), current_session)
            final_SLEAP.(current_animal).(current_session) = struct;
        else
            continue
        end

        idx = ~cellfun('isempty',strfind({files.name},'.csv')); %find the instances of .xlsx in the file list.
        %This command converts the name field into a cell array and searches
        %the cell array with strfind
        csvFiles = files(idx); %build a mat file index
        clear idx
        csv_names = {csvFiles.name};

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

            % if ~contains(csv_names{mm}, 'ABET') | ~contains(lower(csv_names{mm}), '_gpio')
            %     disp('BLANK folder was not analyzed due to missing files! Check contents and try again');
            % 
            %     % Skip the rest of the loop for this folder
            %     continue;
            % end
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
                if exist('SLEAP_csv','var') && exist('ABET_file','var')
                    disp('Required files found, analyzing...');                           
                    if strcmp(current_session, 'SHOCK_TEST')
                        alignment_event = 'choiceTime';
                        [BehavData,ABETfile]=ABET2TableFn_ShockTest(ABET_file);
                        ABET_removeheader = ABETfile(2:end,:);
                        tbl_ABET = cell2table(ABET_removeheader);
                        tbl_ABET.Properties.VariableNames = ABETfile(1,:);
                        % gpio_tbl = readtable(GPIO_file);
                        % shk_times = tbl_ABET.Evnt_Time(strcmp(tbl_ABET.Item_Name, 'shock_on_off') & tbl_ABET.Arg1_Value == 1);
                        % stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Time_s_ > 0);
                        % frames = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName,'BNC Sync Output') & gpio_tbl.Value == 1);
                        %check GPIO file to extract each TTL, since the TTL is 1000ms and is
                        %sampled repeatedly. This will only extract events that are separated by >
                        %8sec, so be sure to change this if the TTL or task structure changes
                        %dramatically!
                        % pp = 2;
                        % ttl_filtered = stTime(1);
                        % for kk = 1:size(stTime,1)-1
                        %     if abs(stTime(kk)-stTime(kk+1)) > 8
                        %         ttl_filtered(pp) = stTime(kk+1);
                        %         pp=pp+1;
                        %     end
                        % end
                        % ttl_filtered = ttl_filtered';

                        BehavData.TrialPossible(:)=BehavData.TrialPossible(:);
                        BehavData.choiceTime(:)=BehavData.choiceTime(:); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+7.39500000000000;
                        % BehavData.collectionTime(:)=BehavData.collectionTime(:)+stTime(1);
                        % BehavData.stTime(:)=BehavData.stTime(:)+stTime(1);
                        % shk_times(:)=shk_times(:)+stTime(1);

                        BehavData.choTime2 = BehavData.choiceTime-BehavData.TrialPossible;
                        % BehavData.choTime3 = BehavData.Insc_TTL+BehavData.choTime2;

                        %filter based on TrialFilter inputs (see TrialFilter.m for full list of
                        %possibilities)
                        BehavData=TrialFilter(BehavData,'SHK',1);

                        % frames3 = frames(1:2:end-2);  %frames3 = frames(1:2:end-1) %frames3 = frames(1:2:end-2); the number of samples to skip (:#:) corresponds to the degree of temporal downsampling that the video underwent

                        % use this if you have specifi sessions within the loop that have
                        % the wrong acquisition rate (10 Hz vs 20 Hz)
                        % if strcmp(current_session,'RM_D1')
                        %     frames3 = frames(1:2:end-2);
                        % else
                        %     frames3 = frames(1:4:end-2);
                        %
                        % end


                        %%
                        
                        eTS = BehavData.(alignment_event); %get time stamps

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

                        SLEAP_time = uv.dt:uv.dt:height(SLEAP_data)*uv.dt; %generate time trace

                        %adjust  time to account for the fact that Inscopix recording
                        %starts first (stTime(1);
                        % SLEAP_time = SLEAP_time + stTime(1);

                        SLEAP_data.idx_time = SLEAP_time';

                        % if ~isempty(SLEAP_time_range_adjustment)
                        %     time_ranges = time_ranges-SLEAP_time_range_adjustment;
                        % end

                        SLEAP_data_vel_filtered_session = SLEAP_data.vel_filtered_2';

                        zscored_SLEAP_data_vel_filtered_session =(SLEAP_data_vel_filtered_session-mean(SLEAP_data_vel_filtered_session)./std(SLEAP_data_vel_filtered_session));
                        %%
                        tic
                        for i = 1 %could loop through multiple mice like this if you had it
                            % eTS = BehavData.(uv.behav); %get time stamps
                            velocity_trace =  SLEAP_data_vel_filtered_session;
                            %     ca = neuron.S; %get binarized calcium
                            velocity_time = SLEAP_data.idx_time';

                            % velocity_time = velocity_time + stTime(1);

                            %calculate time windows for each event
                            evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
                            numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
                            %%
                            tic
                            for u = 1:size(velocity_trace,1)
                                %% initialize trial matrices
                                velocity_trace_trials = NaN(size(eTS,1),numMeasurements); %
                                velocity_unitTrace = velocity_trace(u,:); %get trace
                                %             %%
                                for t = 1:size(eTS,1)
                                    %% set each trial's temporal boundaries
                                    timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
                                    BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
                                    if min(timeWin) > min(velocity_time) & max(timeWin) < max(velocity_time)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                                        %% get unit event counts in trials
                                        %% get unit ca traces in trials
                                        idx = velocity_time > min(timeWin) & velocity_time < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                                        bl_idx = velocity_time > min(BL_win) & velocity_time < max(BL_win);
                                        %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                                        velocity_trace_trials(t,1:sum(idx)) = velocity_unitTrace(idx);
                                        % zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
                                        % zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
                                        velocity_zb(t,:) = mean(velocity_trace_trials(t,:)); %baseline mean
                                        velocity_zsd(t,:) = std(velocity_trace_trials(t,:)); %baseline std


                                        tmp = 0;
                                        for j = 1:size(velocity_trace_trials,2)
                                            tmp = tmp+1;
                                            zall_motion(t,tmp) = (velocity_trace_trials(t,j) - velocity_zb(t))/velocity_zsd(t);
                                        end
                                        clear j;



                                    end
                                end
                                clear idx timeWin BL_win bl_idx
                                %%
                                unitXTrials(u).velocity_trace_trials = velocity_trace_trials;
                                unitXTrials(u).velocity_zb = velocity_zb;
                                unitXTrials(u).velocity_zsd = velocity_zsd;
                                unitXTrials(u).zall_motion = zall_motion;

                                %% store unit averaged data
                                unitAVG.velocity_traces(u,:) = nanmean(velocity_trace_trials);           %store trial averaged calcium traces
                                unitSEM.velocity_traces(u,:) = std(velocity_trace_trials,'omitnan')/sqrt(size(velocity_trace_trials,1));
                                clear caEvtCtTrials caTraceTrials caEvtRateTrials unitTrace idx
                            end
                        end


                        toc
                        %     final(i).name = mouseData(i).mouseID;
                        %     final(i).day = i;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).time = SLEAP_time; %final(i).time = caTime;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).unitAVG = unitAVG;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).unitXTrials = unitXTrials;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).uv = uv;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).unitSEM = unitSEM;
                        final_SLEAP.(current_animal).(current_session).BehavData = BehavData;
                        final_SLEAP.(current_animal).(current_session).SLEAP_data = SLEAP_data;
                        %Because the "neuron" data type is a pain in the ass to
                        %work with, we will instead save some of the variables
                        final_SLEAP.(current_animal).(current_session).SLEAP_data_velocity = SLEAP_data_vel_filtered_session;
                        final_SLEAP.(current_animal).(current_session).zscored_SLEAP_data_velocity = zscored_SLEAP_data_vel_filtered_session;


                        clear unitTS unitTrace unitXTrials unitAVG unitSEM i zall zb zsd zb_window zsd_window zall_window zb_session zsd_session zall_session neuron

                    else
                        for i = 1:length(uv.behav)


                            alignment_event = char(uv.behav(i));
                            [BehavData,ABETfile,Descriptives, block_end]=ABET2TableFn_Chamber_A_v6(ABET_file,[]);
                            ABET_removeheader = ABETfile(2:end,:);
                            tbl_ABET = cell2table(ABET_removeheader);
                            tbl_ABET.Properties.VariableNames = ABETfile(1,:);
                            % gpio_tbl = readtable(GPIO_file);
                            % shk_times = tbl_ABET.Evnt_Time(strcmp(tbl_ABET.Item_Name, 'shock_on_off') & tbl_ABET.Arg1_Value == 1);
                            % stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Time_s_ > 0);
                            % frames = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName,'BNC Sync Output') & gpio_tbl.Value == 1);

                            %check GPIO file to extract each TTL, since the TTL is 1000ms and is
                            %sampled repeatedly. This will only extract events that are separated by >
                            %8sec, so be sure to change this if the TTL or task structure changes
                            %dramatically!
                            % pp = 2;
                            % ttl_filtered = stTime(1);
                            % for kk = 1:size(stTime,1)-1
                            %     if abs(stTime(kk)-stTime(kk+1)) > 8
                            %         ttl_filtered(pp) = stTime(kk+1);
                            %         pp=pp+1;
                            %     end
                            % end
                            % ttl_filtered = ttl_filtered';


                            BehavData.TrialPossible(:)=BehavData.TrialPossible(:);
                            BehavData.choiceTime(:)=BehavData.choiceTime(:); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+7.39500000000000;
                            BehavData.collectionTime(:)=BehavData.collectionTime(:);
                            BehavData.stTime(:)=BehavData.stTime(:);

                            BehavData.choTime2 = BehavData.choiceTime-BehavData.TrialPossible;
                            % BehavData.choTime3 = BehavData.Insc_TTL+BehavData.choTime2;

                            %filter based on TrialFilter inputs (see TrialFilter.m for full list of
                            %possibilities)
                            BehavData=TrialFilter(BehavData,'ALL',1);

                            % frames3 = frames(1:2:end-2);  %frames3 = frames(1:2:end-1) %frames3 = frames(1:2:end-2); the number of samples to skip (:#:) corresponds to the degree of temporal downsampling that the video underwent

                            eTS = BehavData.(alignment_event); %get time stamps

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

                            SLEAP_time = uv.dt:uv.dt:height(SLEAP_data)*uv.dt; %generate time trace

                            %adjust  time to account for the fact that Inscopix recording
                            %starts first (stTime(1);
                            % SLEAP_time = SLEAP_time + stTime(1);

                            SLEAP_data.idx_time = SLEAP_time';

                            % if ~isempty(SLEAP_time_range_adjustment)
                            %     time_ranges = time_ranges-SLEAP_time_range_adjustment;
                            % end

                            SLEAP_data_vel_filtered_session = SLEAP_data.vel_filtered_2';

                            zscored_SLEAP_data_vel_filtered_session =(SLEAP_data_vel_filtered_session-mean(SLEAP_data_vel_filtered_session)./std(SLEAP_data_vel_filtered_session));
                            %%
                            tic
                            for i = 1 %could loop through multiple mice like this if you had it
                                % eTS = BehavData.(uv.behav); %get time stamps
                                velocity_trace =  SLEAP_data_vel_filtered_session;
                                %     ca = neuron.S; %get binarized calcium
                                velocity_time = SLEAP_data.idx_time';

                                % velocity_time = velocity_time + stTime(1);

                                %calculate time windows for each event
                                evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
                                numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
                                %%
                                tic
                                for u = 1:size(velocity_trace,1)
                                    %% initialize trial matrices
                                    velocity_trace_trials = NaN(size(eTS,1),numMeasurements); %
                                    velocity_unitTrace = velocity_trace(u,:); %get trace
                                    %             %%
                                    for t = 1:size(eTS,1)
                                        %% set each trial's temporal boundaries
                                        timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
                                        BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
                                        if min(timeWin) > min(velocity_time) & max(timeWin) < max(velocity_time)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                                            %% get unit event counts in trials
                                            %% get unit ca traces in trials
                                            idx = velocity_time > min(timeWin) & velocity_time < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                                            bl_idx = velocity_time > min(BL_win) & velocity_time < max(BL_win);
                                            %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                                            velocity_trace_trials(t,1:sum(idx)) = velocity_unitTrace(idx);
                                            % zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
                                            % zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
                                            velocity_zb(t,:) = mean(velocity_trace_trials(t,:)); %baseline mean
                                            velocity_zsd(t,:) = std(velocity_trace_trials(t,:)); %baseline std


                                            tmp = 0;
                                            for j = 1:size(velocity_trace_trials,2)
                                                tmp = tmp+1;
                                                zall_motion(t,tmp) = (velocity_trace_trials(t,j) - velocity_zb(t))/velocity_zsd(t);
                                            end
                                            clear j;



                                        end
                                    end
                                    clear idx timeWin BL_win bl_idx
                                    %%
                                    unitXTrials(u).velocity_trace_trials = velocity_trace_trials;
                                    unitXTrials(u).velocity_zb = velocity_zb;
                                    unitXTrials(u).velocity_zsd = velocity_zsd;
                                    unitXTrials(u).zall_motion = zall_motion;

                                    %% store unit averaged data
                                    unitAVG.velocity_traces(u,:) = nanmean(velocity_trace_trials);           %store trial averaged calcium traces
                                    unitSEM.velocity_traces(u,:) = std(velocity_trace_trials,'omitnan')/sqrt(size(velocity_trace_trials,1));
                                    clear caEvtCtTrials caTraceTrials caEvtRateTrials unitTrace idx
                                end
                            end


                        toc
                        %     final(i).name = mouseData(i).mouseID;
                        %     final(i).day = i;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).time = SLEAP_time; %final(i).time = caTime;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).unitAVG = unitAVG;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).unitXTrials = unitXTrials;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).uv = uv;
                        final_SLEAP.(current_animal).(current_session).(alignment_event).unitSEM = unitSEM;
                        final_SLEAP.(current_animal).(current_session).BehavData = BehavData;
                        final_SLEAP.(current_animal).(current_session).SLEAP_data = SLEAP_data;
                        %Because the "neuron" data type is a pain in the ass to
                        %work with, we will instead save some of the variables
                        final_SLEAP.(current_animal).(current_session).SLEAP_data_velocity = SLEAP_data_vel_filtered_session;
                        final_SLEAP.(current_animal).(current_session).zscored_SLEAP_data_velocity = zscored_SLEAP_data_vel_filtered_session;


                        clear unitTS unitTrace unitXTrials unitAVG unitSEM i zall zb zsd zb_window zsd_window zall_window zb_session zsd_session zall_session neuron
                        end
                    end
                else
                    disp('Remember: folder was not analyzed due to missing files! Check contents and try again');

                    % Skip the rest of the loop for this folder
                    continue;
                end
                
            end
            clear ABET_file SLEAP_csv GPIO_file
            % clearvars -except final
    end
end

