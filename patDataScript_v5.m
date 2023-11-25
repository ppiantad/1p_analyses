%% Edit these uservariables with what you want to look at
uv.evtWin = [-10 10]; %what time do you want to look at around each event
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at
ts1 = (-10:.1:10-0.1);
% Define the directory path you want to start with
startDirectory = 'I:\MATLAB\Sean CNMFe\BLA-NAcSh\BLA-Insc-19';

% Use the dir function to get a list of subfolders
subfolders = dir(startDirectory);


% Initialize folder_list as an empty cell array
folder_list = {};


% Loop through the list of subfolders
for i = 1:length(subfolders)
    % Check if the item in subfolders is a directory (not "." or "..") or
    % one of the sets of files that I haven't analyzed yet (PR currently)
    if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..') && ~contains(subfolders(i).name, 'PR')
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
    if isempty(matFiles) || isempty(csvFiles)
        disp('Missing .mat file, skipping folder');
    else
        currentMatFile = strcat(folder_list{ii}, '\', matFiles.name);
        filesFound = true; % Set the flag to true since .mat files were found

        % Add your code to process .mat files here

    end
    % Check the filesFound flag and print the final message
    if filesFound
        disp('Folder analyzed successfully');



        folder_strings = strsplit(folder_list{ii}, '\');
        %     session_strings = strsplit(folder_strings{end}, {'-', ' '});
        %     mat_strings = strsplit(char(matFiles.name),'_');
        %     date_strings = strsplit(mat_strings{4}, '-');
        csv_names = {csvFiles.name};
        current_animal = folder_strings{5}; % Would have to change this depending on your folder structure, but there should be an animal name folder given our current workflow.
        current_session = folder_strings{6};
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
        end


        %     ABET_file = strcat(folder_list{ii}, '\', csvFiles(1).name);
        %     GPIO_file = strcat(folder_list{ii}, '\', csvFiles(2).name);
        if strcmp(current_session, 'SHOCK TEST')
            alignment_event = 'choiceTime';
            [BehavData,ABETfile]=ABET2TableFn_ShockTest(ABET_file);
            ABET_removeheader = ABETfile(2:end,:);
            tbl_ABET = cell2table(ABET_removeheader);
            tbl_ABET.Properties.VariableNames = ABETfile(1,:);
            gpio_tbl = readtable(GPIO_file);
            shk_times = tbl_ABET.Evnt_Time(strcmp(tbl_ABET.Item_Name, 'shock_on_off') & tbl_ABET.Arg1_Value == 1);
            stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Time_s_ > 0);
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


            load(currentMatFile);
            % load(currentMatFile);
            % ts1 = uv.dt:uv.dt:length(neuron.C_raw)*uv.dt;

            %create array of FRAMES for aligning
            %NEED TO FIGURE OUT HOW TO MAKE THE FRAMES AUTOMATICALLY = THE SAME SIZE AS
            %THE neuron.C_ ARRAY LENGTH
            length_ca_trace = size(neuron.C,2);
            trim_frames = size(frames(1:2:end),1)-length_ca_trace;

            frames3 = frames(1:2:end-2);  %frames3 = frames(1:2:end-1) %frames3 = frames(1:2:end-2); the number of samples to skip (:#:) corresponds to the degree of temporal downsampling that the video underwent

            % use this if you have specifi sessions within the loop that have
            % the wrong acquisition rate (10 Hz vs 20 Hz)
            % if strcmp(current_session,'RM_D1')
            %     frames3 = frames(1:2:end-2);
            % else
            %     frames3 = frames(1:4:end-2);
            %
            % end


            %for testing frames vs. ts1
            % figure; plot(ts1, neuron.C_raw(2,:))
            % figure; plot(frames3, neuron.C_raw(2,:))
            % BehavData.choiceTime = BehavData.choiceTime-1;
            %%

            eTS = BehavData.(alignment_event); %get time stamps
            ca = neuron.C_raw; %get calcium
            zb_session = mean(ca,2);
            zsd_session = std(ca,[],2);
            caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


            %calculate time windows for each event
            evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
            numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
            %%
            tic
            for u = 1:size(ca,1)
                %% initialize trial matrices
                caTraceTrials = NaN(size(eTS,1),numMeasurements); %
                unitTrace = ca(u,:); %get trace
                %         current_animal = char(upper(mat_strings(1)));
                current_animal = matlab.lang.makeValidName(current_animal);
                current_session = char(folder_strings(end));
                current_session = regexprep(current_session,{' ', '-'}, '_');
                %         current_session = matlab.lang.makeValidName(current_session);

                %             %%
                for t = 1:size(eTS,1)
                    %% set each trial's temporal boundaries
                    timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
                    BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
                    if min(timeWin) > min(frames3) & max(timeWin) < max(frames3)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                        %% get unit event counts in trials
                        %% get unit ca traces in trials
                        idx = frames3 > min(timeWin) & frames3 < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                        bl_idx = frames3 > min(BL_win) & frames3 < max(BL_win);
                        %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                        caTraceTrials(t,1:sum(idx)) = unitTrace(idx);
                        zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
                        zb_window(t,:) = mean(caTraceTrials(t,:));
                        zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
                        zsd_window(t,:) = std(caTraceTrials(t,:));
                        tmp = 0;
                        for j = 1:size(caTraceTrials,2)
                            tmp = tmp+1;
                            zall(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
                            zall_window(t,tmp) = (caTraceTrials(t,j) - zb_window(t))/zsd_window(t);
                            zall_session(t,tmp) = (caTraceTrials(t,j) - zb_session(u))/zsd_session(u);
                        end
                        clear j;



                    end
                end
                clear idx timeWin BL_win bl_idx
                %%
                unitXTrials(u).caTraces = caTraceTrials;
                unitXTrials(u).zb = zb;
                unitXTrials(u).zb_window = zb_window;
                unitXTrials(u).zb_session = zb_session;
                unitXTrials(u).zsd = zsd;
                unitXTrials(u).zsd_window = zsd_window;
                unitXTrials(u).zsd_session = zsd_session;
                unitXTrials(u).zall = zall;
                unitXTrials(u).zall_window = zall_window;
                unitXTrials(u).zall_session = zall_session;


                %% store unit averaged data
                unitAVG.caTraces(u,:) = nanmean(caTraceTrials);           %store trial averaged calcium traces
                unitSEM.caTraces(u,:) = std(caTraceTrials,'omitnan')/sqrt(size(caTraceTrials,1));
                clear caEvtCtTrials caTraceTrials caEvtRateTrials unitTrace idx
            end


            toc
            %     final(i).name = mouseData(i).mouseID;
            %     final(i).day = i;
            final.(current_animal).(current_session).(alignment_event).time = frames3; %final(i).time = caTime;
            final.(current_animal).(current_session).(alignment_event).unitAVG = unitAVG;
            final.(current_animal).(current_session).(alignment_event).unitXTrials = unitXTrials;
            final.(current_animal).(current_session).(alignment_event).uv = uv;
            final.(current_animal).(current_session).(alignment_event).unitSEM = unitSEM;
            final.(current_animal).(current_session).(alignment_event).uv.BehavData = BehavData;
            final.(current_animal).(current_session).neuron = neuron;

            clear unitTS unitTrace unitXTrials unitAVG unitSEM i zall zb zsd zb_window zsd_window zall_window zb_session zsd_session zall_session neuron

        else
            for i = 1:length(uv.behav)


                alignment_event = char(uv.behav(i));
                [BehavData,ABETfile,Descriptives, block_end]=ABET2TableFn_Chamber_A_v6(ABET_file,[]);
                ABET_removeheader = ABETfile(2:end,:);
                tbl_ABET = cell2table(ABET_removeheader);
                tbl_ABET.Properties.VariableNames = ABETfile(1,:);
                gpio_tbl = readtable(GPIO_file);
                shk_times = tbl_ABET.Evnt_Time(strcmp(tbl_ABET.Item_Name, 'shock_on_off') & tbl_ABET.Arg1_Value == 1);
                stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Time_s_ > 0);
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

                BehavData.TrialPossible(:)=BehavData.TrialPossible(:)+stTime(1);
                BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+7.39500000000000;
                BehavData.collectionTime(:)=BehavData.collectionTime(:)+stTime(1);
                BehavData.stTime(:)=BehavData.stTime(:)+stTime(1);



                % shk_times(:)=shk_times(:)+stTime(1);

                BehavData.choTime2 = BehavData.choiceTime-BehavData.TrialPossible;
                % BehavData.choTime3 = BehavData.Insc_TTL+BehavData.choTime2;

                %filter based on TrialFilter inputs (see TrialFilter.m for full list of
                %possibilities)
                BehavData=TrialFilter(BehavData,'ALL',1);


                load(currentMatFile);
                % load(currentMatFile);
                % ts1 = uv.dt:uv.dt:length(neuron.C_raw)*uv.dt;

                %create array of FRAMES for aligning
                %NEED TO FIGURE OUT HOW TO MAKE THE FRAMES AUTOMATICALLY = THE SAME SIZE AS
                %THE neuron.C_ ARRAY LENGTH
                length_ca_trace = size(neuron.C,2);
                trim_frames = size(frames(1:2:end),1)-length_ca_trace;

                frames3 = frames(1:2:end-2);  %frames3 = frames(1:2:end-1) %frames3 = frames(1:2:end-2); the number of samples to skip (:#:) corresponds to the degree of temporal downsampling that the video underwent

                % use this if you have specifi sessions within the loop that have
                % the wrong acquisition rate (10 Hz vs 20 Hz)
                % if strcmp(current_session,'RM_D1')
                %     frames3 = frames(1:2:end-2);
                % else
                %     frames3 = frames(1:4:end-2);
                %
                % end


                %for testing frames vs. ts1
                % figure; plot(ts1, neuron.C_raw(2,:))
                % figure; plot(frames3, neuron.C_raw(2,:))
                % BehavData.choiceTime = BehavData.choiceTime-1;

                %%

                eTS = BehavData.(alignment_event); %get time stamps
                ca = neuron.C_raw; %get calcium
                zb_session = mean(ca,2);
                zsd_session = std(ca,[],2);
                caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


                %calculate time windows for each event
                evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
                numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
                %%
                tic
                for u = 1:size(ca,1)
                    %% initialize trial matrices
                    caTraceTrials = NaN(size(eTS,1),numMeasurements); %
                    unitTrace = ca(u,:); %get trace
                    %         current_animal = char(upper(mat_strings(1)));
                    current_animal = matlab.lang.makeValidName(current_animal);
                    current_session = char(folder_strings(end));
                    current_session = regexprep(current_session,{' ', '-'}, '_');
                    %         current_session = matlab.lang.makeValidName(current_session);

                    %             %%
                    for t = 1:size(eTS,1)
                        %% set each trial's temporal boundaries
                        timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
                        BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
                        if min(timeWin) > min(frames3) & max(timeWin) < max(frames3)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                            %% get unit event counts in trials
                            %% get unit ca traces in trials
                            idx = frames3 > min(timeWin) & frames3 < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                            bl_idx = frames3 > min(BL_win) & frames3 < max(BL_win);
                            %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                            caTraceTrials(t,1:sum(idx)) = unitTrace(idx);
                            zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
                            zb_window(t,:) = mean(caTraceTrials(t,:));
                            zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
                            zsd_window(t,:) = std(caTraceTrials(t,:));
                            tmp = 0;
                            for j = 1:size(caTraceTrials,2)
                                tmp = tmp+1;
                                zall(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
                                zall_window(t,tmp) = (caTraceTrials(t,j) - zb_window(t))/zsd_window(t);
                                zall_session(t,tmp) = (caTraceTrials(t,j) - zb_session(u))/zsd_session(u);
                            end
                            clear j;



                        end
                    end
                    clear idx timeWin BL_win bl_idx
                    %%
                    unitXTrials(u).caTraces = caTraceTrials;
                    unitXTrials(u).zb = zb;
                    unitXTrials(u).zb_window = zb_window;
                    unitXTrials(u).zb_session = zb_session;
                    unitXTrials(u).zsd = zsd;
                    unitXTrials(u).zsd_window = zsd_window;
                    unitXTrials(u).zsd_session = zsd_session;
                    unitXTrials(u).zall = zall;
                    unitXTrials(u).zall_window = zall_window;
                    unitXTrials(u).zall_session = zall_session;


                    %% store unit averaged data
                    unitAVG.caTraces(u,:) = nanmean(caTraceTrials);           %store trial averaged calcium traces
                    unitSEM.caTraces(u,:) = std(caTraceTrials,'omitnan')/sqrt(size(caTraceTrials,1));
                    clear caEvtCtTrials caTraceTrials caEvtRateTrials unitTrace idx
                end


                toc
                %     final(i).name = mouseData(i).mouseID;
                %     final(i).day = i;
                final.(current_animal).(current_session).(alignment_event).time = frames3; %final(i).time = caTime;
                final.(current_animal).(current_session).(alignment_event).unitAVG = unitAVG;
                final.(current_animal).(current_session).(alignment_event).unitXTrials = unitXTrials;
                final.(current_animal).(current_session).(alignment_event).uv = uv;
                final.(current_animal).(current_session).(alignment_event).unitSEM = unitSEM;
                final.(current_animal).(current_session).(alignment_event).uv.BehavData = BehavData;
                final.(current_animal).(current_session).neuron = neuron;

                clear unitTS unitTrace unitXTrials unitAVG unitSEM i zall zb zsd zb_window zsd_window zall_window zb_session zsd_session zall_session neuron
            end
        end
    else
        disp('Remember: BLANK folder was not analyzed due to missing files! Check contents and try again');

        % Skip the rest of the loop for this folder
        continue;
    end
end
% clearvars -except final

%%
% for c = 1:size(ca,2)
%     subplot(121)
%     imagesc(final.unitXTrials(c).caTraces)
%     xticklabels = (final.uv.evtWin(1,1)):5:(final.uv.evtWin(1,2));
%     xticks = linspace(1, length(final.unitAVG.caTraces), numel(xticklabels));
%     set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
%     set(gca,'FontName','Arial','FontSize',16)
%     xlabel('Time from(s)','FontSize',22)
%     ylabel('Trial number','FontSize',22)
%     subplot(122)
%     trialTime = uv.evtWin(1,1):uv.dt:uv.evtWin(1,2)-uv.dt;
%     plot(trialTime,final.unitAVG.caTraces(c,:))
%     set(gca,'FontName','Arial','FontSize',16)
%     xlabel('Time from  (s)','FontSize',22)
%     xlim([uv.evtWin(1,1) uv.evtWin(1,2)])
%     xline(0,'--')
%     pause
% end
