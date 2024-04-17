%% Edit these uservariables with what you want to look at
%this script has been verified as of 7/30/2022 to be well-aligned based on
%shock-responsive neurons from multiple sessions for multiple mice
uv.evtWin = [-10 10]; %what time do you want to look at around each event [-10 10] [-5 35]
uv.BLper = [-10 -5]; %what baseline period do you want for z-score [-10 -5] [-5 0]
uv.dt = 0.1; %what is your frame rate (check neuron.Fs to be sure) 0.2 0.1
uv.behav = 'choiceTime'; %which behavior/timestamp to look at choiceTime stTime


load('BLA-Insc-9_Session-20211025-100456_bla-insc-9-RM-D12.CNMF_final.mat');
[BehavData,ABETfile,Descriptives, block_end, largeRewSide, smallRewSide]=ABET2TableFn_Chamber_A_v6('BLA-Insc-9 10252021 ABET.csv',[]);
gpio_tbl = readtable('Session-20211025-100456_bla-insc-9-RM-D12_GPIO.csv');

SLEAP_data = readtable('BLA-Insc-9_RM D12_body_sleap_data.csv');
%EDIT FOR EACH MOUSE AS NECESSARY
SLEAP_time_range_adjustment = []; %16.2733; %15.3983; %[]; %-16.5448; %[]; %[]16.2733; -1.23;

boris_file = []; %'BLA-Insc-39_RDT_D1_2023-05-11T14_32_33.avi. BORIS.csv'; %'BLA-Insc-27_RDT_D1.csv';

[BehavData, boris_Extract_tbl] = boris_to_table(boris_file, BehavData, block_end, largeRewSide, smallRewSide, SLEAP_time_range_adjustment);



% [BehavData,ABETfile]=ABET2TableFn_ShockTest('BLA-Insc-3 01272021.csv');


ABET_removeheader = ABETfile(2:end,:);

tbl_ABET = cell2table(ABET_removeheader);
tbl_ABET.Properties.VariableNames = ABETfile(1,:);



shk_times = tbl_ABET.Evnt_Time(strcmp(tbl_ABET.Item_Name, 'shock_on_off') & tbl_ABET.Arg1_Value == 1);

% this filter needs to be set based on the parameters of the GPIO file
% written by Inscopix. occasionally the change things, so make sure this
% works regularly! 
stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Value > 5000); 

frames = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName,'BNC Sync Output') & gpio_tbl.Value == 1);
frames_test = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName,'BNC Sync Output'));
frames_test_2 = frames_test(2:2:end);

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
BehavData.Insc_TTL = zeros(length(BehavData.TrialPossible),1);
dd = 2;
for cc = 1:size(BehavData, 1)
    if BehavData.TrialPossible(cc) > stTime(1)
        BehavData.Insc_TTL(cc) = ttl_filtered(dd);
        dd = dd+1;
    elseif BehavData.TrialPossible(cc) <= stTime(1)
        BehavData.Insc_TTL(cc) = 0;
    end
end

BehavData.TrialPossible(:)=BehavData.TrialPossible(:)+stTime(1);
BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+7.39500000000000;
BehavData.collectionTime(:)=BehavData.collectionTime(:)+stTime(1);
BehavData.stTime(:)=BehavData.stTime(:)+stTime(1);



% shk_times(:)=shk_times(:)+stTime(1);

BehavData.choTime2 = BehavData.choiceTime-BehavData.TrialPossible;
BehavData.choTime3 = BehavData.Insc_TTL+BehavData.choTime2;

%filter based on TrialFilter inputs (see TrialFilter.m for full list of
%possibilities)
BehavData=TrialFilter(BehavData,'ALL',1); %BehavData=TrialFilter(BehavData,'OMITALL',0, 'BLANK_TOUCH', 0);


% BehavData_for_SLEAP = BehavData;
% BehavData_for_SLEAP.TrialPossible(:)=BehavData_for_SLEAP.TrialPossible(:)-stTime(1);
% BehavData_for_SLEAP.choiceTime(:)=BehavData_for_SLEAP.choiceTime(:)-stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1); %BehavData.choiceTime(:)=BehavData.choiceTime(:)+7.39500000000000;
% BehavData_for_SLEAP.collectionTime(:)=BehavData_for_SLEAP.collectionTime(:)-stTime(1);
% BehavData_for_SLEAP.stTime(:)=BehavData_for_SLEAP.stTime(:)-stTime(1);
% adjust SLEAP timestamps because Inscopix is started first (with the start
% time indicated in stTime(1)
% SLEAP_data.idx_time(:) = SLEAP_data.idx_time(:) + stTime(1);

% % Assuming 'Timestamp' is the timestamp variable in your table
% timestamps = SLEAP_data.idx_time;
% 
% % Specify the original and target sampling rates
% originalSamplingRate = 30; % Hz
% targetSamplingRate = 10;    % Hz
% 
% % Calculate the downsampled indices
% downsampledIndices = round(linspace(1, height(SLEAP_data), height(SLEAP_data) / (originalSamplingRate / targetSamplingRate)));
% 
% % Downsample the table
% SLEAP_downsampled_data = SLEAP_data(downsampledIndices, :);
% 
% SLEAP_data = SLEAP_downsampled_data;
% 




% ts1 = uv.dt:uv.dt:length(neuron.C_raw)*uv.dt;

%create array of FRAMES for aligning
%NEED TO FIGURE OUT HOW TO MAKE THE FRAMES AUTOMATICALLY = THE SAME SIZE AS
%THE neuron.C_ ARRAY LENGTH
length_ca_trace = size(neuron.C,2);
trim_frames = size(frames(1:2:end),1)-length_ca_trace;

frames3 = frames(1:2:end-2); % frames3 = frames(1:2:end-2);  %frames3 = frames(1:2:end-1) %frames3 = frames(1:2:end-2); the number of samples to skip (:#:) corresponds to the degree of temporal downsampling that the video underwent
% frames3 = frames_test_2(1:4:end);



%for testing frames vs. ts1
% figure; plot(ts1, neuron.C_raw(2,:))
% figure; plot(frames3, neuron.C_raw(2,:))
% BehavData.choiceTime = BehavData.choiceTime-1;

%%
for i = 1 %could loop through multiple mice like this if you had it
    eTS = BehavData.(uv.behav); %periods_with_low_velocity; %BehavData.(uv.behav); %get time stamps
    ca = neuron.C_raw; %get calcium
%     ca = neuron.S; %get binarized calcium
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
        %             %%
        for t = 1:size(eTS,1)
            %% set each trial's temporal boundaries
            timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
            BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
            if min(timeWin) > min(frames3) & max(timeWin) < max(frames3)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                %% get unit event counts in trials
                %% get unit ca traces in trials
                idx = frames3 > min(timeWin) & frames3 < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
                sum(idx)
                bl_idx = frames3 > min(BL_win) & frames3 < max(BL_win);
                %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                caTraceTrials(t,1:sum(idx)) = unitTrace(idx);
                % zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
                % zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
                zb(t,:) = mean(caTraceTrials(t,:)); %baseline mean
                zsd(t,:) = std(caTraceTrials(t,:)); %baseline std


                tmp = 0;
                for j = 1:size(caTraceTrials,2)
                    tmp = tmp+1;
                    zall(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t); 
                end
                clear j;
            
             
                
            end
        end
        clear idx timeWin BL_win bl_idx
        %% 
        unitXTrials(u).caTraces = caTraceTrials;
        unitXTrials(u).zb = zb;
        unitXTrials(u).zsd = zsd;
        unitXTrials(u).zall = zall;

        %% store unit averaged data
        unitAVG.caTraces(u,:) = nanmean(caTraceTrials);           %store trial averaged calcium traces
        unitSEM.caTraces(u,:) = std(caTraceTrials,'omitnan')/sqrt(size(caTraceTrials,1));
        clear caEvtCtTrials caTraceTrials caEvtRateTrials unitTrace idx
    end
end

toc
%     final(i).name = mouseData(i).mouseID;
%     final(i).day = i;
final(i).time = frames3; %final(i).time = caTime;
final(i).unitAVG = unitAVG;
final(i).unitXTrials = unitXTrials;
final(i).uv = uv;
final(i).unitSEM = unitSEM;
final(i).uv.BehavData = BehavData;


clear unitTS unitTrace unitXTrials unitAVG unitSEM i zall zb zsd



%%
time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);

[numTrials,~]=size(BehavData.collectionTime(:));
Tris=[1:numTrials]';

% window_ts = uv.dt:uv.dt:numMeasurements*uv.dt;

window_ts2 = uv.evtWin(1):uv.dt:uv.evtWin(2)-uv.dt;
window_ts3 = uv.evtWin(1):uv.dt:uv.evtWin(2);


for ii = 1:size(final.unitXTrials,2)
    for jj = 1:size(final.unitXTrials(ii).zall)
        final.unitAVG.zscored_caTraces(ii,:) = mean(final.unitXTrials(ii).zall);
    end
end
    
 figure
imagesc(window_ts3, 1, final.unitAVG.zscored_caTraces);hold on;   
    
% B = sortrows(final.unitAVG.zscored_caTraces,[50:80]);
% 
%  figure
% imagesc(window_ts3, 1, B);hold on; 


%% rank order traces for heatmap based on some subwindow of larger window
% set an ind to get a particular subwindow
ind_zero_to_two = window_ts3 > 0 & window_ts3 < 2;
% get mean of this subwindow
mean_zero_to_two = mean(final.unitAVG.zscored_caTraces(:,ind_zero_to_two),2);
% add mean to the zscored data (temporarily)
add_mean_to_zscored_data = [mean_zero_to_two final.unitAVG.zscored_caTraces];
% rank order by the mean column, which is column 1
rank_by_c1 = sortrows(add_mean_to_zscored_data, 1);
% get rid of mean column
rank_ordered_mean_zscore = rank_by_c1(:,2:end);

figure
imagesc(window_ts3, 1, rank_ordered_mean_zscore);hold on;



%% Use this code to plot heatmaps for each individual cell, across trials for all levels of iter
% **most useful for plotting matched cells within the same experiment, e.g., pan-neuronal matched Pre-RDT RM vs. RDT D1**

% for ii = 1:size(final.unitXTrials, 2)
%     figure;
%     % Initialize variables to store global max and min for heatmap and line graph
%     globalMaxHeatmap = -inf;
%     globalMinHeatmap = inf;
%     globalMaxYLine = -inf;
%     globalMinYLine = inf;
% 
%     % Create subplot with 1 row and 2 columns
%     num_columns_plot = 1;
%     subplot(2, num_columns_plot, 1);
% 
%     % Generate the heatmap
%     imagesc(ts1, 1, final.unitXTrials(ii).zall);
%     hold on;
%     % title({"Cell from " + strrep(mouse_cells{qq, ii}, '_', '-'), "Classified as " + respClass_all_array{qq}(ii), "(Overall cell number = " + (ii) + ")"}, 'FontSize', 9)
% 
%     % Update global max and min for heatmap
%     localMaxHeatmap = max(final.unitXTrials(ii).zall);
%     localMinHeatmap = min(final.unitXTrials(ii).zall);
%     globalMaxHeatmap = max(globalMaxHeatmap, localMaxHeatmap);
%     globalMinHeatmap = min(globalMinHeatmap, localMinHeatmap);
% 
% 
%     time2Collect = BehavData.collectionTime(:) - BehavData.choiceTime(:);
%     trialStartTime = BehavData.stTime(:) - BehavData.choiceTime(:);
%     [numTrials, ~] = size(BehavData.collectionTime(:));
%     Tris = [1:numTrials]';
%     % scatter(time2Collect, Tris               , 'Marker', 'p', 'MarkerFaceColor', 'w', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.7)
%     % scatter(trialStartTime, Tris, 'Marker', 's', 'MarkerFaceColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.7)
%     plot(zeros(numTrials, 1), Tris, 'LineWidth', 3, 'LineStyle', "--", 'Color', 'w')
% 
% colorbar;
% 
% hold off;
% clear time2Collect trialStartTime numTrials Tris behav_tbl
% 
% % Create subplot for the mean and raw data
% subplot(2, num_columns_plot, 1);
% 
% % Plot the mean as a thick black line
% meanData = mean(final.unitXTrials(ii).zall);
% plot(ts1, meanData(:, 1:end-1), 'k', 'LineWidth', 2);
% hold on;
% 
% % Plot the raw data in grey with transparency
% for trial = 1:size(final.unitXTrials(ii).zall, 1)
%     plot(ts1, final.unitXTrials(ii).zall(trial, :), 'Color', [0.1, 0.1, 0.1, 0.1]);
%     hold on;
% end
% 
% title("Mean and Raw Data", 'FontSize', 9)
% % Update global max and min for line graph
% localMaxYLine = max(final.unitXTrials(ii).zall, [], "all");
% localMinYLine = min(final.unitXTrials(ii).zall, [], "all");
% globalMaxYLine = max(globalMaxYLine, localMaxYLine);
% globalMinYLine = min(globalMinYLine, localMinYLine);
% 
% hold off;
% 
% % Set the same colorbar scale for all heatmap subplots
% for qq = 1:iter
%     subplot(2, num_columns_plot, 1);
%     clim([globalMinHeatmap, globalMaxHeatmap]);
%     colorbar;
% end
% 
% % Set the same Y-axis scale for all line graph subplots
% for qq = 1:iter
%     subplot(2, num_columns_plot, 1+1);
%     ylim([globalMinYLine, globalMaxYLine]);
% 
% 
%     pause
%     hold
%     close
% end



%%



cell_select = 10


figure
imagesc(window_ts2, 1, final.unitXTrials(cell_select).zall);hold on;
scatter(time2Collect,Tris,'Marker','p','MarkerFaceColor','w')
plot(zeros(numTrials,1),Tris)
% xticklabels = (final.uv.evtWin(1,1)):5:(final.uv.evtWin(1,2));
% xticks = linspace(1, length(final.unitAVG.caTraces), numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% set(gca,'FontName','Arial','FontSize',16)
xlabel('Time from Large Rew Choice (s)','FontSize',22)
ylabel('Trial number','FontSize',22)
title('Z-scored Ca Traces (normalized)')


figure;
plot(window_ts3, final.unitAVG.zscored_caTraces(cell_select,:));
title('Z-scored Ca Traces (normalized)')


% figure
% imagesc(window_ts2, 1, final.unitXTrials(animalID_select).caTraces);hold on;
% scatter(time2Collect,Tris,'Marker','p','MarkerFaceColor','w')
% plot(zeros(numTrials,1),Tris)
% % xticklabels = (final.uv.evtWin(1,1)):5:(final.uv.evtWin(1,2));
% % xticks = linspace(1, length(final.unitAVG.caTraces), numel(xticklabels));
% % set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% % set(gca,'FontName','Arial','FontSize',16)
% xlabel('Time from Large Rew Choice (s)','FontSize',22)
% ylabel('Trial number','FontSize',22)
% title('Ca Traces (not normalized)')
% 
% 
% figure;
% plot(window_ts3, final.unitAVG.caTraces(animalID_select,:));
% title('Ca Traces (not normalized)')





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


%% organizing data for PCA

large_rew_ind = BehavData.bigSmall == 1.2;

B1_Large_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 1;
B2_Large_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 2;
B3_Large_ind = BehavData.bigSmall == 1.2 & BehavData.Block == 3;

for ii = 1:length(final.unitXTrials)
    for kk = 1:length(final.unitXTrials(ii).zall)
        B1_large_rew_mean(ii,:) = mean(final.unitXTrials(ii).zall(B1_Large_ind, :));
        B1_small_rew_mean(ii,:) = mean(final.unitXTrials(ii).zall(~B1_Large_ind, :));
        B2_large_rew_mean(ii,:) = mean(final.unitXTrials(ii).zall(B2_Large_ind, :));
        B2_small_rew_mean(ii,:) = mean(final.unitXTrials(ii).zall(~B2_Large_ind, :));
        B3_large_rew_mean(ii,:) = mean(final.unitXTrials(ii).zall(B3_Large_ind, :));
        B3_small_rew_mean(ii,:) = mean(final.unitXTrials(ii).zall(~B3_Large_ind, :));
        
        
        
    end
    ii = ii+1;
end



B1_large_mean = mean(B1_large_rew_mean);
B2_large_mean = mean(B2_large_rew_mean);
B3_large_mean = mean(B3_large_rew_mean);

B1_large_SEM = nanstd(B1_large_rew_mean,1)/(sqrt(size(B1_large_rew_mean, 1)));
B2_large_SEM = nanstd(B2_large_rew_mean,1)/(sqrt(size(B2_large_rew_mean, 1)));
B3_large_SEM = nanstd(B3_large_rew_mean,1)/(sqrt(size(B3_large_rew_mean, 1)));




 figure
imagesc(window_ts3, 1, B1_large_rew_mean);hold on;  


 figure
imagesc(window_ts3, 1, B2_large_rew_mean);hold on;  


clims = [-4 4]
 figure
imagesc(window_ts3, 1, B3_large_rew_mean, clims);hold on;  

ind_zero_to_minus_3 = window_ts3 > 0 & window_ts3 < 5;
window_ind = window_ts3 > -5 & window_ts3 < 5;
x_window = window_ts3(:,window_ind);


B1_mean_zero_to_minus_3 = mean(B1_large_rew_mean(:,ind_zero_to_minus_3),2);
B1_add_mean_to_zscored_data = [B1_mean_zero_to_minus_3 B1_large_rew_mean];
B1_rank_by_c1 = sortrows(B1_add_mean_to_zscored_data, 1);
B1_rank_ordered_mean_zscore = B1_rank_by_c1(:,2:end);


 figure
imagesc(window_ts3, 1, B1_rank_ordered_mean_zscore,clims);hold on;  

B2_mean_zero_to_minus_3 = mean(B2_large_rew_mean(:,ind_zero_to_minus_3),2);
B2_add_mean_to_zscored_data = [B2_mean_zero_to_minus_3 B2_large_rew_mean];
B2_rank_by_c1 = sortrows(B2_add_mean_to_zscored_data, 1);
B2_rank_ordered_mean_zscore = B2_rank_by_c1(:,2:end);

 figure
imagesc(window_ts3, 1, B2_rank_ordered_mean_zscore,clims);hold on;  


B3_mean_zero_to_minus_3 = mean(B3_large_rew_mean(:,ind_zero_to_minus_3),2);
B3_add_mean_to_zscored_data = [B3_mean_zero_to_minus_3 B3_large_rew_mean];
B3_rank_by_c1 = sortrows(B3_add_mean_to_zscored_data, 1);
B3_rank_ordered_mean_zscore = B3_rank_by_c1(:,2:end);

 figure
imagesc(window_ts3, 1, B3_rank_ordered_mean_zscore,clims);hold on;  

B1_large_mean = mean(B1_large_rew_mean);
B2_large_mean = mean(B2_large_rew_mean);
B3_large_mean = mean(B3_large_rew_mean);

 figure
imagesc(window_ts3(:,window_ind), 1, B1_rank_ordered_mean_zscore(:,window_ind),clims);hold on;  

 figure
imagesc(window_ts3(:,window_ind), 1, B2_rank_ordered_mean_zscore(:,window_ind),clims);hold on;  

 figure
imagesc(window_ts3(:,window_ind), 1, B3_rank_ordered_mean_zscore(:,window_ind),clims);hold on;  




figure; 
shadedErrorBar(window_ts3(1,1:end-1), B1_large_mean(1,1:end-1), B1_large_SEM(1,1:end-1));
hold on;
shadedErrorBar(window_ts3(1,1:end-1), B2_large_mean(1,1:end-1), B2_large_SEM(1,1:end-1),'lineProps', '-b');
hold on;
shadedErrorBar(window_ts3(1,1:end-1), B3_large_mean(1,1:end-1), B3_large_SEM(1,1:end-1),'lineProps', '-r');
legend('block 1 (0%)','block 2 (50%)', 'block 3 (75%)')


%% PCA attempts
 large_concat = vertcat(B1_large_rew_mean(:,1:end-1), B2_large_rew_mean(:,1:end-1), B3_large_rew_mean(:,1:end-1));
 large_concat_horizontal = [B1_large_rew_mean(:,1:end-1), B2_large_rew_mean(:,1:end-1), B3_large_rew_mean(:,1:end-1)];
 all_concat = [B1_large_rew_mean(:,1:end-1), B2_small_rew_mean(:,1:end-1), B2_large_rew_mean(:,1:end-1), B2_small_rew_mean(:,1:end-1), B3_large_rew_mean(:,1:end-1), B3_small_rew_mean(:,1:end-1)];

[W,q,s] = svd(B2_large_rew_mean, 'econ');
W = W(:,1:20);

[W,~,~] = svd(all_concat, 'econ');
W = W(:,1:20);

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

[W,~,~] = svd(large_concat, 'econ');
W = W(:,1:20);

% minimal plotting
dpca_plot(large_concat, W, W, @dpca_plot_default);



% https://www.mathworks.com/matlabcentral/answers/270329-how-to-select-the-components-that-show-the-most-variance-in-pca

large_concat = large_concat - mean(large_concat);


[coeff,score,latent,~,explained] = pca(large_concat);


covarianceMatrix = cov(large_concat);
[V,D] = eig(covarianceMatrix);

dataInPrincipalComponentSpace = large_concat*coeff;

figure
bar(latent)

figure
bar(explained)

x_comp = 1; % Principal component for x-axis
y_comp = 2; % Principal component for y-axis
figure
hold on
for nc = 1:N
    h = plot([0 coeff(nc,x_comp)],[0 coeff(nc,y_comp)]);
    set(h,'Color','b')
end

% Scree plot
figure
h = plot(explained,'.-');
set(h,'LineWidth',3,'MarkerSize',36)
ylim([0 100])
set(gca,'XTick',1:N)
title('Explained variance by principal components')
xlabel('Principal component number')
ylabel('Fraction of variation explained [%]')
%%
% https://www.mathworks.com/help/bioinfo/ref/mapcaplot.html
mapcaplot(all_concat(:,1:101), all_concat(:,102:203));

