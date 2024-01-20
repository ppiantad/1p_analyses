% 
% SLEAP_data = readtable('BLA-Insc-32_Pre-RDT RM_body_sleap_data.csv');
% %EDIT FOR EACH MOUSE AS NECESSARY
% SLEAP_time_range_adjustment =  []; %16.2733; %15.3983; %[]; %-16.5448; %[]; %[]16.2733; 

% adjust SLEAP timestamps because Inscopix is started first (with the start
% time indicated in stTime(1)
% SLEAP_data.idx_time(:) = SLEAP_data.idx_time(:) + stTime(1);

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
SLEAP_time = SLEAP_time + stTime(1);

SLEAP_data.idx_time = SLEAP_time';

if ~isempty(SLEAP_time_range_adjustment)
    time_ranges = time_ranges-SLEAP_time_range_adjustment;
end

SLEAP_data_vel_filtered_session = SLEAP_data.vel_filtered_2';

zscored_SLEAP_data_vel_filtered_session =(SLEAP_data_vel_filtered_session-mean(SLEAP_data_vel_filtered_session)./std(SLEAP_data_vel_filtered_session));

%%
%%
for i = 1 %could loop through multiple mice like this if you had it
    eTS = BehavData.(uv.behav); %get time stamps
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


time2Collect = BehavData.collectionTime(:) - BehavData.choiceTime(:);
trialStartTime = BehavData.stTime(:) - BehavData.choiceTime(:);

figure
imagesc(window_ts3, 1, zall_motion);hold on;  

[numTrials, ~] = size(BehavData.collectionTime(:));
Tris = [1:numTrials]';
scatter(time2Collect, Tris               , 'Marker', 'p', 'MarkerFaceColor', 'w', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha',1)
scatter(trialStartTime, Tris, 'Marker', 's', 'MarkerFaceColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 1)
plot(zeros(numTrials, 1), Tris, 'LineWidth', 3, 'LineStyle', "--", 'Color', 'w')
yline(30.5,'--',{'Block 2', 'begins'},'LineWidth',3);
yline(60.5,'--',{'Block 3', 'begins'},'LineWidth',3);
colorbar;
hold off;

figure;
plot(window_ts2, velocity_trace_trials);


for qq = 1:3
    block_mean(qq,:) = mean(zall_motion(BehavData.Block == qq, :));
    block_sem(qq,:) = nansem(zall_motion(BehavData.Block == qq, :));
    block_time2Collect_median(qq,:) = median(time2Collect(BehavData.Block == qq, :));
    block_trialStartTime_median(qq,:) = median(trialStartTime(BehavData.Block == qq, :));

end


%%
figure;
subplot('Position', [0.1, 0.45, 0.8, 0.50]);
imagesc(window_ts3, 1, zall_motion);hold on;  
time2Collect = BehavData.collectionTime(:) - BehavData.choiceTime(:);
trialStartTime = BehavData.stTime(:) - BehavData.choiceTime(:);
[numTrials, ~] = size(BehavData.collectionTime(:));
Tris = [1:numTrials]';
scatter(time2Collect, Tris               , 'Marker', 'p', 'MarkerFaceColor', 'w', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha',1)
scatter(trialStartTime, Tris, 'Marker', 's', 'MarkerFaceColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 1)
plot(zeros(numTrials, 1), Tris, 'LineWidth', 3, 'LineStyle', "-", 'Color', 'w')
yline(30.5,'-',{'Block 2', 'begins'},'LineWidth',3);
yline(60.5,'-',{'Block 3', 'begins'},'LineWidth',3);
colorbar;


subplot('Position', [0.1, 0.1, 0.69, 0.25]); %[left of figure, bottom, right of figure, top]
shadedErrorBar(window_ts2, block_mean(1,:), block_sem(1,:));
hold on;
shadedErrorBar(window_ts2, block_mean(2,:), block_sem(2,:),'lineProps', '-b');
hold on;
shadedErrorBar(window_ts2, block_mean(3,:), block_sem(3,:),'lineProps', '-r');

xline(0, 'Color', 'k', 'LineStyle','-', 'LineWidth', 2);
xline(block_trialStartTime_median, 'Color', 'b', 'LineStyle','-', 'LineWidth', 1);
xline(block_time2Collect_median, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 1);
legend('block 1 (0%)','block 2 (50%)', 'block 3 (75%)')