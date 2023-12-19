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
    velocity_trace =  zscored_SLEAP_data_vel_filtered_session;
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

figure
imagesc(window_ts3, 1, velocity_trace_trials);hold on;   
figure;
plot(window_ts2, velocity_trace_trials);



% toc
% %     final(i).name = mouseData(i).mouseID;
% %     final(i).day = i;
% final(i).time = frames3; %final(i).time = caTime;
% final(i).unitAVG = unitAVG;
% final(i).unitXTrials = unitXTrials;
% final(i).uv = uv;
% final(i).unitSEM = unitSEM;
% final(i).uv.BehavData = BehavData;
% 
% 
% clear unitTS unitTrace unitXTrials unitAVG unitSEM i zall zb zsd



