





ca = final.BLA_Insc_24.Pre_RDT_RM.CNMFe_data.C_raw  ;
BehavData = final.BLA_Insc_24.Pre_RDT_RM.uv.BehavData;
[BehavData,trials,varargin]=TrialFilter(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0);
trials = cell2mat(trials);

% shuffle data for comparison
[num_cells, num_samples] = size(ca);
shuffled_data = zeros(num_cells, num_samples); % Preallocate matrix for efficiency
shift_val = randi(num_samples); % Generate a random shift value for each signal RUAIRI RECOMMENDED KEEPING THE SAME SHIFT VAL, rather than randomizing per neuron. this is because then you keep the overall correlation b/w the neurons, but disrupt the relationship to the event timestamps
for i = 1:num_cells
    % shift_val = randi(num_signals); % Generate a random shift value for each signal
    shuffled_data(i,:) = circshift(ca(i,:), shift_val,2); % Perform the circular shuffle
end
ca = shuffled_data;

num_samples = size(ca, 2);
sampling_frequency = (final.BLA_Insc_24.Pre_RDT_RM.uv.dt)*100;
time_array = final.BLA_Insc_24.Pre_RDT_RM.time;
eTS = BehavData.choiceTime; %get time stamps

zb_session = mean(ca,2);
zsd_session = std(ca,[],2);
% caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace


%calculate time windows for each event
evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate



for u = 1:size(ca,1)
    neuron_num = neuron_num+1;
    % initialize trial matrices
    caTraceTrials = NaN(size(eTS,1),numMeasurements); %
    unitTrace = ca(u,:); %get trace
    for t = 1:size(eTS,1)
        unitTrace_zscored = zscore(unitTrace);
        timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
        BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
        if min(timeWin) > min(time_array) && max(timeWin) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
            % get unit event counts in trials
            % get unit ca traces in trials
            idx = time_array > min(timeWin) & time_array < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
            bl_idx = time_array > min(BL_win) & time_array < max(BL_win);
            %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
            caTraceTrials(t,1:sum(idx)) = unitTrace(idx);
            zscored_caTraceTrials(t, 1:sum(idx)) = unitTrace_zscored(idx);
            zb(t,:) = nanmean(unitTrace(bl_idx)); %baseline mean
            zb_window(t,:) = nanmean(caTraceTrials(t,:));
            zsd(t,:) = nanstd(unitTrace(bl_idx)); %baseline std
            zsd_window(t,:) = nanstd(caTraceTrials(t,:));
            tmp = 0;
            for j = 1:size(caTraceTrials,2)
                tmp = tmp+1;
                zall_baselined(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
                zall_window(t,tmp) = (caTraceTrials(t,j) - zb_window(t))/zsd_window(t);
                zall_session(t,tmp) = (caTraceTrials(t,j) - zb_session(u))/zsd_session(u);
                zall_luthi(t, tmp) = (caTraceTrials(t,j) - zb_window(t))/zb_session(u);
            end
            clear j;

        end
    end
    caTraceTrials = caTraceTrials(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
    zall = zall_window(:, 1:size(ts1, 2)); %added to make sure dimensions are the same as ts1
    zall= sgolayfilt(zall, uv.smoothing.params(1), uv.smoothing.params(2));
    zall_array_test{u} = zall;
    zall_mouse_test{ii, 1}(u) = {zall};
    caTraceTrials_mouse{ii, 1}(u) = {caTraceTrials};
    zall_mean(u,:) = mean(zall);
    trials_per_mouse{ii, 1} = trials;
    clear zall caTraceTrials zb zsd;
end

% Define figure width and height (adjust as needed)
fig_width = 200; % Width in pixels
fig_height = 300; % Height in pixels

% Create figure with custom size
figure('Position', [100, 100, fig_width, fig_height]); 
plot(ts1, zall_array{1, 10}(14, :))
hold on; plot(ts1, zall_array{1, 10}(24, :))
hold on; plot(ts1, zall_array{1, 10}(15, :))
hold on; plot(ts1, zall_array{1, 10}(16, :))
ylim([-2 2])
xlim([-8 8])


% Define figure width and height (adjust as needed)
fig_width = 200; % Width in pixels
fig_height = 300; % Height in pixels

% Create figure with custom size
figure('Position', [100, 100, fig_width, fig_height]); 
plot(ts1, zall_mouse_test{10, 1}{1, 10}  (14, :))
hold on; plot(ts1, zall_mouse_test{10, 1}{1, 10}  (24, :))
hold on; plot(ts1, zall_mouse_test{10, 1}{1, 10}  (15, :))
hold on; plot(ts1, zall_mouse_test{10, 1}{1, 10}  (16, :))
ylim([-2 2])
xlim([-8 8])