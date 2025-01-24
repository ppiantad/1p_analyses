function [zall_baselined, zall_window, zall_session, caTraceTrials, trial_ca, StartChoiceCollect_times, zscored_caTraceTrials, zall_luthi] = align_and_zscore(BehavData, unitTrace, eTS, uv, time_array, zb_session, zsd_session, u, use_normalized_time)

for t = 1:size(eTS,1)
    % set each trial's temporal boundaries
    % get trial data from 5s before start to 5s after collection
    trial_based_time{t} = BehavData.stTime(t)-5:uv.dt:BehavData.collectionTime(t)+5;
    StartChoiceCollect_times{t} = [BehavData.stTime(t) BehavData.choiceTime(t) BehavData.collectionTime(t)];
    trial_ca_idx = time_array > min(trial_based_time{t}) & time_array < max(trial_based_time{t});
    trial_ca{t} = unitTrace(trial_ca_idx);
    
    timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
    BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
    unitTrace_zscored = zscore(unitTrace);
    if use_normalized_time == 0
    
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
    elseif use_normalized_time == 1
        [normalized_trial_ca, concatenated_normalized_trial_ca] = normalize_trials_in_time_fn_single_trial(trial_ca);
        if min(timeWin) > min(time_array) && max(timeWin) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
            % get unit event counts in trials
            % get unit ca traces in trials
            idx = time_array > min(timeWin) & time_array < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
            bl_idx = time_array > min(BL_win) & time_array < max(BL_win);
            %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
            caTraceTrials(t,1:sum(idx)) = unitTrace(idx);
            zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
            zb_window(t,:) = mean(caTraceTrials(t,:));
            zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
            zsd_window(t,:) = std(caTraceTrials(t,:));
            tmp = 0;
            for j = 1:size(caTraceTrials,2)
                tmp = tmp+1;
                zall_baselined(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
                zall_window(t,tmp) = (caTraceTrials(t,j) - zb_window(t))/zsd_window(t);
                zall_session(t,tmp) = (caTraceTrials(t,j) - zb_session(u))/zsd_session(u);
            end
            clear j;

        end

    end
    
end