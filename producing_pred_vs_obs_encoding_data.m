 trial_start_pred_activity = X(:,2:26)*(encoding_output.coefs(1,2:26))';

 %%


%% Edit these uservariables with what you want to look at
%this script has been verified as of 7/30/2022 to be well-aligned based on
%shock-responsive neurons from multiple sessions for multiple mice
uv.evtWin = [-10 10]; %what time do you want to look at around each event [-10 10] [-5 35]
uv.BLper = [-10 -5]; %what baseline period do you want for z-score [-10 -5] [-5 0]
uv.dt = 0.1; %what is your frame rate (check neuron.Fs to be sure) 0.2 0.1
uv.behav = 'choiceTime'; %which behavior/timestamp to look at choiceTime stTime



 eTS = BehavData.stTime;
 
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
         caTraceTrials(t,1:sum(idx)) = y(idx);
         % zb(t,:) = mean(unitTrace(bl_idx)); %baseline mean
         % zsd(t,:) = std(unitTrace(bl_idx)); %baseline std
         zb(t,:) = mean(caTraceTrials(t,:)); %baseline mean
         zsd(t,:) = std(caTraceTrials(t,:)); %baseline std


         tmp = 0;
         for j = 1:size(caTraceTrials,2)
             tmp = tmp+1;
             zall_real(t,tmp) = (caTraceTrials(t,j) - zb(t))/zsd(t);
         end
         clear j;



     end
 end


 figure; plot(zall);
 ts1 = (-10:.1:10-0.1);
 figure; plot(ts1, zall(:, 1:200));
 figure; plot(ts1, mean(zall(:, 1:200)));
