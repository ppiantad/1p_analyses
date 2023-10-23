%% Edit these uservariables with what you want to look at
uv.evtWin = [-10 10]; %what time do you want to look at around each event
uv.dt = 0.1; %what is your frame rate
uv.behav = 'choiceTime'; %which behavior/timestamp to look at

[BehavData,ABETfile,Descriptives, block_end]=ABET2TableFn_Chamber_A_v5('BLA-Insc-1 01112021.csv',[]);

% [BehavData,ABETfile]=ABET2TableFn_ShockTest('BLA-Insc-3 01272021.csv');


ABET_removeheader = ABETfile(2:end,:);

tbl_ABET = cell2table(ABET_removeheader);
tbl_ABET.Properties.VariableNames = ABETfile(1,:);

% gpio_tbl = readtable('BLA-Insc-1_2021-01-20_RDT_D1_GPIO.csv');

% shk_times = tbl_ABET.Evnt_Time(strcmp(tbl_ABET.Item_Name, 'shock_on_off') & tbl_ABET.Arg1_Value == 1);

% 
% stTime = gpio_tbl.Time_s_(strcmp(gpio_tbl.ChannelName, 'GPIO-2') & gpio_tbl.Time_s_ > 0);
% 
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

% BehavData.TrialPossible(:)=BehavData.TrialPossible(:)+stTime(1);
% BehavData.choiceTime(:)=BehavData.choiceTime(:)+stTime(1);
% BehavData.collectionTime(:)=BehavData.collectionTime(:)+stTime(1);
% BehavData.stTime(:)=BehavData.stTime(:)+stTime(1);



% shk_times(:)=shk_times(:)+stTime(1);

% BehavData.choTime2 = BehavData.choiceTime-BehavData.TrialPossible;
% BehavData.choTime3 = BehavData.Insc_TTL+BehavData.choTime2;


BehavData=TrialFilter(BehavData,'REW',1.2);



load('BLA-Insc-1_2021-01-11_RM_D7.CNMF_final.mat');
ts1 = uv.dt:uv.dt:length(neuron.C_raw)*uv.dt;
% BehavData.choiceTime = BehavData.choiceTime-1;

%%
for i = 1 %could loop through multiple mice like this if you had it
    eTS = BehavData.(uv.behav); %get time stamps
    ca = neuron.C_raw; %get calcium
    caTime = uv.dt:uv.dt:length(ca)*uv.dt; %generate time trace
    %calculate time windows for each event
    evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
    numMeasurements = round(evtWinSpan/uv.dt); %need to round due to odd frame rate
    %%
    tic
    for u = 1:size(ca,1)
        %% initialize trial matrices
        caTraceTrials = NaN(size(eTS,1),numMeasurements);
        unitTrace = ca(u,:); %get trace
        %             %%
        for t = 1:size(eTS,1)
            %% set each trial's temporal boundaries
            timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
            if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range
                %% get unit event counts in trials
                %% get unit ca traces in trials
                idx = caTime > min(timeWin) & caTime < max(timeWin);      %logical index of time window around each behavioral event time
                %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
                caTraceTrials(t,1:sum(idx)) = unitTrace(idx);
                
            end
        end
        clear idx timeWin
        %% 
        unitXTrials(u).caTraces = caTraceTrials;
        
        %% store unit averaged data
        unitAVG.caTraces(u,:) = nanmean(caTraceTrials);           %store trial averaged calcium traces
        
        clear caEvtCtTrials caTraceTrials caEvtRateTrials unitTrace idx
    end
end

toc
%     final(i).name = mouseData(i).mouseID;
%     final(i).day = i;
final(i).time = caTime;
final(i).unitAVG = unitAVG;
final(i).unitXTrials = unitXTrials;
final(i).uv = uv;
clear unitTS unitTrace unitXTrials unitAVG

%%

% BASELINE_PER = [-10 -5];
% tmp=0;
% 
% for i = 1:size(final.unitXTrials,2)
% %     BL_shifted(pp,:)=[BASELINE_PER(1)+time2EPOC(i) BASELINE_PER(2)+time2EPOC(i)]; %BL_shifted(pp,:)=[BASELINE_PER(1)+(-1*time2Collect(i)) BASELINE_PER(2)+time2Collect(i)];
% %     ind = ts2(1,:) < BL_shifted(pp,2) & ts2(1,:) > BL_shifted(pp,1);
%     ind = window_ts2(1,:) < BASELINE_PER(2) & window_ts2(1,:) >= BASELINE_PER(1);
%     for kk = 1:size(final.unitXTrials(i).caTraces,1)
%         zb(kk,1) = mean(final.unitXTrials(i).caTraces(kk,ind)); % baseline period mean
%         zsd(kk,1) = std(final.unitXTrials(i).caTraces(kk,ind)); % baseline period stdev
%         
%         for j = 1:size(window_ts2,2)
%             tmp=tmp+1;
%             zall(kk,tmp) = (final.unitXTrials(i).caTraces(kk,j)) - zb(kk)/zsd(kk);
%         
%         end
%         tmp=0;
%         meanBaseline(i).caTraces = zb;
%         stdBaseline(i).caTraces = zsd;
%         zscoredALL(i).caTraces = zall;
%         
% 
% 
% 
%     end
% 
% end
% 
%     for j = 1:size(meanBaseline,2) % Z score per bin
%         tmp = tmp + 1;
%         zall(kk,1)=(Y_dF_all(i,j) - zb)/zsd;
%         dfALL(i,tmp)=(Y_dF_all(i,j) - zbmedian)/zbmedian;
%         BL_mean(i) = mean(Y_dF_all(i,:));
%     end
%     tmp=0;
% end
% end

%%
time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);

[numTrials,~]=size(BehavData.collectionTime(:));
Tris=[1:numTrials]';

% window_ts = uv.dt:uv.dt:numMeasurements*uv.dt;

window_ts2 = uv.evtWin(1):uv.dt:uv.evtWin(2)-uv.dt;


figure
imagesc(window_ts2, 1, final.unitXTrials(14).caTraces);hold on;
% scatter(time2Collect,Tris,'Marker','p','MarkerFaceColor','w')
plot(zeros(numTrials,1),Tris)
% xticklabels = (final.uv.evtWin(1,1)):5:(final.uv.evtWin(1,2));
% xticks = linspace(1, length(final.unitAVG.caTraces), numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% set(gca,'FontName','Arial','FontSize',16)
xlabel('Time from Large Rew Choice (s)','FontSize',22)
ylabel('Trial number','FontSize',22)

%%
for c = 1:size(ca,2)
    subplot(121)
    imagesc(final.unitXTrials(c).caTraces)
    xticklabels = (final.uv.evtWin(1,1)):5:(final.uv.evtWin(1,2));
    xticks = linspace(1, length(final.unitAVG.caTraces), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    set(gca,'FontName','Arial','FontSize',16)
    xlabel('Time from(s)','FontSize',22)
    ylabel('Trial number','FontSize',22)
    subplot(122)
    trialTime = uv.evtWin(1,1):uv.dt:uv.evtWin(1,2)-uv.dt;
    plot(trialTime,final.unitAVG.caTraces(c,:))
    set(gca,'FontName','Arial','FontSize',16)
    xlabel('Time from  (s)','FontSize',22)
    xlim([uv.evtWin(1,1) uv.evtWin(1,2)])
    xline(0,'--')
    pause
end
%% Generate a continuous behavior trace of
behavTrace = zeros(1,length(caTime));
for b = 1:length(eTS)
    startTS(b,:) = eTS(b,1);
    [idxPoint idxPoint] = min(abs(caTime-startTS(b))); %grab the index for when the current grooming bout starts
    behavTrace(idxPoint) = 1;
    clear idxPoint
end
