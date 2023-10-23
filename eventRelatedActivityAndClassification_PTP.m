% load('BLA_Risk_Data_struct_01222023.mat')


session_to_analyze = 'RDT_D1';
epoc_to_align = 'choiceTime';
event_to_analyze = {'BLOCK',1,'REW',1.2};

window_sz = (0:.1:20);
ts1 = (-10:.1:10);

iter = 0;
%% %user selected variables
uv.chooseFluoresenceOrRate = 1;                                             %set to 1 to classify fluoresence response; set to 2 to classify firing rate responses
uv.sigma = 1;                                                               %this parameter controls the number of standard deviations that the response must exceed to be classified as a responder. try 1 as a starting value and increase or decrease as necessary.
uv.evtWin = [-10 10];                                                       %time window around each event in sec relative to event times (use long windows here to see more data)
uv.evtSigWin.outcome = [-4 0];                                       %period within time window that response is classified on (sec relative to event)
% uv.evtSigWin.outcome = [0 4]; 
% uv.evtSigWin.groomingStop = [-.5 3];
% uv.evtSigWin.faceGroomingStart = [-.5 2];
% uv.evtSigWin.faceGroomingStop = [-.5 2];
uv.dt = 0.1;                                                                %time step size (sec)
uv.resamples = 1000                                                         %number of resamples to use in shuffle analysis

sub_window_idx = ts1 >= uv.evtSigWin.outcome(1) & ts1 <= uv.evtSigWin.outcome(2);


identity_classification_win = 'Outcome';
identity_classification_str = join(string(uv.evtSigWin.outcome), 'to');

if contains(identity_classification_str, '-')
    identity_classification_str = strrep(identity_classification_str, '-', 'Minus_');
end



identity_classification_str = join([identity_classification_win, identity_classification_str],'_');

iter = iter +1;


evtWinSpan = max(uv.evtWin) - min(uv.evtWin);                               %calculate length of each period to examine neural activity in
numMeasurements = evtWinSpan/uv.dt;                                         %calculate the number of measurements in each time window

animalIDs = (fieldnames(final));
neuron_num = 0;

for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        [data,trials, varargin_identity_class] = TrialFilter(final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.BehavData, 'REW', 1.2, 'BLOCK', 3);
        varargin_strings = string(varargin_identity_class);
        varargin_strings = strrep(varargin_strings, '0.3', 'Small');
        varargin_strings = strrep(varargin_strings, '1.2', 'Large');
        filter_args = strjoin(varargin_strings,'_');

        trials = cell2mat(trials);
        
        for qq = 1:size(final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials,2)
            neuron_num = neuron_num+1; 
           
%           caTraceTrials = final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).caTraces(trials,:);
            %zall is currently the hardcoded z-score which is taken from
            %-10 to -5; should update this so that the zscore used is maybe
            %to the entire session?
            
            % IF Z-SCORE TO WINDOW, UNCOMMENT BELOW
%             caTraceTrials = final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).zall_window(trials,:);
            %IF Z-SCORE TO A BASELINE, UNCOMMENT BELOW
            caTraceTrials = final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).zall(trials,:);

            % for some events, the mice have no trials, therefore there are
            % no traces. this line basically skips those neurons (adding 0
            % to the respClass struct), to maintain the same total # of
            % cells throughout (for easy filtering if wanting to check the
            % mean of the respClass.activated = 1 neurons, for example
            if isempty(caTraceTrials)
                respClass.(identity_classification_str).(filter_args).activated(neuron_num,1) = 0;     
                respClass.(identity_classification_str).(filter_args).inhibited(neuron_num,1) = 0;
                respClass.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).activated(qq,1) = 0;     
                respClass.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).inhibited(qq,1) = 0;
                respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).neutral(qq,1) = 0;
            elseif ~isempty(caTraceTrials)

                trialCt = size(caTraceTrials,1);                                    %number of trials for currently analyzed event

                for e = 1:size(data,1)                                                      %for each event
                    evtWinIdx = ts1 >= uv.evtSigWin.outcome(1,1) &...          %calculate logical index for each event period
                        ts1 <= uv.evtSigWin.outcome(1,2);
                end
                clear evtWinSpan e
                for g = 1:uv.resamples                                              %for each resampling of the data
                    [~,shuffledIDX] = sort(randi...                                 %generate a matrix of random integers from 1 to the number of measurements in each time window (not each number is generated/some are repeated)
                        (numMeasurements,trialCt,numMeasurements),2);               %sort the data index to create a new list of indices
                    for t = 1:trialCt                                               %for each trial
                        shuffledTrace(t,:) = caTraceTrials(t,shuffledIDX(t,:));     %shuffle the calcium trace
                        %         shuffledEvtRate(t,:) = caEvtRateTrials(t,shuffledIDX(t,:)); %shuffle the event rate
                    end
                    nullDistTrace(g,:) = nanmean(shuffledTrace);                    %calculate the NaN mean of the shuffled traces
                    %     nullDistEvtRate(g,:) = nanmean(shuffledEvtRate);                %calculate the NaN mean of the shuffled event rates
                end
                clear shuffled* g t trialCt

                %% choose to classify fluoresence or event rates
                if uv.chooseFluoresenceOrRate == 1                                  %if user selected to classify the fluoresence
                    nullDist = nullDistTrace;                                       %direct transfer
                    empiricalTrialWin = nanmean...
                        (caTraceTrials(:,evtWinIdx));                   %NaN mean of the fluorescent response across trials, within the time window. this gets the within trial mean.
                    empiricalWinAvg = nanmean(empiricalTrialWin);                   %across trial mean
                elseif uv.chooseFluoresenceOrRate == 2                              %if user selected to classify the event rates
                    nullDist = nullDistEvtRate;                                     %direct transfer
                    empiricalTrialWin = nanmean...
                        (caEvtRateTrials(:,evtWinIdx.(evts{e})),2);                 %within trial average
                    empiricalWinAvg = nanmean(empiricalTrialWin);                   %across trial mean
                end
                clear empiricalTrialWin
                %%
                sdNull = nanstd(nullDist(:));                                          %calculate the standard deviation of the null distribution
                upperSD = nanmean(nullDist(:)) + (uv.sigma*sdNull);                    %calculate upper limit of interval around the mean
                lowerSD = nanmean(nullDist(:)) - (uv.sigma*sdNull);                    %calculate lower limit of interval around the mean
                clear sdNull nullDist
                %%
                respClass.(session_to_analyze).(identity_classification_str).(filter_args).activated(neuron_num,1) = empiricalWinAvg > upperSD;     %classify as activated if empirical response exceeds upper limit
                respClass.(session_to_analyze).(identity_classification_str).(filter_args).inhibited(neuron_num,1) = empiricalWinAvg < lowerSD;     %classify as inhibited if empirical response exceeds lower limit
                respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).activated(qq,1) = empiricalWinAvg > upperSD;     %classify as activated if empirical response exceeds upper limit
                respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).inhibited(qq,1) = empiricalWinAvg < lowerSD;     %classify as inhibited if empirical response exceeds lower limit
                respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).neutral(qq,1) = respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).activated(qq,1) == 0 & respClass_mouse.(currentanimal).(session_to_analyze).(epoc_to_align).(identity_classification_str).(filter_args).inhibited(qq,1) == 0;
                if empiricalWinAvg > upperSD
                    respClass_all(iter, neuron_num) = 1;
                elseif empiricalWinAvg < lowerSD
                    respClass_all(iter, neuron_num) = 2;
                else 
                    respClass_all(iter, neuron_num) = 3;
                end
                clear upperSD lowerSD empiricalWinAvg
                %             %% store trial by trial data
                %             unitXTrials(u).(evts{e}).caEvtCts = caEvtCtTrials;                  %store evoked calcium event counts over all trials
                %             unitXTrials(u).(evts{e}).caEvtRate = caEvtRateTrials;               %store evoked calcium event rates over all trials
                %             unitXTrials(u).(evts{e}).caTraces = caTraceTrials;                  %store evoked calcium traces over all trials
                %             %% store unit averaged data
                %             unitAVG.(evts{e}).caEvtCts(u,:) = nanmean(caEvtCtTrials);           %store trial averaged event counts
                %             unitAVG.(evts{e}).caEvtRates(u,:) = nanmean(caEvtRateTrials);       %store trial averaged event rates
                %             unitAVG.(evts{e}).caTraces(u,:) = nanmean(caTraceTrials);           %store trial averaged calcium traces
                %             clear caEvtCtTrials caTraceTrials caEvtRateTrials
            end
        end
    end


end
sum_activated(iter) = sum(respClass_all(iter,:) == 1);
sum_inhibited(iter) = sum(respClass_all(iter,:) == 2);
sum_neutral(iter) = sum(respClass_all(iter,:) == 3);
sum_inhibited_percent(iter) = (sum_inhibited(iter)/neuron_num)*100;
sum_activated_percent(iter) = (sum_activated(iter)/neuron_num)*100;
sum_neutral_percent(iter) = (sum_neutral(iter)/neuron_num)*100;
labels = {'activated', 'inhibited', 'neutral'};
figure; pie([sum_activated(iter) sum_inhibited(iter) sum_neutral(iter)])


%% plot activated neurons

load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
colormap(batlowW); % c1 = colorbar; 


figure; plot(ts1, neuron_mean(respClass_all(iter,:) == 1,:))
hold on;
figure;
% figure; plot(ts1, mean(neuron_mean(respClass_all(iter,:) == 1,:)))
shadedErrorBar(ts1, mean(neuron_mean(respClass_all(iter,:) == 1,:)), mean(neuron_sem(respClass_all(iter,:) == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
hold off; 
% hold on;
% figure; plot(ts1, mean(neuron_mean(respClass_all(iter,:) == 1,:)))
% errorplot3(mean(neuron_mean(respClass_all(iter,:) == 1,:))-mean(neuron_sem(respClass_all(iter,:) == 1,:)),mean(neuron_mean(respClass_all(iter,:) == 1,:))+mean(neuron_sem(respClass_all(iter,:) == 1,:)),[-10 10],batlowW(iter,:),.15);
% hold off;
% plot inhibited neurons
figure; plot(ts1, neuron_mean(respClass_all(iter,:) == 2,:))
figure; shadedErrorBar(ts1, mean(neuron_mean(respClass_all(iter,:) == 2,:)), mean(neuron_sem(respClass_all(iter,:) == 2,:)), 'lineProps', {'color', batlowW(iter,:)});
hold off; 
% plot neutral neurons
figure; plot(ts1, neuron_mean(respClass_all(iter,:) == 3,:))
figure; shadedErrorBar(ts1, mean(neuron_mean(respClass_all(iter,:) == 3,:)), mean(neuron_sem(respClass_all(iter,:) == 3,:)), 'lineProps', {'color', batlowW(iter,:)});



%%
co_activated = respClass_all(1,:) == 1 & respClass_all(2,:) == 1;
co_inhibited = respClass_all(1,:) == 2 & respClass_all(2,:) == 2;
co_activated_sum = sum(co_activated);
co_inhibited_sum = sum(co_inhibited);
exclusive_activated = respClass_all(1,:) ~= 1 & respClass_all(2,:) == 1;
exclusive_inhibited = respClass_all(1,:) ~= 2 & respClass_all(2,:) == 2;
exclusive_activated_sum = sum(exclusive_activated);
exclusive_inhibited_sum = sum(exclusive_inhibited);


% this is the start of checking if neurons are MORE active than during
% other events, i.e. if you wanted to check if REW activated neurons are
% significantly differentially activated by SHK. preliminary poking around
% seems to suggest that few large reward active neurons have their activity
% increase in response to SHK
co_activated_indices = find(co_activated(1,:) == 1);
[h,p] = ttest(neuron_mean_LARGE(co_activated_indices(10),sub_window_idx),neuron_mean(co_activated_indices(10),sub_window_idx));


total_modulated = [sum_activated + sum_inhibited];
total_co_modulated = [co_activated_sum + co_inhibited_sum];

A = total_modulated;
I = total_co_modulated;
figure; venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black')


A = sum_inhibited;
I = co_inhibited_sum;
figure; venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black')

A = sum_activated;
I = co_activated_sum;
figure; venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black')



%%
figure;
hold on; 
shadedErrorBar(ts1, mean(neuron_mean(respClass_all(iter,:) == 1,:)), mean(neuron_sem(respClass_all(iter,:) == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean(respClass_all(iter,:) == 2,:)), mean(neuron_sem(respClass_all(iter,:) == 2,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean(respClass_all(iter,:) == 3,:)), mean(neuron_sem(respClass_all(iter,:) == 3,:)), 'lineProps', {'color', batlowW(iter,:)});
