%% Run me first


iter = 0;


load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)




%% Load the session you want to examine


% load('BLA-NAcShell_Risk_2023_09_15.mat')

% load('BLA_panneuronal_Risk_2023_07_06.mat')

load('BLA_panneuronal_Risk_2023_11_15.mat')

% load('NAcSh_D2_Cre-OFF_GCAMP_all.mat')

% load('BLA_panneuronal_Risk_matched_RM_D1_vs_Pre_RDT_RM.mat')

% load('BLA_panneuronal_Risk_matched_PreRDTRM_RDT_D1.mat')

% load('BLA_panneuronal_Risk_matched_RDT_D1_vs_RDT_D2.mat')

% load('BLA_panneuronal_Risk_matched_RDT_D2_vs_RDT_D3.mat')

% load('BLA_panneuronal_Risk_matched_RDT_D1_vs_SHOCK_TEST.mat')

% load('BLA_NAcSh_Risk_matched_RM_D1_vs_Pre_RDT_RM.mat');

% load('BLA_NAcSh_Risk_matched_Pre_RDT_RM_vs_RDT_D1.mat')

% load('BLA_NAcSh_Risk_matched_RDT_D1_vs_RDT_D2.mat')

%%

session_to_analyze = 'Pre_RDT_RM';
epoc_to_align = 'choiceTime';
event_to_analyze = {'BLOCK',1,'REW',1.2};

window_sz = (0:.1:20-0.1);
ts1 = (-10:.1:10-0.1);


%% %user selected variables
clear neuron_mean neuron_sem neuron_num trials
uv.chooseFluoresenceOrRate = 1;                                             %set to 1 to classify fluoresence response; set to 2 to classify firing rate responses
uv.sigma = 2;                                                               %this parameter controls the number of standard deviations that the response must exceed to be classified as a responder. try 1 as a starting value and increase or decrease as necessary.
uv.evtWin = [-10 10];                                                       %time window around each event in sec relative to event times (use long windows here to see more data)
% uv.evtSigWin.outcome = [-3 0]; %for trial start
uv.evtSigWin.outcome = [-4 0]; %for pre-choice                                     %period within time window that response is classified on (sec relative to event)
% uv.evtSigWin.outcome = [1 3]; %for REW

% uv.evtSigWin.outcome = [0 1]; %for SHK
% uv.evtSigWin.groomingStop = [-.5 3];
% uv.evtSigWin.faceGroomingStart = [-.5 2];
% uv.evtSigWin.faceGroomingStop = [-.5 2];
uv.dt = 0.1;                                                                %time step size (sec)
uv.resamples = 100                                                         %number of resamples to use in shuffle analysis 1000

sub_window_idx = ts1 >= uv.evtSigWin.outcome(1) & ts1 <= uv.evtSigWin.outcome(2);
% neuron_sem = zeros(1, size(ts1, 2));
% neuron_mean = [];

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
        [data,trials, varargin_identity_class] = TrialFilter(final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.BehavData, 'REW', 1.2);
        
        if ~strcmp('stTime',data.Properties.VariableNames)
            data.stTime = data.TrialPossible - 5;
        end
        if ~strcmp('collectionTime',data.Properties.VariableNames)
            data.collectionTime = data.choiceTime + 5;
        end
        behav_tbl_temp{ii,:} = data;
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
            caTraceTrials = final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).caTraces(trials,:);
            

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
                neuron_mean(neuron_num,:) = nan;
                neuron_sem(neuron_num,:) = nan; 
            elseif ~isempty(caTraceTrials)
                for h = 1:size(caTraceTrials,1)
                    zb(h) = mean(caTraceTrials(h,:)); %baseline mean
                    zsd(h) = std(caTraceTrials(h,:)); %baseline std
                    tmp = 0;
                    for j = 1:size(caTraceTrials(1:length(window_sz)),2)
                        tmp = tmp+1;
                        zall(h,tmp) = (caTraceTrials(h,tmp) - zb(h))/zsd(h);
                        
                    end
                    % zall(h,:) = sgolayfilt(zall(h,:), 9, 21);
%                     zall_to_BL_array(iter, neuron_num) = {final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).zall(trials,:)};
                end

                % Loop through each row of zall
                for z = 1:size(zall, 1)
                    % Apply Savitzky-Golay filter to each row
                    zall(z, :) = sgolayfilt(zall(z, :), 9, 21);
                end
                mouse_cells(iter, neuron_num) = {currentanimal};
                zall_array(iter, neuron_num) = {zall};
                % neuron_mean(neuron_num,:) = sgolayfilt((mean(zall,1)), 9, 21);
                neuron_mean(neuron_num,:) = (mean(zall,1));
                if size(zall, 1) == 1
                    neuron_sem(neuron_num,:) = zeros(1, size(neuron_sem, 2));
                else
                    neuron_sem(neuron_num,:) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                end
                zsd_array(iter, neuron_num) = {zsd};
                zall_mouse{ii, iter+1}(qq) = {zall};
                caTraceTrials_mouse{ii, iter+1}(qq) = {caTraceTrials};
                neuron_mean_mouse{ii, iter+1}(qq,: ) = mean(zall, 1);
                neuron_sem_mouse{ii, iter+1}(qq,: ) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                caTraceTrials = zall;
                clear zall zb zsd;
                
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
                    respClass_all(neuron_num) = 1;
                elseif empiricalWinAvg < lowerSD
                    respClass_all(neuron_num) = 2;
                else 
                    respClass_all(neuron_num) = 3;
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
            clear caTraceTrials;
        end
    end


end
varargin_list{iter,:} = varargin_identity_class;
behav_tbl_iter{iter, :} = behav_tbl_temp;
clear behav_tbl_temp

sum_activated(iter) = sum(respClass_all == 1);
sum_inhibited(iter) = sum(respClass_all == 2);
sum_neutral(iter) = sum(respClass_all == 3);
sum_inhibited_percent(iter) = (sum_inhibited(iter)/neuron_num)*100;
sum_activated_percent(iter) = (sum_activated(iter)/neuron_num)*100;
sum_neutral_percent(iter) = (sum_neutral(iter)/neuron_num)*100;
labels = {'activated', 'inhibited', 'neutral'};
figure; pie([sum_activated(iter) sum_inhibited(iter) sum_neutral(iter)])
neuron_mean_array(iter) = {neuron_mean};
neuron_sem_array(iter) = {neuron_sem};
respClass_all_array(:,iter) = {respClass_all}';

clear respClass_all 


%% plot activated neurons

figure; plot(ts1, neuron_mean(respClass_all_array{:,iter} == 1,:))
hold on;
figure;
% figure; plot(ts1, mean(neuron_mean(respClass_all(iter,:) == 1,:)))
shadedErrorBar(ts1, nanmean(neuron_mean(respClass_all_array{:,iter} == 1,:)), nanmean(neuron_sem(respClass_all_array{:,iter} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
hold off; 
% hold on;
% figure; plot(ts1, mean(neuron_mean(respClass_all(iter,:) == 1,:)))
% errorplot3(mean(neuron_mean(respClass_all(iter,:) == 1,:))-mean(neuron_sem(respClass_all(iter,:) == 1,:)),mean(neuron_mean(respClass_all(iter,:) == 1,:))+mean(neuron_sem(respClass_all(iter,:) == 1,:)),[-10 10],batlowW(iter,:),.15);
% hold off;
% plot inhibited neurons
figure; plot(ts1, neuron_mean(respClass_all_array{:,iter} == 2,:))
figure; shadedErrorBar(ts1, nanmean(neuron_mean(respClass_all_array{:,iter} == 2,:)), nanmean(neuron_sem(respClass_all_array{:,iter} == 2,:)), 'lineProps', {'color', batlowW(iter,:)});
hold off; 
% plot neutral neurons
figure; plot(ts1, neuron_mean(respClass_all_array{:,iter} == 3,:))
figure; shadedErrorBar(ts1, nanmean(neuron_mean(respClass_all_array{:,iter} == 3,:)), nanmean(neuron_sem(respClass_all_array{:,iter} == 3,:)), 'lineProps', {'color', batlowW(iter,:)});


%% Use this code to plot heatmaps for each individual cell, across trials for all levels of iter
% **most useful for plotting matched cells within the same experiment, e.g., pan-neuronal matched Pre-RDT RM vs. RDT D1**

for ii = 1:size(zall_array, 2)
    figure;
    % Initialize variables to store global max and min for heatmap and line graph
    globalMaxHeatmap = -inf;
    globalMinHeatmap = inf;
    globalMaxYLine = -inf;
    globalMinYLine = inf;
    for qq = 1:iter
        
        % Create subplot with 1 row and 2 columns
        num_columns_plot = iter;
        subplot(2, num_columns_plot, qq);

        % Generate the heatmap
        imagesc(ts1, 1, zall_array{qq, ii});
        hold on;
        title({"Cell from " + strrep(mouse_cells{qq, ii}, '_', '-'), "Classified as " + respClass_all_array{qq}(ii), "(Overall cell number = " + (ii) + ")"}, 'FontSize', 9)
        
        % Update global max and min for heatmap
        localMaxHeatmap = max(zall_array{qq, ii}(:));
        localMinHeatmap = min(zall_array{qq, ii}(:));
        globalMaxHeatmap = max(globalMaxHeatmap, localMaxHeatmap);
        globalMinHeatmap = min(globalMinHeatmap, localMinHeatmap);
  
        % Find the row index in animalIDs that matches mouse_cells{qq, ii}
        isMatch = find(ismember(animalIDs,       mouse_cells{qq, ii}));

        if ~isempty(isMatch)  
            % Access beha     v_tbl_iter using the first index (assuming there's only one match)
            behav_tbl = behav_tbl_iter{qq, 1}{isMatch, 1};
            time2Collect = behav_tbl.collectionTime(:) - behav_tbl.choiceTime(:);
            trialStartTime = behav_tbl.stTime(:) - behav_tbl.choiceTime(:);
            [numTrials, ~] = size(behav_tbl.collectionTime(:));
            Tris = [1:numTrials]';
            scatter(time2Collect, Tris               , 'Marker', 'p', 'MarkerFaceColor', 'w')
            scatter(trialStartTime, Tris, 'Marker', 's', 'MarkerFaceColor', 'k')
            plot(zeros(numTrials, 1), Tris, 'LineWidth', 3, 'LineStyle', "--", 'Color', 'w')
        end     
        colorbar;
        
        hold off;
        clear time2Collect trialStartTime numTrials Tris behav_tbl

        % Create subplot for the mean and raw data
        subplot(2, num_columns_plot, iter + qq);

        % Plot the mean as a thick black line
        meanData = mean(zall_array{qq, ii});
        plot(ts1, meanData, 'k', 'LineWidth', 2);
        hold on;

        % Plot the raw data in grey with transparency
        for trial = 1:size(zall_array{qq, ii}, 1)
            plot(ts1, zall_array{qq, ii}(trial, :), 'Color', [0.1, 0.1, 0.1, 0.1]);
            hold on;
        end

        title("Mean and Raw Data", 'FontSize', 9)
        % Update global max and min for line graph
        localMaxYLine = max(zall_array{qq, ii}, [], "all");
        localMinYLine = min(zall_array{qq, ii}, [], "all");
        globalMaxYLine = max(globalMaxYLine, localMaxYLine);
        globalMinYLine = min(globalMinYLine, localMinYLine);
        
        hold off;        
    end
    % Set the same colorbar scale for all heatmap subplots
    for qq = 1:iter
        subplot(2, num_columns_plot, qq);
        clim([globalMinHeatmap, globalMaxHeatmap]);
        colorbar;
    end

    % Set the same Y-axis scale for all line graph subplots
    for qq = 1:iter
        subplot(2, num_columns_plot, iter + qq);
        ylim([globalMinYLine, globalMaxYLine]);
    end

    pause
    hold
    close
end

%%
excited_to_excited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1;
inhibited_to_inhibited = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 2;
excited_to_excited_sum = sum(excited_to_excited);
indices_excited_to_excited = find(excited_to_excited == 1);
inhibited_to_inhibited_sum = sum(inhibited_to_inhibited);

excited_to_inhibited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 2;
excited_to_inhibited_sum = sum(excited_to_inhibited);
excited_to_neutral = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3;
excited_to_neutral_sum = sum(excited_to_neutral);

inhibited_to_excited = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 1;
inhibited_to_excited_sum = sum(inhibited_to_excited);
inhibited_to_neutral = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 3;
inhibited_to_neutral_sum = sum(inhibited_to_neutral);

% above this line are good for single session stacked plots (e.g., how do
% neurons on session 1 that are activated (or inhibited) change to session
% 2?) 

neutral_to_excited = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1;
neutral_to_excited_sum = sum(neutral_to_excited);

neutral_to_inhibited = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 2;
neutral_to_inhibited_sum = sum(neutral_to_inhibited);

% need to add excited_to_excited twice because they are present in BOTH
% sessions
total_activated_possible = excited_to_excited_sum+ excited_to_excited_sum + excited_to_inhibited_sum+ excited_to_neutral_sum+ neutral_to_excited_sum + inhibited_to_excited_sum;
% need to add inhibited_to_inhibited twice because they are present in BOTH
% sessions
total_inhibited_possible = inhibited_to_inhibited_sum + inhibited_to_inhibited_sum + inhibited_to_excited_sum + inhibited_to_neutral_sum + neutral_to_inhibited_sum + excited_to_inhibited_sum;






neutral_to_neutral = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3;
neutral_to_neutral_sum = sum(neutral_to_neutral);

exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} ~= 1;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1);
exclusive_activated_session_2 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} == 1;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2);



exclusive_inhibited = respClass_all_array{1,1} ~= 2 & respClass_all_array{1,2} == 2;
% exclusive_activated_sum = sum(exclusive_activated);
% exclusive_inhibited_sum = sum(exclusive_inhibited);





% exclusive_modulated = exclusive_activated_sum + exclusive_inhibited_sum;


%%
%CREATE A STACKED BAR PLOT TO SHOW THE PROPORTION OF MODULATED (CROSS
%SESSION) NEURONS THAT CO OR EXCLUSIVELY MODULATED

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the total height of the stacked bar
total_height = sum(total_activated_possible)-excited_to_excited_sum;
% total_height = sum(total_inhibited_possible)-inhibited_to_inhibited_sum;

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Create a matrix for the bar plot data
stacked_plot_data = [excited_to_excited_sum, excited_to_inhibited_sum, excited_to_neutral_sum, neutral_to_excited_sum, inhibited_to_excited_sum];
% stacked_plot_data = [inhibited_to_inhibited_sum, inhibited_to_excited_sum, inhibited_to_neutral_sum, neutral_to_inhibited_sum, excited_to_inhibited_sum];


% Calculate the remaining portion
% remaining_portion = total_height - sum(stacked_plot_data);
figure;
% Create a vector for the x-axis values
x = 1;

% Create a stacked bar plot
bar(x, [stacked_plot_data], 'stacked');

% Set the colors for each section
colormap([0.3, 0.6, 0.9; 0.5, 0.5, 0.5; 0.2, 0.2, 0.2]); % Blue, Gray, and Dark Gray

% % Set the x-axis tick labels
% xticklabels({'Activation Types'});

% Set the y-axis limits to match the total height
ylim([0, total_height]);

% Add labels and a legend
% ylabel('Total Activation');
legend('Co-Excited', 'Excited to Inhibited', 'Excited to Neutral', 'Neutral to Excited', 'Inhibited to Excited');
% legend('Co-Inhibited', 'Inhibited to Excited', 'Inhibited to Neutral', 'Neutral to Inhibited', 'Excited to Inhibited');

%%
%CREATE A STACKED BAR PLOT TO SHOW THE PROPORTION OF MODULATED neurons
%within a single session

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the total height of the stacked bar
total_height = sum(sum_activated(1));
% total_height = sum(sum_inhibited);

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Create a matrix for the bar plot data
stacked_plot_data = [excited_to_excited_sum, excited_to_inhibited_sum, excited_to_neutral_sum];
% stacked_plot_data = [co_inhibited_sum, exclusive_inhibited_sum, inhibited_to_excited_sum, inhibited_to_neutral_sum];


% Calculate the remaining portion
% remaining_portion = total_height - sum(stacked_plot_data);
figure;
% Create a vector for the x-axis values
x = 1;

% Create a stacked bar plot
bar(x, [stacked_plot_data], 'stacked');

% Set the colors for each section
colormap([0.3, 0.6, 0.9; 0.5, 0.5, 0.5; 0.2, 0.2, 0.2]); % Blue, Gray, and Dark Gray

% % Set the x-axis tick labels
% xticklabels({'Activation Types'});

% Set the y-axis limits to match the total height
ylim([0, total_height]);

legend_string = strcat('Co-Excited =', num2str(excited_to_excited_sum));

% Add labels and a legend
% ylabel('Total Activation');
legend(strcat('Co-Excited = ', num2str(excited_to_excited_sum)), strcat('Excited to Inhibited = ', num2str(excited_to_inhibited_sum)), strcat('Excited to Neutral = ', num2str(excited_to_neutral_sum)));
% legend('Co-Inhibited', 'Exclusive Inhibited', 'Inhib to Excite', 'Inhib to Neutral');

%%
%CREATE A STACKED BAR PLOT TO SHOW THE PROPORTION OF MODULATED (CROSS
%SESSION) NEURONS THAT CO OR EXCLUSIVELY MODULATED

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the total height (total number of neurons)
% total_height = sum(sum_activated);
total_height = sum(sum_inhibited);


% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the percentages for each section
% percentage_co_activated = (co_activated_sum / total_height) * 100;
% percentage_exclusive_activated = (exclusive_activated_sum / total_height) * 100;
% percentage_remaining = 100 - (percentage_co_activated + percentage_exclusive_activated);

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Calculate the percentages for each section
percentage_co_inhibited = (co_inhibited_sum / total_height) * 100;
percentage_exclusive_inhibited = (exclusive_inhibited_sum / total_height) * 100;
percentage_remaining = 100 - (percentage_co_inhibited + percentage_exclusive_inhibited);

% UN-COMMENT THE DATA THAT YOU WANT TO PLOT
% Create a matrix for the bar plot data
% stacked_plot_data = [percentage_co_activated, percentage_exclusive_activated, percentage_remaining];
stacked_plot_data = [percentage_co_inhibited, percentage_exclusive_inhibited, percentage_remaining];

figure;
% Create a vector for the x-axis values
x = 1;

% Create a stacked bar plot
bar(x, stacked_plot_data, 'stacked');

% Set the colors for each section
colormap([0.3, 0.6, 0.9; 0.5, 0.5, 0.5; 0.2, 0.2, 0.2]); % Blue, Gray, and Dark Gray

% Set the x-axis tick labels
xticklabels({'Activation Types'});

% Set the y-axis limits to be percentages (0 to 100)
ylim([0, 100]);


legend('Co', 'Exclusive', 'Remaining');



%%
Block_1_activated = respClass_all_array{1,1} == 1;
Block_1_activated_sum = sum(Block_1_activated);
Block_2_consistent_activated = Block_1_activated & respClass_all_array{1,2} == 1;
disp(sum(Block_2_consistent_activated));
Block_3_consistent_activated = Block_1_activated & respClass_all_array{1,3} == 1;
disp(sum(Block_3_consistent_activated));

Block_1_responsive = respClass_all_array{1,1} == 1 | respClass_all_array{1,2} == 2;
disp(sum(Block_1_responsive));
Block_2_consistent_responsive = (respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1) | (respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 2);
disp(sum(Block_2_consistent_responsive));
Block_3_consistent_responsive = (respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1) & respClass_all_array{1,3} == 1 | (respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 2) & respClass_all_array{1,3} == 2;
disp(sum(Block_3_consistent_responsive));


Block_2_large_to_small = (respClass_all_array{1,2} ~= Block_2_consistent_responsive) & (respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 1) | (respClass_all_array{1,2} == 2 & respClass_all_array{1,3} == 2);
disp(sum(Block_2_large_to_small));
Block_3_large_to_small = (respClass_all_array{1,2} ~= Block_2_consistent_responsive) & (respClass_all_array{1,1} == 1 & respClass_all_array{1,4} == 1) | (respClass_all_array{1,2} == 2 & respClass_all_array{1,4} == 2);
disp(sum(Block_3_large_to_small));

% this is the start of checking if neurons are MORE active than during
% other events, i.e. if you wanted to check if REW activated neurons are
% significantly differentially activated by SHK. preliminary poking around
% seems to suggest that few large reward active neurons have their activity
% increase in response to SHK
co_activated_indices = find(excited_to_excited(1,:) == 1);
co_activated_indices_sum = numel(co_activated_indices);
for qq = 1:size(co_activated_indices, 2)
    [h(qq),p(qq),ci{qq},stats{qq}] = ttest(neuron_mean_array{1, 1}(co_activated_indices(qq),sub_window_idx),neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx));
    mean_diff(qq) = mean(neuron_mean_array{1, 1}(co_activated_indices(qq),sub_window_idx) - mean(neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx)));
end

sig_increase_shk_from_large = co_activated_indices(p < 0.05 & mean_diff < 0);
sig_increase_shk_from_large_sum = numel(sig_increase_shk_from_large);
sig_increase_shk_from_large_ind = zeros(1, size(respClass_all_array{1,2}, 2));

sig_increase_shk_from_large_ind(:, sig_increase_shk_from_large) = 1;

shk_activated = respClass_all_array{1,2} == exclusive_activated_session_2 |  respClass_all_array{1,2} == sig_increase_shk_from_large_ind;
shk_activated_sum = sum(shk_activated);

figure; plot(ts1, mean(neuron_mean_array{1, 2}(shk_activated,:))); hold on; plot(ts1,  mean(neuron_mean_array{1, 1}(Block_1_activated,:)));



figure;
hold on; 
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(shk_activated,:)), mean(neuron_sem_array{1, 2}(shk_activated,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 1}(Block_1_activated,:)), mean(neuron_sem_array{1, 1}(Block_1_activated,:)), 'lineProps', {'color', batlowW(iter,:)});


% total_modulated = [sum_activated + sum_inhibited];
% total_co_modulated = [co_activated_sum + co_inhibited_sum];


total_modulated = [(Block_1_activated_sum/neuron_num)*100 (shk_activated_sum/neuron_num)*100];

A = total_modulated;
I = (sig_increase_shk_from_large_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end

%%
total_modulated = [sum_activated + sum_inhibited];
total_co_modulated = [co_activated_sum + co_inhibited_sum];

A = (total_modulated/neuron_num)*100;
I = (total_co_modulated/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
  %Now label each zone 
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end




A = (sum_activated/neuron_num)*100;
I = (excited_to_excited_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end


A = (sum_activated(iter)/sum_activated(iter-1))*100;
I = (excited_to_excited_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end


A = (sum_activated/neuron_num)*100;
I = (exclusive_activated_sum/co_activated_sum)*100
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end



A = (sum_inhibited/neuron_num)*100;
I = (co_inhibited_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end



large_pre_and_shock_post = (respClass_all_array{1,1} == 1 & shk_activated == 1);
large_pre_and_shock_post_sum = sum(large_pre_and_shock_post)
large_pre_to_small_pre = (respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 1);
large_pre_to_small_pre_sum = sum(large_pre_to_small_pre)

overlap = large_pre_and_shock_post == 1 & large_pre_to_small_pre == 1; 
overlap_sum = sum(overlap)

large_responsive_sum = sum(respClass_all_array{1,1} == 1);


total_modulated = [(large_pre_and_shock_post_sum/large_responsive_sum)*100 (large_pre_to_small_pre_sum/large_responsive_sum)*100];
overlap_sum = (sum(overlap)/large_responsive_sum)*100;

A = total_modulated;
I = overlap_sum;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end



%% calculates RESPONSIVE neuron percents for a series of filtering steps. dependent on what filters are used, so beware! 

large_pre_and_shock_post = (respClass_all_array{1,1} == 1 | rrespClass_all_array{1,1} == 2) & (shk_activated == 1);
large_pre_and_shock_post_sum = sum(large_pre_and_shock_post)
large_pre_to_small_pre = (respClass_all_array{1,1} == 1 | respClass_all_array{1,1} == 2) & (respClass_all_array{1,3} == 1 | respClass_all_array{1,3} == 2);
large_pre_to_small_pre_sum = sum(large_pre_to_small_pre)
large_pre_still_large_pre = (respClass_all_array{1,1} == 1) & (respClass_all_array{1,4} == 1);
large_pre_still_large_sum = sum(large_pre_still_large_pre)


overlap_shock_small = large_pre_and_shock_post == 1 & large_pre_to_small_pre == 1; 
overlap_shock_small_sum = sum(overlap_shock_small)

overlap_shock_large = large_pre_and_shock_post == 1 & large_pre_still_large_pre == 1; 
overlap_shock_large_sum = sum(overlap_shock_large)

overlap_small_large = large_pre_to_small_pre == 1 & large_pre_still_large_pre == 1; 
overlap_small_large_sum = sum(overlap_small_large)

overlap_all = large_pre_to_small_pre == 1 & large_pre_still_large_pre == 1 & large_pre_and_shock_post == 1; 
overlap_all_sum = sum(overlap_all)

large_responsive_sum = sum(respClass_all_array{1,1} == 1 | respClass_all_array{1,1} == 2)


total_modulated = [(large_pre_and_shock_post_sum/large_responsive_sum)*100 (large_pre_to_small_pre_sum/large_responsive_sum)*100 (large_pre_still_large_sum/large_responsive_sum)*100];
overlap_sum = [(overlap_shock_small_sum/large_responsive_sum)*100 (overlap_shock_large_sum/large_responsive_sum)*100 (overlap_small_large_sum/large_responsive_sum)*100 (overlap_all_sum/large_responsive_sum)*100]

A = total_modulated;
I = overlap_sum;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
[H, S] = venn(A,I);
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end



%%
figure;
hold on; 
shadedErrorBar(ts1, mean(neuron_mean(respClass_all_array{:,iter} == 1,:)), mean(neuron_sem(respClass_all_array{:,iter} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean(respClass_all_array{:,iter} == 2,:)), mean(neuron_sem(respClass_all_array{:,iter} == 2,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean(respClass_all_array{:,iter} == 3,:)), mean(neuron_sem(respClass_all_array{:,iter} == 3,:)), 'lineProps', {'color', batlowW(iter,:)});

%%

pun_responsive_B1_large_to_B2_null = (respClass_all_array{1,1} == 1 & respClass_all_array{1,1} == shk_activated) & respClass_all_array{1,3} == 2;
sum(pun_responsive_B1_large_to_B2_null)

