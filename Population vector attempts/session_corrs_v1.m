% load dataset that contains data where events have been identified
% run raster_with_representative_neurons_v2.m for the mouse for whom you
% will analyze below
% also you should first run block_wise_changes_v1.m to get the blocktimes


select_mouse = 'BLA_Insc_32';

% for RDT D1 BLA_Insc_25:
%prechoice neuron num 46
%postchoice rew num 38
%consumption num 39
%shock num 11

select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';


neuronTypes = respClass_all_array_mouse{select_mouse_index, 3};
neuralActivity = final.(select_mouse).(first_session).CNMFe_data.C_raw;
BehavData = final.(select_mouse).(first_session).uv.BehavData;



%%


% Assume neuronTypes is a vector with neuron types (1, 2, 3)
% and neuralActivity is the matrix (neurons x time points)

uniqueTypes = unique(neuronTypes);
for i = 1:length(uniqueTypes)
    typeIdx = (neuronTypes == uniqueTypes(i));
    activitySubset = neuralActivity(typeIdx, :);
    
    % Calculate correlation matrix
    similarityMatrix = corr(activitySubset');
    
    % (Optional) Summarize the similarity, e.g., by taking the mean of the upper triangle
    similarityScore(i) = mean(similarityMatrix(triu(true(size(similarityMatrix)), 1)));
end


%%
windowSize = 100; % Number of time points in the window
uniqueTypes = unique(neuronTypes);
for t = 1:(size(neuralActivity, 2) - windowSize + 1)
    for i = 1:length(uniqueTypes)
        typeIdx = (neuronTypes == uniqueTypes(i));
        % use specific subset of neurons
        activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
        % uncomment if you want to use all neurons
        %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
        similarityMatrix = corr(activitySubset');
        similarityOverTime(i, t) = mean(similarityMatrix(triu(true(size(similarityMatrix)), 1)));
    end
end





% Plot the similarity over time
figure;
plot(time_array, similarityOverTime(1,1:trim_length));
xlabel('Time');
ylabel('Mean Similarity');
% legend({'Type 1', 'Type 2', 'Type 3'});
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')

%%
% Example block variables
% block_1_mouse = [start_time_1, stop_time_1];
% block_2_mouse = [start_time_2, stop_time_2];
% block_3_mouse = [start_time_3, stop_time_3];

% Find indices in time_array that correspond to the start and stop times
[~, startIdx1] = min(abs(time_array - block_1_mouse(1, 1)));
[~, stopIdx1] = min(abs(time_array - block_1_mouse(1, 2)));

[~, startIdx2] = min(abs(time_array - block_2_mouse(1, 1)));
[~, stopIdx2] = min(abs(time_array - block_2_mouse(1, 2)));

[~, startIdx3] = min(abs(time_array - block_3_mouse(1, 1)));
[~, stopIdx3] = min(abs(time_array - block_3_mouse(1, 2)));

% Extract similarity scores for each block
similarityStage1 = similarityOverTime(1, startIdx1:stopIdx1);
similarityStage2 = similarityOverTime(1, startIdx2:stopIdx2);
similarityStage3 = similarityOverTime(1, startIdx3:stopIdx3);

% Calculate the mean similarity for each stage across the specified time points
meanSimilarityStage1 = mean(similarityStage1, 2); % Mean across time for each neuron type
meanSimilarityStage2 = mean(similarityStage2, 2);
meanSimilarityStage3 = mean(similarityStage3, 2);

% Display mean similarity for each stage
disp('Mean Similarity for Stage 1:');
disp(meanSimilarityStage1);

disp('Mean Similarity for Stage 2:');
disp(meanSimilarityStage2);

disp('Mean Similarity for Stage 3:');
disp(meanSimilarityStage3);

% Optionally, plot the results
figure;
bar([meanSimilarityStage1, meanSimilarityStage2, meanSimilarityStage3]);
xlabel('Neuron Type');
ylabel('Mean Similarity');
legend({'Stage 1', 'Stage 2', 'Stage 3'});
title('Mean Similarity Across Stages');


%%

next_neuron = 0;

PV_consumption_mouse = [];
for ff = 1:size(respClass_all_array_mouse{select_mouse_index, 3}, 2)
    if respClass_all_array_mouse{select_mouse_index, 3}(ff) == 1
        next_neuron = next_neuron+1;
        PV_consumption_mouse (next_neuron) = neuron_mean_mouse{select_mouse_index, 3}(ff, 91);



    end

end




block_1_data_for_consumption_neurons = block_1_ca_mouse{1, 1}(respClass_all_array_mouse{select_mouse_index, 3} == 1, :);

num_samples = size(block_1_data_for_consumption_neurons, 2)
% sampling_frequency = (final.(select_mouse).(first_session).choiceTime.uv.dt)*100;
sampling_frequency = (final.(select_mouse).(first_session).uv.dt)*100;
% Create a time array
time_array = (0:(num_samples-1)) / sampling_frequency;


windowSize = 100; % Number of time points in the window
for t = 1:size(block_1_data_for_consumption_neurons, 2)
    % for i = 1:length(uniqueTypes)
        % typeIdx = (neuronTypes == uniqueTypes(i));
        % use specific subset of neurons
        activitySubset = block_1_data_for_consumption_neurons(:, t);
        % uncomment if you want to use all neurons
        %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
        similarityMatrix = corrcoef(activitySubset', PV_consumption_mouse);
        similarityOverTime(t) = similarityMatrix(2);
    % end
end

figure; plot(time_array(1:size(similarityOverTime, 2)), similarityOverTime)
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')



%%
PV_prechoice_all_mouse = [];
for ff = 1:size(zall_mouse{select_mouse_index, 1}, 2)
    for cc = 1:size(zall_mouse{select_mouse_index,1}{1, ff}, 1)
        PV_prechoice_all_mouse(cc, ff) = mean(zall_mouse{select_mouse_index,1}{1, ff}(cc, ts1 >= -4 & ts1 <= 0));

    end

end


PV_prechoice_all_mouse_just_prechoice = PV_prechoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, 1} == 1); 
PV_prechoice_all_mouse_just_prechoice = [PV_prechoice_all_mouse_just_prechoice PV_prechoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, 1} ~= 1)];
mean_PV_prechoice_all_mouse = mean(PV_prechoice_all_mouse)';

hold on;
figure; imagesc(1:size(PV_prechoice_all_mouse_just_prechoice, 1), [], PV_prechoice_all_mouse_just_prechoice')
colormap gray;
colorbar;
clim([0 1.5]);
hold off;

% Define the position for the first green line
y1_start = 1;  % Start at the first neuron
y1_end = sum(respClass_all_array_mouse{select_mouse_index, 1} == 1);  % End at the number of neurons in the first group

% Define the position for the second line
y2_start = y1_end + 1;  % Start right after the first group
y2_end = size(PV_prechoice_all_mouse_just_prechoice, 2);  % End at the last neuron

% Draw the first vertical green line
hold on;
line([1 1], [y1_start y1_end], 'Color', 'g', 'LineWidth', 2);

% Draw the second vertical line (can be any color, e.g., black here)
line([1 1], [y2_start y2_end], 'Color', 'r', 'LineWidth', 2);

hold off;


ca = final.(select_mouse).(first_session).CNMFe_data.(ca_data_type);
ca_zscored = zscore(ca, [], 2);
time_array = final.(select_mouse).(first_session).time;

for t = 1:size(ca_zscored, 2)
    % for i = 1:length(uniqueTypes)
        % typeIdx = (neuronTypes == uniqueTypes(i));
        % use specific subset of neurons
        activitySubset = ca_zscored(:, t);
        % uncomment if you want to use all neurons
        %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
        similarityMatrix = corrcoef(activitySubset, mean_PV_prechoice_all_mouse);
        prechoice_similarityOverTime(t) = similarityMatrix(2);
    % end
end

figure; plot(time_array, prechoice_similarityOverTime(1, 1:size(time_array, 1)))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')


b1_prechoice_corr = mean(prechoice_similarityOverTime(1, time_array > block_1_mouse(select_mouse_index, 1) & time_array <= block_1_mouse(select_mouse_index, 2)));
b2_prechoice_corr = mean(prechoice_similarityOverTime(1, time_array > block_2_mouse(select_mouse_index, 1) & time_array <= block_2_mouse(select_mouse_index, 2)));
b3_prechoice_corr = mean(prechoice_similarityOverTime(1, time_array > block_3_mouse(select_mouse_index, 1) & time_array <= block_3_mouse(select_mouse_index, 2)));


%%

PV_postchoice_all_mouse = [];
for ff = 1:size(zall_mouse{select_mouse_index, 2}, 2)
    for cc = 1:size(zall_mouse{select_mouse_index,2}{1, ff}, 1)
        PV_postchoice_all_mouse(cc, ff) = mean(zall_mouse{select_mouse_index,2}{1, ff}(cc, ts1 >= 0 & ts1 <= 2));

    end

end


PV_postchoice_all_mouse_just_consumption = PV_postchoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, 2} == 1); 
PV_postchoice_all_mouse_just_consumption = [PV_postchoice_all_mouse_just_consumption PV_postchoice_all_mouse(:, respClass_all_array_mouse{select_mouse_index, 2} ~= 1)];
mean_PV_postchoice_all_mouse = mean(PV_postchoice_all_mouse)';
hold on;
figure; imagesc(1:size(PV_postchoice_all_mouse_just_consumption, 1), [], PV_postchoice_all_mouse_just_consumption')
colormap gray;
colorbar;
clim([0 2]);
hold off;

hold off;

% Define the position for the first green line
y1_start = 1;  % Start at the first neuron
y1_end = sum(respClass_all_array_mouse{select_mouse_index, 2} == 1);  % End at the number of neurons in the first group

% Define the position for the second line
y2_start = y1_end + 1;  % Start right after the first group
y2_end = size(PV_postchoice_all_mouse_just_consumption, 2);  % End at the last neuron

% Draw the first vertical green line
hold on;
line([1 1], [y1_start y1_end], 'Color', 'g', 'LineWidth', 2);

% Draw the second vertical line (can be any color, e.g., black here)
line([1 1], [y2_start y2_end], 'Color', 'r', 'LineWidth', 2);

hold off;


ca = final.(select_mouse).(session_to_analyze).CNMFe_data.(ca_data_type);
ca_zscored = zscore(ca, [], 2);
time_array = final.(select_mouse).(first_session).time;

for t = 1:size(ca_zscored, 2)
    % for i = 1:length(uniqueTypes)
        % typeIdx = (neuronTypes == uniqueTypes(i));
        % use specific subset of neurons
        activitySubset = ca_zscored(:, t);
        % uncomment if you want to use all neurons
        %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
        similarityMatrix = corrcoef(activitySubset, mean_PV_postchoice_all_mouse);
        postchoice_similarityOverTime(t) = similarityMatrix(2);
    % end
end

figure; plot(time_array, postchoice_similarityOverTime(1, 1:size(time_array, 1)))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')



b1_postchoice_corr = mean(postchoice_similarityOverTime(1, time_array > block_1_mouse(select_mouse_index, 1) & time_array <= block_1_mouse(select_mouse_index, 2)));
b2_postchoice_corr = mean(postchoice_similarityOverTime(1, time_array > block_2_mouse(select_mouse_index, 1) & time_array <= block_2_mouse(select_mouse_index, 2)));
b3_postchoice_corr = mean(postchoice_similarityOverTime(1, time_array > block_3_mouse(select_mouse_index, 1) & time_array <= block_3_mouse(select_mouse_index, 2)));




%%

PV_consumption_all_mouse = [];
for ff = 1:size(zall_mouse{select_mouse_index, 3}, 2)
    for cc = 1:size(zall_mouse{select_mouse_index,3}{1, ff}, 1)
        PV_consumption_all_mouse(cc, ff) = mean(zall_mouse{select_mouse_index,3}{1, ff}(cc, ts1 >= 1 & ts1 <= 3));

    end

end


PV_consumption_all_mouse_just_consumption = PV_consumption_all_mouse(:, respClass_all_array_mouse{select_mouse_index, 3} == 1); 
PV_consumption_all_mouse_just_consumption = [PV_consumption_all_mouse_just_consumption PV_consumption_all_mouse(:, respClass_all_array_mouse{select_mouse_index, 3} ~= 1)];
mean_PV_consumption_all_mouse = mean(PV_consumption_all_mouse)';
hold on;
figure; imagesc(1:size(PV_consumption_all_mouse_just_consumption, 1), [], PV_consumption_all_mouse_just_consumption')
colormap gray;
colorbar;
clim([0 2]);
hold off;

hold off;

% Define the position for the first green line
y1_start = 1;  % Start at the first neuron
y1_end = sum(respClass_all_array_mouse{select_mouse_index, 3} == 1);  % End at the number of neurons in the first group

% Define the position for the second line
y2_start = y1_end + 1;  % Start right after the first group
y2_end = size(PV_consumption_all_mouse_just_consumption, 2);  % End at the last neuron

% Draw the first vertical green line
hold on;
line([1 1], [y1_start y1_end], 'Color', 'g', 'LineWidth', 2);

% Draw the second vertical line (can be any color, e.g., black here)
line([1 1], [y2_start y2_end], 'Color', 'r', 'LineWidth', 2);

hold off;


ca = final.(select_mouse).(session_to_analyze).CNMFe_data.(ca_data_type);
ca_zscored = zscore(ca, [], 2);
time_array = final.(select_mouse).(first_session).time;

for t = 1:size(ca_zscored, 2)
    % for i = 1:length(uniqueTypes)
        % typeIdx = (neuronTypes == uniqueTypes(i));
        % use specific subset of neurons
        activitySubset = ca_zscored(:, t);
        % uncomment if you want to use all neurons
        %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
        similarityMatrix = corrcoef(activitySubset, mean_PV_consumption_all_mouse);
        consumption_similarityOverTime(t) = similarityMatrix(2);
    % end
end

figure; plot(time_array, consumption_similarityOverTime(1, 1:size(time_array, 1)))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')


b1_consumption_corr = mean(consumption_similarityOverTime(1, time_array > block_1_mouse(select_mouse_index, 1) & time_array <= block_1_mouse(select_mouse_index, 2)));
b2_consumption_corr = mean(consumption_similarityOverTime(1, time_array > block_2_mouse(select_mouse_index, 1) & time_array <= block_2_mouse(select_mouse_index, 2)));
b3_consumption_corr = mean(consumption_similarityOverTime(1, time_array > block_3_mouse(select_mouse_index, 1) & time_array <= block_3_mouse(select_mouse_index, 2)));





%%
figure; plot(time_array, prechoice_similarityOverTime(1, 1:size(time_array, 1)))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')
hold on; plot(time_array, postchoice_similarityOverTime(1, 1:size(time_array, 1)))
hold on; plot(time_array, consumption_similarityOverTime(1, 1:size(time_array, 1)))

%% attempting to get something similar to Courtin 5e
eTS = BehavData.choiceTime;

% Define time windows
pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [1 3];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

% Find indices corresponding to each time window
pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

for t = 1:size(eTS,1)

    timeWin = [eTS(t)+uv.evtWin(1,1):uv.dt:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
    BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];



    if min(timeWin) > min(time_array) && max(timeWin) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
        % get unit event counts in trials
        % get unit ca traces in trials
        idx = time_array > min(timeWin) & time_array < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);
        bl_idx = time_array > min(BL_win) & time_array < max(BL_win);
        %caTraceTrials(t,1:sum(idx)) = unitTrace(idx);               %store the evoked calcium trace around each event   (see below, comment out if dont want normalized to whole trace)
        prechoice_similarityOverTimeTrials(t,1:sum(idx)) = prechoice_similarityOverTime(idx);

        % prechoice_similarityOverTimeTrials_zb(t,:) = nanmean(unitTrace(bl_idx)); %baseline mean
        prechoice_similarityOverTimeTrials_zb_window(t,:) = nanmean(prechoice_similarityOverTimeTrials(t,:));
        % prechoice_similarityOverTimeTrials_zsd(t,:) = nanstd(unitTrace(bl_idx)); %baseline std
        prechoice_similarityOverTimeTrials_zsd_window(t,:) = nanstd(prechoice_similarityOverTimeTrials(t,:));
        tmp = 0;
        for j = 1:size(prechoice_similarityOverTimeTrials,2)
            tmp = tmp+1;
            % prechoice_similarityOverTimeTrials_zall_baselined(t,tmp) = (prechoice_similarityOverTimeTrials(t,j) - prechoice_similarityOverTimeTrials_zb(t))/prechoice_similarityOverTimeTrials_zsd(t);
            prechoice_similarityOverTimeTrials_raw(t, tmp) = prechoice_similarityOverTimeTrials(t,j);
            prechoice_similarityOverTimeTrials_zall_window(t,tmp) = (prechoice_similarityOverTimeTrials(t,j) - prechoice_similarityOverTimeTrials_zb_window(t))/prechoice_similarityOverTimeTrials_zsd_window(t);
            % prechoice_similarityOverTimeTrials_zall_session(t,tmp) = (prechoice_similarityOverTimeTrials(t,j) - zb_session(u))/zsd_session(u);
        end
        clear j;

    end


end

block_1_inds = BehavData.Block == 1 & BehavData.Blank_Touch == 0 & BehavData.omissionALL == 0; 
block_2_inds = BehavData.Block == 2 & BehavData.Blank_Touch == 0 & BehavData.omissionALL == 0; 
block_3_inds = BehavData.Block == 3 & BehavData.Blank_Touch == 0 & BehavData.omissionALL == 0; 

prechoice_period_b1 = nanmean(prechoice_similarityOverTimeTrials_raw(block_1_inds, pre_choice_indices), 2)
prechoice_period_b2 = nanmean(prechoice_similarityOverTimeTrials_raw(block_2_inds, pre_choice_indices), 2)
prechoice_period_b3 = nanmean(prechoice_similarityOverTimeTrials_raw(block_3_inds, pre_choice_indices), 2)

pre_choice_over_time = [prechoice_period_b1; prechoice_period_b2; prechoice_period_b3];
figure; plot(pre_choice_over_time)