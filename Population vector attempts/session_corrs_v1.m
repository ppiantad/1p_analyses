% load dataset that contains data where events have been identified
% run raster_with_representative_neurons_v2.m for the mouse for whom you
% will analyze below
% also you should first run block_wise_changes_v1.m to get the blocktimes


select_mouse = 'BLA_Insc_35';

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
        similarityOverTime(t) = similarityMatrix(2);
    % end
end

figure; plot(time_array, similarityOverTime(1, 1:size(time_array, 1)))
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
clim([0 2]);
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
        similarityOverTime(t) = similarityMatrix(2);
    % end
end

figure; plot(time_array, similarityOverTime(1, 1:size(time_array, 1)))
xline(BehavData.stTime(BehavData.bigSmall == 1.2), '--b')
xline(BehavData.stTime(BehavData.bigSmall == 0.3), '--g')
xline(BehavData.choiceTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--r')
xline(BehavData.collectionTime(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3), '--k')
xline(BehavData.choiceTime(BehavData.shock == 1), '--y')


b1_prechoice_corr = mean(similarityOverTime(1, time_array > block_1_mouse(select_mouse_index, 1) & time_array <= block_1_mouse(select_mouse_index, 2)));
b2_prechoice_corr = mean(similarityOverTime(1, time_array > block_2_mouse(select_mouse_index, 1) & time_array <= block_2_mouse(select_mouse_index, 2)));
b3_prechoice_corr = mean(similarityOverTime(1, time_array > block_3_mouse(select_mouse_index, 1) & time_array <= block_3_mouse(select_mouse_index, 2)));

