% load dataset that contains data where events have been identified
% run raster_with_representative_neurons_v2.m for the mouse for whom you
% will analyze below
% also you should first run block_wise_changes_v1.m to get the blocktimes


neuronTypes = respClass_all_array_mouse{1, 10};
neuralActivity = final.BLA_Insc_24.RDT_D1.CNMFe_data.C_raw;



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

