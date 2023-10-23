%%THIS CODE MIGHT BE REPLACED BECAUSE IT IS EASIER TO SEPARATE THIS
%%ORIGINAL "MATCHED" DATASET INTO "final" DATASETS WITH JUST THE 2 SESSIONS
%%& THE MATCHED CELLS, THAT WAY YOU CAN USE THE REGULAR
%%ACCESS_RISK_INSCOPIX CODE


load('BLA_panneuronal_Risk_MATCHED_07142023.mat')



%%

paired_sessions = {'Pre_RDT_RM', 'RDT_D1'};
paired_sessions_struct_name = [paired_sessions{1}, '_vs_', paired_sessions{2}];
epoc_to_align = 'choiceTime';
event_to_analyze = {'BLOCK',1,'REW',1.2};
sampling_rate = 0.1;
window_sz = (0:sampling_rate:20-sampling_rate);
ts1 = (-10:sampling_rate:10-sampling_rate);
numMeasurements = size(ts1, 2);
clear neuron_mean neuron_sem neuron_num zall_array zall_to_BL_array zsd_array trials neuron_num

%%

animalIDs = (fieldnames(cellreg_struct.(paired_sessions_struct_name)));
% neuron_num = 0;
neuron_num_first = 0;
neuron_num_second = 0;
% filter_names_idx = cellfun(@ischar,event_to_analyze);
% filter_strings = string(event_to_analyze(filter_names_idx));
neuron_num{1} = 0;
neuron_num{2} = 0;
neuron_mean = cell(size(animalIDs, 1), size(paired_sessions(1,:), 2));
neuron_sem = cell(size(animalIDs, 1), size(paired_sessions(1,:), 2));
% neuron_sem = zeros(1, size(ts1, 2));
for ii = 1:size(animalIDs,1)
    current_animal = char(animalIDs(ii));
    for zz = 1:size(paired_sessions(1,:), 2)
        paired_sessions_current = paired_sessions{1, zz};
        neuron_num{ii, zz} = 0;
        if isfield(cellreg_struct.(paired_sessions_struct_name).(current_animal), paired_sessions_current)
            [data,trials,varargin] = TrialFilter(cellreg_struct.(paired_sessions_struct_name).(current_animal).(paired_sessions_current).(epoc_to_align).uv.BehavData, 'REW', 1.2, 'BLOCK', 3);
            trials = cell2mat(trials);

            for qq = 1:size(cellreg_struct.(paired_sessions_struct_name).(current_animal).(paired_sessions_current).(epoc_to_align).unitXTrials,2)
                neuron_num{ii, zz} = neuron_num{ii, zz}+1;
                neuron_mean{ii, zz}(neuron_num{ii, zz},:) = mean(cellreg_struct.(paired_sessions_struct_name).(current_animal).(paired_sessions_current).(epoc_to_align).unitXTrials(qq).zall(trials,1:numMeasurements)); 
                neuron_sem{ii, zz}(neuron_num{ii, zz},:) = nansem(cellreg_struct.(paired_sessions_struct_name).(current_animal).(paired_sessions_current).(epoc_to_align).unitXTrials(qq).zall(trials,1:numMeasurements));
                
                %uncomment if you want to save any of the data to the
                %unitXTrials directory for each mouse / cell. Would only
                %recommend doing this if filtering behavior by 'ALL,1, so that
                %you capture all trials
                %             final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).zall_window = zall;




                %             kk = 1;
                %             for kk = 1:size(final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(kk).zall,1)
                %                 zall_cell{ii,:} =
                %             end
            end
        elseif ~isfield(final.(current_animal), session_to_analyze)


        end
    end
end
%%
% Initialize an empty array to store the concatenated values
concatenated_values_session_1 = [];
concatenated_values_session_2 = [];

% Loop through each row of neuron_mean
for ii = 1:size(neuron_mean, 1)
    % Get the cell in column 1 for the current row
    cell_data_1 = neuron_mean{ii, 1};
    cell_data_2 = neuron_mean{ii, 2};
    % Concatenate the values from the cell to the overall array
    concatenated_values_session_1 = [concatenated_values_session_1; cell_data_1];
    concatenated_values_session_2 = [concatenated_values_session_2; cell_data_2];
end
%% get correlation coefficients across all cells
% Assuming you have two arrays: concatenated_values_session_1 and concatenated_values_session_2

% Initialize an array to store the correlations
correlations = zeros(size(concatenated_values_session_1, 1), 1);

% Loop through each row and calculate the correlation
for i = 1:size(concatenated_values_session_1, 1)
    row1 = concatenated_values_session_1(i, :);
    row2 = concatenated_values_session_2(i, :);
    
    % Calculate the correlation coefficient
    correlations(i) = corr(row1', row2');
end

% 'correlations' now contains the correlation coefficients for each pair of rows

mean_corr = mean(correlations(:,1));


%%
% Define the time windows
baseline_start = -10;
baseline_end = -7;
test_start = 0;
test_end = 3;

% Find the indices corresponding to the time windows in ts1
baseline_indices = find(ts1 >= baseline_start & ts1 <= baseline_end);
test_indices = find(ts1 >= test_start & ts1 <= test_end);

% Initialize arrays to store modulation classification for both sessions
modulation_classification_session_1 = zeros(size(concatenated_values_session_1, 1), 1)';
modulation_classification_session_2 = zeros(size(concatenated_values_session_1, 1), 1)';

% Loop through each row and perform the Wilcoxon signed-rank test
for i = 1:size(concatenated_values_session_1, 1)
    % Extract data for the baseline and test periods for the current row
    baseline_data_session_1 = concatenated_values_session_1(i, baseline_indices);
    test_data_session_1 = concatenated_values_session_1(i, test_indices);
    
    % Perform the Mann-Whitney U test
    [p_value_session_1, ~, stats_session_1] = ranksum(baseline_data_session_1, test_data_session_1);

    % Calculate the mean activity for both windows in session 1
    mean_baseline_session_1 = mean(baseline_data_session_1);
    mean_test_session_1 = mean(test_data_session_1);
    
    % Set a significance threshold (e.g., 0.05) for the Mann-Whitney U test
    significance_threshold = 0.05;
    
    % Classify based on modulation for session 1
    if p_value_session_1 < significance_threshold
        % Responsive
        if mean_test_session_1 > mean_baseline_session_1
            modulation_classification_session_1(i) = 1;  % Positively modulated
        elseif mean_test_session_1 < mean_baseline_session_1
            modulation_classification_session_1(i) = 2;  % Negatively modulated
        else
            modulation_classification_session_1(i) = 3;  % Not significantly modulated
        end
    else
        % Not responsive
        modulation_classification_session_1(i) = 3;  % Not significantly modulated
    end
end

% Loop through each row and perform the Mann-Whitney U test for session 2
for i = 1:size(concatenated_values_session_2, 1)
    % Extract data for the baseline and test periods for the current row in session 2
    baseline_data_session_2 = concatenated_values_session_2(i, baseline_indices);
    test_data_session_2 = concatenated_values_session_2(i, test_indices);
    
    % Perform the Mann-Whitney U test
    [p_value_session_2, ~, stats_session_2] = ranksum(baseline_data_session_2, test_data_session_2);
    
    % Calculate the mean activity for both windows in session 2
    mean_baseline_session_2 = mean(baseline_data_session_2);
    mean_test_session_2 = mean(test_data_session_2);
    
    % Set a significance threshold (e.g., 0.05) for the Mann-Whitney U test
    significance_threshold = 0.05;
    
    % Classify based on modulation for session 2
    if p_value_session_2 < significance_threshold
        % Responsive
        if mean_test_session_2 > mean_baseline_session_2
            modulation_classification_session_2(i) = 1;  % Positively modulated
        elseif mean_test_session_2 < mean_baseline_session_2
            modulation_classification_session_2(i) = 2;  % Negatively modulated
        else
            modulation_classification_session_2(i) = 3;  % Not significantly modulated
        end
    else
        % Not responsive
        modulation_classification_session_2(i) = 3;  % Not significantly modulated
    end
end

% 'modulation_classification' now contains the classification (1 for positively modulated, 2 for negatively modulated, 3 for not significantly modulated) for each row

figure; plot(ts1, mean(concatenated_values_session_1(modulation_classification_session_1 == 1,:)));
hold on; plot(ts1, mean(concatenated_values_session_2(modulation_classification_session_2 == 1,:)));
figure; plot(ts1, mean(concatenated_values_session_1(modulation_classification_session_1 == 2,:)));
hold on; plot(ts1, mean(concatenated_values_session_2(modulation_classification_session_2 == 2,:)));
figure; plot(ts1, mean(concatenated_values_session_1(modulation_classification_session_1 == 3,:)));
hold on; plot(ts1, mean(concatenated_values_session_2(modulation_classification_session_2 == 3,:)));

%%


start_time = 0; % sub-window start time
end_time = 4; % sub-window end time

% Find the indices in ts1 that correspond to the sub-window
sub_window_idx = ts1 >= start_time & ts1 <= end_time;

% Extract the corresponding columns from neuron_mean
sub_window_activity = neuron_mean(:, sub_window_idx);
mean_activity = mean(sub_window_activity, 2);

[~, sorted_idx] = sort(mean_activity, 'descend');

% Sort neuron_mean based on the sorted indices
ranked_neurons = neuron_mean(sorted_idx, :);

figure;
% Generate the heatmap
imagesc(ts1, 1, ranked_neurons);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

% Reverse the y-axis so that the highest mean activity is at the top
set(gca, 'YDir', 'reverse');

%% This code was generated by ChatGPT
% perform PCA
[coeff, score, ~, ~, explained] = pca(neuron_mean);

% plot the mean activity of the main principle components
num_pcs_to_plot = 3; % choose the number of principle components to plot
pcs_to_plot = 1:num_pcs_to_plot;
figure;
hold on;
for i = pcs_to_plot
    plot(ts1, coeff(:,i));
end
xlabel('Time (s)');
ylabel('PCA weight');
legend(strcat('PC', string(pcs_to_plot)));
title('Mean activity of main principle components');

% plot the percentage of variance explained by each principle component
figure;
pareto(explained);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Variance Explained by Principal Components');

% determine which neurons are assigned to each principle component
[~, max_scores_idx] = max(abs(score), [], 2); % find the index of the max score for each neuron
neurons_per_pc = accumarray(max_scores_idx, (1:size(neuron_mean,1))', [], @(x) {x}); % group the neurons by principle component
disp('Neurons assigned to each principle component:');
for i = 1:num_pcs_to_plot
    fprintf('PC %d: Neurons %s\n', i, num2str(neurons_per_pc{i}));
end

% calculate the correlation matrix between principle components
score_norm = zscore(score); % normalize each principle component to have mean 0 and std 1
pc_corr = corrcoef(score_norm);
figure;
imagesc(pc_corr);
colorbar;
xlabel('Principal Component');
ylabel('Principal Component');
title('Correlation Matrix between Principal Components');


%% from Sean to do kmeans on traces
[idx,C,sumdist3] = kmeans(neuron_mean,4,'Distance','correlation','Display','final', 'Replicates', 200,'Start','uniform');

figure; plot(ts1, neuron_mean(idx == 1, :));
figure; plot(ts1, neuron_mean(idx == 2, :));
figure; plot(ts1, neuron_mean(idx == 3, :));
figure; plot(ts1, neuron_mean(idx == 4, :));


figure; plot(ts1, mean(neuron_mean(idx == 1, :)));
hold on; 
plot(ts1, mean(neuron_mean(idx == 2, :)));
hold on; 
plot(ts1, mean(neuron_mean(idx == 3, :)));
hold on; 
plot(ts1, mean(neuron_mean(idx == 4, :)));

%%
% perform k-means clustering on the principle components
num_clusters = 4; % choose the number of clusters
[idx, centroids] = kmeans(score(:,1:num_pcs_to_plot), num_clusters); % cluster the neurons based on the principle components

% assign each neuron to a cluster
neurons_per_cluster = accumarray(idx, (1:size(neuron_mean,1))', [], @(x) {x}); % group the neurons by cluster
disp('Neurons assigned to each cluster:');
for i = 1:num_clusters
    fprintf('Cluster %d: Neurons %s\n', i, num2str(neurons_per_cluster{i}));
end

figure; plot(ts1, mean(neuron_mean(neurons_per_cluster{1,1}, :)));
hold on;
plot(ts1, mean(neuron_mean(neurons_per_cluster{2,1}, :)));
hold on;
plot(ts1, mean(neuron_mean(neurons_per_cluster{3,1}, :)));
hold on;
plot(ts1, mean(neuron_mean(neurons_per_cluster{4,1}, :)));