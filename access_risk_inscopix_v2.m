
iter = 0;
%%
% load('BLA-NAcShell_Risk_2023_09_15.mat')

% load('BLA_panneuronal_Risk_2024_01_04.mat')

% load('NAcSh_D2_Cre-OFF_GCAMP_all.mat')

load('BLA_panneuronal_matched_Pre_RDT_RM_vs_RDT_D1_01042024.mat')

% load('BLA_panneuronal_Risk_matched_RM_D1_vs_Pre_RDT_RM.mat')

% load('BLA_NAcSh_Risk_matched_Pre_RDT_RM_vs_RDT_D1.mat')




%%

session_to_analyze = 'RDT_D1';
epoc_to_align = 'choiceTime';
event_to_analyze = {'BLOCK',1,'REW',1.2};

window_sz = (0:.1:20-0.1);
ts1 = (-10:.1:10-0.1);

neuron_mean_concat = [];

neuron_mean_unnorm_concat = [];
neuron_sem_concat = [];


if exist('iter', 'var') == 1
   
elseif exist('iter', 'var') == 0
    iter = 0;
end

clear neuron_mean neuron_sem neuron_num zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized

%%

animalIDs = (fieldnames(final));
neuron_num = 0;

sum_trials_per_iter = 0;
filter_names_idx = cellfun(@ischar,event_to_analyze);
filter_strings = string(event_to_analyze(filter_names_idx));

% neuron_sem = zeros(1, size(ts1, 2));
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    
   
    if isfield(final.(currentanimal), session_to_analyze)
        [data,trials,varargin] = TrialFilter(final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.BehavData, 'REW', 1.2, 'BLOCK', 1);
        behav_tbl_temp{ii,:} = data;
        trials = cell2mat(trials);
        % trials_per_mouse{ii, iter+1} = trials;
        sum_trials_per_iter = sum_trials_per_iter+size(trials, 1);
        
        evtWinSpan = max(final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.evtWin) - min(final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.evtWin);
        numMeasurements = round(evtWinSpan/final.(currentanimal).(session_to_analyze).(epoc_to_align).uv.dt); %need to round due to odd frame rate
        
        for qq = 1:size(final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials,2)
            neuron_num = neuron_num+1; 

                caTraceTrials = final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).caTraces(trials,1:numMeasurements);
                for h = 1:size(caTraceTrials,1)
                    zb(h) = mean(caTraceTrials(h,:)); %baseline mean
                    zsd(h) = std(caTraceTrials(h,:)); %baseline std
                    tmp = 0;
                    for j = 1:size(caTraceTrials(1:length(window_sz)),2)
                        tmp = tmp+1;
                        zall(h,tmp) = (caTraceTrials(h,tmp) - zb(h))/zsd(h);
                    end

                end
                caTraceTrials_mean(neuron_num,:) = mean(caTraceTrials,1);
                zall_array(neuron_num) = {zall};
                neuron_mean(neuron_num,:) = mean(zall,1);
                neuron_mean_unnormalized(neuron_num,:) = mean(caTraceTrials,1);
                neuron_sem(neuron_num,:) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                zsd_array(neuron_num) = {zsd};
                zall_to_BL_array(neuron_num) = {final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).zall(trials,:)};
                zall_mouse{ii, iter+1}(qq) = {zall};
                caTraceTrials_mouse{ii, iter+1}(qq) = {caTraceTrials};
                neuron_mean_mouse{ii, iter+1}(qq,: ) = mean(zall, 1);
                neuron_sem_mouse{ii, iter+1}(qq,: ) = nanstd(zall,1)/(sqrt(size(zall, 1)));
                trials_per_mouse{ii, iter+1} = trials;
                %uncomment if you want to save any of the data to the
                %unitXTrials directory for each mouse / cell. Would only
                %recommend doing this if filtering behavior by 'ALL,1, so that
                %you capture all trials
                %             final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(qq).zall_window = zall;
           
        clear zall caTraceTrials zb zsd;
            processed_data = struct;
            
            
%             kk = 1;
%             for kk = 1:size(final.(currentanimal).(session_to_analyze).(epoc_to_align).unitXTrials(kk).zall,1)
%                 zall_cell{ii,:} =  
%             end
        end
    elseif ~isfield(final.(currentanimal), session_to_analyze)
        

    end
    
end

iter = iter+1;
varargin_list{iter,:} = varargin;
behav_tbl_iter{iter, :} = behav_tbl_temp;
sum_trials_per_iter_array{iter} = sum_trials_per_iter;
clear behav_tbl_temp

if iter == 1
    neuron_mean_concat_PCA = [neuron_mean];
    neuron_mean_unnorm_concat_PCA = [neuron_mean_unnormalized];
    neuron_sem_concat_PCA = [neuron_sem];

    neuron_mean_concat_DECODING = [neuron_mean];
    neuron_mean_unnorm_concat_DECODING = [neuron_mean_unnormalized];
    neuron_sem_concat_DECODING = [neuron_sem];



elseif iter > 1
    neuron_mean_concat_PCA = [neuron_mean_concat_PCA, neuron_mean];
    neuron_mean_unnorm_concat_PCA = [neuron_mean_unnorm_concat_PCA, neuron_mean_unnormalized];
    neuron_sem_concat_PCA = [neuron_sem_concat_PCA, neuron_sem];

    neuron_mean_concat_DECODING = [neuron_mean_concat_DECODING; neuron_mean];
    neuron_mean_unnorm_concat_DECODING = [neuron_mean_unnorm_concat_DECODING; neuron_mean_unnormalized];
    neuron_sem_concat_DECODING = [neuron_sem_concat_DECODING ; neuron_sem];

end




%USE THIS VERSION FOR PCA FOR NOW - STILL NEED TO FIGURE OUT HOW TO GET THE
%VERSION BELOW TO WORK WITH PCA
% if iter == 1
%     neuron_mean_concat = [neuron_mean];
%     neuron_mean_unnorm_concat = [neuron_mean_unnormalized];
%     neuron_sem_concat = [neuron_sem];
% elseif iter > 1
%     neuron_mean_concat = [neuron_mean_concat, neuron_mean];
%     neuron_mean_unnorm_concat = [neuron_mean_unnorm_concat, neuron_mean_unnormalized];
%     neuron_sem_concat = [neuron_sem_concat, neuron_sem];
% end


% if iter == 1
%     neuron_mean_concat = [neuron_mean];
%     neuron_mean_unnorm_concat = [neuron_mean_unnormalized];
%     neuron_sem_concat = [neuron_sem];
% elseif iter > 1
%     neuron_mean_concat = [neuron_mean_concat; neuron_mean];
%     neuron_mean_unnorm_concat = [neuron_mean_unnorm_concat; neuron_mean_unnormalized];
%     neuron_sem_concat = [neuron_sem_concat; neuron_sem];
% end
% 
%%
% get the behavioral data corresponding to each block & stack them on top
% of each other. These data can be used to plot the median or mean choice
% time on a PCA graph, for example


% Initialize the concatenated table
concatenatedTable = table();

% Iterate through the 3x1 cell array
for i = 1:numel(behav_tbl_iter)
    % Assuming each cell contains a 12x1 cell array of tables
    twelveByOneCellArray = behav_tbl_iter{i};
    
    % Initialize a temporary table to store the concatenated tables for this cell
    tempTable = table();
    
    % Iterate through the 12x1 cell array
    for j = 1:numel(twelveByOneCellArray)
        % Assuming each cell in the 12x1 cell array contains a table
        currentTable = twelveByOneCellArray{j};
        
        % Concatenate the current table to the temporary table vertically
        tempTable = vertcat(tempTable, currentTable);
    end
    
    % Concatenate the temporary table to the overall concatenated table vertically
    concatenatedTable = vertcat(concatenatedTable, tempTable);
end

median_choice_time_block_1 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_choice_time_block_2 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_choice_time_block_3 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));

median_collect_time_block_1 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_collect_time_block_2 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_collect_time_block_3 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));



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
    plot(ts1(:, 1:numMeasurements), coeff(:,i));
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