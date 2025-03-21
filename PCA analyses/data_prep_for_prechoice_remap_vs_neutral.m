%% This code sets up PCA for Prechoice Remapped vs. Neutral comparison on Block 1

% load 10x dataset - run block_wise_changes_v1
% run data_loop to get Block 1 data
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1

remapped_pre_on_block_1_means = zall_mean_all_array{1, 11}(remapped_prechoice == 1, :);

true_neutral_on_block_1_means = zall_mean_all_array{1, 11}(true_neutral == 1, :);


remapped_pre_on_block_1_sems = sem_all_array{1, 11}(remapped_prechoice == 1, :);

true_neutral_on_block_1_sems = sem_all_array{1, 11}(true_neutral == 1, :);

figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, nanmean(remapped_pre_on_block_1_means), nanmean(remapped_pre_on_block_1_sems), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, nanmean(true_neutral_on_block_1_means), nanmean(true_neutral_on_block_1_sems), 'lineProps', {'color', 'b'});
legend([h(1).mainLine h(2).mainLine], 'new (safe block)', 'new (risky blocks)')
xlim([-4 4]);
ylim([-0.6 0.7])


%%

prechoice_and_neutral_array = {remapped_pre_on_block_1_means, true_neutral_on_block_1_means}

% match the number of neurons in each group (if necessary - data will stay
% the same if sizes are equal
clear idx temp
for i = 1:length(prechoice_and_neutral_array)

    temp(i) = (size(prechoice_and_neutral_array{1, i},1));
end

[minSize,minIdx] = min(temp);

for i = 1:length(prechoice_and_neutral_array)
    if i ~= minIdx
        idx_for_subsample = randperm(temp(i),minSize);
        idx_for_subsample = sort(idx_for_subsample);
        zdataTemp{i} = prechoice_and_neutral_array{1, i}(idx_for_subsample,:);
    else
        zdataTemp{i} = prechoice_and_neutral_array{1, i};

    end
end

% MY APPROACH
% % use zall array if you want to check how trials compare across block
neuron_mean_concat = horzcat(zdataTemp{:});
% neuron_mean_concat = horzcat(zall_mean_all_array{1, 11}, zall_mean_all_array{1, 12});
%12/14/2024 MAYBE I DONT NEED THIS LINE BELOW? MAYBE IT ARTIFICALLY DEFLATES
%DIFFERENCES?

neuron_mean_concat = zscore(neuron_mean_concat, 0 , 2);

%%


% assuming column 11 are your data if you've loaded the 10x dataset
for zz = 1:size(behav_tbl_iter{11, 1}  , 1)
    current_behav = behav_tbl_iter{11, 1}{zz, 1};

    median_start_time(1, zz) = median(current_behav.stTime - current_behav.choiceTime)
    median_collect_time(1, zz) = median(current_behav.collectionTime - current_behav.choiceTime)
    [~, closest_index_start(zz)] = min(abs(ts1 - median_start_time(1, zz)));
    [~, closest_index_collect(zz)] = min(abs(ts1 - median_collect_time(1, zz)));
    [~, closest_index_zero(zz)] = min(abs(ts1 - 0));
end

% the PCA code assumes you are plotting 2 lines, so we need to get the
% median times for these events, and duplicate it so that it grabs the same
% values for the prechoice remapped line and the true_neutral line
closest_index_start = [median(closest_index_start) median(closest_index_start)]
closest_index_collect = [median(closest_index_collect) median(closest_index_collect)]
closest_index_zero = [median(closest_index_zero) median(closest_index_zero)]


% from here - run PCA_1p_analyses, starting a few lines in