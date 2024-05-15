block_1_ca = ca(1, time_array > block_1(1) & time_array < block_1(2));
figure;
% Find peaks
[peaks, peak_locs] = findpeaks(ca(1,:), 'MinPeakDistance',4);

% Plot the peaks
plot(ca(1,:));
hold on;
plot(peak_locs, peaks, 'r*');
hold off;
xlabel('Index');
ylabel('Value');
title('Peaks in block_1_ca');

%%
iter = 0

%%
% load('BLA-NAcShell_Risk_2024_01_04.mat')

% load('BLA_panneuronal_Risk_2024_01_04.mat')

load('BLA_panneuronal_Risk_2024_04_19_just_CNMFe_and_BehavData.mat')

% load('NAcSh_D2_Cre-OFF_GCAMP_all.mat')

% load('BLA_panneuronal_matched_Pre_RDT_RM_vs_RDT_D1_01042024.mat')

% load('BLA_panneuronal_Risk_matched_RM_D1_vs_Pre_RDT_RM.mat')

% load('BLA_NAcSh_Risk_matched_Pre_RDT_RM_vs_RDT_D1.mat')

%% Edit these uservariables with what you want to look at
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
% uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at

ca_data_type = "S"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate

session_to_analyze = 'RDT_D1';
epoc_to_align = 'choiceTime';
ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
animalIDs = (fieldnames(final));
neuron_num = 0;

clear neuron_mean neuron_sem neuron_num zall_mean zall_array zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all 


%% FILTER TO GET UN-SHUFFLED DATA
iter = iter+1;
neuron_num = 0;
for ii = 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
        % only use rewarded trials for this, otherwise things get wonky
        [BehavData,trials,varargin]=TrialFilter(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0); 
        block_1 = [BehavData.stTime(BehavData.Block == 1) BehavData.collectionTime(BehavData.Block == 1)]; 
        block_1_mouse(ii,:) = [block_1(1, 1) block_1(end, 2)];
        block_2 = [BehavData.stTime(BehavData.Block == 2) BehavData.collectionTime(BehavData.Block == 2)]; 
        block_2_mouse(ii,:) = [block_2(1, 1) block_2(end, 2)];
        block_3 = [BehavData.stTime(BehavData.Block == 3) BehavData.collectionTime(BehavData.Block == 3)];
        block_3_mouse(ii,:) = [block_3(1, 1) block_3(end, 2)];
        % [BehavData,trials,varargin]=TrialFilter(BehavData,'REW', 1.2, 'BLOCK', 3);
        % trials = cell2mat(trials);

        % % BehavData = BehavData(BehavData.shockIntensity >= 0.08 & BehavData.shockIntensity <= 0.13, :);
        % % trials = trials(BehavData.shockIntensity >= 0.08 & BehavData.shockIntensity <= 0.13, :);
        %
        % %Create a logical index array based on your conditions
        % logical_index = BehavData.stTime - BehavData.TrialPossible >= 10 & BehavData.stTime - BehavData.TrialPossible <= 50;
        %
        % % Use the logical index array to subset BehavData
        % BehavData = BehavData(logical_index,: );
        % trials = trials(logical_index);


        ca = final.(currentanimal).(session_to_analyze).CNMFe_data.(ca_data_type);
        if strcmp(ca_data_type, 'S')
            ca = full(ca);

        end


        time_array = final.(currentanimal).(session_to_analyze).time;

        for dd = 1:size(ca, 1)
            neuron_num = neuron_num + 1;
            block_1_ca = ca(dd, time_array > block_1_mouse(ii, 1) & time_array < block_1_mouse(ii, 2));
            [peaks, peak_locs] = findpeaks(block_1_ca, 'MinPeakDistance',4);
            block_1_peaks_sum(neuron_num) = sum(peaks);
            block_1_peaks_sum_mouse{ii, iter}(dd) = sum(peaks);
            block_1_length(neuron_num) = block_1_mouse(ii, 2)-block_1_mouse(ii, 1);
            block_1_peaks_per_s_mouse{ii, iter}(dd) = block_1_peaks_sum_mouse{ii, iter}(dd)/block_1_length(neuron_num);
            block_2_ca = ca(dd, time_array > block_2_mouse(ii, 1) & time_array < block_2_mouse(ii, 2));
            [peaks, peak_locs] = findpeaks(block_2_ca, 'MinPeakDistance',4);
            block_2_peaks_sum(neuron_num) = sum(peaks);
            block_2_peaks_sum_mouse{ii, iter}(dd) = sum(peaks);
            block_2_length(neuron_num) = block_2_mouse(ii, 2)-block_2_mouse(ii, 1);
            block_2_peaks_per_s_mouse{ii, iter}(dd) = block_2_peaks_sum_mouse{ii, iter}(dd)/block_2_length(neuron_num);
            block_3_ca = ca(dd, time_array > block_3_mouse(ii, 1) & time_array < block_3_mouse(ii, 2));
            [peaks, peak_locs] = findpeaks(block_3_ca, 'MinPeakDistance',4);
            block_3_peaks_sum(neuron_num) = sum(peaks);
            block_3_peaks_sum_mouse{ii, iter}(dd) = sum(peaks);
            block_3_length(neuron_num) = block_3_mouse(ii, 2)-block_3_mouse(ii, 1);
            block_3_peaks_per_s_mouse{ii, iter}(dd) = block_3_peaks_sum_mouse{ii, iter}(dd)/block_3_length(neuron_num);
        end
    end
end

block_1_peaks_per_s = block_1_peaks_sum./block_1_length;
block_2_peaks_per_s = block_2_peaks_sum./block_2_length;
block_3_peaks_per_s = block_3_peaks_sum./block_3_length;

peaks_per_blocks_for_ANOVA = [block_1_peaks_per_s', block_2_peaks_per_s', block_3_peaks_per_s'];

% [p,t,stats] = anova1(MPG,Origin);

[p,t,stats] = anova1(peaks_per_blocks_for_ANOVA)
[c,m,h,gnames] = multcompare(stats);

block_1_mean_SD = [nanmean(block_1_peaks_per_s); std(block_1_peaks_per_s)/sqrt(length(block_1_peaks_per_s))]
block_2_mean_SD = [nanmean(block_2_peaks_per_s); std(block_2_peaks_per_s)/sqrt(length(block_2_peaks_per_s))]
block_3_mean_SD = [nanmean(block_3_peaks_per_s); std(block_3_peaks_per_s)/sqrt(length(block_3_peaks_per_s))]

means = [block_1_mean_SD(1) block_2_mean_SD(1) block_3_mean_SD(1)];
sems = [block_1_mean_SD(2) block_2_mean_SD(2) block_3_mean_SD(2)];

figure;
bar(means);
hold on;
er = errorbar(1:3, means, sems,'k.', 'LineWidth', 1);    
                         

figure;
h = dabarplot(peaks_per_blocks_for_ANOVA,...
    'scatter',1,'scattersize',15,'scatteralpha',0.5);

figure;
h = daviolinplot(peaks_per_blocks_for_ANOVA,...
    'scatter',1,'scattersize',15,'scatteralpha',0.5);

block_1_peaks_per_s_nonzero = nonzeros(block_1_peaks_per_s)';
block_2_peaks_per_s_nonzero = nonzeros(block_2_peaks_per_s)';
block_3_peaks_per_s_nonzero = nonzeros(block_3_peaks_per_s)';

mean(block_1_peaks_per_s_nonzero)
mean(block_2_peaks_per_s_nonzero)
mean(block_3_peaks_per_s_nonzero)

% non_zero_peaks_per_blocks_for_ANOVA = [block_1_peaks_per_s_nonzero', block_2_peaks_per_s_nonzero', block_3_peaks_per_s_nonzero'];

for qq = 1:size(block_1_peaks_per_s_mouse, 1)
    block_1_peaks_per_s_mouse_mean(qq) = mean(block_1_peaks_per_s_mouse{qq, 1});
    block_2_peaks_per_s_mouse_mean(qq) = mean(block_2_peaks_per_s_mouse{qq, 1});
    block_3_peaks_per_s_mouse_mean(qq) = mean(block_3_peaks_per_s_mouse{qq, 1});
end

block_1_to_3_diff = block_1_peaks_per_s_mouse_mean - block_3_peaks_per_s_mouse_mean;
figure; scatter(riskiness, block_1_to_3_diff')


%%
% Calculate the number of neurons
num_neurons = size(block_1_peaks_per_s, 2);

% Initialize counters for each category
increase_count = 0;
decrease_count = 0;
mixed_count = 0;
increase_then_decrease_count = 0;
decrease_then_increase_count = 0;

% Loop through each neuron
for neuron = 1:num_neurons
    % Extract data for the neuron from each block
    neuron_data_block1 = block_1_peaks_per_s(:, neuron);
    neuron_data_block2 = block_2_peaks_per_s(:, neuron);
    neuron_data_block3 = block_3_peaks_per_s(:, neuron);
    
    % Check if values increase across the blocks
    if all(diff([neuron_data_block1, neuron_data_block2, neuron_data_block3]) > 0)
        increase_count = increase_count + 1;
    % Check if values decrease across the blocks
    elseif all(diff([neuron_data_block1, neuron_data_block2, neuron_data_block3]) < 0)
        decrease_count = decrease_count + 1;
    % Otherwise, it's mixed
    elseif all(diff([neuron_data_block1, neuron_data_block2]) > 0 & diff([neuron_data_block2, neuron_data_block3]) < 0)
        increase_then_decrease_count = increase_then_decrease_count + 1;
    elseif all(diff([neuron_data_block1, neuron_data_block2]) < 0 & diff([neuron_data_block2, neuron_data_block3]) > 0)
        decrease_then_increase_count = decrease_then_increase_count + 1;
    end
end

% Calculate proportions
total_neurons = num_neurons;
proportion_increase = increase_count / total_neurons;
proportion_decrease = decrease_count / total_neurons;
% proportion_mixed = mixed_count / total_neurons;
proportion_increase_then_decrease = increase_then_decrease_count / total_neurons;
proportion_decrease_then_increase = decrease_then_increase_count / total_neurons;

% Display results
disp(['Proportion of neurons with increasing trend: ', num2str(proportion_increase)]);
disp(['Proportion of neurons with decreasing trend: ', num2str(proportion_decrease)]);
% disp(['Proportion of neurons with mixed trend: ', num2str(proportion_mixed)]);
disp(['Proportion of neurons with increasing trend: ', num2str(proportion_increase_then_decrease)]);
disp(['Proportion of neurons with decreasing trend: ', num2str(proportion_decrease_then_increase)]);

% Create a pie chart
categories = {'Increasing', 'Decreasing', 'Increase then Decrease', 'Decrease then Increase'};
proportions = [proportion_increase, proportion_decrease, proportion_increase_then_decrease, proportion_decrease_then_increase];
% explode = [0.1, 0.1, 0]; % Explode the "Increasing" slice for better visibility

figure;
pie(proportions, categories);
title('Proportion of Neurons with Different Trends');