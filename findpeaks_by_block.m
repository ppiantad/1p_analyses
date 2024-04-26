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

load('BLA_panneuronal_Risk_2024_03_07_just_CNMFe_and_BehavData.mat')

% load('NAcSh_D2_Cre-OFF_GCAMP_all.mat')

% load('BLA_panneuronal_matched_Pre_RDT_RM_vs_RDT_D1_01042024.mat')

% load('BLA_panneuronal_Risk_matched_RM_D1_vs_Pre_RDT_RM.mat')

% load('BLA_NAcSh_Risk_matched_Pre_RDT_RM_vs_RDT_D1.mat')

%% Edit these uservariables with what you want to look at
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5]
uv.BLper = [-10 -5];
uv.dt = 0.1; %what is your frame rate
% uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at

ca_data_type = "C"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes

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
            block_1_length(neuron_num) = block_1_mouse(ii, 2)-block_1_mouse(ii, 1);
            block_2_ca = ca(dd, time_array > block_2_mouse(ii, 1) & time_array < block_2_mouse(ii, 2));
            [peaks, peak_locs] = findpeaks(block_2_ca, 'MinPeakDistance',4);
            block_2_peaks_sum(neuron_num) = sum(peaks);
            block_2_length(neuron_num) = block_2_mouse(ii, 2)-block_2_mouse(ii, 1);
            block_3_ca = ca(dd, time_array > block_3_mouse(ii, 1) & time_array < block_3_mouse(ii, 2));
            [peaks, peak_locs] = findpeaks(block_3_ca, 'MinPeakDistance',4);
            block_3_peaks_sum(neuron_num) = sum(peaks);
            block_3_length(neuron_num) = block_3_mouse(ii, 2)-block_3_mouse(ii, 1);
        end
    end
end

block_1_peaks_per_s = block_1_peaks_sum./block_1_length;
block_2_peaks_per_s = block_2_peaks_sum./block_2_length;
block_3_peaks_per_s = block_3_peaks_sum./block_3_length;

mean(block_1_peaks_per_s)
mean(block_2_peaks_per_s)
mean(block_3_peaks_per_s)

block_1_peaks_per_s_nonzero = nonzeros(block_1_peaks_per_s)';
block_2_peaks_per_s_nonzero = nonzeros(block_2_peaks_per_s)';
block_3_peaks_per_s_nonzero = nonzeros(block_3_peaks_per_s)';

mean(block_1_peaks_per_s_nonzero)
mean(block_2_peaks_per_s_nonzero)
mean(block_3_peaks_per_s_nonzero)