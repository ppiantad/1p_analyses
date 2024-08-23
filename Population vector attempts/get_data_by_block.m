function [block_1_ca_mouse] = get_data_by_block(animalIDs, session_to_analyze, iter, final, ca_data_type)


%% FILTER TO GET UN-SHUFFLED DATA
iter = iter+1;
neuron_num = 0;
for ii = 1:size(animalIDs,1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;
        % only use rewarded trials for this, otherwise things get wonky
        [BehavData,trials,varargin]=TrialFilter_test(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0); 
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
            session_ca = ca(dd, :);
            [peaks, peak_locs] = findpeaks(session_ca, 'MinPeakDistance',4);
            session_peaks_sum(neuron_num) = sum(peaks);
            session_peaks_sum_mouse{ii}(dd) = sum(peaks);
            session_length(neuron_num) = final.BLA_Insc_24.RDT_D1.time(end)  - final.BLA_Insc_24.RDT_D1.time(1) ;
            session_peaks_per_s_mouse{ii}(dd) = session_peaks_sum_mouse{ii}(dd)/session_length(neuron_num);
            session_peaks_per_min_mouse{ii}(dd) = (session_peaks_sum_mouse{ii}(dd)/session_length(neuron_num))*60;
            block_1_ca = ca(dd, time_array > block_1_mouse(ii, 1) & time_array < block_1_mouse(ii, 2));
            [peaks, peak_locs] = findpeaks(block_1_ca, 'MinPeakDistance',4);
            block_1_peaks_sum(neuron_num) = sum(peaks);
            block_1_peaks_sum_mouse{ii}(dd) = sum(peaks);
            block_1_length(neuron_num) = block_1_mouse(ii, 2)-block_1_mouse(ii, 1);
            block_1_peaks_per_s_mouse{ii}(dd) = block_1_peaks_sum_mouse{ii}(dd)/block_1_length(neuron_num);
            block_1_ca_mouse{ii}(dd,:) = block_1_ca;
            block_2_ca = ca(dd, time_array > block_2_mouse(ii, 1) & time_array < block_2_mouse(ii, 2));
            [peaks, peak_locs] = findpeaks(block_2_ca, 'MinPeakDistance',4);
            block_2_peaks_sum(neuron_num) = sum(peaks);
            block_2_peaks_sum_mouse{ii}(dd) = sum(peaks);
            block_2_length(neuron_num) = block_2_mouse(ii, 2)-block_2_mouse(ii, 1);
            block_2_peaks_per_s_mouse{ii}(dd) = block_2_peaks_sum_mouse{ii}(dd)/block_2_length(neuron_num);
            block_2_ca_mouse{ii}(dd,:) = block_2_ca;
            block_3_ca = ca(dd, time_array > block_3_mouse(ii, 1) & time_array < block_3_mouse(ii, 2));
            [peaks, peak_locs] = findpeaks(block_3_ca, 'MinPeakDistance',4);
            block_3_peaks_sum(neuron_num) = sum(peaks);
            block_3_peaks_sum_mouse{ii}(dd) = sum(peaks);
            block_3_length(neuron_num) = block_3_mouse(ii, 2)-block_3_mouse(ii, 1);
            block_3_peaks_per_s_mouse{ii}(dd) = block_3_peaks_sum_mouse{ii}(dd)/block_3_length(neuron_num);
            block_3_ca_mouse{ii}(dd,:) = block_3_ca;
        end
    end
end


%%
% Assume block_2_ca_mouse{1, 2} is a 112x5212 double array
% data = block_2_ca_mouse{1, 2};
% 
% % Initialize a cell array to store correlation matrices for each column
% correlation_results = cell(1, size(data, 2));
% 
% % Loop over each column
% for col = 1:size(data, 2)
%     % Extract the current column
%     column_data = data(:, col);
% 
%     % Compute the correlation matrix for the current column
%     % Since the column data is a single vector, corr(column_data) will return NaN
%     % Therefore, we need to compute the correlation for each pair manually
%     correlation_matrix = zeros(size(column_data, 1));
%     for i = 1:size(column_data, 1)
%         for j = 1:size(column_data, 1)
%             correlation_matrix(i, j) = corr(column_data(i), column_data(j));
%         end
%     end
% 
%     % Store the correlation matrix
%     correlation_results{col} = correlation_matrix;
% end

% Now correlation_results contains the correlation matrix for each column

