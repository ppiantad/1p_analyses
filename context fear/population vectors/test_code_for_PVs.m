

animalIDs = (fieldnames(final));

session_to_analyze = 'D3';

for ii = 1:size(fieldnames(final),1)


    select_mouse = animalIDs{ii};
    select_mouse_index = find(strcmp(animalIDs, select_mouse));

    ca = final.(select_mouse).(session_to_analyze).CNMFe_data.(ca_data_type);



    % ca = final.B51618.D3.CNMFe_data.S;
    % ca = ca(pre_choice_cells_mouse == 1, :);
    % filter based on variable above
    % ca = ca(selected_rows, :);
    % ca = ca(prechoice_conserved_mouse{gg, 1} ~= 1 & prechoice_lost_mouse{gg, 1} ~= 1 & prechoice_remapped_mouse{gg, 1} ~= 1, :);
    % ca_zscored = zscore(ca, [], 2);
    ca_zscored = full(ca);

    % Bin size (number of samples per bin)
    bin_size = 3;

    % Number of bins
    n_bins = floor(size(ca_zscored, 2) / bin_size);

    % Trim data to fit evenly into bins
    ca_trimmed = ca_zscored(:, 1:n_bins*bin_size);

    % Reshape and average within bins for each neuron
    ca_binned = reshape(ca_trimmed, size(ca_zscored,1), bin_size, n_bins);
    ca_binned = squeeze(mean(ca_binned, 2));


    % ca_zscored = ca(prechoice_lost_mouse{gg, 1} == 1, :);
    % time_array = final.(select_mouse).(first_session).time;


    spikes_to_use = calcium_events.(select_mouse).stimulus_type_2;  

    % calcium_spikes_mean = mean(calcium_events.(select_mouse).stimulus_type_2, 2);
    calcium_spikes_mean = calcium_events.(select_mouse).mean_events_stim2;  

    similarityOverTime = [];

    

    for t = 1:size(ca_binned, 2)
        % for i = 1:length(uniqueTypes)
        % typeIdx = (neuronTypes == uniqueTypes(i));
        % use specific subset of neurons
        activitySubset = ca_binned(:, t);
        % uncomment if you want to use all neurons
        %activitySubset = neuralActivity(typeIdx, t:(t + windowSize - 1));
        similarityMatrix = corrcoef(activitySubset, calcium_spikes_mean);
        % similarityMatrix = corrcoef(activitySubset, spikes_to_use(:, 1));
        similarityOverTime(t) = similarityMatrix(2);
        % end
    end
    similarityOverTime_array{ii} = similarityOverTime;
    
end

%%
% Ensure all arrays are exactly 2000 columns
for i = 1:length(similarityOverTime_array)
    if size(similarityOverTime_array{i}, 2) > 2000
        similarityOverTime_array{i} = similarityOverTime_array{i}(:, 1:2000);
    end
end

% Get unique treatment groups
unique_treatments = unique(current_animal_treatment);

% Initialize structure to store results
treatment_means = struct();

% Calculate nanmean for each treatment group
for t = 1:length(unique_treatments)
    treatment_name = unique_treatments{t};
    
    % Find indices for this treatment
    treatment_idx = strcmp(current_animal_treatment, treatment_name);
    
    % Get data for this treatment
    treatment_data = similarityOverTime_array(treatment_idx);
    
    % Stack all arrays for this treatment (each row is one animal)
    stacked_data = vertcat(treatment_data{:});
    
    % Calculate nanmean across animals (column-wise)
    treatment_means.(matlab.lang.makeValidName(treatment_name)) = nanmean(stacked_data, 1);
end

% Access results like:
% treatment_means.Experimental
% treatment_means.NoShock
% treatment_means.OneContext