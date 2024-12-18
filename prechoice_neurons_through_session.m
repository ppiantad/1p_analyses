

mouse_to_check = 9; 

prechoice_block_2_3_column_1_data = zall_mouse{mouse_to_check, 1}(prechoice_remapped_mouse{mouse_to_check, 1} == 1);
prechoice_block_2_3_column_2_data = zall_mouse{mouse_to_check, 5} (prechoice_remapped_mouse{mouse_to_check, 1} == 1);

behav_part_1 = behav_tbl_iter{1, 1}{mouse_to_check, 1}  
behav_part_2 = behav_tbl_iter{5, 1}{mouse_to_check, 1} 
concat_behav = vertcat(behav_part_1, behav_part_2 )

% Initialize the concatenated cell array
concatenated_columns = cell(1, size(prechoice_block_2_3_column_1_data, 2));


% Iterate through each column and concatenate the data
for i = 1:size(prechoice_block_2_3_column_1_data, 2)
    concatenated_columns{i} = vertcat(prechoice_block_2_3_column_1_data{i}, ...
                                      prechoice_block_2_3_column_2_data{i});
end

% Initialize an array to store the resulting mean values
[numRows, numCols] = size(concatenated_columns{1}); % Assume all cells have the same dimensions
mean_array = zeros(numRows, numCols); % Resulting double array

% Loop through each cell array
for i = 1:length(concatenated_columns)
    % Add the values from the current cell array to the mean calculation
    mean_array = mean_array + concatenated_columns{i};
end

% Divide by the number of cell arrays to calculate the mean
mean_array = mean_array / length(concatenated_columns);

figure; imagesc(ts1, [], mean_array)

% figure; plot(concat_behav.bigSmall)

%%
if size(respClass_all_array, 2) == 10
    comparison_arrays_full = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays_full = [1 2 3; 4 5 6]
elseif size(respClass_all_array, 2) == 7
    comparison_arrays_full = [1 2 3; 5 6 7]
end

% ca_data_type = uv.ca_data_type

% Assume 'ca' is your calcium imaging data, sampled every 100 ms
% bin_size is the desired bin size in ms
bin_size = 100;  % Set this as needed
bin_factor = bin_size / (uv.dt*1000);  % Determine how many samples to bin together

first_session = session_to_analyze;

ts1 = (uv.evtWin(1):bin_factor/10:uv.evtWin(2)-0.1);
%%
% for aa = 1:size(comparison_arrays_full, 1) %size(comparison_arrays_full, 1)
% 
%     comparison_arrays = comparison_arrays_full(aa, :)
%     first_session = session_to_analyze;

for gg = 1:size(animalIDs, 1) %size(animalIDs, 1)
    % if gg ~= 2 && gg ~= 3 && gg ~= 6
    select_mouse = animalIDs{gg};

    select_mouse_index = find(strcmp(animalIDs, select_mouse));
    BehavData = final_behavior.(select_mouse).(first_session).uv.BehavData;


    ca = final.(select_mouse).(first_session).CNMFe_data.(ca_data_type);
    % ca = ca(prechoice_indices_for_PV{aa, select_mouse_index} == 1, :);
    % ca_zscored = zscore(ca, [], 2);
    ca_zscored = normalize(ca, 2);
    time_array = final.(select_mouse).(first_session).time;
    % time_array_all{aa, gg} = time_array;
    % Reshape the data and take the mean of each bin
    [n_neurons, n_samples] = size(ca_zscored);
    ca_binned = squeeze(mean(reshape(ca_zscored(:, 1:floor(n_samples/bin_factor)*bin_factor), n_neurons, bin_factor, []), 2));

    % Squeeze the data to remove the singleton dimension
    % ca_binned = squeeze(ca_binned);
    ca = ca_binned;
    % Assuming time_array is a 1D array of time points corresponding to each sample
    time_array_binned = mean(reshape(time_array(1:floor(length(time_array)/bin_factor)*bin_factor), bin_factor, []), 1);
    time_array = time_array_binned;





    %%
    PV_prechoice_all_mouse = [];
    prechoice_block_1_column_1_data = {};
    prechoice_block_1_column_2_data  = {};
    prechoice_mean_mouse = [];
    prechoice_block_1_column_1_data = zall_mouse{select_mouse_index, 1}(prechoice_block_1_mouse{select_mouse_index, 1} == 1);
    prechoice_block_1_column_2_data = zall_mouse{select_mouse_index, 5}(prechoice_block_1_mouse{select_mouse_index, 1} == 1);
    concatenated_columns = {};
    for i = 1:size(prechoice_block_1_column_1_data, 2)
        concatenated_columns{i} = vertcat(prechoice_block_1_column_1_data{i}, ...
            prechoice_block_1_column_2_data{i});
    end
    for cc = 1:size(concatenated_columns, 2)
        prechoice_data = concatenated_columns{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);


    end
    prechoice_block_1_mean_mouse_array{gg} = prechoice_mean_mouse;

        %%
    PV_prechoice_all_mouse = [];
    prechoice_block_2_3_column_1_data = {};
    prechoice_block_2_3_column_2_data  = {};
    prechoice_mean_mouse = [];
    prechoice_block_2_3_column_1_data = zall_mouse{select_mouse_index, 1}(prechoice_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);
    prechoice_block_2_3_column_2_data = zall_mouse{select_mouse_index, 5}(prechoice_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);
    concatenated_columns = {};
    for i = 1:size(prechoice_block_2_3_column_1_data, 2)
        concatenated_columns{i} = vertcat(prechoice_block_2_3_column_1_data{i}, ...
            prechoice_block_2_3_column_2_data{i});
    end
    for cc = 1:size(concatenated_columns, 2)
        prechoice_data = concatenated_columns{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);


    end
    prechoice_block_2_3_mean_mouse_array{gg} = prechoice_mean_mouse;

    PV_prechoice_all_mouse = [];
    prechoice_block_all_column_1_data = {};
    prechoice_block_all_column_2_data  = {};
    prechoice_mean_mouse = [];
    prechoice_block_all_column_1_data = zall_mouse{select_mouse_index, 1};
    prechoice_block_all_column_2_data = zall_mouse{select_mouse_index, 5};
    concatenated_columns = {};
    for i = 1:size(prechoice_block_all_column_1_data, 2)
        concatenated_columns{i} = vertcat(prechoice_block_all_column_1_data{i}, ...
            prechoice_block_all_column_2_data{i});
    end
    for cc = 1:size(concatenated_columns, 2)
        prechoice_data = concatenated_columns{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);


    end
    prechoice_all_mean_mouse_array{gg} = prechoice_mean_mouse;


end



%% TO RUN CODE BELOW, load 10x variables then run eventRelatedActivityAndClassification_PTP_v4:
% 1. choiceTime 'OMITALL', 0, 'BLANK_TOUCH', 0
% 2. collectionTime'OMITALL', 0, 'BLANK_TOUCH', 0
%%
if size(respClass_all_array, 2) == 10
    comparison_arrays_full = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays_full = [1 2 3; 4 5 6]
elseif size(respClass_all_array, 2) == 7
    comparison_arrays_full = [1 2 3; 5 6 7]
end

% ca_data_type = uv.ca_data_type

% Assume 'ca' is your calcium imaging data, sampled every 100 ms
% bin_size is the desired bin size in ms
bin_size = 100;  % Set this as needed
bin_factor = bin_size / (uv.dt*1000);  % Determine how many samples to bin together

first_session = session_to_analyze;

ts1 = (uv.evtWin(1):bin_factor/10:uv.evtWin(2)-0.1);
%%
% for aa = 1:size(comparison_arrays_full, 1) %size(comparison_arrays_full, 1)
% 
%     comparison_arrays = comparison_arrays_full(aa, :)
%     first_session = session_to_analyze;

for gg = 1:size(animalIDs, 1) %size(animalIDs, 1)
    % if gg ~= 2 && gg ~= 3 && gg ~= 6
    select_mouse = animalIDs{gg};

    select_mouse_index = find(strcmp(animalIDs, select_mouse));
    BehavData = final_behavior.(select_mouse).(first_session).uv.BehavData;


    ca = final.(select_mouse).(first_session).CNMFe_data.(ca_data_type);
    % ca = ca(prechoice_indices_for_PV{aa, select_mouse_index} == 1, :);
    % ca_zscored = zscore(ca, [], 2);
    ca_zscored = normalize(ca, 2);
    time_array = final.(select_mouse).(first_session).time;
    % time_array_all{aa, gg} = time_array;
    % Reshape the data and take the mean of each bin
    [n_neurons, n_samples] = size(ca_zscored);
    ca_binned = squeeze(mean(reshape(ca_zscored(:, 1:floor(n_samples/bin_factor)*bin_factor), n_neurons, bin_factor, []), 2));

    % Squeeze the data to remove the singleton dimension
    % ca_binned = squeeze(ca_binned);
    ca = ca_binned;
    % Assuming time_array is a 1D array of time points corresponding to each sample
    time_array_binned = mean(reshape(time_array(1:floor(length(time_array)/bin_factor)*bin_factor), bin_factor, []), 1);
    time_array = time_array_binned;





    %%
    prechoice_block_1_column_1_data = {};
    prechoice_mean_mouse = [];
    prechoice_block_1_column_1_data = zall_mouse{select_mouse_index, 11}(prechoice_block_1_mouse{select_mouse_index, 1} == 1);
    for cc = 1:size(prechoice_block_1_column_1_data, 2)
        prechoice_data = prechoice_block_1_column_1_data{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    end
    prechoice_block_1_mean_mouse_array{gg} = prechoice_mean_mouse;

    prechoice_block_2_3_column_1_data = {};
    prechoice_mean_mouse = [];
    prechoice_block_2_3_column_1_data = zall_mouse{select_mouse_index, 11}(prechoice_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);

    for cc = 1:size(prechoice_block_2_3_column_1_data, 2)
        prechoice_data = prechoice_block_2_3_column_1_data{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    end
    prechoice_block_2_3_mean_mouse_array{gg} = prechoice_mean_mouse;

    prechoice_block_all_column_1_data = {};
    prechoice_mean_mouse = [];
    prechoice_block_all_column_1_data = zall_mouse{select_mouse_index, 11};

    for cc = 1:size(prechoice_block_all_column_1_data, 2)
        prechoice_data = prechoice_block_all_column_1_data{1, cc};
        prechoice_mean_mouse(:, cc) = mean(prechoice_data (:, ts1 >= -4 & ts1 <= 0), 2);
    end
    prechoice_all_mean_mouse_array{gg} = prechoice_mean_mouse;
    %%
    postchoice_block_1_column_1_data = {};
    postchoice_mean_mouse = [];
    postchoice_block_1_column_1_data = zall_mouse{select_mouse_index, 11}(postchoice_reward_block_1_mouse{select_mouse_index, 1} == 1);
    for cc = 1:size(postchoice_block_1_column_1_data, 2)
        postchoice_data = postchoice_block_1_column_1_data{1, cc};
        postchoice_mean_mouse(:, cc) = mean(postchoice_data (:, ts1 >= 0 & ts1 <= 2), 2);
    end
    postchoice_block_1_mean_mouse_array{gg} = postchoice_mean_mouse;

    postchoice_block_2_3_column_1_data = {};
    postchoice_mean_mouse = [];
    postchoice_block_2_3_column_1_data = zall_mouse{select_mouse_index, 11}(postchoice_reward_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);

    for cc = 1:size(postchoice_block_2_3_column_1_data, 2)
        postchoice_data = postchoice_block_2_3_column_1_data{1, cc};
        postchoice_mean_mouse(:, cc) = mean(postchoice_data (:, ts1 >= 0 & ts1 <= 2), 2);
    end
    postchoice_block_2_3_mean_mouse_array{gg} = postchoice_mean_mouse;

    postchoice_block_all_column_1_data = {};
    postchoice_mean_mouse = [];
    postchoice_block_all_column_1_data = zall_mouse{select_mouse_index, 11};

    for cc = 1:size(postchoice_block_all_column_1_data, 2)
        postchoice_data = postchoice_block_all_column_1_data{1, cc};
        postchoice_mean_mouse(:, cc) = mean(postchoice_data (:, ts1 >= 0 & ts1 <= 2), 2);
    end
    postchoice_all_mean_mouse_array{gg} = postchoice_mean_mouse;

    %%
    collect_block_1_column_1_data = {};
    collect_mean_mouse = [];
    collect_block_1_column_1_data = zall_mouse{select_mouse_index, 12}(collect_block_1_mouse{select_mouse_index, 1} == 1);
    for cc = 1:size(collect_block_1_column_1_data, 2)
        collect_data = collect_block_1_column_1_data{1, cc};
        collect_mean_mouse(:, cc) = mean(collect_data (:, ts1 >= 1 & ts1 <= 3), 2);
    end
    collect_block_1_mean_mouse_array{gg} = collect_mean_mouse;

    collect_block_2_3_column_1_data = {};
    collect_mean_mouse = [];
    collect_block_2_3_column_1_data = zall_mouse{select_mouse_index, 12}(collect_blocks_2_and_3_mouse{select_mouse_index, 1} == 1);

    for cc = 1:size(collect_block_2_3_column_1_data, 2)
        collect_data = collect_block_2_3_column_1_data{1, cc};
        collect_mean_mouse(:, cc) = mean(collect_data (:, ts1 >= 1 & ts1 <= 3), 2);
    end
    collect_block_2_3_mean_mouse_array{gg} = collect_mean_mouse;

    collecte_block_all_column_1_data = {};
    collect_mean_mouse = [];
    collect_block_all_column_1_data = zall_mouse{select_mouse_index, 12};

    for cc = 1:size(collect_block_all_column_1_data, 2)
        collect_data = collect_block_all_column_1_data{1, cc};
        collect_mean_mouse(:, cc) = mean(collect_data (:, ts1 >= 1 & ts1 <= 3), 2);
    end
    collect_all_mean_mouse_array{gg} = collect_mean_mouse;

end
% end

%%
for hh = 1:size(prechoice_all_mean_mouse_array, 2)

    extracted = prechoice_all_mean_mouse_array{hh};
    session_long_prechoice(:, hh) = mean(extracted, 2);

end

mean_session_long_prechoice = mean(session_long_prechoice, 2)



for hh = 1:size(prechoice_block_1_mean_mouse_array, 2)

    extracted = prechoice_block_1_mean_mouse_array{hh};
    session_long_prechoice_block_1(:, hh) = mean(extracted, 2);

end

mean_session_long_prechoice = mean(session_long_prechoice_block_1, 2)

for hh = 1:size(prechoice_block_2_3_mean_mouse_array, 2)

    extracted = prechoice_block_2_3_mean_mouse_array{hh};
    session_long_prechoice_block_2_3(:, hh) = mean(extracted, 2);

end

mean_session_long_prechoice = mean(session_long_prechoice_block_2_3, 2)

%%
% Create the figure
figure;

% Plot the first dataset on the left Y-axis
yyaxis left;
plot(mean_session_long_prechoice, '-b', 'LineWidth', 1.5); % Blue line for the first dataset
ylabel('Pre-choice ensemble PV'); % Label for the left Y-axis
xlabel('Index'); % Common X-axis label

% Plot the second dataset on the right Y-axis
yyaxis right;
plot(meanData_large, '-r', 'LineWidth', 1.5); % Red line for the second dataset
ylabel('% large rew choice'); % Label for the right Y-axis





