%% This code takes as input filtered, un-normalized data from access_risk_inscopix. For example, large rew choice and small rew choice, with each being associated with 1 column of caTraceTrials_mouse

select_mouse = 'BLA_Insc_26';

select_mouse_index = find(strcmp(animalIDs, select_mouse));

caTraceTrials_mouse_decoding = caTraceTrials_mouse(select_mouse_index,:);

trials_per_mouse_decoding = trials_per_mouse(select_mouse_index,:);


%%


columnArray = [];
columnArray_all = []; 
columnArray_mouse = {}; 
trialArray_mouse = {}; 
trialArray_all = [];
time_seriesArray_all = [];
time_seriesArray_mouse = {};

for i = 1:size(caTraceTrials_mouse_decoding, 1)
    cells = caTraceTrials_mouse_decoding(i,:);
    for z = 1:size(cells, 2)
        nestedCellArray = cells{1, z};
        for k = 1:size(nestedCellArray, 2)
            columnArray = [];
            trialArray = [];
            time_series = [];
            % Extract the double array and transpose it to a column vector
            for qq = 1:size(nestedCellArray{1, k}, 1)
                trialArray = vertcat(trialArray, ones(size(nestedCellArray{1, k}, 2), 1)*qq);
                columnArray = vertcat(columnArray, nestedCellArray{1, k}(qq,:)');  % Assuming the data is in the first column of the nested cell array
                time_series = vertcat(time_series, ones(size(nestedCellArray{1, k}, 2), 1) .* ts1');
                % trial_column(qq) = vertcat(trial_column, trial_num); % flattenedColumnData{k} = vertcat(flattenedColumnData(k), columnArray(:));
            end
            columnArray_mouse(z, k) = {columnArray};
            trialArray_mouse(z, k) = {trialArray};
            time_seriesArray_mouse(z, k) = {time_series};
            clear columnArray trialArray time_series
        end
    end
    columnArray_all = horzcat(columnArray_all, columnArray_mouse);
    trialArray_all = horzcat(trialArray_all, trialArray_mouse); 
    time_seriesArray_all = horzcat(time_seriesArray_all, time_seriesArray_mouse);
    clear columnArray_mouse trialArray_mouse time_seriesArray_mouse
end

% Initialize an empty cell array to store the concatenated columns
concatenatedColumns = cell(1, size(columnArray_all, 2));

% get the events (based on the number of iters - which corresponds to the
% trial filtering done in access_risk_inscopix) into an array

for col = 1:size(columnArray_all, 2)
    for z = 1:iter
        event{z} = ones(size(trialArray_all{z, col}, 1), 1)*z;
    end
    concatenatedEvents{col} = vertcat(event{:});
    
end 

% Loop over columns
for col = 1:size(columnArray_all, 2)
    % Extract the current column
    currentColumn = columnArray_all(:, col);
    currentColumn_trials = trialArray_all(:,col); 
    currentColumn_time = time_seriesArray_all(:,col);
    % event_1 = ones(size(currentColumn{1, 1}, 1), 1);
    % event_2 = ones(size(currentColumn{2, 1}, 1), 1)*2;
    % Vertically concatenate the rows of the current column
    concatenatedColumns{col} = vertcat(currentColumn{:});
    concatenatedColumns_trials{col} = vertcat(currentColumn_trials{:});
    concatenatedColumns_time{col} = vertcat(currentColumn_time{:});
    
end



%% get minimum event numbers across mice

%this needs to keep track of max for each event and choose the LOWEST VALUE
% and keep that conistent across events! edit me!!
% Loop through each column
for col = 1:size(concatenatedColumns_trials, 2)

    for iters = 1:iter
        % Extract cell arrays for the current column
        values = concatenatedColumns_trials{1, col}(concatenatedEvents{1, col} == iters)';
        max_values(iters, col) = max(values);
        values_array{iters} = values; 
        % Calculate the maximum values for each case
        

    end
    % max_values(:, col) = max(values);
    
end

% Find the minimum value present in either array
min_value = min(min(max_values));

figure;

for qq = 1:size(caTraceTrials_mouse_decoding, 1)
    for z = 1:size(cells, 2)
        nestedCellArray = cells{1, z};
        for k = 1:size(nestedCellArray, 2)
            ca_data{z, k} = nestedCellArray{1, k}(1:min_value, :);
            ca_data_means{z}(k,:) = mean(nestedCellArray{1, k}(1:min_value, :));
        end
        event_mean(z,:) = mean(ca_data_means{1, z});
    end
end

figure;
for means = 1:size(event_mean, 1)
    hold on;
    plot(ts1, event_mean(means,:))

end
%% trim arrays to only include data from the mouse with the fewest of a given trial type.
% given that we have so many more large reward trials, this will be limited
% by the small rew trials in nearly every case. this could get a bit
% challenging when conducting things across trial block etc. It could be
% better to randomly select large rew. trials up to the minimum # of small
% reward trials, but for now im just selecting the first bunch of them


% Loop through each column
for col = 1:size(concatenatedColumns_trials, 2)
    % Find indices where concatenatedEvents == 1 and the value is less than or equal to min_value
    indices_to_keep = find(concatenatedColumns_trials{1, col} <= min_value);
    trimmed_concatenatedColumns{1, col}  = concatenatedColumns{1, col}(indices_to_keep);
    trimmed_concatenatedColumns_time{1, col}  = concatenatedColumns_time{1, col}(indices_to_keep);
    trimmed_concatenatedColumns_trials{1, col}  = concatenatedColumns_trials{1, col}(indices_to_keep);
    trimmed_concatenatedEvents{1, col}  = concatenatedEvents{1, col}(indices_to_keep);
    clear indices_to_keep
    % Trim the data based on the indices
    % trimmed_data{1, col} = concatenatedColumns_trials{1, col}(indices_to_keep);
end


%% Get cell arrays that correspond to each time bin

for col = 1:size(concatenatedColumns_time, 2)

    for zz = 1:size(ts1, 2)
        current_step = ts1(zz);
        indices_to_keep = find(trimmed_concatenatedColumns_time{1, col} == current_step);
        trimmed_concatenatedColumns_offsets{zz, col} = trimmed_concatenatedColumns{1, col}(indices_to_keep);
        trimmed_concatenatedColumns_time_offsets{zz, col} = trimmed_concatenatedColumns_time{1, col}(indices_to_keep);
        trimmed_concatenatedColumns_trials_offsets{zz, col} = trimmed_concatenatedColumns_trials{1, col}(indices_to_keep);
        trimmed_concatenatedEvents_offsets{zz, col} = trimmed_concatenatedEvents{1, col}(indices_to_keep);

    end
end


%% Uncomment this section if you do not want to use the "trimmed" array to decode, but beware that depending on the event, the # of trials may differ dramatically!  
% Get cell arrays that correspond to each time bin 

% for col = 1:size(concatenatedColumns_time, 2)
% 
%     for zz = 1:size(ts1, 2)
%         current_step = ts1(zz);
%         indices_to_keep = find(concatenatedColumns_time{1, col} == current_step);
%         trimmed_concatenatedColumns_offsets{zz, col} = concatenatedColumns{1, col}(indices_to_keep);
%         trimmed_concatenatedColumns_time_offsets{zz, col} = concatenatedColumns_time{1, col}(indices_to_keep);
%         trimmed_concatenatedColumns_trials_offsets{zz, col} = concatenatedColumns_trials{1, col}(indices_to_keep);
%         trimmed_concatenatedEvents_offsets{zz, col} = concatenatedEvents{1, col}(indices_to_keep);
% 
%     end
% end