%% This code takes as input filtered, un-normalized data from access_risk_inscopix. For example, large rew choice and small rew choice, with each being associated with 1 column of caTraceTrials_mouse



columnArray = [];
columnArray_all = []; 
columnArray_mouse = {}; 
trialArray_mouse = {}; 
trialArray_all = [];
time_seriesArray_all = [];
time_seriesArray_mouse = {};

for i = 1:size(caTraceTrials_mouse, 1)
    cells = caTraceTrials_mouse(i,:);
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

% Loop over columns
for col = 1:size(columnArray_all, 2)
    % Extract the current column
    currentColumn = columnArray_all(:, col);
    currentColumn_trials = trialArray_all(:,col); 
    currentColumn_time = time_seriesArray_all(:,col);
    event_1 = ones(size(currentColumn{1, 1}, 1), 1);
    event_2 = ones(size(currentColumn{2, 1}, 1), 1)*2;
    % Vertically concatenate the rows of the current column
    concatenatedColumns{col} = vertcat(currentColumn{:});
    concatenatedColumns_trials{col} = vertcat(currentColumn_trials{:});
    concatenatedColumns_time{col} = vertcat(currentColumn_time{:});
    concatenatedEvents{col} = vertcat(event_1, event_2);
end



%% get minimum event numbers across mice

% Assuming concatenatedColumns_trials and concatenatedEvents are already defined

% Initialize arrays to store maximum values for each case
max_values_event1 = zeros(1, size(concatenatedColumns_trials, 2));
max_values_event2 = zeros(1, size(concatenatedColumns_trials, 2));

% Loop through each column
for col = 1:size(concatenatedColumns_trials, 2)
    % Extract cell arrays for the current column
    values_event1 = concatenatedColumns_trials{1, col}(concatenatedEvents{1, col} == 1);
    values_event2 = concatenatedColumns_trials{1, col}(concatenatedEvents{1, col} == 2);
    
    % Calculate the maximum values for each case
    max_values_event1(col) = max(values_event1);
    max_values_event2(col) = max(values_event2);
end

% Find the minimum value present in either array
min_value = min(min([max_values_event1; max_values_event2]));

%% trim arrays to only include data from the mouse with the fewest of a given trial type.
% given that we have so many more large reward trials, this will be limited
% by the small rew trials in nearly every case. this could get a bit
% challenging when conducting things across trial block etc. It could be
% better to randomly select large rew. trials up to the minimum # of small
% reward trials, but for now im just selecting the first bunch of them


% Initialize a new cell array to store the trimmed data
trimmed_data = cell(size(concatenatedColumns_trials));

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