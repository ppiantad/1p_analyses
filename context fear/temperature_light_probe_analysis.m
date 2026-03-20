

metaDirectory = 'd:\MATLAB\my_repo\context fear\temperature & light probe';
metaDirectory_subfolders = dir(metaDirectory );
metafolder_list = {};


% Loop through the list of subfolders
for i = 1:length(metaDirectory_subfolders)
    % Check if the item in subfolders is a directory (not "." or "..") or
    % one of the sets of files that I haven't analyzed yet (PR currently)
    if metaDirectory_subfolders(i).isdir && ~strcmp(metaDirectory_subfolders(i).name, '.') && ~strcmp(metaDirectory_subfolders(i).name, '..') && ~contains(metaDirectory_subfolders(i).name, 'not in final dataset')
  % if metaDirectory_subfolders(i).isdir && ~strcmp(metaDirectory_subfolders(i).name, '.') && ~strcmp(metaDirectory_subfolders(i).name, '..') && ~contains(metaDirectory_subfolders(i).name, 'PR') && ~contains(metaDirectory_subfolders(i).name, 'not in final dataset')
        % if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..') && ~contains(lower(subfolders(i).name), 'shock')
        % if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..')
        % Get the full path of the subfolder
        metasubfolderPath = fullfile(metaDirectory, metaDirectory_subfolders(i).name);

        % Create a cell array for the subfolder path and append it
        % vertically to folder_list
        metafolder_list = vertcat(metafolder_list, {metasubfolderPath});



        % Add your analysis code here

    end
end



folderMask = ~[metaDirectory_subfolders.isdir]; %find all of the folders in the directory and remove them from the list
files = {metaDirectory_subfolders.name}
files = files(folderMask);

for ii = 1:size(files, 2)

    current_file = files{ii};
    file_to_load = readtable(current_file);

    if contains(current_file, 'Monday')
        session_name = 'D1_Conditioning';
        light_temp_struct.(session_name) = file_to_load;
    elseif contains(current_file, 'Tuesday')
        session_name = 'D2_Conditioning';
        light_temp_struct.(session_name) = file_to_load;
    elseif contains(current_file, 'Thursday')
        session_name = 'D3_Discrimination';
        light_temp_struct.(session_name) = file_to_load;
    elseif contains(current_file, 'Friday')
        session_name = 'D4_Discrimination';
        light_temp_struct.(session_name) = file_to_load;
    elseif contains(current_file, 'context_fear_times')
        opts = detectImportOptions(current_file);
        opts.VariableOptions
        opts = setvartype(opts,{'Time', 'Time_1', 'Time_2', 'Time_3'},{'datetime', 'datetime', 'datetime', 'datetime'});
        mouse_session_times = readtimetable(current_file, opts);

    end



end

%%

session_to_analyze = 'D3_Discrimination';


session_start_inds =  find(strcmp(mouse_session_times.Var1, 'Mouse 1 - start'));
session_end_inds =  find(contains(mouse_session_times.Var1, 'Mouse 4'));



start_time = timeofday(mouse_session_times.Time(session_start_inds(1)));
end_time = timeofday(mouse_session_times.Time(session_end_inds(2)))

% Extract the relevant columns
time_column = timeofday(light_temp_struct.(session_to_analyze){:, 2}); % Convert table column to array
temperature_column = light_temp_struct.(session_to_analyze){:, 3}; % Extract temperature data
brightness_column = light_temp_struct.(session_to_analyze){:, 4}; % Extract temperature data

% % Define the time range
% start_datetime = time_column(1) + start_time;
% end_datetime = time_column(1) + end_time;

% Find indices within the time range
valid_indices = (time_column >= start_time) & (time_column <= end_time);

% Extract the corresponding temperature data
filtered_temperature = temperature_column(valid_indices);
filtered_time = time_column(valid_indices); 
filtered_brightness = brightness_column(valid_indices);




%%

max_temp = max(filtered_temperature)
min_temp = min(filtered_temperature)

temp_diff = max_temp-min_temp;

time_elapsed = end_time - start_time;

time_elapsed_seconds = seconds(time_elapsed);

time_elapsed_minutes = minutes(time_elapsed);

temp_change_in_farenheit_per_min = temp_diff/time_elapsed_minutes

approx_temp_change_per_hour = 60*temp_change_in_farenheit_per_min


% temp_change_in_farenheit_per_min = (temp_change_in_celcius_per_min * 9/5);


% Sample data: temperature readings taken once per second
% temperature = randn(1, 600) + linspace(20, 25, 600); % Example data
time_seconds = 1:length(filtered_temperature); % Time vector in seconds

% Calculate the rate of change per second
rate_per_second = diff(filtered_temperature); % First-order difference

% Calculate the rate of change per minute (60-second intervals)
window_size = 60;
rate_per_minute = (filtered_temperature(window_size+1:end) - filtered_temperature(1:end-window_size)) / window_size;

% Time vectors for plotting
time_rate_second = time_seconds(2:end); % Adjusted time for per-second rate
time_rate_minute = time_seconds(window_size+1:end); % Adjusted time for per-minute rate

% Plot the results
figure;
subplot(2,1,1);
plot(time_rate_second, rate_per_second, '-b');
xlabel('Time (seconds)');
ylabel('Rate of Change (°C/s)');
title('Rate of Temperature Change per Second');

subplot(2,1,2);
plot(time_rate_minute, rate_per_minute, '-r');
xlabel('Time (seconds)');
ylabel('Rate of Change (°C/min)');
title('Rate of Temperature Change per Minute');


figure; plot(time_seconds, filtered_temperature)
figure; plot(time_seconds, filtered_brightness)

discrimination_bins = [0 2; 2 4; 4 6; 6 8; 8 10; 10 12]
discrimination_bins_seconds = discrimination_bins*60;

total_bins = size(discrimination_bins, 1); % Get number of bins

binned_brightness = zeros(total_bins, 1); % Preallocate for efficiency

for ff = 1:total_bins
    % Find indices where time falls within the bin
    valid_indices = (time_seconds > discrimination_bins_seconds(ff,1)) & ...
                    (time_seconds < discrimination_bins_seconds(ff,2));

    % Take the mean of the filtered brightness values within the time bin
    binned_brightness(ff) = mean(filtered_brightness(valid_indices));
end

figure; plot(1:6, binned_brightness)

%%

session_to_analyze = 'D3_Discrimination';

if strcmp(session_to_analyze, 'D3_Discrimination')
    timing_info = mouse_session_times(:, 6:7);
    timing_info.Properties.VariableNames = {'ID', 'Time_session'};

elseif strcmp(session_to_analyze, 'D4_Discrimination')
    timing_info = mouse_session_times(:, 9:10);
    timing_info.Properties.VariableNames = {'ID', 'Time_session'};


end


% "pseudo" mice used for the discrete experiment
mouse_number = [1 2 3 4]




for ii = 1:size(mouse_number, 2)
    current_mouse = sprintf('Mouse %d', mouse_number(ii));
    session_start_inds =  find(strcmp(timing_info.ID, strcat(current_mouse, ' - start')));
    session_end_inds =  find(strcmp(timing_info.ID, strcat(current_mouse, ' - end')));


    start_time = timeofday(timing_info{session_start_inds(1), 2}); 


    end_time = timeofday(timing_info{session_end_inds(1), 2}); 

    % Extract the relevant columns
    time_column = timeofday(light_temp_struct.(session_to_analyze){:, 2}); % Convert table column to array
    temperature_column = light_temp_struct.(session_to_analyze){:, 3} % Extract temperature data
    brightness_column = light_temp_struct.(session_to_analyze){:, 4}; % Extract temperature data

    % % Define the time range
    % start_datetime = time_column(1) + start_time;
    % end_datetime = time_column(1) + end_time;

    % Find indices within the time range
    valid_indices = (time_column >= start_time) & (time_column <= end_time);

    % Extract the corresponding temperature data
    filtered_temperature = temperature_column(valid_indices);
    filtered_temp_mouse{ii} = filtered_temperature';
    filtered_time = time_column(valid_indices);
    filtered_brightness = brightness_column(valid_indices);

    time_seconds = 1:length(filtered_temperature); % Time vector in seconds

    discrimination_bins = [0 2; 2 4; 4 6; 6 8; 8 10; 10 12]
    discrimination_bins_seconds = discrimination_bins*60;

    total_bins = size(discrimination_bins, 1); % Get number of bins

    binned_brightness = zeros(total_bins, 1); % Preallocate for efficiency

    for ff = 1:total_bins
        % Find indices where time falls within the bin
        valid_indices = (time_seconds > discrimination_bins_seconds(ff,1)) & ...
            (time_seconds < discrimination_bins_seconds(ff,2));

        % Take the mean of the filtered brightness values within the time bin
        binned_brightness(ff) = mean(filtered_brightness(valid_indices));
        binned_brightness_sem(ff) = std(filtered_brightness(valid_indices))/sqrt(size(filtered_brightness, 1));
    end

   binned_brightness_mouse(:, ii) = binned_brightness;
   binned_brightness_sem_mouse(:, ii) = binned_brightness_sem;


   binned_temp = zeros(total_bins, 1); % Preallocate for efficiency

   for ff = 1:total_bins
       % Find indices where time falls within the bin
       valid_indices = (time_seconds > discrimination_bins_seconds(ff,1)) & ...
           (time_seconds < discrimination_bins_seconds(ff,2));

       % Take the mean of the filtered brightness values within the time bin
       binned_temp(ff) = mean(filtered_temperature(valid_indices));
       binned_temp_sem(ff) = std(filtered_temperature(valid_indices))/sqrt(size(filtered_temperature, 1));
   end

   binned_temp_mouse(:, ii) = binned_temp;
   binned_temp_sem_mouse(:, ii) = binned_temp_sem;

end

% figure;
% shadedErrorBar(1:6, mean(binned_brightness_mouse, 2), mean(binned_brightness_sem_mouse, 2))
% ylim([0 100])
% 
% figure;
% shadedErrorBar(1:6, mean(binned_temp_mouse, 2), mean(binned_temp_sem_mouse, 2))
% ylim([76 78])


figure;
shadedErrorBar(1:6, mean(binned_brightness_mouse, 2), std(binned_brightness_mouse, [], 2))
ylim([0 100])

figure;
shadedErrorBar(1:6, mean(binned_temp_mouse, 2), std(binned_temp_mouse, [], 2))
ylim([74 80])
yticks([74:1:80])
%%
% Assuming filtered_temp_mouse is your 1x4 cell array
% where each cell contains a 1xN array of temperature data

% Initialize cell array to store temperature changes
temp_changes = cell(1, 4);

% Loop through each cell
for i = 1:4
    % Extract the temperature data from current cell
    temp_data = filtered_temp_mouse{i};
    
    % Calculate baseline mean from first 20 columns
    baseline_mean = mean(temp_data(1:20));
    
    % Calculate change from baseline for each column
    % (positive values = increase, negative = decrease)
    temp_changes{i} = temp_data - baseline_mean;
    
    % Optional: Display results for this cell
    fprintf('Cell %d:\n', i);
    fprintf('  Baseline mean (first 20 samples): %.2f\n', baseline_mean);
    fprintf('  Max temperature change: %.2f\n', max(temp_changes{i}));
    fprintf('  Min temperature change: %.2f\n', min(temp_changes{i}));
    fprintf('  Mean temperature change: %.2f\n', mean(temp_changes{i}));
    fprintf('\n');
end

% Optional: Plot the temperature changes for visualization
figure;
for i = 1:4
    subplot(2, 2, i);
    plot(temp_changes{i});
    hold on;
    yline(0, 'k--', 'LineWidth', 1); % Add zero reference line
    xlabel('Time (sec)');
    % ylabel('Temperature change from baseline');
    title(sprintf('Mouse %d: Temperature Change', i));
    grid on;
    ylim([0 1])
    yticks([0 1])
end
% sgtitle('Temperature Changes Relative to Baseline (First 20 Samples)');

%%
% Find the minimum length across all cells
min_length = inf;
for i = 1:4
    min_length = min(min_length, length(temp_changes{i}));
end

% Trim all arrays to the minimum length
temp_changes_trimmed = zeros(4, min_length);  % Create matrix for easier mean/std calculation
for i = 1:4
    temp_changes_trimmed(i, :) = temp_changes{i}(1:min_length);
end

% Calculate mean and standard deviation across mice (rows)
mean_temp_change = mean(temp_changes_trimmed, 1);  % Mean across rows
std_temp_change = std(temp_changes_trimmed, 0, 1);  % Std across rows

% Plot using shadedErrorBar
figure;
shadedErrorBar(1:min_length, mean_temp_change, std_temp_change)
xlabel('Time (samples)');
ylabel('Temperature change from baseline (°C)');
title('Average Temperature Change Across Mice');
grid on;

% Add a reference line at zero
hold on;
yline(0, 'k--', 'LineWidth', 1, 'Alpha', 0.5);

% Optional: Set y-axis limits based on your data range
% You might want to adjust these based on your actual temperature changes
ylim([-2 2]);  % Adjust as needed for your data
yticks([-2:0.5:2]);  % Adjust as needed

% Display info about trimming
fprintf('Arrays trimmed to uniform length: %d samples\n', min_length);
for i = 1:4
    original_length = length(temp_changes{i});
    if original_length > min_length
        fprintf('  Cell %d: trimmed from %d to %d samples (%d samples removed)\n', ...
            i, original_length, min_length, original_length - min_length);
    else
        fprintf('  Cell %d: no trimming needed (%d samples)\n', i, original_length);
    end
end