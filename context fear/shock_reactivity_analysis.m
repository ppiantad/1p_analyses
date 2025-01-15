%%
uv.evtWin = [-2 4]; %what time do y
frame_rate = 30; % Frames per second
% Parameters
shock_start_time = 4 * 60; % First shock in seconds
shock_interval = 60; % Interval between shocks in seconds
shock_duration = 2; % Duration of each shock in seconds
num_shocks = 6;

% Initialize variables
shk_on = zeros(1, num_shocks);
shk_off = zeros(1, num_shocks);

% Calculate on and off times for each shock
for i = 0:(num_shocks-1)
    shk_on(i+1) = shock_start_time + i * shock_interval; % Shock start time in seconds
    shk_off(i+1) = shk_on(i+1) + shock_duration; % Shock end time in seconds
end

ts1 = (uv.evtWin(1):1/frame_rate:uv.evtWin(2)-1/frame_rate);
animalIDs = (fieldnames(final_DLC));


sessions_to_analyze = {'D1_Afternoon', 'D2_Afternoon'};
for gg = 1:size(sessions_to_analyze, 2)
    session_to_analyze = sessions_to_analyze{gg};
    mouse_count = 0;
    for ii = 1:size(animalIDs,1)
        currentanimal = char(animalIDs(ii));
        if isfield(final_DLC.(currentanimal), session_to_analyze)
            mouse_count = mouse_count+1;

            % Find the row index where the 'mouse' column matches 'current_mouse'
            row_idx = strcmp(experimental_grps.mouse, currentanimal);
            experimental_grps_updated(mouse_count, :) = experimental_grps(row_idx, :);

            body_velocity = [];
            labels = [];
            % Get the body_velocity column
            body_velocity = final_DLC.(currentanimal).(session_to_analyze).movement_data.body_velocity;

            eTS = shk_on'; %get time stamps

            %calculate time windows for each event
            evtWinSpan = max(uv.evtWin) - min(uv.evtWin);
            % Define parameters


            % Calculate the total number of samples from body_velocity
            num_samples = height(body_velocity); % Assuming body_velocity is a table column

            % Generate time array in seconds
            time_array = (0:num_samples-1) / frame_rate;

            % Convert time array to minutes if needed
            time_array_minutes = time_array / 60;
            for t = 1:size(eTS,1)
                % set each trial's temporal boundaries
                timeWin = [eTS(t)+uv.evtWin(1,1):1/frame_rate:eTS(t)+uv.evtWin(1,2)];  %calculate time window around each event
                % BL_win = [eTS(t)+uv.BLper(1,1):uv.dt:eTS(t)+uv.BLper(1,2)];
                % unitTrace_zscored = zscore(unitTrace);


                if min(timeWin) > min(time_array) && max(timeWin) < max(time_array)    %if the beginning and end of the time window around the event occurred during the recording period. if not, the time window is out of range %if min(timeWin) > min(caTime) & max(timeWin) < max(caTime)
                    % get unit event counts in trials
                    % get unit ca traces in trials
                    idx = time_array >= min(timeWin) & time_array < max(timeWin);      %logical index of time window around each behavioral event time  %idx = caTime > min(timeWin) & caTime < max(timeWin);


                    velocity_for_shocks(t,:) = body_velocity(idx);


                end

            end
            velocity_for_shocks_mouse{gg}{ii} = velocity_for_shocks;
            mean_velocity_for_shocks{gg}(ii, :) = mean(velocity_for_shocks);
            sem_velocity_for_shocks{gg}(ii, :) = std(velocity_for_shocks)/sqrt(size(velocity_for_shocks, 1));
            clear velocity_for_shocks

        end

    end
end

for ff = 1:size(velocity_for_shocks_mouse, 2)
    current_velocity_array = velocity_for_shocks_mouse{1, ff}
        for hh = 1:size(current_velocity_array, 2)
            mouse_velocity_array = current_velocity_array{1, hh};
            if isempty(mouse_velocity_array)
                shock_epoch_means_mouse{ff}(:, hh) = NaN;
            else
                shock_epoch_means_mouse{ff}(:, hh) = mean(mouse_velocity_array(:, ts1 >= 0 & ts1 < 1), 2);
            end

        end


end



experimental_data_shks_D1 = shock_epoch_means_mouse{1, 1}(:, strcmp(experimental_grps_updated.group, 'Experimental'));
experimental_data_shks_D2 = shock_epoch_means_mouse{1, 2}(:, strcmp(experimental_grps_updated.group, 'Experimental'));

one_context_data_shks_D1 = shock_epoch_means_mouse{1, 1}(:, strcmp(experimental_grps_updated.group, 'One Context'));
one_context_data_shks_D2 = shock_epoch_means_mouse{1, 2}(:, strcmp(experimental_grps_updated.group, 'One Context'));

no_shock_data_shks_D1 = shock_epoch_means_mouse{1, 1}(:, strcmp(experimental_grps_updated.group, 'No Shock'));
no_shock_data_shks_D2 = shock_epoch_means_mouse{1, 2}(:, strcmp(experimental_grps_updated.group, 'No Shock'));


mean_experimental_data_shks_D1 = nanmean(experimental_data_shks_D1, 1)
mean_experimental_data_shks_D2 = nanmean(experimental_data_shks_D2, 1)

[h p] = ttest(mean_experimental_data_shks_D1, mean_experimental_data_shks_D2);

mean_one_context_data_shks_D1 = nanmean(one_context_data_shks_D1, 1)
mean_one_context_data_shks_D2 = nanmean(one_context_data_shks_D2, 1)

[h p] = ttest(mean_one_context_data_shks_D1, mean_one_context_data_shks_D2);

mean_no_shock_data_shks_D1 = nanmean(no_shock_data_shks_D1, 1)
mean_no_shock_data_shks_D2 = nanmean(no_shock_data_shks_D2, 1)

[h p] = ttest(mean_no_shock_data_shks_D1, mean_no_shock_data_shks_D2);