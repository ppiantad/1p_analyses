
fs_cam = 30; %set sampling rate according to camera, this is hard coded for now

animalIDs = (fieldnames(final_SLEAP));
session_to_analyze = 'RDT_D1';
% IDs_from_list = hM4Di_IDs;
% IDs_from_list = stGtACR_IDs;
% IDs_from_list = PdCO_IDs;
IDs_from_list = ChrimsonR_IDs;

reward_cup_time = [];
right_screen_time = [];
left_screen_time = [];


animals_with_sessions = {};

for dd = 1:size(animalIDs)

    select_mouse = animalIDs{dd};
    if isfield(final_SLEAP.(select_mouse), session_to_analyze)
        [~, loc] = ismember(select_mouse, IDs_from_list);
        shapeData = final_SLEAP.(select_mouse).(session_to_analyze).shapeData;
        animals_with_sessions{dd} = select_mouse;
        large_rew_side = large_screen_side{loc};
        SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;
        % X_data = SLEAP_data.corrected_x_pix;
        % Y_data = SLEAP_data.corrected_y_pix;
        X_data = SLEAP_data.x_pix;
        Y_data = SLEAP_data.y_pix;
        % Apply Savitzky-Golay filter to each row. this mostly removes high
        % frequency changes that occur due to slight jumps in keypoint location
        X_data = sgolayfilt(X_data, 9, 33);
        Y_data = sgolayfilt(Y_data, 9, 33);
        BehavData = final_behavior.(select_mouse).(session_to_analyze).uv.BehavData;
        onset_trials = BehavData.stTime';
        choice_trials = BehavData.choiceTime';
        offset_trials = BehavData.collectionTime';

        time_ranges_trials = [onset_trials; choice_trials; offset_trials];

        SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;

        % velocity_data = final_SLEAP.(select_mouse).(session_to_analyze).zscored_SLEAP_data_velocity';

        velocity_data = zscore(SLEAP_data.vel_cm_s)';

        
        % adjusted_start_time = BehavData.TrialPossible(1)-60;
        % SLEAP_data.idx_time = SLEAP_data.idx_time+adjusted_start_time;

        %%
        % FILTER ALL EXISTING DATA ON THESE TIME RANGES
        % filter streams
        if ~isempty(SLEAP_data)
            filtered_motion = [];
            max_ind = SLEAP_data.idx_frame(end);
            idx_frame_redo = 1:1:size(SLEAP_data, 1);
            good_index = 1;
            for j = 1:size(time_ranges_trials,2)
                onset = round(time_ranges_trials(1,j)*fs_cam)+1;
                choice = round(time_ranges_trials(2,j)*fs_cam)+1;
                offset = round(time_ranges_trials(3,j)*fs_cam)+1;

                % throw it away if onset or offset extends beyond recording window
                if isinf(offset)
                    if onset <= max_ind && onset > 0
                        filtered_motion{good_index} = SLEAP_data(onset:end);
                        break %return
                    end
                else
                    % if offset <= max_ind && offset > 0 && onset <= max_ind && onset > 0
                        % buffering this by adding +1 to the end time for now, for
                        % some reason the array seems too short without?
                        % after some extensive checking, it seems like the
                        % strangeness where the body @ start and @ end does not
                        % overlap often comes from the fact that the mouse's tail
                        % can trigger the IR beam in the food cup on a non-trivial
                        % # of trials
                        filtered_motion{j}= [X_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j))'; Y_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j))']; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
                        choice_times{j} = [X_data(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'; Y_data(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'];
                        filtered_velocity{j}= velocity_data(SLEAP_data.idx_time > time_ranges_trials(1,j) & SLEAP_data.idx_time < time_ranges_trials(3,j));
                        % filtered_gcamp{j}= Y_dF_all_session(gcamp_samples(onset:offset));
                        % filtered{j} = Y_data_filtered(SLEAP_data.idx_frame(onset:offset))'; %SLEAP_data.vel_cm_s(SLEAP_data.idx_frame(onset:offset));
                        good_index = good_index + 1;
                    % end
                end
            end
            % if KEEPDATA
            %     data.streams.Motion.filtered = filtered;
            % else
            %     data.streams.Motion.data = filtered;
            %     data.streams.Motion.filtered = [];
            % end
        end

        %%
        for col = 1:size(filtered_motion, 2)
            coordinates = filtered_motion{col};
            x = coordinates(1,:);
            y = coordinates(2,:);

            colResult = zeros(1, length(x)); % Initialize result for this column

            for hh = 1:length(x)
                % Check if coordinates are within the circle
                if sqrt((x(hh) - shapeData{1}.Center(1)).^2 + (y(hh) - shapeData{1}.Center(2)).^2) <= shapeData{1}.Radius
                    colResult(hh) = 1;
                    % Check if coordinates are within the first rectangle
                elseif x(hh) >= shapeData{2}.Center(1) - shapeData{2}.Size(1)/2 && x(hh) <= shapeData{2}.Center(1) + shapeData{2}.Size(1)/2 && ...
                        y(hh) >= shapeData{2}.Center(2) - shapeData{2}.Size(2)/2 && y(hh) <= shapeData{2}.Center(2) + shapeData{2}.Size(2)/2
                    colResult(hh) = 2;
                    % Check if coordinates are within the second rectangle
                elseif x(hh) >= shapeData{3}.Center(1) - shapeData{3}.Size(1)/2 && x(hh) <= shapeData{3}.Center(1) + shapeData{3}.Size(1)/2 && ...
                        y(hh) >= shapeData{3}.Center(2) - shapeData{3}.Size(2)/2 && y(hh) <= shapeData{3}.Center(2) + shapeData{3}.Size(2)/2
                    colResult(hh) = 3;
                end
            end

            resultArray{col} = colResult; % Store result for this column in resultArray
        end

        for m = 1:size(resultArray, 2)
            results = resultArray{1, m};
            % reward_cup_time{dd}(m) = (sum(results == 1)/fs_cam;
            % left_screen_time{dd}(m) = sum(results == 2)/fs_cam;
            % right_screen_time{dd}(m) = sum(results == 3)/fs_cam;
            reward_cup_time{dd}(m) = sum(results == 1)/size(resultArray{1, m}, 2);

            if strcmp(large_rew_side, 'left')
                large_rew_screen_time{dd}(m) = sum(results == 2)/size(resultArray{1, m}, 2);
                small_rew_screen_time{dd}(m) = sum(results == 3)/size(resultArray{1, m}, 2);

            elseif strcmp(large_rew_side, 'right')
                small_rew_screen_time{dd}(m) = sum(results == 2)/size(resultArray{1, m}, 2);
                large_rew_screen_time{dd}(m) = sum(results == 3)/size(resultArray{1, m}, 2);

            end

            left_screen_time{dd}(m) = sum(results == 2)/size(resultArray{1, m}, 2);
            right_screen_time{dd}(m) = sum(results == 3)/size(resultArray{1, m}, 2);
            other_zone_time{dd}(m) = sum(results == 0)/size(resultArray{1, m}, 2);

        end
        clear resultArray

        mean_reward_cup_large_B1(dd) = mean(reward_cup_time{1, dd}(BehavData.bigSmall == 1.2 & BehavData.Block == 1));
        mean_reward_cup_large_B2(dd) = mean(reward_cup_time{1, dd}(BehavData.bigSmall == 1.2 & BehavData.Block == 2));
        mean_reward_cup_large_B3(dd)= mean(reward_cup_time{1, dd}(BehavData.bigSmall == 1.2 & BehavData.Block == 3));

        mean_reward_cup_small_B1(dd) = mean(reward_cup_time{1, dd}(BehavData.bigSmall == 0.3 & BehavData.Block == 1));
        mean_reward_cup_small_B2(dd) = mean(reward_cup_time{1, dd}(BehavData.bigSmall == 0.3 & BehavData.Block == 2));
        mean_reward_cup_small_B3(dd) = mean(reward_cup_time{1, dd}(BehavData.bigSmall == 0.3 & BehavData.Block == 3));

        mean_left_screen_time_large_B1(dd) = mean(left_screen_time{1, dd}(BehavData.bigSmall == 1.2 & BehavData.Block == 1));
        mean_left_screen_time_large_B2(dd) = mean(left_screen_time{1, dd}(BehavData.bigSmall == 1.2 & BehavData.Block == 2));
        mean_left_screen_time_large_B3(dd) = mean(left_screen_time{1, dd}(BehavData.bigSmall == 1.2 & BehavData.Block == 3));

        mean_left_screen_time_small_B1(dd) = mean(left_screen_time{1, dd}(BehavData.bigSmall == 0.3 & BehavData.Block == 1));
        mean_left_screen_time_small_B2(dd) = mean(left_screen_time{1, dd}(BehavData.bigSmall == 0.3 & BehavData.Block == 2));
        mean_left_screen_time_small_B3(dd) = mean(left_screen_time{1, dd}(BehavData.bigSmall == 0.3 & BehavData.Block == 3));

        mean_right_screen_time_large_B1(dd) = mean(right_screen_time{1, dd}(BehavData.bigSmall == 1.2 & BehavData.Block == 1));
        mean_right_screen_time_large_B2(dd) = mean(right_screen_time{1, dd}(BehavData.bigSmall == 1.2 & BehavData.Block == 2));
        mean_right_screen_time_large_B3(dd) = mean(right_screen_time{1, dd}(BehavData.bigSmall == 1.2 & BehavData.Block == 3));

        mean_right_screen_time_small_B1(dd) = mean(right_screen_time{1, dd}(BehavData.bigSmall == 0.3 & BehavData.Block == 1));
        mean_right_screen_time_small_B2(dd) = mean(right_screen_time{1, dd}(BehavData.bigSmall == 0.3 & BehavData.Block == 2));
        mean_right_screen_time_small_B3(dd) = mean(right_screen_time{1, dd}(BehavData.bigSmall == 0.3 & BehavData.Block == 3));
        

        mean_reward_cup_B1(dd) = mean(reward_cup_time{1, dd}(BehavData.Block == 1 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_right_screen_time_B1(dd) = mean(right_screen_time{1, dd}(BehavData.Block == 1 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_left_screen_time_B1(dd) = mean(left_screen_time{1, dd}(BehavData.Block == 1 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        other_zone_time_B1(dd) = mean(other_zone_time{1, dd}(BehavData.Block == 1 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));

        mean_reward_cup_B2(dd) = mean(reward_cup_time{1, dd}(BehavData.Block == 2 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_right_screen_time_B2(dd) = mean(right_screen_time{1, dd}(BehavData.Block == 2 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_left_screen_time_B2(dd) = mean(left_screen_time{1, dd}(BehavData.Block == 2 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        other_zone_time_B2(dd) = mean(other_zone_time{1, dd}(BehavData.Block == 2 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));

        mean_reward_cup_B3(dd) = mean(reward_cup_time{1, dd}(BehavData.Block == 3 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_right_screen_time_B3(dd) = mean(right_screen_time{1, dd}(BehavData.Block == 3 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_left_screen_time_B3(dd) = mean(left_screen_time{1, dd}(BehavData.Block == 3 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        other_zone_time_B3(dd) = mean(other_zone_time{1, dd}(BehavData.Block == 3 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));

        mean_large_rew_screen_time_B1(dd) = mean(large_rew_screen_time{1, dd}(BehavData.Block == 1 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_large_rew_screen_time_B2(dd) = mean(large_rew_screen_time{1, dd}(BehavData.Block == 2 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_large_rew_screen_time_B3(dd) = mean(large_rew_screen_time{1, dd}(BehavData.Block == 3 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));


        mean_small_rew_screen_time_B1(dd) = mean(small_rew_screen_time{1, dd}(BehavData.Block == 1 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_small_rew_screen_time_B2(dd) = mean(small_rew_screen_time{1, dd}(BehavData.Block == 2 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));
        mean_small_rew_screen_time_B3(dd) = mean(small_rew_screen_time{1, dd}(BehavData.Block == 3 & BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0));

    end
end


dwell_times_table = table;

dwell_times_table.Animals = animalIDs;
dwell_times_table.mean_large_screen_time_B1 = mean_large_rew_screen_time_B1';
dwell_times_table.mean_large_screen_time_B2 = mean_large_rew_screen_time_B2';
dwell_times_table.mean_large_screen_time_B3 = mean_large_rew_screen_time_B3';

dwell_times_table.mean_small_screen_time_B1 = mean_small_rew_screen_time_B1';
dwell_times_table.mean_small_screen_time_B2 = mean_small_rew_screen_time_B2';
dwell_times_table.mean_small_screen_time_B3 = mean_small_rew_screen_time_B3';

dwell_times_table.mean_reward_cup_B1 = mean_reward_cup_B1';
dwell_times_table.mean_reward_cup_B2 = mean_reward_cup_B2';
dwell_times_table.mean_reward_cup_B3 = mean_reward_cup_B3';

dwell_times_table.other_zone_time_B1 = other_zone_time_B1';
dwell_times_table.other_zone_time_B2 = other_zone_time_B2';
dwell_times_table.other_zone_time_B3 = other_zone_time_B3';

%%
[~, loc] = ismember(IDs_from_list, animalIDs);

dwell_times_table_reorg = dwell_times_table(loc, :);

% dwell_times_table_reorg.treatment = hM4Di_treatment_groups;
% dwell_times_table_reorg.treatment = stGtACR_treatment_groups;
% dwell_times_table_reorg.treatment = PdCO_treatment_groups;
dwell_times_table_reorg.treatment = ChrimsonR_treatment_groups;

% remove any mice where all values are 0 - this means these mice had no
% session. should be double checked to make sure 0s aren't from a bug or
% similar
rows_to_remove = all(dwell_times_table_reorg{:, 2:end-1} == 0, 2);

dwell_times_table_reorg(rows_to_remove, :) = [];

% 
% control_means = mean(dwell_times_table_reorg(strcmp(dwell_times_table_reorg.treatment, 'mCherry'), 2:7))
% stgtacr_means = mean(dwell_times_table_reorg(strcmp(dwell_times_table_reorg.treatment, 'hM4Di'), 2:7))

%% hM4Di 

large_choice_mCherry = [dwell_times_table_reorg.mean_large_screen_time_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_large_screen_time_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_large_screen_time_B1(strcmp('hM4Di', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B2(strcmp('hM4Di', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B3(strcmp('hM4Di', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% hM4Di 

large_choice_mCherry = [dwell_times_table_reorg.mean_small_screen_time_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_small_screen_time_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_small_screen_time_B1(strcmp('hM4Di', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B2(strcmp('hM4Di', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B3(strcmp('hM4Di', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% hM4Di 

large_choice_mCherry = [dwell_times_table_reorg.mean_reward_cup_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_reward_cup_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_reward_cup_B1(strcmp('hM4Di', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B2(strcmp('hM4Di', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B3(strcmp('hM4Di', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));


% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% hM4Di 

large_choice_mCherry = [dwell_times_table_reorg.other_zone_time_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.other_zone_time_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.other_zone_time_B1(strcmp('hM4Di', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B2(strcmp('hM4Di', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B3(strcmp('hM4Di', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));


% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', hM4Di_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, hM4Di_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', hM4Di_color, 'MarkerFaceColor', hM4Di_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% stGtACR 

large_choice_mCherry = [dwell_times_table_reorg.mean_large_screen_time_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_large_screen_time_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_large_screen_time_B1(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B2(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B3(strcmp('stGtACR', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% stGtACR 

large_choice_mCherry = [dwell_times_table_reorg.mean_small_screen_time_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_small_screen_time_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_small_screen_time_B1(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B2(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B3(strcmp('stGtACR', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% stGtACR 

large_choice_mCherry = [dwell_times_table_reorg.mean_reward_cup_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_reward_cup_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_reward_cup_B1(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B2(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B3(strcmp('stGtACR', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));


% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% stGtACR 

large_choice_mCherry = [dwell_times_table_reorg.other_zone_time_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.other_zone_time_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.other_zone_time_B1(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B2(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B3(strcmp('stGtACR', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));


% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% PdCO 

large_choice_mCherry = [dwell_times_table_reorg.mean_large_screen_time_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_large_screen_time_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_large_screen_time_B1(strcmp('PdCO', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B2(strcmp('PdCO', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B3(strcmp('PdCO', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% PdCO 

large_choice_mCherry = [dwell_times_table_reorg.mean_small_screen_time_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_small_screen_time_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_small_screen_time_B1(strcmp('PdCO', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B2(strcmp('PdCO', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B3(strcmp('PdCO', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', PdCO_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, PdCO_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', PdCO_color, 'MarkerFaceColor', PdCO_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% stGtACR 

large_choice_mCherry = [dwell_times_table_reorg.mean_reward_cup_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_reward_cup_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_reward_cup_B1(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B2(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B3(strcmp('stGtACR', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));


% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% stGtACR 

large_choice_mCherry = [dwell_times_table_reorg.other_zone_time_B1(strcmp('mCherry', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.other_zone_time_B2(strcmp('mCherry', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B3(strcmp('mCherry', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.other_zone_time_B1(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B2(strcmp('stGtACR', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B3(strcmp('stGtACR', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));


% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', stGtACR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, stGtACR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', stGtACR_color, 'MarkerFaceColor', stGtACR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%



















%%
%% Females vs Males 

large_choice_mCherry = [dwell_times_table_reorg.mean_large_screen_time_B1(strcmp('Female', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_large_screen_time_B2(strcmp('Female', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B3(strcmp('Female', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_large_screen_time_B1(strcmp('Male', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B2(strcmp('Male', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_large_screen_time_B3(strcmp('Male', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% Females vs Males 

large_choice_mCherry = [dwell_times_table_reorg.mean_small_screen_time_B1(strcmp('Female', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_small_screen_time_B2(strcmp('Female', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B3(strcmp('Female', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_small_screen_time_B1(strcmp('Male', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B2(strcmp('Male', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_small_screen_time_B3(strcmp('Male', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));







% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% hM4Di 

large_choice_mCherry = [dwell_times_table_reorg.mean_reward_cup_B1(strcmp('Female', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.mean_reward_cup_B2(strcmp('Female', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B3(strcmp('Female', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.mean_reward_cup_B1(strcmp('Male', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B2(strcmp('Male', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.mean_reward_cup_B3(strcmp('Male', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));


% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% hM4Di 

large_choice_mCherry = [dwell_times_table_reorg.other_zone_time_B1(strcmp('Female', dwell_times_table_reorg.treatment)),dwell_times_table_reorg.other_zone_time_B2(strcmp('Female', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B3(strcmp('Female', dwell_times_table_reorg.treatment))]*100;
large_choice_hM4Di = [dwell_times_table_reorg.other_zone_time_B1(strcmp('Male', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B2(strcmp('Male', dwell_times_table_reorg.treatment)), dwell_times_table_reorg.other_zone_time_B3(strcmp('Male', dwell_times_table_reorg.treatment))]*100;

mean_large = nanmean(large_choice_mCherry, 1);
mean_small = nanmean(large_choice_hM4Di, 1);
sem_large = nanstd(large_choice_mCherry, 0, 1) ./ sqrt(size(large_choice_mCherry, 1));
sem_small = nanstd(large_choice_hM4Di, 0, 1) ./ sqrt(size(large_choice_hM4Di, 1));


% X-axis points
x_points = 1:size(large_choice_mCherry, 2);


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_choice_mCherry, 1)
    plot(x_points, large_choice_mCherry(i, :), '-', ...
        'Color', mCherry_color, ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(large_choice_hM4Di, 1)
    plot(x_points, large_choice_hM4Di(i, :), '-', ...
        'Color', ChrimsonR_color, ... % Red with 60% opacity
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large, sem_large, mCherry_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', mCherry_color, 'MarkerFaceColor', mCherry_color, ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small, sem_small, ChrimsonR_symbol, ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', ChrimsonR_color, 'MarkerFaceColor', ChrimsonR_color, ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50', '75'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 50]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:50);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

























