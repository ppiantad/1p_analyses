animalIDs = (fieldnames(final_SLEAP));
session_to_analyze = 'RDT_OPTO_CHOICE';

b1_large_path_length = [];
b2_large_path_length = [];
b3_large_path_length = [];

b1_small_path_length = [];
b2_small_path_length = [];
b3_small_path_length = [];

animals_with_sessions = {}; 

for dd = 1:size(animalIDs)

    select_mouse = animalIDs{dd};
    if isfield(final_SLEAP.(select_mouse), session_to_analyze)
        animals_with_sessions{dd} = select_mouse; 





        SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;
        X_data = SLEAP_data.x_pix;
        Y_data = SLEAP_data.y_pix;


        % Apply Savitzky-Golay filter to each row. this mostly removes high
        % frequency changes that occur due to slight jumps in keypoint location
        X_data = sgolayfilt(X_data, 9, 33);
        Y_data = sgolayfilt(Y_data, 9, 33);
        % [X_data, Y_data] = correct_XY_outliers_v1(X_data, Y_data);

        % onset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.stTime';
        % choice_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.choiceTime';
        % offset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.collectionTime';
        % fs_cam = 30; %set sampling rate according to camera, this is hard coded for now
        % time_ranges_trials = [onset_trials; choice_trials; offset_trials];


        % gcamp_samples = 1:1:size(Y_dF_all_session, 2);

        % gcamp_time = (0:length(F405_downsampled_data)-1)/fs_cam;

        % velocity_data = final_SLEAP.(select_mouse).(session_to_analyze).zscored_SLEAP_data_velocity';

        

        % SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data;


        BehavData = final_behavior.(select_mouse).(session_to_analyze).uv.BehavData;
        % adjusted_start_time = BehavData.TrialPossible(1)-60;
        % adjusted_start_time = final_behavior.(select_mouse).(session_to_analyze).uv.session_start_adjustment;
        % SLEAP_data.idx_time = SLEAP_data.idx_time+adjusted_start_time;
        velocity_data = zscore(SLEAP_data.vel_cm_s)';

        onset_trials = BehavData.stTime';
        choice_trials = BehavData.choiceTime';
        offset_trials = BehavData.collectionTime';
        fs_cam = 30; %set sampling rate according to camera, this is hard coded for now
        time_ranges_trials = [onset_trials; choice_trials; offset_trials];


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

        figure;
        hold on
        % Adding labels and title (customize as needed)
        xlabel('X-axis Label');
        ylabel('Y-axis Label');
        title('Plotting the Line');

        % Display grid if desired
        grid on;

        for ii = 1:size(filtered_motion, 2)
            % Assuming filtered{1, 1} contains your X and Y coordinates
            X = filtered_motion{1, ii}(1, :);
            Y = filtered_motion{1, ii}(2, :);

            % Plotting the line
            plot(X, Y, '-x');  % You can use different markers or line styles if needed


        end
        hold off

        path_length_array = [];
        
        for qq = 1:size(filtered_motion, 2)
            coordinates = filtered_motion{1, qq}';
            distances = pdist2(coordinates, coordinates);

            % Sum up the distances to get the total path length
            path_length = sum(diag(distances, 1));

            % disp(['Path Length: ', num2str(path_length)]);

            distances_matrix{qq} = distances;
            path_length_array(qq) = path_length;
            clear coordinates distances path_length
        end


        b1_large_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 1));
        b2_large_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 2));
        b3_large_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 3));

        b1_small_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 1));
        b2_small_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 2));
        b3_small_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 3));
    end
end
path_length_table = table;
path_length_table.Animals = animalIDs;
path_length_table.b1_large_path_length = b1_large_path_length';
path_length_table.b2_large_path_length = b2_large_path_length';
path_length_table.b3_large_path_length = b3_large_path_length';

path_length_table.b1_small_path_length = b1_small_path_length';
path_length_table.b2_small_path_length = b2_small_path_length';
path_length_table.b3_small_path_length = b3_small_path_length';

%% need to run generate_behav_figs first

%for hM4Di vs mCherry

% Find indices where Animals match valid_animalIDs
valid_idx = ismember(path_length_table.Animals, risk_table.valid_animalIDs);

% Filter path_length_table to only include valid animals
filtered_path_length_table = path_length_table(valid_idx, :);

% Extract TreatmentCondition for the matched valid_animalIDs
[~, loc] = ismember(filtered_path_length_table.Animals, risk_table.valid_animalIDs);

% Get the corresponding TreatmentCondition from risk_table
filtered_treatment_conditions = risk_table.TreatmentCondition(loc);

% Find indices where TreatmentCondition is 'mCherry'
mCherry_idx = strcmp(filtered_treatment_conditions, 'mCherry');

% Extract the relevant data
large_choice_mCherry = [filtered_path_length_table.b1_large_path_length(mCherry_idx), ...
                         filtered_path_length_table.b2_large_path_length(mCherry_idx), ...
                         filtered_path_length_table.b3_large_path_length(mCherry_idx)];



% Find indices where TreatmentCondition is 'mCherry'
hM4Di_idx = strcmp(filtered_treatment_conditions, 'stGtACR');

% Extract the relevant data
large_choice_hM4Di = [filtered_path_length_table.b1_large_path_length(hM4Di_idx), ...
                         filtered_path_length_table.b2_large_path_length(hM4Di_idx), ...
                         filtered_path_length_table.b3_large_path_length(hM4Di_idx)];


% large_choice_mCherry = [path_length_table.b1_large_path_length(strcmp('mCherry', risk_table.TreatmentCondition)), path_length_table.b2_large_path_length(strcmp('mCherry', risk_table.TreatmentCondition)), path_length_table.b3_large_path_length(strcmp('mCherry', risk_table.TreatmentCondition))];
% large_choice_hM4Di = [path_length_table.b1_large_path_length(strcmp('hM4Di', risk_table.TreatmentCondition)), path_length_table.b2_large_path_length(strcmp('hM4Di', risk_table.TreatmentCondition)), path_length_table.b3_large_path_length(strcmp('hM4Di', risk_table.TreatmentCondition))];

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
% ylim([0 1.1 * max([mean_large + sem_large, ...
%                    mean_small + sem_small])]); % Adjust ylim dynamically
% set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%% for hM4Di vs mCherry

% Find indices where Animals match valid_animalIDs
valid_idx = ismember(path_length_table.Animals, risk_table.valid_animalIDs);

% Filter path_length_table to only include valid animals
filtered_path_length_table = path_length_table(valid_idx, :);

% Extract TreatmentCondition for the matched valid_animalIDs
[~, loc] = ismember(filtered_path_length_table.Animals, risk_table.valid_animalIDs);

% Get the corresponding TreatmentCondition from risk_table
filtered_treatment_conditions = risk_table.TreatmentCondition(loc);

% Find indices where TreatmentCondition is 'mCherry'
mCherry_idx = strcmp(filtered_treatment_conditions, 'mCherry');

% Extract the relevant data
large_choice_mCherry = [filtered_path_length_table.b1_small_path_length(mCherry_idx), ...
                         filtered_path_length_table.b2_small_path_length(mCherry_idx), ...
                         filtered_path_length_table.b3_small_path_length(mCherry_idx)];



% Find indices where TreatmentCondition is 'mCherry'
hM4Di_idx = strcmp(filtered_treatment_conditions, 'stGtACR');

% Extract the relevant data
large_choice_hM4Di = [filtered_path_length_table.b1_small_path_length(hM4Di_idx), ...
                         filtered_path_length_table.b2_small_path_length(hM4Di_idx), ...
                         filtered_path_length_table.b3_small_path_length(hM4Di_idx)];

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
% ylim([0 1.1 * max([mean_large + sem_large, ...
%                    mean_small + sem_small])]); % Adjust ylim dynamically
% set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;