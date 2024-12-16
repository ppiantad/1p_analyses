% animalIDs = (fieldnames(final_SLEAP));
session_to_analyze = 'Pre_RDT_RM';

b1_large_path_length = [];
b2_large_path_length = [];
b3_large_path_length = [];

b1_small_path_length = [];
b2_small_path_length = [];
b3_small_path_length = [];

animals_with_sessions = {}; 

for dd = 1:size(animalIDs)

    select_mouse = animalIDs{dd};
    path_length_array = [];
    if isfield(final_SLEAP.(select_mouse), session_to_analyze)
        animals_with_sessions{dd} = select_mouse; 

        



        SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;

        X_data = SLEAP_data.corrected_x_pix;
        Y_data = SLEAP_data.corrected_y_pix;
        [X_data, Y_data] = correct_XY_outliers_v1(X_data, Y_data);

        % onset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.stTime';
        % choice_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.choiceTime';
        % offset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.collectionTime';
        fs_cam = 30; %set sampling rate according to camera, this is hard coded for now
        % time_ranges_trials = [onset_trials; choice_trials; offset_trials];


        % gcamp_samples = 1:1:size(Y_dF_all_session, 2);

        % gcamp_time = (0:length(F405_downsampled_data)-1)/fs_cam;

        SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;

        % velocity_data = final_SLEAP.(select_mouse).(session_to_analyze).zscored_SLEAP_data_velocity';

        velocity_data = zscore(SLEAP_data.vel_cm_s)';

        % SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data;

        
        BehavData = final_SLEAP.(select_mouse).(session_to_analyze).BehavData;

        % uncomment below if you want to run this on specific trials, eg
        % using the 10x variable data 
        % BehavData = behav_tbl_iter{1, 1}{dd};

        onset_trials = BehavData.stTime';
        choice_trials = BehavData.choiceTime';
        offset_trials = BehavData.collectionTime';
        time_ranges_trials = [onset_trials; choice_trials; offset_trials];
        adjusted_start_time = BehavData.TrialPossible(1)-60;
        SLEAP_data.idx_time = SLEAP_data.idx_time+adjusted_start_time;



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
                        % choice_times{j} = [X_data(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'; Y_data(interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), time_ranges_trials(2,j), 'nearest'))'];
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
        % filtered_motion_array{}
      
        b1_large_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 1));
        b2_large_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 2));
        b3_large_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 3));

        b1_small_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 1));
        b2_small_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 2));
        b3_small_path_length(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 3));

        large_path_length_all_blocks(dd) = mean(path_length_array(1, BehavData.bigSmall == 1.2));
        small_path_length_all_blocks(dd) = mean(path_length_array(1, BehavData.bigSmall == 0.3));
    end

    b1_path_length_mouse{dd} = path_length_array(1, BehavData.Block == 1 &  BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0);
    b2_path_length_mouse{dd} = path_length_array(1, BehavData.Block == 2 &  BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0);
    b3_path_length_mouse{dd} = path_length_array(1, BehavData.Block == 3 &  BehavData.omissionALL == 0 &  BehavData.Blank_Touch == 0);
    path_length_array_mouse{dd} = path_length_array; 
end

path_length_concat = cat(2, path_length_array_mouse{:})';

%%

large_path_length = ([b1_large_path_length; b2_large_path_length; b3_large_path_length])';
small_path_length = ([b1_small_path_length; b2_small_path_length; b3_small_path_length])';



mean_large_path = nanmean(large_path_length, 1);
mean_small_path = nanmean(small_path_length, 1);

sem_large = nanstd(large_path_length, 0, 1) ./ sqrt(size(large_path_length, 1));
sem_small = nanstd(small_path_length, 0, 1) ./ sqrt(size(small_path_length, 1));


% X-axis points
x_points = 1:3;


% Plotting
figure;
hold on;

% Set figure size
width = 200; % Width of the figure
height = 450; % Height of the figure
set(gcf, 'Position', [50, 25, width, height]); % Set position and size

% Plot individual lines for "Large" data
for i = 1:size(large_path_length, 1)
    plot(x_points, large_path_length(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end


% Plot individual lines for "Small" data
for i = 1:size(small_path_length, 1)
    plot(x_points, small_path_length(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % 
        'LineWidth', 1.2);
end

% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large_path, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small_path, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50%', '75%'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% % Set axis limits, labels, and legend
% ylim([0 1.1 * max([mean_large_path + sem_large, ...
%                    mean_small_path + sem_small])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:200:1000);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;



%%
% Extract path length array (Y-axis)
path_length_array = path_length_array_mouse{1, 10};


% Extract the variable to correlate (X-axis)
variable_x = variable_to_correlate{1, 10};

% Extract the bigSmall column
bigSmall = behav_tbl_iter{1, 1}{10, 1}.bigSmall;

% Find indices for bigSmall == 0.3 and bigSmall == 1.2
indices_red = bigSmall == 0.3;
indices_blue = bigSmall == 1.2;
figure;
% Plot red points (bigSmall == 0.3)
scatter(variable_x(indices_red), path_length_array(indices_red), 'r', 'filled'); 
hold on;

% Plot blue points (bigSmall == 1.2)
scatter(variable_x(indices_blue), path_length_array(indices_blue), 'b', 'filled');

% Add labels and legend
xlabel('Variable to Correlate');
ylabel('Path Length');
legend({'bigSmall = 0.3', 'bigSmall = 1.2'}, 'Location', 'Best');
title('Scatter Plot of Path Length with bigSmall Grouping');
hold off;
