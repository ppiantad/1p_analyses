animalIDs = fieldnames(final_DLC);

session_to_analyze = 'D3';

recorded_fps = 30;

% Parameters
stimulus_duration = 2 * 60; % 2 minutes in seconds
num_repeats = 3;
total_stimuli = 2;

% Initialize variables
stimulus_times = cell(total_stimuli, 1);
current_time = 0;

% Loop through each stimulus alternately and calculate start and end times
for j = 1:num_repeats
    for i = 1:total_stimuli
        start_time = current_time;
        end_time = start_time + stimulus_duration;
        
        if isempty(stimulus_times{i})
            stimulus_times{i} = [start_time, end_time];
        else
            stimulus_times{i} = [stimulus_times{i}; start_time, end_time];
        end
        
        current_time = end_time;
    end
end


% Multiply every value in stimulus_times by the FPS
for i = 1:total_stimuli
    stimulus_frames{i} = stimulus_times{i} * recorded_fps;
end



%%

mouse_count = 0;
for gg = 1:size(animalIDs, 1)
    current_mouse = animalIDs{gg};
    if strcmp(final_DLC.(current_mouse).experimental_grp, 'One Context')
        mouse_count = mouse_count+1;
        freeze_data = final_DLC.(current_mouse).(session_to_analyze).DLC_data_raw.freeze; 
        frames_data = final_DLC.(current_mouse).(session_to_analyze).DLC_data_raw.frame;
        
        for qq = 1:num_repeats
            safe_context(:, qq) = freeze_data(stimulus_frames{1, 1}(qq,1)+1:stimulus_frames{1, 1}(qq,2));
        end
        
        
        for pp = 1:size(safe_context, 2)
            mean_safe_context(mouse_count, pp) = sum(freeze_data(stimulus_frames{1, 1}(pp,1)+1:stimulus_frames{1, 1}(pp,2)))/size(safe_context, 1);

        end


        for qq = 1:num_repeats
            if qq <= 2
                aversive_context{qq} = freeze_data(stimulus_frames{1, 2}(qq,1)+1:stimulus_frames{1, 2}(qq,2));
            elseif qq > 2
                aversive_context{qq} = freeze_data(stimulus_frames{1, 2}(qq,1)+1:end);
            end
        end


        for pp = 1:size(aversive_context,2 )

            mean_aversive_context(mouse_count, pp) = sum(aversive_context{1,pp})/size(aversive_context{1,pp}, 1);

        end
    end




end


% Calculate the mean of each column
mean_safe = mean(mean_safe_context);
mean_aversive = mean(mean_aversive_context);

% Initialize the interleaved array
interleaved_means = zeros(1, 6);

% Interleave the means
interleaved_means(1:2:end) = mean_safe;
interleaved_means(2:2:end) = mean_aversive;

% Display the interleaved means
disp('Interleaved Means:');
disp(interleaved_means);

figure; plot(interleaved_means);

%%
% Plotting the data
figure;
plot(final_DLC.B46859.D3.DLC_data_raw  .frame, final_DLC.B46859.D3.DLC_data_raw  .freeze);
hold on;

% Define colors and transparency
color_1 = [0, 0, 1]; % Blue for stimulus 1
color_2 = [1, 0, 0]; % Red for stimulus 2
transparency = 0.3; % Transparency level

% Plot rectangles for each stimulus presentation
for i = 1:total_stimuli
    for j = 1:num_repeats
        % Get start and end frames
        start_frame = stimulus_frames{i}(j, 1);
        end_frame = stimulus_frames{i}(j, 2);
        
        % Determine the Y-axis range
        y_limits = ylim;
        
        % Plot the rectangle
        rectangle('Position', [start_frame, y_limits(1), end_frame-start_frame, y_limits(2)-y_limits(1)], ...
                  'FaceColor', [color_1 transparency] * (i == 1) + [color_2 transparency] * (i == 2), ...
                  'EdgeColor', 'none');
    end
end

% Finalize the plot
hold off;
