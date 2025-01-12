%% discrimination sessions
clear aversive_context safe_context mean_aversive mean_safe mean_aversive_context mean_safe_context

recorded_fps = 30;



% Parameters
stimulus_duration = 2 * 60; % 2 minutes in seconds
num_repeats = 3;
total_stimuli = 2;

% Initialize variables
stimulus_times = cell(total_stimuli, 1);
current_time = 0;

session_length_in_min = 12;

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


% experimental_grps = readtable('I:\MATLAB\my_repo\context fear\organize_DLC_data\PFC mice.xlsx');
% experimental_grps = readtable('E:\MATLAB\my_repo\context fear\organize_SLEAP_data\full_pilot_mice.xlsx');
experimental_grps = readtable('i:\MATLAB\my_repo\context fear\organize_DLC_data\pilot groups.xlsx');

% Define parameters
threshold = 1; % Velocity threshold
sample_duration = 0.03; % Duration of each sample in seconds
min_duration = 2; % Minimum duration to trigger labeling in seconds
frame_rate = 30; % Frames per second
% Calculate the minimum number of consecutive rows needed
min_samples = min_duration / sample_duration;

animalIDs = fieldnames(final_DLC);

sessions_to_analyze = {'D3', 'D4'};

for ff = 1:size(sessions_to_analyze, 2)
    session_to_analyze = sessions_to_analyze{ff};
    mouse_count = 0;
    for gg = 1:size(animalIDs, 1)
        current_mouse = animalIDs{gg};
        if strcmp(session_to_analyze, 'D1_Morning') & strcmp(current_mouse, 'B57417')
            continue
        elseif strcmp(session_to_analyze, 'D3') & strcmp(current_mouse, 'B51618')
            continue
        else
            if isfield(final_DLC.(current_mouse), session_to_analyze)
                mouse_count = mouse_count+1;
                DLC_data_mouse = final_DLC.(current_mouse).(session_to_analyze).movement_data;

                body_velocity = [];
                labels = [];
                % Get the body_velocity column
                body_velocity = final_DLC.(current_mouse).(session_to_analyze).movement_data.body_velocity;

                % Initialize the new column
                labels = zeros(size(body_velocity));

                % Find consecutive segments where body_velocity < threshold
                below_threshold = body_velocity < threshold;
                start_idx = find(diff([0; below_threshold]) == 1); % Start indices
                end_idx = find(diff([below_threshold; 0]) == -1); % End indices

                % Iterate through each segment and label
                for i = 1:length(start_idx)
                    segment_length = end_idx(i) - start_idx(i) + 1;
                    if segment_length >= min_samples
                        labels(start_idx(i):end_idx(i)) = 1;
                    end
                end

                % Add the labels as a new column to the table
                % final_DLC.B46837.D1_Afternoon.movement_data.freeze_label = labels;

                % Find the row index where the 'mouse' column matches 'current_mouse'
                row_idx = strcmp(experimental_grps.mouse, current_mouse);

                freeze_data{ff}(mouse_count, :) = labels(1:21590)';
                % freeze_data(mouse_count, :) = DLC_data_mouse.was_freezing(1:21590)';
                experimental_grps_updated{ff}(mouse_count, :) = experimental_grps(row_idx, :);
            end
        end
    end
end

%%

for ff = 1:size(sessions_to_analyze, 2)
    session_to_analyze = sessions_to_analyze{ff};
    freeze_data_extract = freeze_data{ff};
    if strcmp(session_to_analyze, 'D3')
        for gg = 1:size(freeze_data_extract, 1)
            freeze_data_session = freeze_data_extract(gg, :);

            for qq = 1:num_repeats
                safe_context{ff}{gg, qq} = freeze_data_session(stimulus_frames{1, 1}(qq,1)+1:stimulus_frames{1, 1}(qq,2));
            end

            for qq = 1:num_repeats
                if qq <= 2
                    aversive_context{ff}{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:stimulus_frames{1, 2}(qq,2));
                elseif qq > 2
                    aversive_context{ff}{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:end);
                end
            end
        end

    elseif strcmp(session_to_analyze, 'D4')

        for gg = 1:size(freeze_data_extract, 1)
            freeze_data_session = freeze_data_extract(gg, :);

            for qq = 1:num_repeats
                aversive_context{ff}{gg, qq} = freeze_data_session(stimulus_frames{1, 1}(qq,1)+1:stimulus_frames{1, 1}(qq,2));
            end

            for qq = 1:num_repeats
                if qq <= 2
                    safe_context{ff}{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:stimulus_frames{1, 2}(qq,2));
                elseif qq > 2
                    safe_context{ff}{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:end);
                end
            end
        end
    end

end


% for gg = 1:size(freeze_data, 1)
%     freeze_data_session = freeze_data(gg, :);
% 
%     for qq = 1:num_repeats
%         safe_context{gg, qq} = freeze_data_session(stimulus_frames{1, 1}(qq,1)+1:stimulus_frames{1, 1}(qq,2));
%     end
% 
%     for qq = 1:num_repeats
%         if qq <= 2
%             aversive_context{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:stimulus_frames{1, 2}(qq,2));
%         elseif qq > 2
%             aversive_context{gg, qq} = freeze_data_session(stimulus_frames{1, 2}(qq,1)+1:end);
%         end
%     end
% end
for gg = 1:size(sessions_to_analyze, 2)

    safe_context_extract = safe_context{gg, 1};
    for pp = 1:size(safe_context_extract, 1)
        for ff = 1:size(safe_context_extract, 2)
            % mean_safe_context(mouse_count, pp) = sum(freeze_data(stimulus_frames{1, 1}(pp,1)+1:stimulus_frames{1, 1}(pp,2)))/size(safe_context, 1);
            % Calculate the proportion of 1s in the column
            b = sum(safe_context_extract{pp,ff}) / size(safe_context_extract{pp,ff}, 2);

            % Store the proportion in mean_aversive_context
            mean_safe_context{gg}(pp, ff) = b;

            sem_safe_context{gg}(pp, ff) = b/sqrt(size(safe_context_extract{pp,ff}, 2));

            % % Calculate the number of rows in the column
            % n = size(safe_context{pp,ff}, 2);
            %
            % % Calculate the standard error of the proportion
            % standard_error_safe(pp, ff) = sqrt(b * (1 - b) / n);
        end
    end
end

for gg = 1:size(sessions_to_analyze, 2)
    session_to_analyze = sessions_to_analyze{ff};
    for pp = 1:size(aversive_context, 1)
        for ff = 1:size(aversive_context, 2)
            % mean_safe_context(mouse_count, pp) = sum(freeze_data(stimulus_frames{1, 1}(pp,1)+1:stimulus_frames{1, 1}(pp,2)))/size(safe_context, 1);
            % Calculate the proportion of 1s in the column
            b = sum(aversive_context{pp,ff}) / size(aversive_context{pp,ff}, 2);

            % Store the proportion in mean_aversive_context
            mean_aversive_context(pp, ff) = b;

            sem_aversive_context(pp, ff) = b/sqrt(size(aversive_context{pp,ff}, 2));

            % % Calculate the number of rows in the column
            % n = size(aversive_context{pp,ff}, 2);
            %
            % % Calculate the standard error of the proportion
            % standard_error_safe(pp, ff) = sqrt(b * (1 - b) / n);
        end
    end
end



if any("sex" == string(experimental_grps.Properties.VariableNames))


    experimental_data_safe_male = mean_safe_context(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_mean_safe_male = mean(experimental_data_safe_male);
    experimental_mean_safe_collapsed_male = mean(experimental_data_safe_male, 2);
    experimental_sem_safe_male = std(experimental_data_safe_male)/sqrt(size(experimental_data_safe_male, 1));
    experimental_mice_safe_male = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);

    experimental_data_aversive_male = mean_aversive_context(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);
    experimental_mean_aversive_male = mean(experimental_data_aversive_male);
    experimental_mean_aversive_collapsed_male = mean(experimental_data_aversive_male, 2);
    experimental_sem_aversive_male = std(experimental_data_aversive_male)/sqrt(size(experimental_data_aversive_male, 1));
    experimental_mice_aversive_male = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'male'), :);

    experimental_data_safe_female = mean_safe_context(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_mean_safe_female = mean(experimental_data_safe_female);
    experimental_mean_safe_collapsed_female = mean(experimental_data_safe_female, 2);
    experimental_sem_safe_female = std(experimental_data_safe_female)/sqrt(size(experimental_data_safe_female, 1));
    experimental_mice_safe_female = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);

    experimental_data_aversive_female = mean_aversive_context(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);
    experimental_mean_aversive_female = mean(experimental_data_aversive_female);
    experimental_mean_aversive_collapsed_female = mean(experimental_data_aversive_female, 2);
    experimental_sem_aversive_female = std(experimental_data_aversive_female)/sqrt(size(experimental_data_aversive_female, 1));
    experimental_mice_aversive_female = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental') & strcmp(experimental_grps_updated.sex, 'female'), :);


else

    experimental_data_safe = mean_safe_context(strcmp(experimental_grps_updated.group, 'Experimental'), :);
    experimental_mean_safe = mean(experimental_data_safe);
    experimental_mean_safe_collapsed = mean(experimental_data_safe, 2);
    experimental_sem_safe = std(experimental_data_safe)/sqrt(size(experimental_data_safe, 1));
    experimental_mice_safe = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental'), :);

    one_context_data_safe = mean_safe_context(strcmp(experimental_grps_updated.group, 'One Context'), :);
    one_context_mean_safe = mean(one_context_data_safe);
    one_context_mean_safe_collapsed = mean(one_context_data_safe, 2);
    one_context_sem_safe = std(one_context_data_safe)/sqrt(size(one_context_data_safe, 1));
    one_context_mice_safe = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'One Context'), :);

    no_shock_data_safe = mean_safe_context(strcmp(experimental_grps_updated.group, 'No Shock'), :);
    no_shock_mean_safe = mean(no_shock_data_safe);
    no_shock_mean_safe_collapsed = mean(no_shock_data_safe, 2);
    no_shock_sem_safe = std(no_shock_data_safe)/sqrt(size(no_shock_data_safe, 1));
    no_shock_mice_safe = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'No Shock'), :);

    experimental_data_aversive = mean_aversive_context(strcmp(experimental_grps_updated.group, 'Experimental'), :);
    experimental_mean_aversive = mean(experimental_data_aversive);
    experimental_mean_aversive_collapsed = mean(experimental_data_aversive, 2);
    experimental_sem_aversive = std(experimental_data_aversive)/sqrt(size(experimental_data_aversive, 1));
    experimental_mice_aversive = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'Experimental'), :);

    one_context_data_aversive = mean_aversive_context(strcmp(experimental_grps_updated.group, 'One Context'), :);
    one_context_mean_aversive = mean(one_context_data_aversive);
    one_context_mean_aversive_collapsed = mean(one_context_data_aversive, 2);
    one_context_sem_aversive = std(one_context_data_aversive)/sqrt(size(one_context_data_aversive, 1));
    one_context_mice_aversive = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'One Context'), :);

    no_shock_data_aversive = mean_aversive_context(strcmp(experimental_grps_updated.group, 'No Shock'), :);
    no_shock_mean_aversive = mean(no_shock_data_aversive);
    no_shock_mean_aversive_collapsed = mean(no_shock_data_aversive, 2);
    no_shock_sem_aversive = std(no_shock_data_aversive)/sqrt(size(no_shock_data_aversive, 1));
    no_shock_mice_aversive = experimental_grps_updated(strcmp(experimental_grps_updated.group, 'No Shock'), :);
end