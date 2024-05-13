
fs_cam = 30; %set sampling rate according to camera, this is hard coded for now

animalIDs = (fieldnames(final_SLEAP));
session_to_analyze = 'Pre_RDT_RM';

reward_cup_time = [];
right_screen_time = [];
left_screen_time = [];


animals_with_sessions = {};

for dd = 1:size(animalIDs)

    select_mouse = animalIDs{dd};
    if isfield(final_SLEAP.(select_mouse), session_to_analyze)
        shapeData = final_SLEAP.(select_mouse).(session_to_analyze).shapeData;
        animals_with_sessions{dd} = select_mouse;

        SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;
        X_data = SLEAP_data.corrected_x_pix;
        Y_data = SLEAP_data.corrected_y_pix;


        onset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.stTime';
        choice_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.choiceTime';
        offset_trials = final_SLEAP.(select_mouse).(session_to_analyze).BehavData.collectionTime';

        time_ranges_trials = [onset_trials; choice_trials; offset_trials];

        SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;

        % velocity_data = final_SLEAP.(select_mouse).(session_to_analyze).zscored_SLEAP_data_velocity';

        velocity_data = zscore(SLEAP_data.vel_cm_s)';

        BehavData = final_SLEAP.(select_mouse).(session_to_analyze).BehavData;
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






    end
end

mean_left_screen_time_all_blocks = mean([mean_left_screen_time_B1; mean_left_screen_time_B2; mean_left_screen_time_B3]);
mean_right_screen_time_all_blocks = mean([mean_right_screen_time_B1; mean_right_screen_time_B2; mean_right_screen_time_B3]);
mean_reward_cup_all_blocks = mean([mean_reward_cup_B1; mean_reward_cup_B2; mean_reward_cup_B3]);
mean_other_zone_time_all_blocks = mean([other_zone_time_B1; other_zone_time_B2; other_zone_time_B3]);