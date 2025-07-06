
session_to_analyze = 'Pre_RDT_RM';


if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_SLEAP, fieldsToRemove{i})
            final_SLEAP = rmfield(final_SLEAP, fieldsToRemove{i});
        end
    end
elseif strcmp('RDT_D2', session_to_analyze)

    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_SLEAP, fieldsToRemove{i})
            final_SLEAP = rmfield(final_SLEAP, fieldsToRemove{i});
        end
    end
end


fs_cam = 30; %set sampling rate according to camera, this is hard coded for now

animalIDs = (fieldnames(final_SLEAP));

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
%%
mean_left_screen_time_all_blocks = mean([mean_left_screen_time_B1; mean_left_screen_time_B2; mean_left_screen_time_B3]);
mean_right_screen_time_all_blocks = mean([mean_right_screen_time_B1; mean_right_screen_time_B2; mean_right_screen_time_B3]);
mean_reward_cup_all_blocks = mean([mean_reward_cup_B1; mean_reward_cup_B2; mean_reward_cup_B3]);
mean_other_zone_time_all_blocks = mean([other_zone_time_B1; other_zone_time_B2; other_zone_time_B3]);

large_rew = [];

large_rew = zeros(1, 10); % Initialize a result array
for i = 1:size(mean_left_screen_time_B1, 2)
    if mean_right_screen_time_B1(i) > mean_left_screen_time_B1(i)
        large_rew(i) = 2;
    elseif mean_right_screen_time_B1(i) < mean_left_screen_time_B1(i)
        large_rew(i) = 1;
    else
        large_rew(i) = 0; % Default for equal values
    end
end


for ii = 1:size(large_rew, 2)
    if large_rew(ii) == 2
        mean_large_screen_time_B1(ii) = mean_right_screen_time_B1(ii)
        mean_large_screen_time_B2(ii) = mean_right_screen_time_B2(ii)
        mean_large_screen_time_B3(ii) = mean_right_screen_time_B3(ii)

        mean_small_screen_time_B1(ii) = mean_left_screen_time_B1(ii)
        mean_small_screen_time_B2(ii) = mean_left_screen_time_B2(ii)
        mean_small_screen_time_B3(ii) = mean_left_screen_time_B3(ii)

    elseif large_rew(ii) == 1 
        mean_large_screen_time_B1(ii) = mean_left_screen_time_B1(ii)
        mean_large_screen_time_B2(ii) = mean_left_screen_time_B2(ii)
        mean_large_screen_time_B3(ii) = mean_left_screen_time_B3(ii)


        mean_small_screen_time_B1(ii) = mean_right_screen_time_B1(ii)
        mean_small_screen_time_B2(ii) = mean_right_screen_time_B2(ii)
        mean_small_screen_time_B3(ii) = mean_right_screen_time_B3(ii)
    end
end



%%

large_zone_time = ([mean_large_screen_time_B1; mean_large_screen_time_B2; mean_large_screen_time_B3]*100)';
small_zone_time = ([mean_small_screen_time_B1; mean_small_screen_time_B2; mean_small_screen_time_B3]*100)';
rew_cup_zone_time = ([mean_reward_cup_B1; mean_reward_cup_B2; mean_reward_cup_B3]*100)';
other_zone_time = ([other_zone_time_B1; other_zone_time_B2; other_zone_time_B3]*100)';

mean_large_zone = nanmean(large_zone_time, 1);
mean_small_zone = nanmean(small_zone_time, 1);
mean_rew_cup_zone = nanmean(rew_cup_zone_time, 1);
mean_other_zone = nanmean(other_zone_time, 1);


sem_large = nanstd(large_zone_time, 0, 1) ./ sqrt(size(large_zone_time, 1));
sem_small = nanstd(small_zone_time, 0, 1) ./ sqrt(size(small_zone_time, 1));
sem_rew = nanstd(rew_cup_zone_time, 0, 1) ./ sqrt(size(rew_cup_zone_time, 1));
sem_other = nanstd(other_zone_time, 0, 1) ./ sqrt(size(other_zone_time, 1));


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
for i = 1:size(large_zone_time, 1)
    plot(x_points, large_zone_time(i, :), '-', ...
        'Color', [0 0 1 0.6], ... % Blue with 60% opacity
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(rew_cup_zone_time, 1)
    plot(x_points, rew_cup_zone_time(i, :), '-', ...
        'Color', [0 1 0 0.6], ... % 
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(small_zone_time, 1)
    plot(x_points, small_zone_time(i, :), '-', ...
        'Color', [1 0 0 0.6], ... % 
        'LineWidth', 1.2);
end

% Plot individual lines for "Small" data
for i = 1:size(other_zone_time, 1)
    plot(x_points, other_zone_time(i, :), '-', ...
        'Color', 'black', ... % 
        'LineWidth', 1.2);
end


% Plot with error bars for "Large" and "Small"
errorbar(x_points, mean_large_zone, sem_large, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, mean_small_zone, sem_small, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

errorbar(x_points, mean_rew_cup_zone, sem_rew, 'square-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'green', 'MarkerFaceColor', 'green', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'

errorbar(x_points, mean_other_zone, sem_other, 'diamond-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'black', 'MarkerFaceColor', 'black', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'


% Format the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'0', '50%', '75%'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 1.1 * max([mean_large_zone + sem_large, ...
                   mean_small_zone + sem_small])]); % Adjust ylim dynamically
set(gca, 'ytick', 0:25:100);
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;

%%
%%  FOR GENERAL analyses of 1p: Tukey's for multiple comparisons
clc
% Combine data for analysis - now with four zone types
data = [large_zone_time, small_zone_time, rew_cup_zone_time, other_zone_time]; % 10x12 matrix (4 zones × 3 blocks each)

num_subjects = size(data, 1); % Number of subjects (mice)

% Create table for analysis
varNames = {'Large_Block1', 'Large_Block2', 'Large_Block3', ...
            'Small_Block1', 'Small_Block2', 'Small_Block3', ...
            'RewardCup_Block1', 'RewardCup_Block2', 'RewardCup_Block3', ...
            'Other_Block1', 'Other_Block2', 'Other_Block3'};
tbl = array2table(data, 'VariableNames', varNames);

% Define within-subject factors
Zone = categorical({'Large', 'Large', 'Large', 'Small', 'Small', 'Small', ...
                   'RewardCup', 'RewardCup', 'RewardCup', 'Other', 'Other', 'Other'}); % Zone factor
TrialBlock = categorical([1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]); % Trial Block factor
WithinDesign = table(Zone', TrialBlock', 'VariableNames', {'Zone', 'TrialBlock'});

% Fit repeated measures model
rm = fitrm(tbl, 'Large_Block1-Other_Block3 ~ 1', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovaResults = ranova(rm, 'WithinModel', 'Zone*TrialBlock');

mauchlyTest = mauchly(rm)
% Display results
disp(ranovaResults);

% Display ANOVA results
fprintf('Zone effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{3,2}, ranovaResults{4,2}, ranovaResults{3,4}, ranovaResults{3,5});

% Display ANOVA results
fprintf('Trial Block effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{5,2}, ranovaResults{6,2}, ranovaResults{5,4}, ranovaResults{5,5});
% Display ANOVA results
fprintf('Interaction effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{7,2}, ranovaResults{8,2}, ranovaResults{7,4}, ranovaResults{7,5});

% Check if interaction is significant (using p < 0.05 threshold)
pValueInteraction = ranovaResults{7, 5}; % p-value for Zone*TrialBlock interaction
disp(['Interaction p-value: ', num2str(pValueInteraction)]);

% If interaction is significant, conduct post-hoc tests
if pValueInteraction < 0.05
    disp('Significant interaction detected. Conducting Tukey''s multiple comparisons...');
    fprintf('\n========== POST-HOC TESTS ==========\n');
    
    % Reshape data for Tukey's multiple comparisons
    % Create long-format data for all Zone-Block combinations
    all_values = [];
    all_groups = {};
    subject_ids = [];
    
    zone_names = {'Large', 'Small', 'RewardCup', 'Other'};
    for zone = 1:4 % 1 = Large, 2 = Small, 3 = RewardCup, 4 = Other
        for block = 1:3
            block_data = data(:, (zone-1)*3 + block); % Extract appropriate zone and block data
            
            all_values = [all_values; block_data];
            all_groups = [all_groups; repmat({sprintf('%s_Block%d', zone_names{zone}, block)}, length(block_data), 1)];
            subject_ids = [subject_ids; (1:num_subjects)'];
        end
    end
    
    % Use anova1 to get stats structure for multcompare
    [~, ~, stats_tukey] = anova1(all_values, all_groups, 'off');
    
    % Perform Tukey's multiple comparisons
    [c, m, h, gnames] = multcompare(stats_tukey, 'CType', 'tukey-kramer', 'Display', 'off');
    
    % Display Tukey's results in a readable format
    fprintf('\nTukey''s HSD Multiple Comparisons Results:\n');
    fprintf('%-20s %-20s %10s %10s %10s %10s\n', 'Group 1', 'Group 2', 'Diff', 'Lower CI', 'Upper CI', 'p-value');
    fprintf('%-20s %-20s %10s %10s %10s %10s\n', repmat('-', 1, 20), repmat('-', 1, 20), repmat('-', 1, 10), repmat('-', 1, 10), repmat('-', 1, 10), repmat('-', 1, 10));
    
    for i = 1:size(c, 1)
        group1_idx = c(i, 1);
        group2_idx = c(i, 2);
        group1_name = gnames{group1_idx};
        group2_name = gnames{group2_idx};
        diff = c(i, 4);
        lower_ci = c(i, 3);
        upper_ci = c(i, 5);
        p_val = c(i, 6);
        
        fprintf('%-20s %-20s %10.3f %10.3f %10.3f %10.3e', ...
            group1_name, group2_name, diff, lower_ci, upper_ci, p_val);
        
        if p_val < 0.05
            fprintf(' *');
        end
        fprintf('\n');
    end
    
    fprintf('\n* indicates significant difference (p < 0.05)\n');
    
    % Display group means for reference
    fprintf('\nGroup Means:\n');
    for i = 1:length(gnames)
        fprintf('%-20s: %.3f ± %.3f (SEM)\n', gnames{i}, m(i, 1), m(i, 2));
    end
    
    % Highlight specific comparisons of interest
    fprintf('\n--- Key Comparisons ---\n');
    fprintf('Zone comparisons within each block:\n');
    for block = 1:3
        fprintf('\n--- Block %d ---\n', block);
        
        % Compare all zone pairs within each block
        zone_pairs = {{'Large', 'Small'}, {'Large', 'RewardCup'}, {'Large', 'Other'}, ...
                     {'Small', 'RewardCup'}, {'Small', 'Other'}, {'RewardCup', 'Other'}};
        
        for pair_idx = 1:length(zone_pairs)
            zone1_name = sprintf('%s_Block%d', zone_pairs{pair_idx}{1}, block);
            zone2_name = sprintf('%s_Block%d', zone_pairs{pair_idx}{2}, block);
            
            % Find the comparison in the results
            for i = 1:size(c, 1)
                group1_name = gnames{c(i, 1)};
                group2_name = gnames{c(i, 2)};
                
                if (strcmp(group1_name, zone1_name) && strcmp(group2_name, zone2_name)) || ...
                   (strcmp(group1_name, zone2_name) && strcmp(group2_name, zone1_name))
                    diff = c(i, 4);
                    lower_ci = c(i, 3);
                    upper_ci = c(i, 5);
                    p_val = c(i, 6);
                    
                    fprintf('%s vs %s: Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.3e', ...
                        zone_pairs{pair_idx}{1}, zone_pairs{pair_idx}{2}, abs(diff), lower_ci, upper_ci, p_val);
                    
                    if p_val < 0.05
                        fprintf(' *');
                    end
                    fprintf('\n');
                    break;
                end
            end
        end
    end
    
    % Block comparisons within each zone
    for zone_idx = 1:4
        fprintf('\nBlock comparisons within %s Zone:\n', zone_names{zone_idx});
        zone_blocks = {sprintf('%s_Block1', zone_names{zone_idx}), ...
                      sprintf('%s_Block2', zone_names{zone_idx}), ...
                      sprintf('%s_Block3', zone_names{zone_idx})};
        
        for i = 1:length(zone_blocks)
            for j = (i+1):length(zone_blocks)
                % Find the comparison in the results
                for k = 1:size(c, 1)
                    group1_name = gnames{c(k, 1)};
                    group2_name = gnames{c(k, 2)};
                    
                    if (strcmp(group1_name, zone_blocks{i}) && strcmp(group2_name, zone_blocks{j})) || ...
                       (strcmp(group1_name, zone_blocks{j}) && strcmp(group2_name, zone_blocks{i}))
                        diff = c(k, 4);
                        lower_ci = c(k, 3);
                        upper_ci = c(k, 5);
                        p_val = c(k, 6);
                        
                        fprintf('Block %d vs Block %d: Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.3e', ...
                            i, j, abs(diff), lower_ci, upper_ci, p_val);
                        
                        if p_val < 0.05
                            fprintf(' *');
                        end
                        fprintf('\n');
                        break;
                    end
                end
            end
        end
    end
    
else
    disp('Interaction not significant. Post-hoc tests for interaction not needed.');
    
    % If main effects are significant, analyze those
    pValueZone = ranovaResults{3, 5}; % p-value for Zone main effect
    pValueBlock = ranovaResults{5, 5};  % p-value for TrialBlock main effect
    
    % Check main effect of Zone
    if pValueZone < 0.05
        fprintf('\n--- Main effect of Zone is significant (p = %.4f) ---\n', pValueZone);
        % Compare overall means across zones
        largeMean = mean(mean(data(:,1:3)));
        smallMean = mean(mean(data(:,4:6)));
        rewardCupMean = mean(mean(data(:,7:9)));
        otherMean = mean(mean(data(:,10:12)));
        
        fprintf('Large Zone mean: %.3f\n', largeMean);
        fprintf('Small Zone mean: %.3f\n', smallMean);
        fprintf('Reward Cup Zone mean: %.3f\n', rewardCupMean);
        fprintf('Other Zone mean: %.3f\n', otherMean);
        
        % Perform post-hoc comparisons for Zone
        zoneComp = multcompare(rm, 'Zone', 'ComparisonType', 'bonferroni');
        
        % Display results
        fprintf('\nPost-hoc comparisons for Zone (across all trial blocks):\n');
        for i = 1:size(zoneComp, 1)
            fprintf('%s vs %s: Mean Diff = %.3f, CI = [%.3f, %.3f], p = %.3e\n', ...
                zoneComp.Zone_1(i), zoneComp.Zone_2(i), ...
                zoneComp.Difference(i), zoneComp.Lower(i), ...
                zoneComp.Upper(i), zoneComp.pValue(i));
        end
    end
    
    % Check main effect of Trial Block
    if pValueBlock < 0.05
        fprintf('\n--- Main effect of Trial Block is significant (p = %.4f) ---\n', pValueBlock);
        % Run pairwise comparisons for Trial Block
        blockComp = multcompare(rm, 'TrialBlock', 'ComparisonType', 'bonferroni');
        
        % Display results
        fprintf('Post-hoc comparisons for Trial Block (across all zones):\n');
        for i = 1:size(blockComp, 1)
            fprintf('Block %s vs Block %s: Mean Diff = %.3f, CI = [%.3f, %.3f], p = %.3e\n', ...
                blockComp.TrialBlock_1(i), blockComp.TrialBlock_2(i), ...
                blockComp.Difference(i), blockComp.Lower(i), ...
                blockComp.Upper(i), blockComp.pValue(i));
        end
    end
end