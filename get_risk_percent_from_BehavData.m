% iter = iter+1;
% neuron_num = 0;
risk_table = table;
for ii = 1:size(animalIDs,1) % 1:size(fieldnames(final),1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        BehavData = final.(currentanimal).(session_to_analyze).uv.BehavData;

        for BehavDataRow = 1:size(BehavData,1)
            if BehavData.shock(BehavDataRow) == 1
                kk = 1;
                while true
                    if (BehavDataRow + kk) > size(BehavData, 1)  % Check if index exceeds the number of rows
                        break;
                    end
                    if ~isnan(BehavData.bigSmall(BehavDataRow + kk)) & BehavData.ForceFree(BehavDataRow + kk) ~= 999
                        BehavData.trial_after_shk(BehavDataRow + kk) = 1;
                        break;
                    else
                        kk = kk + 1;
                    end
                end
            end
            % if BehavDataRow > 1
            %     BehavData.initiation_delay(BehavDataRow+1) = BehavData.stTime(BehavDataRow+1)-BehavData.collectionTime(BehavDataRow); 
            % end

        end

        for BehavDataRow = 1:size(BehavData,1)
            if BehavData.bigSmall(BehavDataRow) ~= 999 & ~isnan(BehavData.bigSmall(BehavDataRow))
                kk = 1;
                while true
                    if (BehavDataRow + kk) > size(BehavData, 1)  % Check if index exceeds the number of rows
                        break;
                    end
                    if ~isnan(BehavData.bigSmall(BehavDataRow + kk)) & BehavData.ForceFree(BehavDataRow + kk) ~= 999
                        BehavData.initiation_delay(BehavDataRow + kk) = BehavData.stTime(BehavDataRow + kk)- BehavData.collectionTime(BehavDataRow);
                        break;
                    else
                        BehavData.initiation_delay(BehavDataRow + kk) = nan;
                        kk = kk + 1;
                    end
                end
            else
                BehavData.initiation_delay(BehavDataRow) = nan;
            end
            % if BehavDataRow > 1
            %     BehavData.initiation_delay(BehavDataRow+1) = BehavData.stTime(BehavDataRow+1)-BehavData.collectionTime(BehavDataRow);
            % end

        end
        
        % only use rewarded trials for this, otherwise things get wonky
        [BehavData,trials,varargin]=TrialFilter_test(BehavData,'ALL', 1); 
        block_1_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 1 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 1); 
        block_1_small_choice_percent = sum(BehavData.bigSmall == 0.3 & BehavData.Block == 1 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 1); 
        % block_1_mouse(ii,:) = [block_1(1, 1) block_1(end, 2)];
        block_2_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 2 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 2);  
        block_2_small_choice_percent = sum(BehavData.bigSmall == 0.3 & BehavData.Block == 2 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 2);  
        % block_2_mouse(ii,:) = [block_2(1, 1) block_2(end, 2)];
        block_3_large_choice_percent = sum(BehavData.bigSmall == 1.2 & BehavData.Block == 3 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);  
        block_3_small_choice_percent = sum(BehavData.bigSmall == 0.3 & BehavData.Block == 3 & BehavData.ForceFree == 0)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3);  
        % block_3_mouse(ii,:) = [block_3(1, 1) block_3(end, 2)];
        lose_shift_percent = sum(BehavData.lose_shift == 1)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3 | BehavData.Block == 2);
        lose_omit_percent = sum(BehavData.lose_omit == 1)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3 | BehavData.Block == 2);
        lose_stay_percent = sum(BehavData.lose_stay == 1)/sum(((BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3) & BehavData.ForceFree == 0) & BehavData.Block == 3 | BehavData.Block == 2);
        win_stay_percent = sum(BehavData.win_stay == 1)/sum(((BehavData.bigSmall == 1.2) & BehavData.ForceFree == 0) & BehavData.Block == 3 | BehavData.Block == 2);
        if ismember('type_binary', BehavData.Properties.VariableNames)
            large_aborts = sum(BehavData.type_binary == 1); %[] sum(BehavData.type_binary == 1)
        else 
            large_aborts = 0;
        end
        trials_completed = sum(BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3);
    end
    risk_table(ii,:) = array2table([block_1_large_choice_percent, block_2_large_choice_percent, block_3_large_choice_percent, block_1_small_choice_percent, block_2_small_choice_percent, block_3_small_choice_percent, large_aborts, lose_shift_percent, lose_omit_percent, lose_stay_percent, win_stay_percent, trials_completed]);
    if ismember('trial_after_shk', BehavData.Properties.VariableNames)
        mean_initiation_latency(ii,:) = [nanmean(BehavData.initiation_delay(BehavData.trial_after_shk == 1)); nanmean(BehavData.initiation_delay(BehavData.trial_after_shk == 0))];
    end
end

% some mice have NaNs if they didn't make it to this trial block. replace
% the NaNs with 0 because not making it to the trial block is basically
% being as risk averse as possible. 
risk_table{:, :}(isnan(risk_table{:, :})) = 0;
row_means = nanmean(risk_table{:, 1:3}, 2);
risk_table.Mean_1_to_3 = row_means;
riskiness = risk_table.Mean_1_to_3;
aborts = risk_table.Var4;
lose_shift = risk_table.Var5;

%%

% Calculate means and SEMs
mean_1_3 = table2array(mean(risk_table(:, 1:3), 1));
sem_1_3 = table2array(std(risk_table(:, 1:3), 0, 1) ./ sqrt(size(risk_table(:, 1:3), 1)));

mean_4_6 = table2array(mean(risk_table(:, 4:6), 1));
sem_4_6 = table2array(std(risk_table(:, 4:6), 0, 1) ./ sqrt(size(risk_table(:, 4:6), 1)));

% X-axis points
x_points = 1:size(mean_1_3, 2);

% Plotting
figure;
hold on;

% Plot lines for risk_table(:, 1:3) and risk_table(:, 4:6)
plot(mean_1_3, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'DisplayName', 'Risk 1:3');
plot(mean_4_6, 's-', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'DisplayName', 'Risk 4:6');

% Add error bars manually using "line"
for i = 1:length(x_points)
    % Error bars for risk_table(:, 1:3)
    line([x_points(i), x_points(i)], [mean_1_3(i) - sem_1_3(i), mean_1_3(i) + sem_1_3(i)], ...
        'Color', 'blue', 'LineWidth', 1.2);

    % Error bars for risk_table(:, 4:6)
    line([x_points(i), x_points(i)], [mean_4_6(i) - sem_4_6(i), mean_4_6(i) + sem_4_6(i)], ...
        'Color', 'red', 'LineWidth', 1.2);
end

% Formatting
ylim([0 1]);
xlabel('Time Points');
ylabel('Mean Value');
legend('Location', 'Best');
title('Mean Risk Values with SEM');
grid on;

%%
mean_large = table2array(mean(risk_table(:, 1:3), 2));
mean_small = table2array(mean(risk_table(:, 4:6), 2));
sem_large = table2array(std(risk_table(:, 1:3), 0, 2) ./ sqrt(size(risk_table(:, 1:3), 2)));
sem_small = table2array(std(risk_table(:, 4:6), 0, 2) ./ sqrt(size(risk_table(:, 4:6), 2)));

% Calculate means for the bar plot
group_means = [mean(mean_large), mean(mean_small)];
group_sems = [mean(sem_large), mean(sem_small)];

% Create a bar plot
figure;
hold on;
bar_handle = bar(group_means, 'FaceAlpha', 0.7, 'BarWidth', 0.5); % Create bar plot with some transparency
errorbar(bar_handle.XEndPoints(1), group_means(1), group_sems(1)); % Create bar plot with some transparency
errorbar(bar_handle.XEndPoints(2), group_means(2), group_sems(2)); % Create bar plot with some transparency


colors = lines(2); % Generate distinct colors for each group

% Add scatter points for each group
x_locations = bar_handle.XEndPoints; % X locations of the bars
scatter_jitter = 0.1; % Jitter width for scatter points

% Scatter points for the first group (mean_large)
scatter(x_locations(1) + (rand(size(mean_large)) - 0.5) * scatter_jitter, mean_large, ...
    50, colors(1, :), 'filled');

% Scatter points for the second group (mean_small)
scatter(x_locations(2) + (rand(size(mean_small)) - 0.5) * scatter_jitter, mean_small, ...
    50, colors(2, :), 'filled');

% Customize plot
set(gca, 'XTick', 1:2, 'XTickLabel', {'Large', 'Small'});
ylabel('Values');
hold off;
%% load Pre_RDT_RM 10 variable dataset, adjust the session_to_analyze = 'RM_D1', run top of script. save risk_table as RM_D1_risk_table.
% set session_to_analyze = 'Pre_RDT_RM', run top of script. then should be
% able to run code below this comment without issue

mean_large_RM_D1 = table2array(mean(RM_D1_risk_table(:, 1:3), 2));
mean_small_RM_D1 = table2array(mean(RM_D1_risk_table(:, 4:6), 2));
sem_large_RM_D1 = table2array(std(RM_D1_risk_table(:, 1:3), 0, 2) ./ sqrt(size(RM_D1_risk_table(:, 1:3), 2)));
sem_small_RM_D1 = table2array(std(RM_D1_risk_table(:, 4:6), 0, 2) ./ sqrt(size(RM_D1_risk_table(:, 4:6), 2)));

mean_large = table2array(mean(risk_table(:, 1:3), 2));
mean_small = table2array(mean(risk_table(:, 4:6), 2));
sem_large = table2array(std(risk_table(:, 1:3), 0, 2) ./ sqrt(size(risk_table(:, 1:3), 2)));
sem_small = table2array(std(risk_table(:, 4:6), 0, 2) ./ sqrt(size(risk_table(:, 4:6), 2)));

% Calculate means for the bar plot
cross_sess_large_means = [mean(mean_large_RM_D1), mean(mean_large)];
cross_sess_large_sems = [mean(sem_large_RM_D1), mean(sem_large)];

% Calculate means for the bar plot
cross_sess_small_means = [mean(mean_small_RM_D1), mean(mean_small)];
cross_sess_small_sems = [mean(sem_small_RM_D1), mean(sem_small)];

% X-axis points
x_points = 1:size(cross_sess_large_means, 2);

% Plotting
figure;
hold on;

% Plot with error bars for "Large" and "Small"
errorbar(x_points, cross_sess_large_means, cross_sess_large_sems, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'blue', 'MarkerFaceColor', 'blue', ...
    'CapSize', 10, 'DisplayName', 'Large'); % Add caps with 'CapSize'

errorbar(x_points, cross_sess_small_means, cross_sess_small_sems, '^-', ...
    'LineWidth', 1.5, 'MarkerSize', 10, 'Color', 'red', 'MarkerFaceColor', 'red', ...
    'CapSize', 10, 'DisplayName', 'Small'); % Add caps with 'CapSize'
% 
% % Add error bars manually using "line"
% for i = 1:length(x_points)
%     % Error bars for "Large"
%     line([x_points(i), x_points(i)], ...
%         [cross_sess_large_means(i) - cross_sess_large_sems(i), cross_sess_large_means(i) + cross_sess_large_sems(i)], ...
%         'Color', 'blue', 'LineWidth', 1.2);
% 
%     % Error bars for "Small"
%     line([x_points(i), x_points(i)], ...
%         [cross_sess_small_means(i) - cross_sess_small_sems(i), cross_sess_small_means(i) + cross_sess_small_sems(i)], ...
%         'Color', 'red', 'LineWidth', 1.2);
% end

% Formatting the X-axis
xticks(x_points); % Set x-ticks at valid x_points
xticklabels({'Early RM', 'Late RM'}); % Provide labels for each x_point
xlim([0.5, length(x_points) + 0.5]); % Add buffer on both sides of x-axis

% Set axis limits, labels, and legend
ylim([0 1.1 * max([cross_sess_large_means + cross_sess_large_sems, ...
                   cross_sess_small_means + cross_sess_small_sems])]); % Adjust ylim dynamically
% xlabel('Condition');
% ylabel('Mean ± SEM');
% legend('Location', 'Best');

% Title and grid for clarity
% title('Cross-Session Risk Analysis');
% grid on;

hold off;