session_to_analyze = 'RDT_D1'

if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39', 'BLA_Insc_13'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_behavior, fieldsToRemove{i})
            final_behavior = rmfield(final_behavior, fieldsToRemove{i});
        end
    end
elseif strcmp('RDT_D2', session_to_analyze)

    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_behavior, fieldsToRemove{i})
            final_behavior = rmfield(final_behavior, fieldsToRemove{i});
        end
    end
end
%%
animalIDs = (fieldnames(final_behavior));

current_behav_data = [];

for ii = 1:size(animalIDs, 1)
    currentanimal = char(animalIDs(ii));
    if isfield(final_behavior.(currentanimal), session_to_analyze)
        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
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





        [BehavData,trials,varargin]=TrialFilter_test(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0);

        trial_after_shock_mouse{ii} = BehavData.trial_after_shk;
        block_data_for_mouse{ii} = BehavData.Block;

        initiation_latency = [];
        for dd = 1:size(BehavData, 1)
            if dd == 1
                initiation_latency(dd) = BehavData.stTime(dd) -BehavData.TrialPossible(dd);
            else
                initiation_latency(dd) = BehavData.stTime(dd) - BehavData.choiceTime(dd-1);
            end

        end
        
        initiation_latency_mouse{ii} =  initiation_latency;
        clear initiation_latency


        collect_latency = [];
        for dd = 1:size(BehavData, 1)
    
            collect_latency(dd) = BehavData.collectionTime(dd) - BehavData.choiceTime(dd);



        end
        collect_latency_mouse{ii} =  collect_latency;
        clear collect_latency

    end
end

initiation_times_concat = cat(1, initiation_latency_mouse{:});

figure; plot(mean(initiation_times_concat))

collect_times_concat = cat(1, collect_latency_mouse{:});

figure; plot(mean(collect_times_concat))

mean_initiation_times_concat = mean(initiation_times_concat);

mean_collect_times_concat = mean(collect_times_concat);


for hh = 1:size(trial_after_shock_mouse, 2)
    current_trial_data = trial_after_shock_mouse{1, hh};
    current_init_data = initiation_latency_mouse{1, hh};
    current_collect_data = collect_latency_mouse{1, hh};
    current_block_data = block_data_for_mouse{1, hh};

    mean_initiation_after_shock(hh) = mean(current_init_data(:, current_trial_data == 1));
    mean_initiation_after_no_shock(hh) = mean(current_init_data(:, current_trial_data ~= 1 & current_block_data == 2 | current_block_data == 3));


    mean_collect_after_shock(hh) = mean(current_collect_data(:, current_trial_data == 1));
    mean_collect_after_no_shock(hh) = mean(current_collect_data(:, current_trial_data ~= 1 & current_block_data == 2 | current_block_data == 3));


    mean_initiation_after_no_shock_all_trials(hh) = mean(current_init_data(:, current_trial_data ~= 1));

    mean_collect_after_no_shock_all_trials(hh) = mean(current_collect_data(:, current_trial_data ~= 1));
end


%%
% Define parameters
num_bootstraps = 10000;  % Number of bootstrap samples
alpha = 0.001;            % Significance level (95% CI)

[num_rows, num_cols] = size(initiation_times_concat);
baseline_range = 1:30;
test_range = 31:num_cols;

% Compute the mean baseline across all rows
baseline_data = initiation_times_concat(:, baseline_range); % 10x30
mean_baseline = mean(baseline_data, 1); % Mean across columns for each row

% Bootstrap the mean baseline distribution
boot_baseline = bootstrp(num_bootstraps, @mean, baseline_data);
baseline_CI = prctile(boot_baseline, [alpha/2 * 100, (1 - alpha/2) * 100]);

% Compute mean across all rows for each column
mean_across_rows = mean(initiation_times_concat, 1);

% Identify time points where the mean is outside the baseline CI
sig_diff_timepoints = (mean_across_rows > baseline_CI(2)) | (mean_across_rows < baseline_CI(1));

% Plot results
figure;
hold on;
plot(1:num_cols, mean_across_rows, 'k', 'LineWidth', 1.5); % Mean across rows
yline(baseline_CI(1), 'r--', 'LineWidth', 1.5); % Lower CI boundary
yline(baseline_CI(2), 'r--', 'LineWidth', 1.5); % Upper CI boundary
scatter(find(sig_diff_timepoints), mean_across_rows(sig_diff_timepoints), 50, 'b', 'filled'); % Highlight significant points
xlabel('Time Points');
ylabel('Mean Across Rows');
title('Significant Deviations from Baseline');
legend('Mean Value', 'Lower CI', 'Upper CI', 'Significant Points');
grid on;
hold off;

% Display significant time points
disp('Significant time points (where mean across rows deviates from baseline CI):');
disp(find(sig_diff_timepoints));






figure;
hold on
% Create a histogram for allCorrelations

width = 700; % Width of the figure
height = 200; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([1 90]);

ylim([0 150]);
% Set X-axis ticks
set(gca, 'XTick', [1, 30, 60, 90], 'YTick', [0 75 150]);
shadedErrorBar(1:90, nanmean(initiation_times_concat), nansem(initiation_times_concat), 'lineProps', {'color', 'r'});
scatter(find(sig_diff_timepoints), mean_initiation_times_concat(sig_diff_timepoints), 50, 'b', 'filled'); % Highlight significant points





%%

% Define parameters
num_bootstraps = 10000;  % Number of bootstrap samples
alpha = 0.001;           % Significance level (95% CI)

[num_rows, num_cols] = size(collect_times_concat);
baseline_range = 1:30;
test_range = 31:num_cols;

% Compute the mean baseline across all rows
baseline_data = collect_times_concat(:, baseline_range); % 10x30
mean_baseline = mean(baseline_data, 1); % Mean across columns for each row

% Bootstrap the mean baseline distribution
boot_baseline = bootstrp(num_bootstraps, @mean, baseline_data);
baseline_CI = prctile(boot_baseline, [alpha/2 * 100, (1 - alpha/2) * 100]);

% Compute mean across all rows for each column
mean_across_rows = mean(collect_times_concat, 1);

% Identify time points where the mean is outside the baseline CI
sig_diff_timepoints = (mean_across_rows > baseline_CI(2)) | (mean_across_rows < baseline_CI(1));

% Plot results
figure;
hold on;
plot(1:num_cols, mean_across_rows, 'k', 'LineWidth', 1.5); % Mean across rows
yline(baseline_CI(1), 'r--', 'LineWidth', 1.5); % Lower CI boundary
yline(baseline_CI(2), 'r--', 'LineWidth', 1.5); % Upper CI boundary
scatter(find(sig_diff_timepoints), mean_across_rows(sig_diff_timepoints), 50, 'b', 'filled'); % Highlight significant points
xlabel('Time Points');
ylabel('Mean Across Rows');
title('Significant Deviations from Baseline');
legend('Mean Value', 'Lower CI', 'Upper CI', 'Significant Points');
grid on;
hold off;

% Display significant time points
disp('Significant time points (where mean across rows deviates from baseline CI):');
disp(find(sig_diff_timepoints));





figure;
hold on
% Create a histogram for allCorrelations

width = 700; % Width of the figure
height = 200; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([1 90]);
ylim([0 8]);
% Set X-axis ticks
set(gca, 'XTick', [1, 30, 60, 90], 'YTick', [0 4 8]);
shadedErrorBar(1:90, nanmean(collect_times_concat), nansem(collect_times_concat), 'lineProps', {'color', 'r'});
scatter(find(sig_diff_timepoints), mean_collect_times_concat(sig_diff_timepoints), 50, 'b', 'filled'); % Highlight significant points

%%
% Define parameters
% num_bootstraps = 1000;  % Number of bootstrap samples
% alpha = 0.05;           % Significance level (95% CI)
% 
% [num_rows, num_cols] = size(initiation_times_concat);
% baseline_range = 1:30;
% test_range = 31:num_cols;
% 
% % Preallocate logical array for significance results
% sig_diff_matrix = false(num_rows, num_cols);
% 
% for row = 1:num_rows
%     % Extract baseline data for this row
%     baseline_data = initiation_times_concat(row, baseline_range);
% 
%     % Bootstrap the baseline mean distribution
%     boot_baseline = bootstrp(num_bootstraps, @mean, baseline_data);
%     baseline_CI = prctile(boot_baseline, [alpha/2 * 100, (1 - alpha/2) * 100]);
% 
%     % Compare each test column against the baseline CI
%     for col = test_range
%         test_value = initiation_times_concat(row, col);
% 
%         % Check if the test value is outside the baseline CI
%         if test_value > baseline_CI(2) || test_value < baseline_CI(1)
%             sig_diff_matrix(row, col) = true;
%         end
%     end
% end
% 
% % Display results
% disp('Matrix of significant differences (rows = trials, columns = timepoints):');
% disp(sig_diff_matrix);
% 
% % Optional: Visualize the results as a heatmap
% figure;
% imagesc(sig_diff_matrix);
% colormap('hot');
% colorbar;
% xlabel('Time Points');
% ylabel('Trials (Rows)');
% title('Significant Differences from Baseline');

%%
% Define parameters
num_bootstraps = 1000;  % Number of bootstrap samples
alpha = 0.05;           % Significance level (95% CI)

[num_rows, num_cols] = size(initiation_times_concat);
baseline_range = 1:30;
test_range = 31:num_cols;

% Compute the mean baseline across all rows
baseline_data = initiation_times_concat(:, baseline_range); % 10x30
mean_baseline = mean(baseline_data, 1); % Mean across columns for each row

% Bootstrap the mean baseline distribution
boot_baseline = bootstrp(num_bootstraps, @mean, baseline_data);
baseline_CI = prctile(boot_baseline, [alpha/2 * 100, (1 - alpha/2) * 100]);

% Compute mean across all rows for each column
mean_across_rows = mean(initiation_times_concat, 1);

% Identify time points where the mean is outside the baseline CI
sig_diff_timepoints = (mean_across_rows > baseline_CI(2)) | (mean_across_rows < baseline_CI(1));

% Plot results
figure;
hold on;
plot(1:num_cols, mean_across_rows, 'k', 'LineWidth', 1.5); % Mean across rows
yline(baseline_CI(1), 'r--', 'LineWidth', 1.5); % Lower CI boundary
yline(baseline_CI(2), 'r--', 'LineWidth', 1.5); % Upper CI boundary
scatter(find(sig_diff_timepoints), mean_across_rows(sig_diff_timepoints), 50, 'b', 'filled'); % Highlight significant points
xlabel('Time Points');
ylabel('Mean Across Rows');
title('Significant Deviations from Baseline');
legend('Mean Value', 'Lower CI', 'Upper CI', 'Significant Points');
grid on;
hold off;

% Display significant time points
disp('Significant time points (where mean across rows deviates from baseline CI):');
disp(find(sig_diff_timepoints));


