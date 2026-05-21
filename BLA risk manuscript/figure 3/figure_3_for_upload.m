


%% Fig. 3B


shk_mean = mean(neuron_mean_array{1, 4}(:, ts1 > 0 & ts1 < 2),  2);

% [peak_values, time_of_peak_activity] = max(neuron_mean_array{1, 1}, [], 2);
[~, sort_indices] = sort(shk_mean);
neuron_mean_sorted = neuron_mean_array{1, 4}(sort_indices, :);


%
% custom_colormap = [
%     1, 1, 1; % white
%     1, 0.9, 0.9;
%     1, 0.8, 0.8;
%     1, 0.7, 0.7;
%     1, 0.6, 0.6;
%     1, 0.5, 0.5;
%     1, 0.4, 0.4;
%     1, 0.3, 0.3;
%     1, 0.2, 0.2;
%     1, 0.1, 0.1;
%     1, 0, 0;   % red
% ];
% 

custom_colormap = [
    1, 1, 1;   % white
    1, 0.95, 0.9;
    1, 0.9, 0.8;
    1, 0.85, 0.7;
    1, 0.75, 0.55;
    1, 0.65, 0.4;
    1, 0.55, 0.3;
    1, 0.45, 0.2;
    1, 0.35, 0.1;
    1, 0.25, 0.05;
    1, 0.15, 0;  % orange
];


n = 256; 
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));


figure('Position', [100, 100, 350, 600]);
hold on

imagesc(ts1, 1, neuron_mean_sorted);


colormap(custom_colormap);

clim([-1 1]);

c = colorbar('eastoutside');
set(c, 'YTick', clim); % 
ylim([1, neuron_num]);
xlim([-8 8]);
set(gca, 'XTick', [-8, 0, 8]);
set(gca, 'YTick', [1, neuron_num]);
xline(0)
% scatter(time2Collect, Tris               , 'Marker', 'p')
% scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

figure;
hold on


width = 300;
height = 600; 
set(gcf, 'Position', [50, 25, width, height]);
xlim([-8 8]);
ylim([-0.5 1.2]);

set(gca, 'XTick', [-8, 0, 8], 'YTick', [-0.5 0 0.5 1]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(respClass_all_array{1, 4}==1, :)), std(neuron_mean_array{1, 4}(respClass_all_array{1, 4}==1, :))/sqrt(size(neuron_mean_array{1, 4}(respClass_all_array{1, 4}==1, :), 1)), 'lineProps', {'color', 'r'});
xline(0);
hold off



%% Fig. 3D


shk_event = respClass_all_array{1,4} == 1;

consumption_event = respClass_all_array{1, 10} == 1;
shk_and_consum_both_excited = respClass_all_array{1,10} == 1 & respClass_all_array{1,4} == 1;
co_activated_indices = find(shk_and_consum_both_excited(1,:) == 1);
co_activated_indices_sum = numel(co_activated_indices);

only_shk = (sum(shk_event)/neuron_num)*100 - (co_activated_indices_sum/neuron_num)*100; 
only_consumption = (sum(consumption_event)/neuron_num)*100 - (co_activated_indices_sum/neuron_num)*100; 
both = (co_activated_indices_sum/neuron_num)*100; 
not_modulated = 100 - (only_shk + only_consumption + both); 

data_for_bar_plot = [not_modulated, only_consumption, both, only_shk];

figure;
bar(1, data_for_bar_plot, 'stacked'); 
colormap([0.7 0.7 0.7; 0 0 1; 0.8 0.8 0; 1 0 0]); 



ylabel('Responsive neurons (%)');
ylim([0 100]); 
legend({'Other', 'Consum', 'Mixed', 'SHK'}, 'Location', 'eastoutside');


y_offset = 0; 
for i = 1:length(data_for_bar_plot)
    if data_for_bar_plot(i) > 0  %
        text(1, y_offset + data_for_bar_plot(i)/2, sprintf('%.1f%%', data_for_bar_plot(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', 'white', 'FontWeight', 'bold');
        y_offset = y_offset + data_for_bar_plot(i); 
    end
end

for zz = 1:size(respClass_all_array_mouse, 1)
    exclusive_shk_activated_mouse{zz} = respClass_all_array_mouse{zz,4} == 1 & respClass_all_array_mouse{zz,1} == 3 & respClass_all_array_mouse{zz,2} == 3 & respClass_all_array_mouse{zz,3} == 3;
    sum_exclusive_shk_activated_mouse(zz) = sum(exclusive_shk_activated_mouse{zz});
    percent_exclusive_shk_activated_mouse(zz) = sum_exclusive_shk_activated_mouse(zz)/size(exclusive_shk_activated_mouse{zz}, 2);
    max_activity_exclusive_shk(zz) = max(mean(neuron_mean_mouse{zz, 4}(exclusive_shk_activated_mouse{1, zz}   == 1, :)));
end


only_shk_sum = sum(shk_event) - co_activated_indices_sum; % SHK only
only_consumption_sum = sum(consumption_event) -co_activated_indices_sum; % Consum
percent_shk_only = (only_shk_sum/neuron_num)*100
percent_consum_only = (only_consumption_sum/neuron_num)*100

%% Fig. 3E
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
%
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

%params
num_bootstraps = 10000;
alpha = 0.001;      % 95% CI

[num_rows, num_cols] = size(initiation_times_concat);
%for block 1 - 30 trials
baseline_range = 1:30;
test_range = 31:num_cols;
baseline_data = initiation_times_concat(:, baseline_range);
mean_baseline = mean(baseline_data, 1); 

boot_baseline = bootstrp(num_bootstraps, @mean, baseline_data);
baseline_CI = prctile(boot_baseline, [alpha/2 * 100, (1 - alpha/2) * 100]);

mean_across_rows = mean(initiation_times_concat, 1);

sig_diff_timepoints = (mean_across_rows > mean(baseline_CI(2, :))) | (mean_across_rows < mean(baseline_CI(1, :)));


figure;
hold on;
plot(1:num_cols, mean_across_rows, 'k', 'LineWidth', 1.5); 
yline(baseline_CI(1), 'r--', 'LineWidth', 1.5); 
yline(baseline_CI(2), 'r--', 'LineWidth', 1.5); 
scatter(find(sig_diff_timepoints), mean_across_rows(sig_diff_timepoints), 50, 'b', 'filled'); 
xlabel('Trials');
ylabel('Latency (s)');
hold off;

disp('Significant time points (where mean across rows deviates from baseline CI):');
disp(find(sig_diff_timepoints));






figure;
hold on
width = 700;
height = 200; 
set(gcf, 'Position', [50, 25, width, height]); 
xlim([1 90]);

ylim([0 150]);
set(gca, 'XTick', [1, 30, 60, 90], 'YTick', [0 75 150]);
shadedErrorBar(1:90, nanmean(initiation_times_concat), nansem(initiation_times_concat), 'lineProps', {'color', 'r'});
scatter(find(sig_diff_timepoints), mean_initiation_times_concat(sig_diff_timepoints), 50, 'b', 'filled'); 

num_bootstraps = 10000;  
alpha = 0.001;       

[num_rows, num_cols] = size(collect_times_concat);
baseline_range = 1:30;
test_range = 31:num_cols;
baseline_data = collect_times_concat(:, baseline_range); 
mean_baseline = mean(baseline_data, 1); 
boot_baseline = bootstrp(num_bootstraps, @mean, baseline_data);
baseline_CI = prctile(boot_baseline, [alpha/2 * 100, (1 - alpha/2) * 100]);


mean_across_rows = mean(collect_times_concat, 1);

sig_diff_timepoints = (mean_across_rows > mean(baseline_CI(2, :))) | (mean_across_rows < mean(baseline_CI(1, :)));

% Plot results
figure;
hold on;
plot(1:num_cols, mean_across_rows, 'k', 'LineWidth', 1.5); 
yline(baseline_CI(1), 'r--', 'LineWidth', 1.5); 
yline(baseline_CI(2), 'r--', 'LineWidth', 1.5); 
scatter(find(sig_diff_timepoints), mean_across_rows(sig_diff_timepoints), 50, 'b', 'filled'); 
xlabel('Trials');
ylabel('Latency (s)');
hold off;

disp('Significant time points (where mean across rows deviates from baseline CI):');
disp(find(sig_diff_timepoints));





figure;
hold on

width = 700;
height = 200; 
set(gcf, 'Position', [50, 25, width, height]);
xlim([1 90]);
ylim([0 8]);
% Set X-axis ticks
set(gca, 'XTick', [1, 30, 60, 90], 'YTick', [0 4 8]);
shadedErrorBar(1:90, nanmean(collect_times_concat), nansem(collect_times_concat), 'lineProps', {'color', 'r'});
scatter(find(sig_diff_timepoints), mean_collect_times_concat(sig_diff_timepoints), 50, 'b', 'filled'); 

num_bootstraps = 1000; 
alpha = 0.05;   

[num_rows, num_cols] = size(initiation_times_concat);
baseline_range = 1:30;
test_range = 31:num_cols;
baseline_data = initiation_times_concat(:, baseline_range); % 
mean_baseline = mean(baseline_data, 1); 
boot_baseline = bootstrp(num_bootstraps, @mean, baseline_data);
baseline_CI = prctile(boot_baseline, [alpha/2 * 100, (1 - alpha/2) * 100]);

mean_across_rows = mean(initiation_times_concat, 1);
sig_diff_timepoints = (mean_across_rows > mean(baseline_CI(2, :))) | (mean_across_rows < mean(baseline_CI(1, :)));


figure;
hold on;
plot(1:num_cols, mean_across_rows, 'k', 'LineWidth', 1.5);
yline(baseline_CI(1), 'r--', 'LineWidth', 1.5); 
yline(baseline_CI(2), 'r--', 'LineWidth', 1.5); %
scatter(find(sig_diff_timepoints), mean_across_rows(sig_diff_timepoints), 50, 'b', 'filled'); 
xlabel('Trial');
ylabel('Latency (s)');
hold off;

disp('Significant time points (where mean across rows deviates from baseline CI):');
disp(find(sig_diff_timepoints));

%% Fig. 3F
%
for q = 1:length (behav_tbl_iter{4, 1})
    nestedCellArray_1 = behav_tbl_iter{4, 1}{q};
    if ~isempty(nestedCellArray_1)
        nestedCellArray_2 = behav_tbl_iter{4, 1}{q};
        if size(nestedCellArray_1, 1) > size(nestedCellArray_2, 1)
            delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime(1:end-1,:);
        else
            delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
        end
        for zz = 1:size(nestedCellArray_1, 1)
            valid_start_times = nestedCellArray_1.stTime(2:end);
            valid_choice_times = nestedCellArray_1.choiceTime(1:end-1);
            trial_types = nestedCellArray_1.bigSmall;
            trial_types_2 = nestedCellArray_2.bigSmall;
        end

        trial_choice_times = nestedCellArray_1.choiceTime - nestedCellArray_1.stTime;
        delay_to_collect_post_shk = nestedCellArray_1.collectionTime - nestedCellArray_1.choiceTime;
        trial_choice_times_by_mouse{q} = trial_choice_times;
        delay_to_initiation_by_mouse{q} = delay_to_initiation;
        delay_to_collect_post_shk_by_mouse{q} = delay_to_collect_post_shk;
        trial_types_by_mouse{q} = trial_types;
        trial_types_second_var_by_mouse{q} = trial_types_2;
        clear trial_choice_times delay_to_initiation delay_to_collect_post_shk trial_types trial_types_2
    end


end

trial_choice_times_concat = cat(1, trial_choice_times_by_mouse{:});
rew_collect_times_concat = cat(1, delay_to_collect_post_shk_by_mouse{:});
delay_to_initiation_concat = cat(1, delay_to_initiation_by_mouse{:});
delay_to_collect_concat = cat(1, delay_to_collect_post_shk_by_mouse{:});

bar_separation_value = 3;

figure;
width = 100; 
height = 500; 
set(gcf, 'Position', [50, 25, width, height]);

swarmchart(ones(1, length(delay_to_collect_concat)), delay_to_collect_concat);
hold on;
yline(mean(delay_to_collect_concat), 'k', 'LineWidth', 2); 

figure;
width = 100; 
height = 500; 
set(gcf, 'Position', [50, 25, width, height]); 

swarmchart(ones(1, length(delay_to_initiation_concat)) * bar_separation_value, delay_to_initiation_concat);
yticks([0 50 100 150 200])
hold on;
yline(mean(delay_to_initiation_concat), 'k', 'LineWidth', 2); 

yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off;


variable_to_correlate = delay_to_collect_post_shk_by_mouse;

% select data
array_for_means = 4; 
meanZallMouse = cell(size(zall_mouse, 2), 1);

% select appropriate time window
% timeRange = (ts1 >= -4) & (ts1 <= 0);
timeRange = (ts1 >= 0) & (ts1 <= 5); % use a longer time window here so that the entirety of the shk activity can be used
% timeRange = (ts1 >= 1) & (ts1 <= 3);

for i = 1:length(zall_mouse)-1
    nestedCellArray_1 = zall_mouse{i, array_for_means};
    nestedCellArray_2 = zall_mouse{i, 2};
    meanNestedCellArray = cell(size(nestedCellArray_1));
    for j = 1:length(nestedCellArray_1)
        currentArray = nestedCellArray_1{j};
        comparisonArray_for_size = nestedCellArray_2{j};
        
        if isequal(variable_to_correlate, delay_to_collect_post_shk_by_mouse)
            currentArray = nestedCellArray_1{j};
        else
            if size(currentArray, 1) > size(comparisonArray_for_size, 1)
                currentArray = currentArray(1:end-1,:);
            else

            end
        end
        % uncomment below if you want to mean center
        % currentArray_mean = mean(currentArray, 2);
        % currentArray = currentArray-currentArray_mean;
        % Compute the mean activity for each row in the time range 0 to 2 seconds
        meanValues = mean(currentArray(:, timeRange), 2);
        % meanValues = max(currentArray(:, timeRange), [], 2);

        meanNestedCellArray{j} = meanValues;
    end
    meanZallMouse{i} = meanNestedCellArray;
end

correlationResults = cell(size(meanZallMouse));
correlationResults_sig = cell(size(meanZallMouse));

for i = 1:length(meanZallMouse)-1
    meanNestedCellArray = meanZallMouse{i};
    
    correlationNestedArray = zeros(size(meanNestedCellArray));
    corr_sig_NestedArray = zeros(size(meanNestedCellArray));
    trialIndex = mod(i-1, length(variable_to_correlate)) + 1;
    trialChoiceTimes = variable_to_correlate{i};
    for j = 1:length(meanNestedCellArray)
        meanValues = meanNestedCellArray{j};
        if length(trialChoiceTimes) == length(meanValues)
            [correlationCoeff, corr_sig_vals] = corr(meanValues, trialChoiceTimes(:));
        elseif length(trialChoiceTimes) < length(meanValues)
            [correlationCoeff, corr_sig_vals] = corr(meanValues(1:end-1), trialChoiceTimes(:));
        else
            [correlationCoeff, corr_sig_vals] = NaN;
        end
        correlationNestedArray(j) = correlationCoeff;
        corr_sig_NestedArray(j) = corr_sig_vals;
    end
    clear meanValues
    correlationResults{i} = correlationNestedArray;
    correlationResults_sig{i} = corr_sig_NestedArray;
end

allCorrelations = [];

for i = 1:length(correlationResults)
    correlationNestedArray = correlationResults{i};
    for j = 1:length(correlationNestedArray)
        correlationCoeff = correlationNestedArray(j);

        if ~isnan(correlationCoeff)
            allCorrelations = [allCorrelations; correlationCoeff];
        end
    end
end

figure;
histogram(allCorrelations);
xlabel('Correlation coefficient');
ylabel('Frequency');

hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;


%
% change as necessary
only_shk_responsive_corrs = allCorrelations(respClass_all_array{1, 4}  ==1);
not_shk_responsive_corrs = allCorrelations(respClass_all_array{1, 4}  ==3);
figure;
histogram(only_shk_responsive_corrs);
xlabel('Correlation coefficient');
ylabel('Frequency');

hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;
%

mean_only_shk = mean(only_shk_responsive_corrs);
mean_not_shk = mean(not_shk_responsive_corrs);


figure;
width = 250; 
height = 250; 
set(gcf, 'Position', [50, 25, width, height]); 
histogram(not_shk_responsive_corrs , 'Normalization', 'probability', 'FaceColor', 'blue','BinWidth', 0.05,'LineStyle','none');
hold on;
histogram(only_shk_responsive_corrs, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05, 'LineStyle','none');
xline(mean_only_shk, 'r')
xline(mean_not_shk, 'g')

xlabel('Correlation coefficient');
ylabel('Probability');

yLimits = ylim;
plot([0 0], yLimits, 'k', 'LineWidth', 2);
xtickformat('%.2f');
ytickformat('%.2f');
hold off;

% stats
[h, p, k] = kstest2(not_shk_responsive_corrs , only_shk_responsive_corrs)

[h,p,ci,stats] = ttest2(not_shk_responsive_corrs , only_shk_responsive_corrs)


bar_separation_value = 3;

figure;
width = 250; 
height = 250; 
set(gcf, 'Position', [50, 25, width, height]); 
swarmchart(ones(1, length(only_shk_responsive_corrs)), only_shk_responsive_corrs)
hold on
swarmchart(ones(1, length(not_shk_responsive_corrs))*bar_separation_value, not_shk_responsive_corrs)

plot([0.5; 1.5], [mean(only_shk_responsive_corrs); mean(only_shk_responsive_corrs)], 'LineWidth',3)
plot([bar_separation_value-.5; bar_separation_value+.5], [mean(not_shk_responsive_corrs); mean(not_shk_responsive_corrs)], 'LineWidth',3)
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off




%shk representative:
find(correlationResults{5, 1} > 0.5)
start_time = 0;% sub-window start time
end_time = 5; % sub-window end time
sub_window_idx = ts1 >= start_time & ts1 <= end_time;
sub_window_activity_session_1 = zall_mouse{5, 4}{1, 80}(:, sub_window_idx);
r_val_for_representative = correlationResults{5, 1}(1, 80)
p_val_for_representative = correlationResults_sig{5, 1}(1, 80)
choice_times_mouse = variable_to_correlate{1, 5};
trial_types = trial_types_second_var_by_mouse{1, 5};

mean_sub_window_activity_session_1 = mean(sub_window_activity_session_1, 2);


x = mean_sub_window_activity_session_1;
y = choice_times_mouse;

colors = repmat([0.5, 0.5, 0.5], length(trial_types), 1); % Default to gray
colors(trial_types == 1.2, :) = repmat([0, 0, 1], sum(trial_types == 1.2), 1); % Blue for trial_types == 1.2
colors(trial_types == 0.3, :) = repmat([1, 0, 0], sum(trial_types == 0.3), 1); % Red for trial_types == 0.3

figure;
set(gcf, 'Position', [100, 100, 200, 200]); % Adjust figure position and size
scatter(x, y, 36, colors, 'filled', 'MarkerEdgeColor', 'k'); % Use 'colors' for MarkerFaceColor

hold on;

coefficients = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);
plot(x_fit, y_fit, 'r');

y_pred = polyval(coefficients, x);
ssr = sum((y_pred - mean(y)).^2);
sst = sum((y - mean(y)).^2);
r_squared = ssr / sst;

text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared)], 'FontSize', 12);

xlabel('Mean dF/F');
ylabel('Latency (s)');
hold off;

%% Fig. 3G

AA_large_data = zall_mean_all_array{1, 12}(respClass_all_array{1, 11}==1, :);
AA_small_data = zall_mean_all_array{1, 13}(respClass_all_array{1, 11}==1, :);

AA_large_data_sems = sem_all_array{1, 12}(respClass_all_array{1, 11}==1, :);
AA_small_data_sems = sem_all_array{1, 13}(respClass_all_array{1, 11}==1, :);


AA_small_no_trials = find(isnan(AA_small_data(:, 1)));

AA_large_data(AA_small_no_trials, :) = [];
AA_small_data(AA_small_no_trials, :) = [];

AA_large_data_sems(AA_small_no_trials, :) = [];
AA_small_data_sems(AA_small_no_trials, :) = [];

mean_data_array = {AA_large_data, AA_small_data};
sem_data_array = {AA_large_data_sems, AA_small_data_sems};

[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-1 4], [-0.2 0.4], 3);

%% Fig. 3H 
% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav


pre_choice_active_block_2_3_ind = find(respClass_all_array{1,8} == 1);
aa_active_ind = find(respClass_all_array{1,11} == 1);
% consum_inhibited_ind = find(all_consum_inhibited == 1);
setListData = {pre_choice_active_block_2_3_ind, aa_active_ind};
setLabels = ["pre_choice_active_ind", "Approach-Abort excited"];
figure;
ve_diagram = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

ve_diagram.ShowIntersectionCounts = true;
ve_diagram.ShowIntersectionAreas = true;


%% Fig. 3I-J see alternate code

%% Fig. K see alternate code

%% Fig. 3L
uv.evtWin = [-8 8]; %what time do you want to look at around each event [-2 8] [-10 5]
% uv.BLper = [-10 -5];
uv.BLper = [uv.evtWin(1) uv.evtWin(1)+3];
uv.dt = 0.1; 


ca_data_type = "C_raw"; % C % C_raw %S
% CNMFe_data.C_raw: CNMFe traces
% CNMFe_data.C: denoised CNMFe traces
% CNMFe_data.S: inferred spikes
% CNMFe_data.spike_prob: CASCADE inferred spikes - multiply x sampling rate
% (10) for spike rate

session_to_analyze = 'RDT_D1';
uv.yoke_data = 0; % set to 1 if you want to be prompted to yoke the number of trials analyzed, set to 0 otherwise

epoc_to_align = 'choiceTime'; % stTime choiceTime collectionTime
period_of_interest = 'postchoice';

if strcmp(epoc_to_align, 'stTime')
    period_of_interest = 'trial_start';
    uv.evtSigWin.outcome = [-1 1]; %for trial start
elseif strcmp(epoc_to_align, 'choiceTime')
    if strcmp(period_of_interest, 'prechoice')
        uv.evtSigWin.outcome = [-4 0]; %for pre-choice   [-4 0]    [-4 1]
    elseif strcmp(period_of_interest, 'postchoice')
        uv.evtSigWin.outcome = [0 2]; %for SHK or immediate post-choice [0 2]
    end
elseif strcmp(epoc_to_align, 'collectionTime')
    period_of_interest = 'reward_collection';
    uv.evtSigWin.outcome = [1 3]; %for REW collection [1 3]
end

ts1 = (uv.evtWin(1):.1:uv.evtWin(2)-0.1);
animalIDs = (fieldnames(final));
neuron_num = 0;
use_normalized_time = 0;
clear neuron_mean neuron_sem neuron_num zall_mean zall_to_BL_array zsd_array trials ii neuron_mean_unnorm_concat neuron_mean_unnormalized sem_all zall_mean_all 

if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final, fieldsToRemove{i})
            final = rmfield(final, fieldsToRemove{i});
        end
    end
elseif strcmp('RDT_D2', session_to_analyze)

    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final, fieldsToRemove{i})
            final = rmfield(final, fieldsToRemove{i});
        end
    end
end

% get trials where shocks were followed by abort (get first abort) or choice (any choice, large/small)

neuron_num = 0;
animalIDs = (fieldnames(final));
for ii = 1:size(animalIDs, 1)
    currentanimal = char(animalIDs(ii));
    if isfield(final.(currentanimal), session_to_analyze)
        data_table = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData  ;

        total_shock_abort_trials = 0;
        total_shock_no_abort_trials = 0;
        rows_for_new_table = 1;
        i = 1;
        while i <= size(data_table, 1) - 1
            if data_table.shock(i) == 1
                match_found = false;
                behav_data_extracted(rows_for_new_table, :) = data_table(i, :);
                rows_for_new_table = rows_for_new_table + 1;
                j = i + 1;
                while j <= size(data_table, 1) && ~match_found
                    if data_table.type_binary(j) == 1 || data_table.type_binary(j) == 2
                        total_shock_abort_trials = total_shock_abort_trials + 1;
                        match_found = true;
                        behav_data_extracted(rows_for_new_table, :) = data_table(j, :);
                        rows_for_new_table = rows_for_new_table + 1;
                    elseif data_table.bigSmall(j) == 1.2 || data_table.bigSmall(j) == 0.3
                        total_shock_no_abort_trials = total_shock_no_abort_trials + 1;
                        match_found = true;
                        behav_data_extracted(rows_for_new_table, :) = data_table(j, :);
                        rows_for_new_table = rows_for_new_table + 1;
                    else
                        j = j + 1;
                    end
                end
                i = j;
            else
                i = i + 1;
            end
        end

        disp(['Total number of trials where a shock was followed by an abort: ', num2str(total_shock_abort_trials)]);
        disp(['Total number of trials where a shock was followed by a bigSmall event (no abort): ', num2str(total_shock_no_abort_trials)]);
        total_shock_abort_trials_array(ii) = total_shock_abort_trials;
        total_shock_no_abort_trials_array(ii) = total_shock_no_abort_trials;
        behav_data_extracted_array{ii} = behav_data_extracted;
        clear behav_data_extracted
    end
end

shk_abort_to_shk_choice_ratio = total_shock_abort_trials_array./total_shock_no_abort_trials_array

% Plot https://stackoverflow.com/questions/54528239/boxplot-for-paired-observations

combined_data = [total_shock_abort_trials_array ; total_shock_no_abort_trials_array]';
figure();
coordLineStyle = 'k.';
boxplot(combined_data, 'Notch', 'off', 'Symbol', coordLineStyle); hold on;
parallelcoords(combined_data, 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
  'Marker', '.', 'MarkerSize', 10);

TF = isoutlier(combined_data, 'grubbs');

[h,p,ci,stats] = ttest(combined_data(:, 1), combined_data(:, 2))

%% Fig. 3M
% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav

shk_ind = find(respClass_all_array{1,4} == 1);
consum_active_block_2_3 = find(respClass_all_array{1,10} == 1);
aa_active_ind = find(respClass_all_array{1,11} == 1);
setListData = {shk_ind, aa_active_ind, consum_active_block_2_3};
setLabels = ["SHK", "Approach-Abort", "Consum"];
figure;
ve_diagram = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

ve_diagram.ShowIntersectionCounts = true;
ve_diagram.ShowIntersectionAreas = true;
% h.SetLabels = [];


% chi-square test of proportions to see if more the proprtion of AA that are also shock cells is greater than consumption that are also shock
% based on https://www.mathworks.com/matlabcentral/answers/96572-how-can-i-perform-a-chi-square-test-to-determine-how-statistically-different-two-proportions-are-in
aa_and_shk = respClass_all_array{1,11} == 1 & respClass_all_array{1,4} == 1 & respClass_all_array{1,10} ~= 1;
aa_and_shk_sum = sum(aa_and_shk)
aa_not_shk = respClass_all_array{1,11} == 1 & respClass_all_array{1,4} ~= 1 & respClass_all_array{1,10} ~= 1;
aa_not_shk_sum = sum(aa_not_shk)

% using block 1 consumption neurons
consumption_and_shk = respClass_all_array{1,10} == 1 & respClass_all_array{1,4} == 1 & respClass_all_array{1,11} ~= 1;
consumption_and_shk_sum = sum(consumption_and_shk)
consumption_not_shk = respClass_all_array{1,10} == 1 & respClass_all_array{1,4} ~= 1 & respClass_all_array{1,11} ~= 1;
consumption_not_shk_sum = sum(consumption_not_shk)


n1 = aa_and_shk_sum; N1 = aa_and_shk_sum+aa_not_shk_sum;

n2 = consumption_and_shk_sum; N2 = consumption_and_shk_sum+consumption_not_shk_sum;


% pooled estimate of proportion

p0 = (n1+n2) / (N1+N2);

% H0 expected(null hypothesis)

n10 = N1 * p0;

n20 = N2 * p0;

% chi-square test

observed = [n1 N1-n1 n2 N2-n2];

expected = [n10 N1-n10 n20 N2-n20];

chi2stat = sum((observed-expected).^2 ./ expected)

p = 1 - chi2cdf(chi2stat,1)


%% Fig. 3N
mean_data_array = {neuron_mean_array{1, 4}(respClass_all_array{1, 11} == 1, :), neuron_mean_array{1, 4}(respClass_all_array{1, 11} == 3, :)};
sem_data_array = {neuron_sem_array{1, 4}(respClass_all_array{1, 11} == 1, :), neuron_sem_array{1, 4}(respClass_all_array{1, 11} == 3, :)};

[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-2 3], [-0.25 0.60], 3);

%% Fig. 3O
mean_data_array = {neuron_mean_array{1, 11}(respClass_all_array{1, 4} == 1, :), neuron_mean_array{1, 11}(respClass_all_array{1, 4} == 3, :)};
sem_data_array = {neuron_sem_array{1, 11}(respClass_all_array{1, 4} == 1, :), neuron_sem_array{1, 11}(respClass_all_array{1, 4} == 3, :)};

[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, [-2 3], [-0.1 0.10], 3);





