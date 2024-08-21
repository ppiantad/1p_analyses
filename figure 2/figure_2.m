load('BLA_C_raw_no_additional_filtering_RDT_D1_only_completed_sessions_zall_window_base_workspace_10_categories.mat')
%% for heatmap, change "plot_num" (what neuron to plot) and array_to_plot (# corresponds to which dataset)

plot_num = 46; 

array_to_plot = 1; % depends on the structure of zall

select_mouse = 'BLA_Insc_25';

% for RDT D1 BLA_Insc_25:
%prechoice neuron num 46
%postchoice rew num 38
%consumption num 39
%shock num 11

select_mouse_index = find(strcmp(animalIDs, select_mouse));

first_session = 'RDT_D1';

second_session = 'RDT_D1';



% in order to trim off excess (because calcium recording starts before
% behavior), you'll need the start time. unfortunately I don't save that
% variable in the "final" struct, but you can get it from the adjustment to
% some of the columns of BehavData. e.g., see below - but make sure to
% update the session etc as necessary! 

% BehavData = final.(select_mouse).(first_session).choiceTime.uv.BehavData;
BehavData = final.(select_mouse).(first_session).uv.BehavData;
% because the first trial possible is ALWAYS 60 seconds after ABET is
% issued, you can determine what adjustment has been made to this column
% (adding time to account for calcium recording starting first) by
% subtracting off 60 from the first element
stTime = BehavData.TrialPossible(1)-60; 



time2Collect = BehavData.collectionTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
trialStartTime = BehavData.stTime(trials_per_mouse{select_mouse_index, array_to_plot}) - BehavData.choiceTime(trials_per_mouse{select_mouse_index, array_to_plot});
median_trialStartTime = median(trialStartTime)
median_time2Collect = median(time2Collect)
xline(median_trialStartTime)
xline(median_time2Collect)
[numTrials, ~] = size(time2Collect);
Tris = [1:numTrials]';

% Define the custom colormap from white to orange
% custom_colormap = [
%     1, 1, 1; % white
%     1, 0.9, 0.8;
%     1, 0.8, 0.6;
%     1, 0.7, 0.4;
%     1, 0.6, 0.2;
%     1, 0.5, 0; % orange
% ];


custom_colormap = [
    1, 1, 1;       % white
    0.9, 0.95, 0.9;
    0.8, 0.9, 0.8;
    0.6, 0.8, 0.6;
    0.4, 0.7, 0.4;
    0.2, 0.6, 0.2;
    0.13, 0.55, 0.13; % forest green
];

% custom_colormap = [
%     1, 1, 1;         % white
%     0.9, 0.95, 0.95;
%     0.8, 0.9, 0.9;
%     0.6, 0.85, 0.85;
%     0.4, 0.8, 0.8;
%     0.2, 0.8, 0.8;
%     0.0, 0.8, 0.8;   % robin's egg blue
% ];

% custom_colormap = [
%     1, 1, 1;         % white
%     0.9, 0.9, 0.95;
%     0.8, 0.8, 0.9;
%     0.6, 0.6, 0.8;
%     0.4, 0.4, 0.7;
%     0.2, 0.2, 0.6;
%     0.0, 0.0, 0.55;   % dark blue
% ];

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


% Generate more intermediate colors for a smoother transition
n = 256; % Number of colors
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 250, 600]); % [left, bottom, width, height]

% Create a tiled layout with 2 rows and 1 column
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
ax1 = nexttile;
hold on;

% Plot the heatmap
imagesc(ts1, 1:size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1), zall_mouse{select_mouse_index, array_to_plot}{1, plot_num});


% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-1 1]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar(ax1, 'eastoutside');
set(c, 'YTick', clim); % 

ylim([0.5, size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1) + 0.5]);

xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', []);
set(gca, 'YTick', [1, size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1)]);
xline(0)
scatter(time2Collect, Tris               , 'Marker', 'p')
scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;

% Second tile (mean and raw data)
ax2 = nexttile;
hold on;


% Plot the raw data in grey with transparency
for trial = 1:size(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}, 1)
    plot(ts1, zall_mouse{select_mouse_index, array_to_plot}{1, plot_num}(trial, :), 'Color', [custom_colormap(end, :), 0.5]);
end

% Plot the mean as a thick black line
meanData = mean(zall_mouse{select_mouse_index, array_to_plot}{1, plot_num});
plot(ts1, meanData, 'r', 'LineWidth', 2, 'Color', 'k');

ylim([-4 4]);
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
xline(0)
yline(0)
fontsize(18, 'points')
hold off;


%% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav
% this outputs a ever so slightly wonky diagram. a few nodes that do not
% actually overlap minimally overlap (but intersections are 0), and 1 node
% that has 1 overlap does not overlap at all. 
pre_choice_active_ind = find(respClass_all_array{1,1} == 1);
consum_active_ind = find(respClass_all_array{1,3} == 1);
post_choice_active_ind = find(respClass_all_array{1,2} == 1);
% consum_inhibited_ind = find(all_consum_inhibited == 1);
setListData = {pre_choice_active_ind, consum_active_ind, post_choice_active_ind};
setLabels = ["Pre-choice excited", "Consumption excited", "Post-choice excited"];

h = vennEulerDiagram(setListData, setLabels, 'drawProportional', false);

h.ShowIntersectionCounts = true;
h.ShowIntersectionAreas = true;
% h.SetLabels = [];

%%

%% use these data for mouse x mouse, which is likely better
% load a RDT dataset with the following variables filtered (in order):
    % "choiceTime.Outcome_Minus_4to0.OMITALL_0_BLANK_TOUCH_0_BLOCK_1"
    % "choiceTime.Outcome_0to2.OMITALL_0_BLANK_TOUCH_0_BLOCK_1"
    % "collectionTime.Outcome_1to3.OMITALL_0_BLANK_TOUCH_0_BLOCK_1"
    % "choiceTime.Outcome_0to2.SHK_1"
    % "choiceTime.Outcome_0to2.LOSS_PLUS_ONE_1"

%%
for q = 1:length (behav_tbl_iter{1, 1})
    nestedCellArray_1 = behav_tbl_iter{1, 1}{q};
    % nestedCellArray_2 = behav_tbl_iter{2, 1}{q};
    % if size(nestedCellArray_1, 1) > size(nestedCellArray_2, 1)
    %     delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime(1:end-1,:);
    % else 
    %     delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
    % end
    % 

    for zz = 1:size(nestedCellArray_1, 1)
        valid_start_times = nestedCellArray_1.stTime(2:end);
        valid_choice_times = nestedCellArray_1.choiceTime(1:end-1);
        delay_to_initiation = valid_start_times - valid_choice_times;
    end

    trial_choice_times = nestedCellArray_1.choiceTime - nestedCellArray_1.stTime;
    % delay_to_initiation = nestedCellArray_2.stTime - nestedCellArray_1.choiceTime;
    delay_to_collect_post_shk = nestedCellArray_1.collectionTime - nestedCellArray_1.choiceTime;
    trial_choice_times_by_mouse{q} = trial_choice_times;
    delay_to_initiation_by_mouse{q} = delay_to_initiation;
    delay_to_collect_post_shk_by_mouse{q} = delay_to_collect_post_shk;
    clear trial_choice_times delay_to_initiation delay_to_collect_post_shk



end


variable_to_correlate = trial_choice_times_by_mouse;


%%
array_for_means = 1; 

% Initialize the new cell array to store the mean values
meanZallMouse = cell(size(zall_mouse, 1), 1);

% Define the time range for 0 to 2 seconds
timeRange = (ts1 >= -4) & (ts1 <= 0);
% timeRange = (ts1 >= 0) & (ts1 <= 2);
% timeRange = (ts1 >= 1) & (ts1 <= 3);


% Iterate through each cell in the zall_mouse array
for i = 1:length(zall_mouse)
    % Get the current nested cell array
    nestedCellArray_1 = zall_mouse{i, array_for_means};
    nestedCellArray_2 = zall_mouse{i, 5};
    % Initialize the nested cell array for storing mean values
    meanNestedCellArray = cell(size(nestedCellArray_1));
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(nestedCellArray_1)
        % Get the current double array
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
        
        % Store the mean values in the corresponding cell of the nested cell array
        meanNestedCellArray{j} = meanValues;
    end
    
    % Store the nested cell array of mean values in the corresponding cell of the main cell array
    meanZallMouse{i} = meanNestedCellArray;
end

% Now, meanZallMouse contains the mean activity for each row in the time period 0 to 2 seconds
% Each cell in meanZallMouse contains a nested cell array with the



%% work in progress to do some correlations w/ block 1 activity.

% prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite;
prechoice_blocks_2_and_3 = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;

% postchoice_reward_block_1 = respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
postchoice_reward_block_1 = respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite;
postchoice_reward_blocks_2_and_3 = respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;

% collect_block_1 = respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
collect_block_1 = respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite;
collect_blocks_2_and_3 = respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite;

% shk_neurons = respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} == inhib_or_excite;
shk_neurons = 0;

rest_of_neurons = neuron_num - [sum(prechoice_block_1)+sum(postchoice_reward_block_1)+sum(collect_block_1)+sum(shk_neurons)];
figure; pie([sum(prechoice_block_1), sum(postchoice_reward_block_1), sum(collect_block_1), sum(shk_neurons), rest_of_neurons])

block_1_pre_and_post = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite;
sum(block_1_pre_and_post)
block_1_post_and_consumption = respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite;
sum(block_1_post_and_consumption)
block_1_pre_and_consumption = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite;
sum(block_1_pre_and_consumption)


block_2_and_3_pre_and_post = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} ~= inhib_or_excite;
sum(block_2_and_3_pre_and_post)
block_2_and_3_post_and_consumption = respClass_all_array{1, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite;
sum(block_2_and_3_post_and_consumption)
block_2_and_3_pre_and_consumption = respClass_all_array{1, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(2, 3)} == inhib_or_excite;
sum(block_2_and_3_pre_and_consumption)


% Assuming the following variables are defined:
% meanZallMouse: 14x1 cell array where each cell contains another cell array with mean values
% trial_choice_times_by_mouse: 1x11 cell array containing values to correlate with

% Initialize the new cell array to store the correlation results
correlationResults = cell(size(meanZallMouse));



% Iterate through each level of meanZallMouse
for i = 1:length(meanZallMouse)
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{i};
    
    % Initialize the nested cell array for storing correlation results
    correlationNestedArray = zeros(size(meanNestedCellArray));
    
    % Determine the corresponding index in trial_choice_times_by_mouse
    % Adjust this logic based on how the indices are mapped
    trialIndex = mod(i-1, length(variable_to_correlate)) + 1;
    
    % Get the corresponding trial choice times
    trialChoiceTimes = variable_to_correlate{i};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        % Get the current mean values array
        meanValues = meanNestedCellArray{j};
        
        % Check if trialChoiceTimes has the same length as meanValues
        if length(trialChoiceTimes) == length(meanValues)
            % Compute the correlation
            correlationCoeff = corr(meanValues, trialChoiceTimes(:));
        else
            % If lengths do not match, handle the mismatch (e.g., set correlation to NaN)
            correlationCoeff = NaN;
        end
        
        % Store the correlation coefficient in the nested cell array
        correlationNestedArray(j) = correlationCoeff;
    end
    clear meanValues
    % Store the nested cell array of correlation coefficients in the main cell array
    correlationResults{i} = correlationNestedArray;
end

% Now, correlationResults contains the correlation coefficients for each nested structure in meanZallMouse
%%
% Assuming correlationResults is defined and contains the correlation coefficients

% Initialize an empty array to collect all correlation coefficients
allCorrelations = [];

% Iterate through each level of correlationResults
for i = 1:length(correlationResults)
    % Get the current nested cell array of correlation coefficients
    correlationNestedArray = correlationResults{i};
    
    % Iterate through each cell in the nested cell array
    for j = 1:length(correlationNestedArray)
        % Get the current correlation coefficient
        correlationCoeff = correlationNestedArray(j);
        
        % Check if the coefficient is not NaN (if applicable)
        if ~isnan(correlationCoeff)
            % Append the coefficient to the allCorrelations array
            allCorrelations = [allCorrelations; correlationCoeff];
        end
    end
end

% Now, allCorrelations contains all the correlation coefficients
% Create a histogram of the correlation coefficients
figure;
histogram(allCorrelations);
xlabel('Correlation Coefficient');
ylabel('Frequency');
title('Histogram of Correlation Coefficients');

% Optionally, you can add a vertical line at 0 for reference
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;

%% SHK responsive neurons assumed to be stored in respClass_all_array{1, 1} for this purpose - change as necessary
only_shk_responsive_corrs = allCorrelations(prechoice_block_1==1);
not_shk_responsive_corrs = allCorrelations(prechoice_block_1~=3);
% Now, allCorrelations contains all the correlation coefficients
% Create a histogram of the correlation coefficients
figure;
histogram(only_shk_responsive_corrs);
xlabel('Correlation Coefficient');
ylabel('Frequency');
title('Histogram of Correlation Coefficients');

% Optionally, you can add a vertical line at 0 for reference
hold on;
yLimits = ylim;
plot([0 0], yLimits, 'r--', 'LineWidth', 2);
hold off;

%%
% Assuming the following variables are defined:
% allCorrelations: array containing correlation coefficients from correlationResults
% only_shk_responsive_corrs: array containing correlation coefficients from a different variable

% Calculate means
mean_only_shk = mean(only_shk_responsive_corrs);
mean_not_shk = mean(not_shk_responsive_corrs);

% Create a histogram for allCorrelations
figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
histogram(not_shk_responsive_corrs , 'Normalization', 'probability', 'FaceColor', 'blue','BinWidth', 0.05);
hold on;

% Create a histogram for only_shk_responsive_corrs on the same figure
histogram(only_shk_responsive_corrs, 'Normalization', 'probability', 'FaceColor', 'red', 'BinWidth', 0.05);
xline(mean_only_shk, 'r')
xline(mean_not_shk, 'g')
% Add labels and title
xlabel('Correlation Coefficient');
ylabel('Probability');
% title('Histograms of Correlation Coefficients');
% legend('All Correlations', 'Only SHK Responsive Correlations');

% Optionally, you can add a vertical line at 0 for reference
yLimits = ylim;
plot([0 0], yLimits, 'k', 'LineWidth', 2);
xtickformat('%.2f');
ytickformat('%.2f');
hold off;

% Perform a Kolmogorov-Smirnov test to compare the two distributions
[h, p] = kstest2(not_shk_responsive_corrs , only_shk_responsive_corrs);

% Display the results of the statistical test
fprintf('Kolmogorov-Smirnov test result:\n');
fprintf('h = %d (0 means the null hypothesis cannot be rejected, 1 means it can be rejected)\n', h);
fprintf('p-value = %.4f\n', p);

[h,p,ci,stats] = ttest2(not_shk_responsive_corrs , only_shk_responsive_corrs)


%%
for hh = 1:length(meanZallMouse)
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{hh};
    only_shk_meanNestedCellArray = meanNestedCellArray(:, respClass_all_array_mouse{hh, 4}==1);

    only_shk_meanNestedCellArray_mat = cell2mat(only_shk_meanNestedCellArray);
    mean_only_shk_meanNestedCellArray_mat = mean(only_shk_meanNestedCellArray_mat, 2);
    mean_mean_only_shk_meanNestedCellArray_mat(hh) = mean(mean_only_shk_meanNestedCellArray_mat);

end

scatter(mean_mean_only_shk_meanNestedCellArray_mat, riskiness)
corr(mean_mean_only_shk_meanNestedCellArray_mat', riskiness)


%%





% % Plot means as bars
% figure;
% hold on;
% bar(1, mean_only_shk, 'FaceColor', 'r'); % Red bar for 'only_shk_responsive_corrs'
% bar(2, mean_not_shk, 'FaceColor', 'b'); % Blue bar for 'not_shk_responsive_corrs'
% 
% swarmchart(1, only_shk_responsive_corrs, 5)
% 
% % Scatter individual data points
% scatter(ones(size(only_shk_responsive_corrs)), only_shk_responsive_corrs, 'r', 'filled', 'jitter', 'on', 'jitterAmount', 0.15); % Red points with jitter
% scatter(2 * ones(size(not_shk_responsive_corrs)), not_shk_responsive_corrs, 'b', 'filled', 'jitter', 'on', 'jitterAmount', 0.15); % Blue points with jitter
% 
% % Customize plot
% xlim([0.5, 2.5]);
% xticks([1 2]);
% xticklabels({'Only Shk Responsive', 'Not Shk Responsive'});
% ylabel('Correlation Values');
% title('Correlation Values and Means');
% grid on;
% 
% % Add legend
% legend({'Mean Only Shk', 'Mean Not Shk', 'Individual Only Shk', 'Individual Not Shk'}, 'Location', 'best');
% 
% hold off;

bar_separation_value = 3;

figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
swarmchart(ones(1, length(only_shk_responsive_corrs)), only_shk_responsive_corrs)
hold on
swarmchart(ones(1, length(not_shk_responsive_corrs))*bar_separation_value, not_shk_responsive_corrs)

% yline(mean(only_shk_responsive_corrs), ones(length(only_shk_responsive_corrs)))
plot([0.5; 1.5], [mean(only_shk_responsive_corrs); mean(only_shk_responsive_corrs)], 'LineWidth',3)
plot([bar_separation_value-.5; bar_separation_value+.5], [mean(not_shk_responsive_corrs); mean(not_shk_responsive_corrs)], 'LineWidth',3)
yline(0);
xtickformat('%.1f');
ytickformat('%.1f');
hold off


%% attempting logistic regression to predict large vs. small vs. omit from magnitude of SHK response


% Initialize arrays to store the data
shockResponses = [];
trialChoices = [];

% Iterate through each level of meanZallMouse
for i = 1:1 %1:length(meanZallMouse)
    % Get the current nested cell array of mean values
    meanNestedCellArray = meanZallMouse{i};
    
    % Get the trial choices for the current mouse
    currentTrialChoices = [behav_tbl_iter{2, 1}{i}.ForceFree behav_tbl_iter{2, 1}{i}.bigSmall];
    currentTrialChoices = currentTrialChoices(currentTrialChoices(:,1) ~= 1, 2);
    % Iterate through each cell in the nested cell array
    for j = 1:length(meanNestedCellArray)
        meanValues = meanNestedCellArray{j};
        % Get the current mean values array
        meanValues = meanNestedCellArray{j};
        meanValues = meanValues(currentTrialChoices(:,1) ~= 1);
        
        % Here we use the meanValues as is, no averaging across trials
        % Flatten the meanValues to a single row vector, if needed
        shockResponses = [shockResponses; meanValues];
        
        % Append the corresponding trial choice
        trialChoices = [trialChoices; currentTrialChoices];
    end
end

% Convert trialChoices to categorical if not already
trialChoices = categorical(trialChoices, {'LargeReward', 'SmallReward', 'NoChoice'});

y = trialChoices;
mdl = fitmnr(shockResponses, y);

% Display the model coefficients
disp('Model Coefficients:');
disp(mdl);

%% minor attempts to correlate with MOTION (velocity after shock)

meanNestedCellArray = meanZallMouse{6};

for j = 1:length(meanNestedCellArray)
    % Get the current mean values array
    meanValues = meanNestedCellArray{j};

    % Check if trialChoiceTimes has the same length as meanValues
    if length(meanValues_motion) == length(meanValues)
        % Compute the correlation
        correlationCoeff = corr(meanValues, meanValues_motion(:));
    else
        % If lengths do not match, handle the mismatch (e.g., set correlation to NaN)
        correlationCoeff = NaN;
    end

    % Store the correlation coefficient in the nested cell array
    correlationNestedArray(j) = correlationCoeff;
end

