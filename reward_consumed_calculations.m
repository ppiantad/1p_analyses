% run generate_behav_figs_RDT.m with Pre_RDT_RM first:


large_choice = [risk_table.block_1_large*22, risk_table.block_2_large*22, risk_table.block_3_large*22];
small_choice = [risk_table.block_1_small*22, risk_table.block_2_small*22, risk_table.block_3_small*22];

large_choice_volume = large_choice * 20;
small_choice_volume = small_choice * 5;

sum_large_choice_volume_Late_RM = sum(large_choice_volume');
sum_small_choice_volume_Late_RM = sum(small_choice_volume');
%%
% then run generate_behav_figs_RDT.m with RDT_D1:

large_choice = [risk_table.block_1_large*22, risk_table.block_2_large*22, risk_table.block_3_large*22];
small_choice = [risk_table.block_1_small*22, risk_table.block_2_small*22, risk_table.block_3_small*22];

large_choice_volume = large_choice * 20;
small_choice_volume = small_choice * 5;

sum_large_choice_volume_RDT_D1 = sum(large_choice_volume');
sum_small_choice_volume_RDT_D1 = sum(small_choice_volume');


%%
% Assuming postchoice_over_time_mean_iter is a 1x2 cell array, each containing a 90x10 double

mean_sum_large_choice_volume_Late_RM = mean(sum_large_choice_volume_Late_RM)
mean_sum_small_choice_volume_Late_RM = mean(sum_small_choice_volume_Late_RM)


mean_Late_RM = [mean_sum_large_choice_volume_Late_RM mean_sum_small_choice_volume_Late_RM]

mean_sum_large_choice_volume_RDT_D1 = mean(sum_large_choice_volume_RDT_D1)
mean_sum_small_choice_volume_RDT_D1 = mean(sum_small_choice_volume_RDT_D1)

mean_RDT_D1 = [mean_sum_large_choice_volume_RDT_D1 mean_sum_small_choice_volume_RDT_D1]

individual_means ={sum_large_choice_volume_Late_RM, sum_small_choice_volume_Late_RM; sum_large_choice_volume_RDT_D1, sum_small_choice_volume_RDT_D1}

mean_values = [mean_Late_RM; mean_RDT_D1]

% Plot bar graph
figure; hold on;
width = 200; % Width of the figure
height = 400; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
bar_handle = bar(mean_values, 'grouped');

% Set bar colors
bar_handle(1).FaceColor = 'b'; % First cell (blue)
bar_handle(2).FaceColor = 'r'; % Second cell (red)

% Get X positions for scatter overlay
x_positions = zeros(2,2); % Store bar center positions
for i = 1:2
    x_positions(:,i) = bar_handle(i).XEndPoints;
end

% Overlay individual data points correctly aligned
for i = 1:2
    scatter(x_positions(i,1) * ones(1,10), individual_means{i,1}, 'k', 'filled'); % First dataset
    scatter(x_positions(i,2) * ones(1,10), individual_means{i,2}, 'k', 'filled'); % Second dataset
end

% Labels and legend
xticks(mean(x_positions,2));
xticklabels({'Late RM', 'RDT'}); % Adjusted labels
ylabel('Reward consumed (ul)');
legend({'Large', 'Small'}, 'Location', 'SouthEast');
ytickformat('%.2f');

hold off;

[h_prechoice_pv_mean, p_pv_mean, ~, stats_prechoice_pv] = ttest(individual_means{1, 1}, individual_means{1, 2})

% 2×2 Repeated Measures ANOVA for behavioral task data 
% Factors: Identification (Safe-identified, Risky-identified) × Block (Safe, Risky)
% Using cell array input structure: individual_means (2×2 cell array)
% - Row 1: Safe block, Row 2: Risky block
% - Column 1: Safe-identified, Column 2: Risky-identified

% Step 1: Extract data from the cell array
late_RM_large = individual_means{1,1};    % Safe block, Safe-identified
late_RM_small = individual_means{1,2};   % Safe block, Risky-identified
RDT_large = individual_means{2,1};   % Risky block, Safe-identified
RDT_small = individual_means{2,2};  % Risky block, Risky-identified

% Verify dimensions to ensure all cells have the same number of subjects
num_subjects = length(late_RM_large);
if length(late_RM_small) ~= num_subjects || ...
   length(RDT_large) ~= num_subjects || ...
   length(RDT_small) ~= num_subjects
    error('All cells must contain the same number of subjects');
end

% Step 2: Prepare data in the format required for repeated measures ANOVA
% Create a table with subject IDs and the 4 condition measurements
subject_id = (1:num_subjects)';
data_wide = table(subject_id, ...
                 late_RM_large', RDT_large', ...
                 late_RM_small', RDT_small', ...
                 'VariableNames', {'Subject', 'Late_RM_Large_vol', 'RDT_Large_vol', 'Late_RM_Small_vol', 'RDT_Small_vol'});

% Create the within-subjects design table
% This defines the structure of our repeated measures design
withinDesign = table(categorical({'Late_RM'; 'RDT'; 'Late_RM'; 'RDT'}), ...
                    categorical({'Large'; 'Large'; 'Small'; 'Small'}), ...
                    'VariableNames', {'Stage', 'Size'});

% Step 3: Conduct the 2×2 repeated measures ANOVA
% Define the repeated measures model - all measures with the same between-subjects model (~ 1)
rm = fitrm(data_wide, 'Late_RM_Large_vol,RDT_Large_vol,Late_RM_Small_vol,RDT_Small_vol ~ 1', 'WithinDesign', withinDesign);

% Run the repeated measures ANOVA
ranova_table = ranova(rm, 'WithinModel', 'Stage*Size');

% Display results
disp('Repeated Measures ANOVA Results:');
disp(ranova_table);

% Display ANOVA results
fprintf('Size effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranova_table{5,2}, ranova_table{6,2}, ranova_table{5,4}, ranova_table{5,5});
% Display ANOVA results
fprintf('Stage (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranova_table{3,2}, ranova_table{4,2}, ranova_table{3,4}, ranova_table{3,5});

fprintf('Interaction (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranova_table{7,2}, ranova_table{8,2}, ranova_table{7,4}, ranova_table{7,5});


% Step 6: Simple effects analysis (if needed)
% If there's a significant interaction, you may want to perform simple effects tests
disp('Simple Effects Analysis:');

% Effect of Block for Safe-identified stimuli
[~, p_safe_id, ~, stats_safe_id] = ttest(late_RM_large, RDT_large);
disp('Effect of Large Across Tasks:');
disp(['t(', num2str(num_subjects-1), ') = ', num2str(stats_safe_id.tstat), ...
      ', p = ', num2str(p_safe_id)]);

% Effect of Block for Risky-identified stimuli
[~, p_risky_id, ~, stats_risky_id] = ttest(late_RM_small, RDT_small);
disp('Effect of Small Across Tasks:');
disp(['t(', num2str(num_subjects-1), ') = ', num2str(stats_risky_id.tstat), ...
      ', p = ', num2str(p_risky_id)]);

% Effect of Identification in Safe Block
[~, p_safe_block, ~, stats_safe_block] = ttest(late_RM_large, late_RM_small);
disp('Effect within RM:');
disp(['t(', num2str(num_subjects-1), ') = ', num2str(stats_safe_block.tstat), ...
      ', p = ', num2str(p_safe_block)]);

% Effect of Identification in Risky Block
[~, p_risky_block, ~, stats_risky_block] = ttest(RDT_large, RDT_small);
disp('Effect of within RDT:');
disp(['t(', num2str(num_subjects-1), ') = ', num2str(stats_risky_block.tstat), ...
      ', p = ', num2str(p_risky_block)]);