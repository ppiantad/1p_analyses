%%  FOR GENERAL analyses of 1p: Tukey's for multiple comparisons
clc
% Combine data for analysis
% data = [large_choice, small_choice]; % Combine large and small choices
data = [large_path_length, small_path_length]; % Combine large and small choices

num_subjects = size(data, 1); % Number of subjects (mice)

% Create table for analysis
varNames = {'Large_Block1', 'Large_Block2', 'Large_Block3', ...
            'Small_Block1', 'Small_Block2', 'Small_Block3'};
tbl = array2table(data, 'VariableNames', varNames);

% Define within-subject factors
Choice = categorical({'Large', 'Large', 'Large', 'Small', 'Small', 'Small'}); % Choice factor
TrialBlock = categorical([1, 2, 3, 1, 2, 3]); % Trial Block factor
WithinDesign = table(Choice', TrialBlock', 'VariableNames', {'Choice', 'TrialBlock'});

% Fit repeated measures model
rm = fitrm(tbl, 'Large_Block1-Small_Block3 ~ 1', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovaResults = ranova(rm, 'WithinModel', 'Choice*TrialBlock');

mauchlyTest = mauchly(rm)
% Display results
disp(ranovaResults);

% Display ANOVA results
fprintf('Reward Size (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{3,2}, ranovaResults{4,2}, ranovaResults{3,4}, ranovaResults{3,5});

% Display ANOVA results
fprintf('Trial Block effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{5,2}, ranovaResults{6,2}, ranovaResults{5,4}, ranovaResults{5,5});
% Display ANOVA results
fprintf('Interaction effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{7,2}, ranovaResults{8,2}, ranovaResults{7,4}, ranovaResults{7,5});

% Check if interaction is significant (using p < 0.05 threshold)
pValueInteraction = ranovaResults{7, 5}; % p-value for Choice*TrialBlock interaction
disp(['Interaction p-value: ', num2str(pValueInteraction)]);

% If interaction is significant, conduct post-hoc tests
if pValueInteraction < 0.05
    disp('Significant interaction detected. Conducting Tukey''s multiple comparisons...');
    fprintf('\n========== POST-HOC TESTS ==========\n');
    
    % Reshape data for Tukey's multiple comparisons
    % Create long-format data for all Choice-Block combinations
    all_values = [];
    all_groups = {};
    subject_ids = [];
    
    for choice = 1:2 % 1 = Large, 2 = Small
        choice_name = {'Large', 'Small'};
        for block = 1:3
            if choice == 1
                block_data = data(:, block); % Large choice data
            else
                block_data = data(:, block + 3); % Small choice data
            end
            
            all_values = [all_values; block_data];
            all_groups = [all_groups; repmat({sprintf('%s_Block%d', choice_name{choice}, block)}, length(block_data), 1)];
            subject_ids = [subject_ids; (1:num_subjects)'];
        end
    end
    
    % Use anova1 to get stats structure for multcompare
    [~, ~, stats_tukey] = anova1(all_values, all_groups, 'off');
    
    % Perform Tukey's multiple comparisons
    [c, m, h, gnames] = multcompare(stats_tukey, 'CType', 'tukey-kramer', 'Display', 'off');
    
    % Display Tukey's results in a readable format
    fprintf('\nTukey''s HSD Multiple Comparisons Results:\n');
    fprintf('%-15s %-15s %10s %10s %10s %10s\n', 'Group 1', 'Group 2', 'Diff', 'Lower CI', 'Upper CI', 'p-value');
    fprintf('%-15s %-15s %10s %10s %10s %10s\n', repmat('-', 1, 15), repmat('-', 1, 15), repmat('-', 1, 10), repmat('-', 1, 10), repmat('-', 1, 10), repmat('-', 1, 10));
    
    for i = 1:size(c, 1)
        group1_idx = c(i, 1);
        group2_idx = c(i, 2);
        group1_name = gnames{group1_idx};
        group2_name = gnames{group2_idx};
        diff = c(i, 4);
        lower_ci = c(i, 3);
        upper_ci = c(i, 5);
        p_val = c(i, 6);
        
        fprintf('%-15s %-15s %10.3f %10.3f %10.3f %10.3e', ...
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
        fprintf('%-15s: %.3f ± %.3f (SEM)\n', gnames{i}, m(i, 1), m(i, 2));
    end
    
    % Highlight specific comparisons of interest
    fprintf('\n--- Key Comparisons ---\n');
    fprintf('Choice comparisons within each block:\n');
    for block = 1:3
        large_name = sprintf('Large_Block%d', block);
        small_name = sprintf('Small_Block%d', block);
        
        % Find the comparison in the results
        for i = 1:size(c, 1)
            group1_name = gnames{c(i, 1)};
            group2_name = gnames{c(i, 2)};
            
            if (strcmp(group1_name, large_name) && strcmp(group2_name, small_name)) || ...
               (strcmp(group1_name, small_name) && strcmp(group2_name, large_name))
                diff = c(i, 4);
                lower_ci = c(i, 3);
                upper_ci = c(i, 5);
                p_val = c(i, 6);
                
                fprintf('Block %d (Large vs Small): Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.3e', ...
                    block, abs(diff), lower_ci, upper_ci, p_val);
                
                if p_val < 0.05
                    fprintf(' *');
                end
                fprintf('\n');
                break;
            end
        end
    end
    
    fprintf('\nBlock comparisons within Large Choice:\n');
    large_blocks = {'Large_Block1', 'Large_Block2', 'Large_Block3'};
    for i = 1:length(large_blocks)
        for j = (i+1):length(large_blocks)
            % Find the comparison in the results
            for k = 1:size(c, 1)
                group1_name = gnames{c(k, 1)};
                group2_name = gnames{c(k, 2)};
                
                if (strcmp(group1_name, large_blocks{i}) && strcmp(group2_name, large_blocks{j})) || ...
                   (strcmp(group1_name, large_blocks{j}) && strcmp(group2_name, large_blocks{i}))
                    diff = c(k, 4);
                    lower_ci = c(k, 3);
                    upper_ci = c(k, 5);
                    p_val = c(k, 6);
                    
                    fprintf('Large Block %d vs Block %d: Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.3e', ...
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
    
    fprintf('\nBlock comparisons within Small Choice:\n');
    small_blocks = {'Small_Block1', 'Small_Block2', 'Small_Block3'};
    for i = 1:length(small_blocks)
        for j = (i+1):length(small_blocks)
            % Find the comparison in the results
            for k = 1:size(c, 1)
                group1_name = gnames{c(k, 1)};
                group2_name = gnames{c(k, 2)};
                
                if (strcmp(group1_name, small_blocks{i}) && strcmp(group2_name, small_blocks{j})) || ...
                   (strcmp(group1_name, small_blocks{j}) && strcmp(group2_name, small_blocks{i}))
                    diff = c(k, 4);
                    lower_ci = c(k, 3);
                    upper_ci = c(k, 5);
                    p_val = c(k, 6);
                    
                    fprintf('Small Block %d vs Block %d: Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.3e', ...
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
    
else
    disp('Interaction not significant. Post-hoc tests for interaction not needed.');
    
    % If main effects are significant, analyze those
    pValueChoice = ranovaResults{3, 5}; % p-value for Choice main effect
    pValueBlock = ranovaResults{5, 5};  % p-value for TrialBlock main effect
    
    % Check main effect of Choice
    if pValueChoice < 0.05
        fprintf('\n--- Main effect of Choice is significant (p = %.4f) ---\n', pValueChoice);
        % Compare overall means (Large vs Small)
        largeMean = mean(mean(data(:,1:3)));
        smallMean = mean(mean(data(:,4:6)));
        fprintf('Large Choice mean: %.3f\n', largeMean);
        fprintf('Small Choice mean: %.3f\n', smallMean);
        fprintf('Difference: %.3f\n', largeMean - smallMean);
    end
    
    % Check main effect of Trial Block
    if pValueBlock < 0.05
        fprintf('\n--- Main effect of Trial Block is significant (p = %.4f) ---\n', pValueBlock);
        % Run pairwise comparisons for Trial Block
        blockComp = multcompare(rm, 'TrialBlock', 'ComparisonType', 'bonferroni');
        
        % Display results
        fprintf('Post-hoc comparisons for Trial Block (across both choice types):\n');
        for i = 1:size(blockComp, 1)
            fprintf('Block %s vs Block %s: Mean Diff = %.3f, CI = [%.3f, %.3f], p = %.3e\n', ...
                blockComp.TrialBlock_1(i), blockComp.TrialBlock_2(i), ...
                blockComp.Difference(i), blockComp.Lower(i), ...
                blockComp.Upper(i), blockComp.pValue(i));
        end
    end
end

%%
clc
% Combine data for analysis
data = [large_choice, small_choice]; % Combine large and small choices
% data = [large_path_length, small_path_length]; % Combine large and small choices

num_subjects = size(data, 1); % Number of subjects (mice)

% Create table for analysis
varNames = {'Large_Block1', 'Large_Block2', 'Large_Block3', ...
            'Small_Block1', 'Small_Block2', 'Small_Block3'};
tbl = array2table(data, 'VariableNames', varNames);

% Define within-subject factors
Choice = categorical({'Large', 'Large', 'Large', 'Small', 'Small', 'Small'}); % Choice factor
TrialBlock = categorical([1, 2, 3, 1, 2, 3]); % Trial Block factor
WithinDesign = table(Choice', TrialBlock', 'VariableNames', {'Choice', 'TrialBlock'});

% Fit repeated measures model
rm = fitrm(tbl, 'Large_Block1-Small_Block3 ~ 1', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovaResults = ranova(rm, 'WithinModel', 'Choice*TrialBlock');


mauchlyTest = mauchly(rm)
% Display results
disp(ranovaResults);

% Display ANOVA results
fprintf('Reward Size (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{3,2}, ranovaResults{4,2}, ranovaResults{3,4}, ranovaResults{3,5});

% Display ANOVA results
fprintf('Trial Block effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{5,2}, ranovaResults{6,2}, ranovaResults{5,4}, ranovaResults{5,5});
% Display ANOVA results
fprintf('Interaction effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{7,2}, ranovaResults{8,2}, ranovaResults{7,4}, ranovaResults{7,5});

% Check if interaction is significant (using p < 0.05 threshold)
pValueInteraction = ranovaResults{7, 5}; % p-value for Choice*TrialBlock interaction
disp(['Interaction p-value: ', num2str(pValueInteraction)]);

% If interaction is significant, conduct post-hoc tests
if pValueInteraction < 0.05
    disp('Significant interaction detected. Conducting post-hoc tests...');
    fprintf('\n========== POST-HOC TESTS ==========\n');
    
    % APPROACH 1: Simple main effects - analyzing one factor at each level of the other factor
    % Compare Choice types at each level of Trial Block
    fprintf('\n--- Simple main effects: Choice at each Trial Block ---\n');
    for block = 1:3
        blockStr = num2str(block);
        % Extract data for this block
        largeData = data(:, block);
        smallData = data(:, block + 3);
        
        % Perform paired t-test for this block
        [h, p, ci, stats] = ttest(largeData, smallData);
        
        fprintf('Block %s (Large vs Small): t(%d) = %.3f, p = %.4f\n', ...
            blockStr, stats.df, stats.tstat, p);
    end
    
    % APPROACH 2: Simple main effects - analyzing Trial Blocks for each Choice type
    fprintf('\n--- Simple main effects: Trial Blocks for Large Choice ---\n');
    % Run one-way repeated measures ANOVA for Large Choice
    largeModel = fitrm(tbl, 'Large_Block1-Large_Block3 ~ 1');
    largeRanova = ranova(largeModel);
    fprintf('Large Choice - Effect of Trial Block: F(%.1f,%.1f) = %.3f, p = %.4f\n', ...
        largeRanova{1,2}, largeRanova{2,2}, largeRanova{1,4}, largeRanova{1,5});
    
    % If significant, run post-hoc pairwise comparisons with Bonferroni correction
    if largeRanova{1,5} < 0.05
        fprintf('  Post-hoc comparisons for Large Choice across blocks:\n');
        % Pairwise comparisons for Large
        largeComparisons = multcompare(largeModel, 'Time', 'ComparisonType', 'bonferroni');
        % Display results
        for i = 1:size(largeComparisons, 1)
            fprintf('  Large Rew Block %d vs Block %d: Mean Diff = %.3f, CI = [%.3f, %.3f], p = %.4f\n', ...
                largeComparisons.Time_1(i), largeComparisons.Time_2(i), ...
                largeComparisons.Difference(i), largeComparisons.Lower(i), ...
                largeComparisons.Upper(i), largeComparisons.pValue(i));
        end
    end
    
    fprintf('\n--- Simple main effects: Trial Blocks for Small Choice ---\n');
    % Run one-way repeated measures ANOVA for Small Choice
    smallModel = fitrm(tbl, 'Small_Block1-Small_Block3 ~ 1');
    smallRanova = ranova(smallModel);
    fprintf('Small Choice - Effect of Trial Block: F(%.1f,%.1f) = %.3f, p = %.4f\n', ...
        smallRanova{1,2}, smallRanova{2,2}, smallRanova{1,4}, smallRanova{1,5});
    
    % If significant, run post-hoc pairwise comparisons with Bonferroni correction
    if smallRanova{1,5} < 0.05
        fprintf('  Post-hoc comparisons for Small Choice across blocks:\n');
        % Pairwise comparisons for Small
        smallComparisons = multcompare(smallModel, 'Time', 'ComparisonType', 'bonferroni');
        % Display results
        for i = 1:size(smallComparisons, 1)
            fprintf('  Small Rew Block %d vs Block %d: Mean Diff = %.3f, CI = [%.3f, %.3f], p = %.4f\n', ...
                smallComparisons.Time_1(i), smallComparisons.Time_2(i), ...
                smallComparisons.Difference(i), smallComparisons.Lower(i), ...
                smallComparisons.Upper(i), smallComparisons.pValue(i));
        end
    end
    
    % % Visualize the interaction
    % fprintf('\n--- Visualizing the Interaction ---\n');
    % 
    % % Calculate means and standard errors
    % largeMeans = mean(data(:,1:3));
    % smallMeans = mean(data(:,4:6));
    % largeSE = std(data(:,1:3)) / sqrt(num_subjects);
    % smallSE = std(data(:,4:6)) / sqrt(num_subjects);
    % 
    % figure;
    % errorbar([1, 2, 3], largeMeans, largeSE, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'auto');
    % hold on;
    % errorbar([1, 2, 3], smallMeans, smallSE, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'auto');
    % 
    % xlabel('Trial Block');
    % ylabel('Response');
    % legend('Large Choice', 'Small Choice', 'Location', 'best');
    % title('Interaction between Choice Type and Trial Block');
    % xticks([1, 2, 3]);
    % xticklabels({'Block 1', 'Block 2', 'Block 3'});
    % grid on;
    % box on;
    
else
    disp('Interaction not significant. Post-hoc tests for interaction not needed.');
    
    % If main effects are significant, analyze those
    pValueChoice = ranovaResults{1, 5}; % p-value for Choice main effect
    pValueBlock = ranovaResults{2, 5};  % p-value for TrialBlock main effect
    
    % Check main effect of Choice
    if pValueChoice < 0.05
        fprintf('\n--- Main effect of Choice is significant (p = %.4f) ---\n', pValueChoice);
        % Compare overall means (Large vs Small)
        largeMean = mean(mean(data(:,1:3)));
        smallMean = mean(mean(data(:,4:6)));
        fprintf('Large Choice mean: %.3f\n', largeMean);
        fprintf('Small Choice mean: %.3f\n', smallMean);
        fprintf('Difference: %.3f\n', largeMean - smallMean);
    end
    
    % Check main effect of Trial Block
    if pValueBlock < 0.05
        fprintf('\n--- Main effect of Trial Block is significant (p = %.4f) ---\n', pValueBlock);
        % Run pairwise comparisons for Trial Block
        blockComp = multcompare(rm, 'TrialBlock', 'ComparisonType', 'bonferroni');
        
        % Display results
        fprintf('Post-hoc comparisons for Trial Block (across both choice types):\n');
        for i = 1:size(blockComp, 1)
            fprintf('Block %s vs Block %s: Mean Diff = %.3f, CI = [%.3f, %.3f], p = %.4f\n', ...
                blockComp.TrialBlock_1{i}, blockComp.TrialBlock_2{i}, ...
                blockComp.Difference(i), blockComp.Lower(i), ...
                blockComp.Upper(i), blockComp.pValue(i));
        end
    end
end
%%
% Combine data for analysis
data = [large_choice]; % Assuming large_choice contains data for the three blocks
num_subjects = size(data, 1); % Number of subjects (mice)

% Create table for analysis
varNames = {'Large_Block1', 'Large_Block2', 'Large_Block3'};
tbl = array2table(data, 'VariableNames', varNames);

% Define within-subject factor: TrialBlock
TrialBlock = categorical([1, 2, 3]); % Trial Block factor
WithinDesign = table(TrialBlock', 'VariableNames', {'TrialBlock'});

% Fit repeated measures model
rm = fitrm(tbl, 'Large_Block1-Large_Block3 ~ 1', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovaResults = ranova(rm, 'WithinModel', 'TrialBlock');

% Display results
disp(ranovaResults);

% Perform post hoc pairwise comparisons
disp('Post-hoc pairwise comparisons (Tukey-corrected):');
mc = multcompare(rm, 'TrialBlock', 'ComparisonType', 'tukey');
disp(mc);
%% between subjects for 2 treatment groups TUKEYS
clc
% Create wide-format data for analysis
num_mice_mCherry = size(large_choice_mCherry, 1);
num_mice_hM4Di = size(large_choice_hM4Di, 1);
% Combine data from both groups into a wide-format table
data = [large_choice_mCherry; large_choice_hM4Di];
Treatment = [repmat({'mCherry'}, num_mice_mCherry, 1); repmat({'hM4Di'}, num_mice_hM4Di, 1)];
Subject = (1:size(data, 1))';
% Convert to table
tbl = array2table(data, 'VariableNames', {'Block1', 'Block2', 'Block3'});
tbl.Treatment = categorical(Treatment);
tbl.Subject = Subject;
% Define within-subject factors
TrialBlock = categorical([1, 2, 3]');
WithinDesign = table(TrialBlock, 'VariableNames', {'TrialBlock'});
% Fit repeated measures model with Treatment as a between-subject factor
rm = fitrm(tbl, 'Block1-Block3 ~ Treatment', 'WithinDesign', WithinDesign);
% Run repeated measures ANOVA
ranovaResults = ranova(rm);
% Display results
disp(ranovaResults);

% Display ANOVA results
fprintf('Trial Block effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{1,2}, ranovaResults{3,2}, ranovaResults{1,4}, ranovaResults{1,5});
% Display ANOVA results
fprintf('Interaction effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{2,2}, ranovaResults{3,2}, ranovaResults{2,4}, ranovaResults{2,5});
% Check if the interaction effect is significant
p_interaction = ranovaResults.pValue(2); % The third row corresponds to the Treatment × TrialBlock interaction
if p_interaction < 0.05
    disp('Significant interaction found! Performing Tukey''s multiple comparisons...');
    
    % Reshape data for post-hoc analysis
    % Create long-format data for Tukey's test
    all_values = [];
    all_groups = {};
    
    for block = 1:3
        % mCherry group for this block
        mCherry_block_data = data(1:num_mice_mCherry, block);
        all_values = [all_values; mCherry_block_data];
        all_groups = [all_groups; repmat({sprintf('mCherry_Block%d', block)}, length(mCherry_block_data), 1)];
        
        % hM4Di group for this block
        hM4Di_block_data = data(num_mice_mCherry+1:end, block);
        all_values = [all_values; hM4Di_block_data];
        all_groups = [all_groups; repmat({sprintf('hM4Di_Block%d', block)}, length(hM4Di_block_data), 1)];
    end
    
    % Perform one-way ANOVA to get stats structure for multcompare
    [~, ~, stats_posthoc] = anova1(all_values, all_groups, 'off');
    
    % Perform Tukey's multiple comparisons
    [c, m, h, gnames] = multcompare(stats_posthoc, 'CType', 'tukey-kramer', 'Display', 'off');
    
    % Display Tukey's results in a more readable format
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
        
        fprintf('%-20s %-20s %10.3f %10.3f %10.3f %10.6f', ...
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
    
    % Highlight comparisons of interest (same blocks across treatments)
    fprintf('\nTreatment comparisons within each block:\n');
    block_comparisons = {};
    for block = 1:3
        mCherry_name = sprintf('mCherry_Block%d', block);
        hM4Di_name = sprintf('hM4Di_Block%d', block);
        
        % Find the comparison in the results
        for i = 1:size(c, 1)
            group1_name = gnames{c(i, 1)};
            group2_name = gnames{c(i, 2)};
            
            if (strcmp(group1_name, mCherry_name) && strcmp(group2_name, hM4Di_name)) || ...
               (strcmp(group1_name, hM4Di_name) && strcmp(group2_name, mCherry_name))
                diff = c(i, 4);
                lower_ci = c(i, 3);
                upper_ci = c(i, 5);
                p_val = c(i, 6);
                
                fprintf('Block %d (mCherry vs hM4Di): Diff = %.3f, 95%% CI [%.3f, %.3f], p = %.6f', ...
                    block, abs(diff), lower_ci, upper_ci, p_val);
                
                if p_val < 0.05
                    fprintf(' *');
                end
                fprintf('\n');
                break;
            end
        end
    end
    
else
    disp('No significant interaction effect found.');
end

% Calculate mean response for each subject across all blocks
mCherry_means = mean(data(1:num_mice_mCherry,:), 2); % Average across blocks
hM4Di_means = mean(data(num_mice_mCherry+1:end,:), 2); % Average across blocks

% Create variables for ANOVA
group_means = [mCherry_means; hM4Di_means];
group_labels = [repmat({'mCherry'}, length(mCherry_means), 1); repmat({'hM4Di'}, length(hM4Di_means), 1)];

% Perform one-way ANOVA
[p, tbl_anova, stats] = anova1(group_means, group_labels, 'off');

% Display ANOVA results
fprintf('Treatment effect (one-way ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    tbl_anova{2,3}, tbl_anova{3,3}, tbl_anova{2,5}, p);

if p < 0.05
    fprintf('Treatment effect is significant (p = %e)\n', p);
    fprintf('mCherry group mean: %.3f\n', mean(mCherry_means));
    fprintf('hM4Di group mean: %.3f\n', mean(hM4Di_means));
    fprintf('Difference: %.3f\n', mean(mCherry_means) - mean(hM4Di_means));
    
    % If you want to add multiple comparison test
    [c,m,h,gnames] = multcompare(stats, 'Display', 'off');
    fprintf('Multiple comparison p-value: %e\n', c(1,6));
else
    fprintf('Treatment effect is not significant (p = %e)\n', p);
end

%% between subjects for 2 treatment groups t-tests
clc
% Create wide-format data for analysis
num_mice_mCherry = size(large_choice_mCherry, 1);
num_mice_hM4Di = size(large_choice_hM4Di, 1);
% Combine data from both groups into a wide-format table
data = [large_choice_mCherry; large_choice_hM4Di];
Treatment = [repmat({'mCherry'}, num_mice_mCherry, 1); repmat({'hM4Di'}, num_mice_hM4Di, 1)];
Subject = (1:size(data, 1))';
% Convert to table
tbl = array2table(data, 'VariableNames', {'Block1', 'Block2', 'Block3'});
tbl.Treatment = categorical(Treatment);
tbl.Subject = Subject;
% Define within-subject factors
TrialBlock = categorical([1, 2, 3]');
WithinDesign = table(TrialBlock, 'VariableNames', {'TrialBlock'});
% Fit repeated measures model with Treatment as a between-subject factor
rm = fitrm(tbl, 'Block1-Block3 ~ Treatment', 'WithinDesign', WithinDesign);
% Run repeated measures ANOVA
ranovaResults = ranova(rm);
% Display results
disp(ranovaResults);

% Display ANOVA results
fprintf('Trial Block effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{1,2}, ranovaResults{3,2}, ranovaResults{1,4}, ranovaResults{1,5});
% Display ANOVA results
fprintf('Interaction effect (RM ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    ranovaResults{2,2}, ranovaResults{3,2}, ranovaResults{2,4}, ranovaResults{2,5});
% Check if the interaction effect is significant
p_interaction = ranovaResults.pValue(2); % The third row corresponds to the Treatment × TrialBlock interaction
if p_interaction < 0.05
    disp('Significant interaction found! Performing post-hoc comparisons at each TrialBlock level...');
    % Extract data for each block
    block1_data = data(:,1);
    block2_data = data(:,2);
    block3_data = data(:,3);
    % Perform independent t-tests at each block
    [~, p_block1, ~, stats_block1] = ttest2(block1_data(1:num_mice_mCherry), block1_data(num_mice_mCherry+1:end));
    [~, p_block2, ~, stats_block2] = ttest2(block2_data(1:num_mice_mCherry), block2_data(num_mice_mCherry+1:end));
    [~, p_block3, ~, stats_block3] = ttest2(block3_data(1:num_mice_mCherry), block3_data(num_mice_mCherry+1:end));
    % Apply Bonferroni correction for multiple comparisons
    p_adjusted = min([p_block1, p_block2, p_block3] * 3, 1); % Ensures values do not exceed 1
    % Create table with t-values, df, p-values, and adjusted p-values
    t_values = [stats_block1.tstat; stats_block2.tstat; stats_block3.tstat];
    df_values = [stats_block1.df; stats_block2.df; stats_block3.df];
    p_values = [p_block1; p_block2; p_block3];
    
    posthoc_results = table(t_values, df_values, p_values, p_adjusted', ...
        'VariableNames', {'t_value', 'df', 'Raw_pValue', 'Bonferroni_pValue'}, ...
        'RowNames', {'Block1', 'Block2', 'Block3'});
    disp('Pairwise comparisons (Treatment effect at each TrialBlock):');
    disp(posthoc_results);
    
    % Also print formatted results
    for i = 1:3
        block_name = sprintf('Block%d', i);
        fprintf('%s: t(%.0f) = %.3f, p = %.4f (adjusted p = %e\n', ...
            block_name, posthoc_results.df(i), posthoc_results.t_value(i), ...
            posthoc_results.Raw_pValue(i), posthoc_results.Bonferroni_pValue(i));
    end
else
    disp('No significant interaction effect found.');
end

% Calculate mean response for each subject across all blocks
mCherry_means = mean(data(1:num_mice_mCherry,:), 2); % Average across blocks
hM4Di_means = mean(data(num_mice_mCherry+1:end,:), 2); % Average across blocks

% Create variables for ANOVA
group_means = [mCherry_means; hM4Di_means];
group_labels = [repmat({'mCherry'}, length(mCherry_means), 1); repmat({'hM4Di'}, length(hM4Di_means), 1)];

% Perform one-way ANOVA
[p, tbl_anova, stats] = anova1(group_means, group_labels, 'off');

% Display ANOVA results
fprintf('Treatment effect (one-way ANOVA): F(%d,%d) = %.3f, p = %e\n', ...
    tbl_anova{2,3}, tbl_anova{3,3}, tbl_anova{2,5}, p);

if p < 0.05
    fprintf('Treatment effect is significant (p = %e)\n', p);
    fprintf('mCherry group mean: %.3f\n', mean(mCherry_means));
    fprintf('hM4Di group mean: %.3f\n', mean(hM4Di_means));
    fprintf('Difference: %.3f\n', mean(mCherry_means) - mean(hM4Di_means));
    
    % If you want to add multiple comparison test
    [c,m,h,gnames] = multcompare(stats, 'Display', 'off');
    fprintf('Multiple comparison p-value: %e\n', c(1,6));
else
    fprintf('Treatment effect is not significant (p = %e)\n', p);
end


%% between subjects for 2 treatment groups DUNNETT'S TEST ONLY
clc
% Create wide-format data for analysis
num_mice_mCherry = size(large_choice_mCherry, 1);
num_mice_hM4Di = size(large_choice_hM4Di, 1);

% Combine data from both groups into a wide-format table
data = [large_choice_mCherry; large_choice_hM4Di];
Treatment = [repmat({'mCherry'}, num_mice_mCherry, 1); 
             repmat({'hM4Di'}, num_mice_hM4Di, 1)];
Subject = (1:size(data, 1))';

% Convert to table
tbl = array2table(data, 'VariableNames', {'Block1', 'Block2', 'Block3'});
tbl.Treatment = categorical(Treatment);
tbl.Subject = Subject;

% Define within-subject factors
TrialBlock = categorical([1, 2, 3]');
WithinDesign = table(TrialBlock, 'VariableNames', {'TrialBlock'});

% Fit repeated measures model with Treatment as a between-subject factor
rm = fitrm(tbl, 'Block1-Block3 ~ Treatment', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovaResults = ranova(rm);

% Get between-subjects effects table using the correct method
betweenSubjectsTable = anova(rm);

% Display results
disp('Mixed ANOVA Results:');
disp(ranovaResults);
disp('Between-Subjects Effects:');
disp(betweenSubjectsTable);

% Correctly access the main effect of Treatment from between-subjects table
p_treatment = betweenSubjectsTable.pValue(1); % First row of between-subjects table is Treatment effect
fprintf('Main effect of Treatment: p = %.4f\n', p_treatment);

% Conduct Dunnett's test for Treatment if main effect is significant
if p_treatment < 0.05
    disp('Significant main effect of Treatment found! Performing Dunnett''s test with mCherry as control...');
    
    % Calculate mean response for each subject across all blocks, grouped by treatment
    mCherry_means = mean(data(1:num_mice_mCherry,:), 2);
    hM4Di_means = mean(data(num_mice_mCherry+1:end,:), 2);
    
    % Prepare data for multcompare with Dunnett's test
    all_means = [mCherry_means; hM4Di_means];
    group_labels = [ones(size(mCherry_means)); 2*ones(size(hM4Di_means))];
    
    % One-way ANOVA for the means
    [p_anova, ~, stats_anova] = anova1(all_means, group_labels, 'off');
    
    % Dunnett's test with mCherry as control (control group index = 1)
    [c, m, h, gnames] = multcompare(stats_anova, 'CType', 'dunnett', 'Control', 1, 'Display', 'off');
    
    % Display results of Dunnett's test
    fprintf('Dunnett''s test results (Control group: mCherry):\n');
    dunnett_results = array2table(c, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});
    
    % Create a mapping from numeric indices to treatment names
    treatment_names = {'mCherry', 'hM4Di'};
    
    % Replace numeric indices with actual treatment names
    dunnett_results.Group1Str = treatment_names(dunnett_results.Group1)';
    dunnett_results.Group2Str = treatment_names(dunnett_results.Group2)';
    
    % Display only the relevant columns with treatment names
    disp(dunnett_results(:, [7 8 4 6])); % Display Group1Str, Group2Str, Difference, and pValue
    
    % Show means for each group
    fprintf('Group means:\n');
    fprintf('mCherry mean: %.3f\n', mean(mCherry_means));
    fprintf('hM4Di mean: %.3f\n', mean(hM4Di_means));
else
    disp('No significant main effect of Treatment found.');
end

% Check if the interaction effect is significant
p_interaction = ranovaResults.pValue(2); % Treatment × TrialBlock interaction
if p_interaction < 0.05
    disp('Significant interaction found! Performing Dunnett''s test at each TrialBlock level...');
    
    % For each block, compare the control to the other treatment group using Dunnett's test
    for block = 1:3
        block_data = data(:,block);
        
        % Extract data for each treatment group in this block
        mCherry_block = block_data(1:num_mice_mCherry);
        hM4Di_block = block_data(num_mice_mCherry+1:end);
        
        % Prepare data for one-way ANOVA and subsequent Dunnett's test
        all_block_data = block_data;
        group_labels = [ones(num_mice_mCherry, 1); 2*ones(num_mice_hM4Di, 1)];
        
        % Run one-way ANOVA for this block
        [p_block_anova, ~, stats_block] = anova1(all_block_data, group_labels, 'off');
        
        fprintf('\nBlock %d: One-way ANOVA p = %.4f\n', block, p_block_anova);
        
        % If ANOVA is significant, run Dunnett's test
        if p_block_anova < 0.05
            % Dunnett's test with mCherry as control (control group index = 1)
            [c_block, m_block, h_block, gnames_block] = multcompare(stats_block, 'CType', 'dunnett', 'Control', 1, 'Display', 'off');
            
            % Display results of Dunnett's test for this block
            fprintf('Dunnett''s test results for Block %d (Control group: mCherry):\n', block);
            block_dunnett_results = array2table(c_block, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});
            
            % Create a mapping from numeric indices to treatment names
            treatment_names = {'mCherry', 'hM4Di'};
            
            % Replace numeric indices with actual treatment names
            block_dunnett_results.Group1Str = treatment_names(block_dunnett_results.Group1)';
            block_dunnett_results.Group2Str = treatment_names(block_dunnett_results.Group2)';
            
            % Display only the relevant columns with treatment names
            disp(block_dunnett_results(:, [7 8 4 6])); % Display Group1Str, Group2Str, Difference, and pValue
            
            % Show means for each group in this block
            fprintf('Block %d means: mCherry = %.3f, hM4Di = %.3f\n', ...
                block, mean(mCherry_block), mean(hM4Di_block));
        else
            fprintf('No significant differences between treatments in Block %d\n', block);
        end
    end
else
    disp('No significant interaction effect found.');
end

% Additional analysis: Test for main effect of Block using repeated measures
p_block = ranovaResults.pValue(1); % Main effect of TrialBlock
if p_block < 0.05
    disp('Significant main effect of TrialBlock found!');
    
    % Perform Mauchly's test for sphericity
    mauchlyTest = rm.mauchly;
    fprintf('Mauchly''s Test of Sphericity: p = %.4f\n', mauchlyTest.pValue);
    
    % Display block means only (no pairwise comparisons needed)
    block_means = [mean(data(:,1)), mean(data(:,2)), mean(data(:,3))];
    fprintf('Block means: Block1 = %.3f, Block2 = %.3f, Block3 = %.3f\n', block_means(1), block_means(2), block_means(3));
else
    disp('No significant main effect of TrialBlock found.');
end

%% between subjects for 3 treatment groups
clc
% Create wide-format data for analysis
num_mice_mCherry = size(large_choice_mCherry, 1);
num_mice_PdCO = size(large_choice_PdCO, 1);
num_mice_ChrimsonR = size(large_choice_ChrimsonR, 1);

% Combine data from all three groups into a wide-format table
data = [large_choice_mCherry; large_choice_PdCO; large_choice_ChrimsonR];
Treatment = [repmat({'mCherry'}, num_mice_mCherry, 1); 
             repmat({'PdCO'}, num_mice_PdCO, 1); 
             repmat({'ChrimsonR'}, num_mice_ChrimsonR, 1)];
Subject = (1:size(data, 1))';

% Convert to table
tbl = array2table(data, 'VariableNames', {'Block1', 'Block2', 'Block3'});
tbl.Treatment = categorical(Treatment);
tbl.Subject = Subject;

% Define within-subject factors
TrialBlock = categorical([1, 2, 3]');
WithinDesign = table(TrialBlock, 'VariableNames', {'TrialBlock'});

% Fit repeated measures model with Treatment as a between-subject factor
rm = fitrm(tbl, 'Block1-Block3 ~ Treatment', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovaResults = ranova(rm);

% Get between-subjects effects table using the correct method
betweenSubjectsTable = anova(rm);

% Display results
disp('Mixed ANOVA Results:');
disp(ranovaResults);
disp('Between-Subjects Effects:');
disp(betweenSubjectsTable);

% Correctly access the main effect of Treatment from between-subjects table
p_treatment = betweenSubjectsTable.pValue(1); % First row of between-subjects table is Treatment effect
fprintf('Main effect of Treatment: p = %.4f\n', p_treatment);

% Conduct post-hoc tests for Treatment if main effect is significant
if p_treatment < 0.05
    disp('Significant main effect of Treatment found! Performing post-hoc comparisons...');
    
    % Calculate mean response for each subject across all blocks, grouped by treatment
    mCherry_means = mean(data(1:num_mice_mCherry,:), 2);
    PdCO_means = mean(data(num_mice_mCherry+1:num_mice_mCherry+num_mice_PdCO,:), 2);
    ChrimsonR_means = mean(data(num_mice_mCherry+num_mice_PdCO+1:end,:), 2);
    
    % Perform pairwise comparisons with Bonferroni correction
    [~, p_mCherry_vs_PdCO, ~, stats_mCherry_vs_PdCO] = ttest2(mCherry_means, PdCO_means);
    [~, p_mCherry_vs_ChrimsonR, ~, stats_mCherry_vs_ChrimsonR] = ttest2(mCherry_means, ChrimsonR_means);
    [~, p_PdCO_vs_ChrimsonR, ~, stats_PdCO_vs_ChrimsonR] = ttest2(PdCO_means, ChrimsonR_means);
    
    % Apply Bonferroni correction (multiply by number of comparisons)
    p_adjusted_mCherry_vs_PdCO = min(p_mCherry_vs_PdCO * 3, 1);
    p_adjusted_mCherry_vs_ChrimsonR = min(p_mCherry_vs_ChrimsonR * 3, 1);
    p_adjusted_PdCO_vs_ChrimsonR = min(p_PdCO_vs_ChrimsonR * 3, 1);
    
    % Display pairwise comparison results
    fprintf('mCherry vs PdCO: p = %.4f (adjusted p = %.4f)\n', p_mCherry_vs_PdCO, p_adjusted_mCherry_vs_PdCO);
    fprintf('mCherry vs ChrimsonR: p = %.4f (adjusted p = %.4f)\n', p_mCherry_vs_ChrimsonR, p_adjusted_mCherry_vs_ChrimsonR);
    fprintf('PdCO vs ChrimsonR: p = %.4f (adjusted p = %.4f)\n', p_PdCO_vs_ChrimsonR, p_adjusted_PdCO_vs_ChrimsonR);
    
    % Show means for each group
    fprintf('Group means:\n');
    fprintf('mCherry mean: %.3f\n', mean(mCherry_means));
    fprintf('PdCO mean: %.3f\n', mean(PdCO_means));
    fprintf('ChrimsonR mean: %.3f\n', mean(ChrimsonR_means));
else
    disp('No significant main effect of Treatment found.');
end

% Check if the interaction effect is significant
p_interaction = ranovaResults.pValue(2); % Treatment × TrialBlock interaction
if p_interaction < 0.05
    disp('Significant interaction found! Performing post-hoc comparisons at each TrialBlock level...');
    
    % For each block, compare the three treatment groups using one-way ANOVA
    for block = 1:3
        block_data = data(:,block);
        
        % Create a table for ANOVA
        block_tbl = table(block_data, categorical(Treatment), 'VariableNames', {'Response', 'Treatment'});
        
        % Run one-way ANOVA for this block
        [p, tbl_anova, stats] = anova1(block_data, Treatment, 'off');
        
        fprintf('\nBlock %d: One-way ANOVA p = %.4f\n', block, p);
        
        % If ANOVA is significant, run multiple comparisons
        if p < 0.05
            % Extract data for each treatment group in this block
            mCherry_block = block_data(1:num_mice_mCherry);
            PdCO_block = block_data(num_mice_mCherry+1:num_mice_mCherry+num_mice_PdCO);
            ChrimsonR_block = block_data(num_mice_mCherry+num_mice_PdCO+1:end);
            
            % Perform pairwise t-tests with Bonferroni correction
            [~, p_mCherry_PdCO] = ttest2(mCherry_block, PdCO_block);
            [~, p_mCherry_ChrimsonR] = ttest2(mCherry_block, ChrimsonR_block);
            [~, p_PdCO_ChrimsonR] = ttest2(PdCO_block, ChrimsonR_block);
            
            % Apply Bonferroni correction
            p_adjusted = [p_mCherry_PdCO, p_mCherry_ChrimsonR, p_PdCO_ChrimsonR] * 3;
            p_adjusted = min(p_adjusted, 1); % Ensure values do not exceed 1
            
            % Display results
            pairwise_results = table([p_mCherry_PdCO; p_mCherry_ChrimsonR; p_PdCO_ChrimsonR], ...
                                     [p_adjusted(1); p_adjusted(2); p_adjusted(3)], ...
                                     'VariableNames', {'Raw_pValue', 'Bonferroni_pValue'}, ...
                                     'RowNames', {'mCherry vs PdCO', 'mCherry vs ChrimsonR', 'PdCO vs ChrimsonR'});
            fprintf('Pairwise comparisons for Block %d:\n', block);
            disp(pairwise_results);
            
            % Show means for each group in this block
            fprintf('Block %d means: mCherry = %.3f, PdCO = %.3f, ChrimsonR = %.3f\n', ...
                block, mean(mCherry_block), mean(PdCO_block), mean(ChrimsonR_block));
        else
            fprintf('No significant differences between treatments in Block %d\n', block);
        end
    end
else
    disp('No significant interaction effect found.');
end

% Additional analysis: Test for main effect of Block using repeated measures
p_block = ranovaResults.pValue(1); % Main effect of TrialBlock
if p_block < 0.05
    disp('Significant main effect of TrialBlock found!');
    
    % Perform Mauchly's test for sphericity
    mauchlyTest = rm.mauchly;
    fprintf('Mauchly''s Test of Sphericity: p = %.4f\n', mauchlyTest.pValue);
    
    % Perform post-hoc pairwise comparisons between blocks
    blockComparisons = multcompare(rm, 'TrialBlock');
    disp('Pairwise comparisons between Trial Blocks:');
    disp(blockComparisons);
    
    % Calculate means for each block
    block_means = [mean(data(:,1)), mean(data(:,2)), mean(data(:,3))];
    fprintf('Block means: Block1 = %.3f, Block2 = %.3f, Block3 = %.3f\n', block_means(1), block_means(2), block_means(3));
else
    disp('No significant main effect of TrialBlock found.');
end

%% between subjects for 3 treatment groups DUNNETT'S TEST ONLY!


clc
% Create wide-format data for analysis
num_mice_mCherry = size(large_choice_mCherry, 1);
num_mice_PdCO = size(large_choice_PdCO, 1);
num_mice_ChrimsonR = size(large_choice_ChrimsonR, 1);

% Combine data from all three groups into a wide-format table
data = [large_choice_mCherry; large_choice_PdCO; large_choice_ChrimsonR];
Treatment = [repmat({'mCherry'}, num_mice_mCherry, 1); 
             repmat({'PdCO'}, num_mice_PdCO, 1); 
             repmat({'ChrimsonR'}, num_mice_ChrimsonR, 1)];
Subject = (1:size(data, 1))';

% Convert to table
tbl = array2table(data, 'VariableNames', {'Block1', 'Block2', 'Block3'});
tbl.Treatment = categorical(Treatment);
tbl.Subject = Subject;

% Define within-subject factors
TrialBlock = categorical([1, 2, 3]');
WithinDesign = table(TrialBlock, 'VariableNames', {'TrialBlock'});

% Fit repeated measures model with Treatment as a between-subject factor
rm = fitrm(tbl, 'Block1-Block3 ~ Treatment', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovaResults = ranova(rm);

% Get between-subjects effects table using the correct method
betweenSubjectsTable = anova(rm);

% Display results
disp('Mixed ANOVA Results:');
disp(ranovaResults);
disp('Between-Subjects Effects:');
disp(betweenSubjectsTable);

% Correctly access the main effect of Treatment from between-subjects table
p_treatment = betweenSubjectsTable.pValue(1); % First row of between-subjects table is Treatment effect
fprintf('Main effect of Treatment: p = %.4f\n', p_treatment);

% Conduct Dunnett's test for Treatment if main effect is significant
if p_treatment < 0.05
    disp('Significant main effect of Treatment found! Performing Dunnett''s test with mCherry as control...');
    
    % Calculate mean response for each subject across all blocks, grouped by treatment
    mCherry_means = mean(data(1:num_mice_mCherry,:), 2);
    PdCO_means = mean(data(num_mice_mCherry+1:num_mice_mCherry+num_mice_PdCO,:), 2);
    ChrimsonR_means = mean(data(num_mice_mCherry+num_mice_PdCO+1:end,:), 2);
    
    % Prepare data for multcompare with Dunnett's test
    all_means = [mCherry_means; PdCO_means; ChrimsonR_means];
    group_labels = [ones(size(mCherry_means)); 2*ones(size(PdCO_means)); 3*ones(size(ChrimsonR_means))];
    
    % One-way ANOVA for the means
    [p_anova, ~, stats_anova] = anova1(all_means, group_labels, 'off');
    
    % Dunnett's test with mCherry as control (control group index = 1)
    [c, m, h, gnames] = multcompare(stats_anova, 'CType', 'dunnett', 'Control', 1, 'Display', 'off');
    
    % Display results of Dunnett's test
    fprintf('Dunnett''s test results (Control group: mCherry):\n');
    dunnett_results = array2table(c, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});
    
    % Create a mapping from numeric indices to treatment names
    treatment_names = {'mCherry', 'PdCO', 'ChrimsonR'};
    
    % Replace numeric indices with actual treatment names
    dunnett_results.Group1Str = treatment_names(dunnett_results.Group1)';
    dunnett_results.Group2Str = treatment_names(dunnett_results.Group2)';
    
    % Display only the relevant columns with treatment names
    disp(dunnett_results(:, [7 8 4 6])); % Display Group1Str, Group2Str, Difference, and pValue
    
    % Show means for each group
    fprintf('Group means:\n');
    fprintf('mCherry mean: %.3f\n', mean(mCherry_means));
    fprintf('PdCO mean: %.3f\n', mean(PdCO_means));
    fprintf('ChrimsonR mean: %.3f\n', mean(ChrimsonR_means));
else
    disp('No significant main effect of Treatment found.');
end

% Check if the interaction effect is significant
p_interaction = ranovaResults.pValue(2); % Treatment × TrialBlock interaction
if p_interaction < 0.05
    disp('Significant interaction found! Performing Dunnett''s test at each TrialBlock level...');
    
    % For each block, compare the control to other treatment groups using Dunnett's test
    for block = 1:3
        block_data = data(:,block);
        
        % Extract data for each treatment group in this block
        mCherry_block = block_data(1:num_mice_mCherry);
        PdCO_block = block_data(num_mice_mCherry+1:num_mice_mCherry+num_mice_PdCO);
        ChrimsonR_block = block_data(num_mice_mCherry+num_mice_PdCO+1:end);
        
        % Prepare data for one-way ANOVA and subsequent Dunnett's test
        all_block_data = block_data;
        group_labels = [ones(num_mice_mCherry, 1); 2*ones(num_mice_PdCO, 1); 3*ones(num_mice_ChrimsonR, 1)];
        
        % Run one-way ANOVA for this block
        [p_block_anova, ~, stats_block] = anova1(all_block_data, group_labels, 'off');
        
        fprintf('\nBlock %d: One-way ANOVA p = %.4f\n', block, p_block_anova);
        
        % If ANOVA is significant, run Dunnett's test
        if p_block_anova < 0.05
            % Dunnett's test with mCherry as control (control group index = 1)
            [c_block, m_block, h_block, gnames_block] = multcompare(stats_block, 'CType', 'dunnett', 'Control', 1, 'Display', 'off');
            
            % Display results of Dunnett's test for this block
            fprintf('Dunnett''s test results for Block %d (Control group: mCherry):\n', block);
            block_dunnett_results = array2table(c_block, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});
            
            % Create a mapping from numeric indices to treatment names
            treatment_names = {'mCherry', 'PdCO', 'ChrimsonR'};
            
            % Replace numeric indices with actual treatment names
            block_dunnett_results.Group1Str = treatment_names(block_dunnett_results.Group1)';
            block_dunnett_results.Group2Str = treatment_names(block_dunnett_results.Group2)';
            
            % Display only the relevant columns with treatment names
            disp(block_dunnett_results(:, [7 8 4 6])); % Display Group1Str, Group2Str, Difference, and pValue
            
            % Show means for each group in this block
            fprintf('Block %d means: mCherry = %.3f, PdCO = %.3f, ChrimsonR = %.3f\n', ...
                block, mean(mCherry_block), mean(PdCO_block), mean(ChrimsonR_block));
        else
            fprintf('No significant differences between treatments in Block %d\n', block);
        end
    end
else
    disp('No significant interaction effect found.');
end

% Additional analysis: Test for main effect of Block using repeated measures
p_block = ranovaResults.pValue(1); % Main effect of TrialBlock
if p_block < 0.05
    disp('Significant main effect of TrialBlock found!');
    
    % Perform Mauchly's test for sphericity
    mauchlyTest = rm.mauchly;
    fprintf('Mauchly''s Test of Sphericity: p = %.4f\n', mauchlyTest.pValue);
    
    % Display block means only (no pairwise comparisons needed)
    block_means = [mean(data(:,1)), mean(data(:,2)), mean(data(:,3))];
    fprintf('Block means: Block1 = %.3f, Block2 = %.3f, Block3 = %.3f\n', block_means(1), block_means(2), block_means(3));
else
    disp('No significant main effect of TrialBlock found.');
end



%% between subjects


% Create wide-format data for analysis
num_mice_mCherry = size(large_choice_risky, 1);
num_mice_hM4Di = size(large_choice_not_risky, 1);

% Combine data from both groups into a wide-format table
data = [large_choice_risky; large_choice_not_risky];
Treatment = [repmat({'risky'}, num_mice_mCherry, 1); repmat({'not_risky'}, num_mice_hM4Di, 1)];
Subject = (1:size(data, 1))';

% Convert to table
tbl = array2table(data, 'VariableNames', {'Block1', 'Block2', 'Block3'});
tbl.Treatment = categorical(Treatment);
tbl.Subject = Subject;

% Define within-subject factors
TrialBlock = categorical([1, 2, 3]');
WithinDesign = table(TrialBlock, 'VariableNames', {'TrialBlock'});

% Fit repeated measures model with Treatment as a between-subject factor
rm = fitrm(tbl, 'Block1-Block3 ~ Treatment', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovaResults = ranova(rm);

% Display results
disp(ranovaResults);


% Extract average values across trial blocks for between-subject analysis
mean_data = mean(data, 2); % Average across the 3 blocks for each subject

% Perform one-way ANOVA for Treatment effect
[p, tbl_between, stats] = anova1(mean_data, tbl.Treatment, 'off'); % Suppress boxplot
disp('Between-Subject Effect (Treatment):');
disp(tbl_between);

% Optional: Perform post-hoc test if needed
c = multcompare(stats, 'Display', 'off'); % Uncomment for post-hoc tests