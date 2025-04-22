

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

% Display results
disp(ranovaResults);

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
        largeRanova{1,2}, largeRanova{1,3}, largeRanova{1,4}, largeRanova{1,5});
    
    % If significant, run post-hoc pairwise comparisons with Bonferroni correction
    if largeRanova{1,5} < 0.05
        fprintf('  Post-hoc comparisons for Large Choice across blocks:\n');
        % Pairwise comparisons for Large
        largeComparisons = multcompare(largeModel, 'Time', 'ComparisonType', 'bonferroni');
        % Display results
        for i = 1:size(largeComparisons, 1)
            fprintf('  Block %d vs Block %d: Mean Diff = %.3f, CI = [%.3f, %.3f], p = %.4f\n', ...
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
        smallRanova{1,2}, smallRanova{1,3}, smallRanova{1,4}, smallRanova{1,5});
    
    % If significant, run post-hoc pairwise comparisons with Bonferroni correction
    if smallRanova{1,5} < 0.05
        fprintf('  Post-hoc comparisons for Small Choice across blocks:\n');
        % Pairwise comparisons for Small
        smallComparisons = multcompare(smallModel, 'Time', 'ComparisonType', 'bonferroni');
        % Display results
        for i = 1:size(smallComparisons, 1)
            fprintf('  Block %d vs Block %d: Mean Diff = %.3f, CI = [%.3f, %.3f], p = %.4f\n', ...
                smallComparisons.Time_1(i), smallComparisons.Time_2(i), ...
                smallComparisons.Difference(i), smallComparisons.Lower(i), ...
                smallComparisons.Upper(i), smallComparisons.pValue(i));
        end
    end
    
    % Visualize the interaction
    fprintf('\n--- Visualizing the Interaction ---\n');
    
    % Calculate means and standard errors
    largeMeans = mean(data(:,1:3));
    smallMeans = mean(data(:,4:6));
    largeSE = std(data(:,1:3)) / sqrt(num_subjects);
    smallSE = std(data(:,4:6)) / sqrt(num_subjects);
    
    figure;
    errorbar([1, 2, 3], largeMeans, largeSE, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'auto');
    hold on;
    errorbar([1, 2, 3], smallMeans, smallSE, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'auto');
    
    xlabel('Trial Block');
    ylabel('Response');
    legend('Large Choice', 'Small Choice', 'Location', 'best');
    title('Interaction between Choice Type and Trial Block');
    xticks([1, 2, 3]);
    xticklabels({'Block 1', 'Block 2', 'Block 3'});
    grid on;
    box on;
    
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
%% between subjects


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


% Check if the interaction effect is significant
p_interaction = ranovaResults.pValue(2); % The third row corresponds to the Treatment × TrialBlock interaction

if p_interaction < 0.05
    disp('Significant interaction found! Performing post-hoc comparisons at each TrialBlock level...');
    
    % Extract data for each block
    block1_data = data(:,1);
    block2_data = data(:,2);
    block3_data = data(:,3);
    
    % Perform independent t-tests at each block
    [~, p_block1] = ttest2(block1_data(1:num_mice_mCherry), block1_data(num_mice_mCherry+1:end));
    [~, p_block2] = ttest2(block2_data(1:num_mice_mCherry), block2_data(num_mice_mCherry+1:end));
    [~, p_block3] = ttest2(block3_data(1:num_mice_mCherry), block3_data(num_mice_mCherry+1:end));
    
    % Apply Bonferroni correction for multiple comparisons
    p_adjusted = min([p_block1, p_block2, p_block3] * 3, 1); % Ensures values do not exceed 1
    
    % Display results
    posthoc_results = table([p_block1; p_block2; p_block3], p_adjusted', ...
        'VariableNames', {'Raw_pValue', 'Bonferroni_pValue'}, ...
        'RowNames', {'Block1', 'Block2', 'Block3'});
    
    disp('Pairwise comparisons (Treatment effect at each TrialBlock):');
    disp(posthoc_results);
else
    disp('No significant interaction effect found.');
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