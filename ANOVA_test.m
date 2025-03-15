

% Combine data for analysis
data = [large_choice, small_choice]; % Combine large and small choices
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


% Extract average values across trial blocks for between-subject analysis
mean_data = mean(data, 2); % Average across the 3 blocks for each subject

% Perform one-way ANOVA for Treatment effect
[p, tbl_between, stats] = anova1(mean_data, tbl.Treatment, 'off'); % Suppress boxplot
disp('Between-Subject Effect (Treatment):');
disp(tbl_between);

% Optional: Perform post-hoc test if needed
c = multcompare(stats, 'Display', 'on'); % Uncomment for post-hoc tests

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