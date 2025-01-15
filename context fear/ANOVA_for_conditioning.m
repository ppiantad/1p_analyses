anova_table = [experimental_mice_no_shk_period; one_context_mice_shk_period; no_shock_mice_no_shk_period]
anova_table.no_shk_period = [experimental_data_no_shk_period; one_context_data_no_shk_period; no_shock_data_no_shk_period]
anova_table.shk_period = [experimental_data_shk_period; one_context_data_shk_period; no_shock_data_shk_period]

% anova_table = [experimental_mice_no_shk_period_females(:, [1,3]); experimental_mice_no_shk_period_males(:, [1,3])]
% anova_table.no_shk_period = [experimental_data_no_shk_period_females; experimental_data_no_shk_period_males]
% anova_table.shk_period = [experimental_data_shk_period_females; experimental_data_shk_period_males]

% Convert 'group' column to categorical if not already
% anova_table.group = categorical(anova_table.sex);
anova_table.group = categorical(anova_table.group);

% Define the within-subjects design (two conditions: 'no_shk_period' and 'shk_period')
within = table({'no_shk_period'; 'shk_period'}, 'VariableNames', {'Condition'});

% Fit the repeated measures model
rm = fitrm(anova_table, 'no_shk_period-shk_period ~ group', 'WithinDesign', within);

% Perform repeated measures ANOVA
ranova_results = ranova(rm, 'WithinModel', 'Condition');

% Display the results
disp(ranova_results);

%%
% Post hoc comparisons for the within-subjects factor
posthoc_within = multcompare(rm, 'Condition');

% Display the results
disp(posthoc_within);

%%
% Post hoc comparisons for the between-subjects factor
posthoc_between = multcompare(rm, 'group');

% Display the results
disp(posthoc_between);

%%
% Interaction post hoc comparisons
posthoc_interaction = multcompare(rm, 'group', 'By', 'Condition');

% Display the results
disp(posthoc_interaction);