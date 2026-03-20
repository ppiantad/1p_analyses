function results = mixed_anova_dunnett(varargin)
% MIXED_ANOVA_DUNNETT Perform mixed ANOVA with Dunnett's post-hoc test
%
% Syntax:
%   results = mixed_anova_dunnett(group1_data, group2_data, group3_data, ...)
%   results = mixed_anova_dunnett(..., 'ControlGroup', control_name)
%   results = mixed_anova_dunnett(..., 'GroupNames', {name1, name2, ...})
%
% Inputs:
%   group_data - Matrices where each row is a subject and columns are blocks
%   'ControlGroup' - Name of control group for Dunnett's test (default: first group)
%   'GroupNames' - Cell array of group names (default: {'Group1', 'Group2', ...})
%
% Output:
%   results - Structure containing all statistical results
%
% Example:
%   results = mixed_anova_dunnett(large_choice_mCherry, large_choice_PdCO, ...
%                                 large_choice_ChrimsonR, ...
%                                 'GroupNames', {'mCherry', 'PdCO', 'ChrimsonR'}, ...
%                                 'ControlGroup', 'mCherry');

    % Parse inputs
    p = inputParser;
    p.KeepUnmatched = false;
    
    % Add parameter name-value pairs
    addParameter(p, 'GroupNames', {}, @iscell);
    addParameter(p, 'ControlGroup', '', @ischar);
    
    % Find where the data ends and parameters begin
    data_end_idx = find(cellfun(@ischar, varargin), 1, 'first') - 1;
    if isempty(data_end_idx)
        data_end_idx = length(varargin);
    end
    
    % Extract data
    group_data = varargin(1:data_end_idx);
    num_groups = length(group_data);
    
    % Parse the remaining parameters
    if data_end_idx < length(varargin)
        parse(p, varargin{data_end_idx+1:end});
    else
        parse(p);
    end
    
    % Set default group names if not provided
    if isempty(p.Results.GroupNames)
        group_names = arrayfun(@(x) sprintf('Group%d', x), 1:num_groups, 'UniformOutput', false);
    else
        group_names = p.Results.GroupNames;
        if length(group_names) ~= num_groups
            error('Number of group names must match number of data groups');
        end
    end
    
    % Set default control group if not provided
    if isempty(p.Results.ControlGroup)
        control_group = group_names{1};
    else
        control_group = p.Results.ControlGroup;
    end
    
    % Find control group index
    control_idx = find(strcmp(group_names, control_group));
    if isempty(control_idx)
        error('Control group "%s" not found in group names', control_group);
    end
    
    % Get sample sizes
    num_mice_per_group = cellfun(@(x) size(x, 1), group_data);
    num_blocks = size(group_data{1}, 2);
    
    % Verify all groups have same number of blocks
    if ~all(cellfun(@(x) size(x, 2), group_data) == num_blocks)
        error('All groups must have the same number of blocks');
    end
    
    % Combine data from all groups
    data = vertcat(group_data{:});
    
    % Create treatment labels
    Treatment = cell(sum(num_mice_per_group), 1);
    start_idx = 1;
    for g = 1:num_groups
        end_idx = start_idx + num_mice_per_group(g) - 1;
        Treatment(start_idx:end_idx) = repmat(group_names(g), num_mice_per_group(g), 1);
        start_idx = end_idx + 1;
    end
    
    Subject = (1:size(data, 1))';
    
    % Convert to table
    block_names = arrayfun(@(x) sprintf('Block%d', x), 1:num_blocks, 'UniformOutput', false);
    tbl = array2table(data, 'VariableNames', block_names);
    tbl.Treatment = categorical(Treatment);
    tbl.Subject = Subject;
    
    % Define within-subject factors
    TrialBlock = categorical(1:num_blocks)';
    WithinDesign = table(TrialBlock, 'VariableNames', {'TrialBlock'});
    
    % Fit repeated measures model
    formula = sprintf('%s-%s ~ Treatment', block_names{1}, block_names{end});
    rm = fitrm(tbl, formula, 'WithinDesign', WithinDesign);
    
    % Run repeated measures ANOVA
    ranovaResults = ranova(rm);
    betweenSubjectsTable = anova(rm);
    
    % Display results
    fprintf('\n=== MIXED ANOVA RESULTS ===\n');
    disp(ranovaResults);
    fprintf('\n=== BETWEEN-SUBJECTS EFFECTS ===\n');
    disp(betweenSubjectsTable);
    
    % Store results
    results.ranova = ranovaResults;
    results.between_subjects = betweenSubjectsTable;
    results.rm_model = rm;
    results.group_names = group_names;
    results.control_group = control_group;
    
    % Main effect of Treatment
    p_treatment = betweenSubjectsTable.pValue(1);
    fprintf('\nMain effect of Treatment: p = %.4f\n', p_treatment);
    results.p_treatment = p_treatment;
    
    % Dunnett's test for Treatment if significant
    if p_treatment < 0.05
        fprintf('\nSignificant main effect of Treatment found!\n');
        fprintf('Performing Dunnett''s test with %s as control...\n', control_group);
        
        % Calculate mean response for each subject across all blocks
        group_means = cell(num_groups, 1);
        start_idx = 1;
        for g = 1:num_groups
            end_idx = start_idx + num_mice_per_group(g) - 1;
            group_means{g} = mean(data(start_idx:end_idx, :), 2);
            start_idx = end_idx + 1;
        end
        
        % Prepare data for Dunnett's test
        all_means = vertcat(group_means{:});
        group_labels = [];
        for g = 1:num_groups
            group_labels = [group_labels; g * ones(length(group_means{g}), 1)];
        end
        
        % One-way ANOVA
        [p_anova, ~, stats_anova] = anova1(all_means, group_labels, 'off');
        
        % Dunnett's test
        [c, m, h, gnames] = multcompare(stats_anova, 'CType', 'dunnett', ...
                                        'Control', control_idx, 'Display', 'off');
        
        % Format results
        fprintf('\n=== DUNNETT''S TEST (Control: %s) ===\n', control_group);
        dunnett_results = array2table(c, 'VariableNames', ...
            {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});
        dunnett_results.Group1Str = group_names(dunnett_results.Group1)';
        dunnett_results.Group2Str = group_names(dunnett_results.Group2)';
        disp(dunnett_results(:, [7 8 4 6]));
        
        fprintf('\nGroup means:\n');
        for g = 1:num_groups
            fprintf('%s mean: %.3f\n', group_names{g}, mean(group_means{g}));
        end
        
        results.dunnett_overall = dunnett_results;
        results.group_means_overall = cellfun(@mean, group_means);
    else
        fprintf('\nNo significant main effect of Treatment found.\n');
        results.dunnett_overall = [];
    end
    
    % Check interaction effect
    p_interaction = ranovaResults.pValue(2);
    fprintf('\nTreatment × TrialBlock interaction: p = %.4f\n', p_interaction);
    results.p_interaction = p_interaction;
    
    if p_interaction < 0.05
        fprintf('\nSignificant interaction found!\n');
        fprintf('Performing Dunnett''s test at each TrialBlock level...\n');
        
        dunnett_by_block = cell(num_blocks, 1);
        
        for block = 1:num_blocks
            block_data = data(:, block);
            
            % Create group labels for this block
            group_labels = [];
            for g = 1:num_groups
                group_labels = [group_labels; g * ones(num_mice_per_group(g), 1)];
            end
            
            % One-way ANOVA for this block
            [p_block_anova, ~, stats_block] = anova1(block_data, group_labels, 'off');
            
            fprintf('\n--- Block %d ---\n', block);
            fprintf('One-way ANOVA p = %.4f\n', p_block_anova);
            
            if p_block_anova < 0.05
                % Dunnett's test
                [c_block, ~, ~, ~] = multcompare(stats_block, 'CType', 'dunnett', ...
                                                  'Control', control_idx, 'Display', 'off');
                
                fprintf('Dunnett''s test (Control: %s):\n', control_group);
                block_dunnett = array2table(c_block, 'VariableNames', ...
                    {'Group1', 'Group2', 'LowerCI', 'Difference', 'UpperCI', 'pValue'});
                block_dunnett.Group1Str = group_names(block_dunnett.Group1)';
                block_dunnett.Group2Str = group_names(block_dunnett.Group2)';
                disp(block_dunnett(:, [7 8 4 6]));
                
                % Block means
                fprintf('Block %d means: ', block);
                start_idx = 1;
                for g = 1:num_groups
                    end_idx = start_idx + num_mice_per_group(g) - 1;
                    fprintf('%s = %.3f', group_names{g}, mean(block_data(start_idx:end_idx)));
                    if g < num_groups
                        fprintf(', ');
                    end
                    start_idx = end_idx + 1;
                end
                fprintf('\n');
                
                dunnett_by_block{block} = block_dunnett;
            else
                fprintf('No significant differences in Block %d\n', block);
                dunnett_by_block{block} = [];
            end
        end
        
        results.dunnett_by_block = dunnett_by_block;
    else
        fprintf('\nNo significant interaction effect found.\n');
        results.dunnett_by_block = [];
    end
    
    % Main effect of Block
    p_block = ranovaResults.pValue(1);
    fprintf('\nMain effect of TrialBlock: p = %.4f\n', p_block);
    results.p_block = p_block;
    
    if p_block < 0.05
        fprintf('Significant main effect of TrialBlock found!\n');
        
        % Mauchly's test
        mauchlyTest = rm.mauchly;
        fprintf('Mauchly''s Test of Sphericity: p = %.4f\n', mauchlyTest.pValue);
        results.mauchly = mauchlyTest;
        
        % Block means
        block_means = mean(data, 1);
        fprintf('Block means: ');
        for b = 1:num_blocks
            fprintf('Block%d = %.3f', b, block_means(b));
            if b < num_blocks
                fprintf(', ');
            end
        end
        fprintf('\n');
        
        results.block_means = block_means;
    else
        fprintf('No significant main effect of TrialBlock found.\n');
        results.block_means = [];
    end
    
    fprintf('\n');
end