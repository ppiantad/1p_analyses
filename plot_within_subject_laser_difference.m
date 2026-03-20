function [diff_data, fig] = plot_within_subject_laser_difference(risk_table, measure_name, varargin)
% PLOT_WITHIN_SUBJECT_LASER_DIFFERENCE Calculate and plot within-subject differences
%
% Syntax:
%   [diff_data, fig] = plot_within_subject_laser_difference(risk_table, measure_name)
%   [diff_data, fig] = plot_within_subject_laser_difference(..., 'YLabel', label)
%   [diff_data, fig] = plot_within_subject_laser_difference(..., 'YLim', [min max])
%   [diff_data, fig] = plot_within_subject_laser_difference(..., 'Title', title_text)
%   [diff_data, fig] = plot_within_subject_laser_difference(..., 'GroupOrder', {'Control', 'A2A', 'D1'})
%   [diff_data, fig] = plot_within_subject_laser_difference(..., 'Colors', struct or cell)
%   [diff_data, fig] = plot_within_subject_laser_difference(..., 'Symbols', struct or cell)
%   [diff_data, fig] = plot_within_subject_laser_difference(..., 'ScatterColor', color or 'match')
%   [diff_data, fig] = plot_within_subject_laser_difference(..., 'FigureWidth', width)
%   [diff_data, fig] = plot_within_subject_laser_difference(..., 'FigureHeight', height)
%
% Inputs:
%   risk_table - Table with columns: TreatmentCondition, LaserTreatment, animalID, and measure
%   measure_name - String name of the column to analyze (e.g., 'final_ratio_attained')
%
% Optional Name-Value Pairs:
%   'YLabel' - Y-axis label (default: 'Laser On - Laser Off')
%   'YLim' - Y-axis limits as [min max] (default: auto)
%   'Title' - Plot title (default: none)
%   'GroupOrder' - Cell array specifying order of groups (default: {'Control', 'A2A', 'D1'})
%   'Colors' - Either:
%              1) Struct with fields matching group names (e.g., colors.Control = [0.3 0.5 0.8] or '#1f77b4')
%              2) Cell array of RGB colors or hex strings in same order as GroupOrder
%              3) Nx3 matrix of RGB colors (default: blue/green/red scheme)
%   'Symbols' - Either:
%              1) Struct with fields matching group names (e.g., symbols.Control = 'o')
%              2) Cell array of symbol strings in same order as GroupOrder
%              (default: 'o' for all groups)
%   'ScatterColor' - Color for individual data points:
%                    'match' to match group colors (default)
%                    'k' or [0 0 0] for black
%                    Any RGB color
%   'ScatterAlpha' - Transparency for scatter points (default: 0.6)
%   'FigureWidth' - Width of figure in pixels (default: 400)
%   'FigureHeight' - Height of figure in pixels (default: 400)
%
% Outputs:
%   diff_data - Structure containing difference scores for each group
%   fig - Figure handle
%
% Examples:
%   % Basic usage with custom figure size
%   [diff_data, fig] = plot_within_subject_laser_difference(risk_table, 'final_ratio_attained', ...
%                                                            'FigureWidth', 600, ...
%                                                            'FigureHeight', 500, ...
%                                                            'YLim', [-5 5]);
%
%   % Using your existing color/symbol variables
%   colors_struct.Control = Control_color;
%   colors_struct.A2A = A2A_color;
%   colors_struct.D1 = D1_color;
%   
%   symbols_struct.Control = Control_symbol;
%   symbols_struct.A2A = A2A_symbol;
%   symbols_struct.D1 = D1_symbol;
%   
%   [diff_data, fig] = plot_within_subject_laser_difference(risk_table, 'final_ratio_attained', ...
%                                                            'Colors', colors_struct, ...
%                                                            'Symbols', symbols_struct, ...
%                                                            'YLabel', 'Δ Final PR Ratio', ...
%                                                            'YLim', [-10 10], ...
%                                                            'FigureWidth', 500, ...
%                                                            'FigureHeight', 600);

    % Helper function to convert hex to RGB
    function rgb = hex2rgb_local(hex_str)
        % Remove '#' if present
        if hex_str(1) == '#'
            hex_str = hex_str(2:end);
        end
        % Convert to RGB
        rgb = [hex2dec(hex_str(1:2)), hex2dec(hex_str(3:4)), hex2dec(hex_str(5:6))] / 255;
    end

    % Helper function to process a single color value
    function rgb = process_color(color_val, group_name)
        if ischar(color_val) || isstring(color_val)
            % Hex string
            rgb = hex2rgb_local(char(color_val));
        elseif isnumeric(color_val)
            % Numeric RGB
            if length(color_val) >= 3
                rgb = color_val(1:3);
            else
                error('Color for %s must have at least 3 elements (RGB)', group_name);
            end
        else
            error('Color for %s must be a hex string or RGB array', group_name);
        end
    end
    
    % Helper function to extract just the marker symbol from errorbar format
    function marker = extract_marker(symbol_str)
        % errorbar symbols can be like '-o', 'o-', 'o', etc.
        % Extract just the marker character
        valid_markers = 'o+*.xsd^v><ph';
        marker = 'o'; % default
        for i = 1:length(symbol_str)
            if ismember(symbol_str(i), valid_markers)
                marker = symbol_str(i);
                return;
            end
        end
    end

    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'YLabel', 'Laser On - Laser Off', @ischar);
    addParameter(p, 'YLim', [], @isnumeric);
    addParameter(p, 'Title', '', @ischar);
    addParameter(p, 'GroupOrder', {'Control', 'A2A', 'D1'}, @iscell);
    addParameter(p, 'Colors', []);
    addParameter(p, 'Symbols', []);
    addParameter(p, 'ScatterColor', 'match');
    addParameter(p, 'ScatterAlpha', 0.6, @isnumeric);
    addParameter(p, 'FigureWidth', 400, @isnumeric);
    addParameter(p, 'FigureHeight', 400, @isnumeric);
    parse(p, varargin{:});
    
    group_order = p.Results.GroupOrder;
    num_groups = length(group_order);
    
    % Process colors input
    if isempty(p.Results.Colors)
        % Default colors
        default_colors = [0.3, 0.5, 0.8;    % Blue for Control
                         0.4, 0.8, 0.4;    % Green for A2A
                         0.9, 0.3, 0.3];   % Red for D1
        colors = default_colors(1:num_groups, :);
    elseif isstruct(p.Results.Colors)
        % Struct format - extract colors in group order
        colors = zeros(num_groups, 3);
        for g = 1:num_groups
            if isfield(p.Results.Colors, group_order{g})
                colors(g, :) = process_color(p.Results.Colors.(group_order{g}), group_order{g});
            else
                error('Colors struct missing field for group: %s', group_order{g});
            end
        end
    elseif iscell(p.Results.Colors)
        % Cell array format
        if length(p.Results.Colors) ~= num_groups
            error('Number of colors must match number of groups');
        end
        colors = zeros(num_groups, 3);
        for g = 1:num_groups
            colors(g, :) = process_color(p.Results.Colors{g}, sprintf('group %d', g));
        end
    else
        % Matrix format
        colors = p.Results.Colors;
        if size(colors, 1) ~= num_groups
            error('Number of colors must match number of groups');
        end
        % Take only first 3 columns if more are provided
        if size(colors, 2) >= 3
            colors = colors(:, 1:3);
        else
            error('Colors matrix must have at least 3 columns (RGB)');
        end
    end
    
    % Process symbols input
    if isempty(p.Results.Symbols)
        % Default symbol
        symbols = repmat({'o'}, 1, num_groups);
    elseif isstruct(p.Results.Symbols)
        % Struct format - extract symbols in group order
        symbols = cell(1, num_groups);
        for g = 1:num_groups
            if isfield(p.Results.Symbols, group_order{g})
                symbols{g} = p.Results.Symbols.(group_order{g});
            else
                error('Symbols struct missing field for group: %s', group_order{g});
            end
        end
    elseif iscell(p.Results.Symbols)
        % Cell array format
        if length(p.Results.Symbols) ~= num_groups
            error('Number of symbols must match number of groups');
        end
        symbols = p.Results.Symbols;
    else
        error('Symbols must be a struct or cell array');
    end
    
    % Process scatter color
    use_group_colors = false;
    if ischar(p.Results.ScatterColor) || isstring(p.Results.ScatterColor)
        if strcmpi(p.Results.ScatterColor, 'match')
            use_group_colors = true;
        elseif strcmpi(p.Results.ScatterColor, 'k')
            scatter_color = [0 0 0];
        else
            % Try to process as hex
            scatter_color = process_color(p.Results.ScatterColor, 'ScatterColor');
        end
    else
        scatter_color = p.Results.ScatterColor;
    end
    
    % Filter out rows with NaN or missing values in the measure column
    valid_rows = ~isnan(risk_table.(measure_name)) & ~ismissing(risk_table.(measure_name));
    
    % Also filter out NaN/missing in critical columns
    valid_rows = valid_rows & ~ismissing(risk_table.TreatmentCondition) & ...
                              ~ismissing(risk_table.LaserTreatment) & ...
                              ~ismissing(risk_table.animalID);
    
    % Create filtered table
    risk_table_clean = risk_table(valid_rows, :);
    
    fprintf('Filtered out %d rows with NaN/missing values\n', sum(~valid_rows));
    fprintf('Analyzing %d valid rows\n\n', sum(valid_rows));
    
    % Initialize structure to store differences
    diff_data = struct();
    
    % Get unique animals
    unique_animals = unique(risk_table_clean.animalID);
    
    % For each group, calculate within-subject differences
    for g = 1:num_groups
        group_name = group_order{g};
        group_diffs = [];
        group_animals = {};
        
        for i = 1:length(unique_animals)
            animal = unique_animals{i};
            
            % Get data for this animal in this group
            animal_data = risk_table_clean(strcmp(risk_table_clean.animalID, animal) & ...
                                          strcmp(risk_table_clean.TreatmentCondition, group_name), :);
            
            % Skip if this animal doesn't have data for this group
            if isempty(animal_data)
                continue;
            end
            
            % Get Laser On and Laser Off values
            laser_on_data = animal_data.(measure_name)(strcmp(animal_data.LaserTreatment, 'Laser On'));
            laser_off_data = animal_data.(measure_name)(strcmp(animal_data.LaserTreatment, 'Laser Off'));
            
            % Additional check: remove any NaN values that might have slipped through
            laser_on_data = laser_on_data(~isnan(laser_on_data));
            laser_off_data = laser_off_data(~isnan(laser_off_data));
            
            % Check that we have both conditions
            if ~isempty(laser_on_data) && ~isempty(laser_off_data)
                % If multiple sessions, take the mean
                laser_on_mean = mean(laser_on_data);
                laser_off_mean = mean(laser_off_data);
                
                % Calculate difference (Laser On - Laser Off)
                diff = laser_on_mean - laser_off_mean;
                
                group_diffs = [group_diffs; diff];
                group_animals{end+1} = animal;
            end
        end
        
        % Store in structure
        diff_data.(group_name).differences = group_diffs;
        diff_data.(group_name).animals = group_animals';
        diff_data.(group_name).mean = mean(group_diffs);
        diff_data.(group_name).sem = std(group_diffs) / sqrt(length(group_diffs));
        diff_data.(group_name).n = length(group_diffs);
    end
    
    % Create figure with specified size
    fig = figure;
    set(fig, 'Position', [100, 100, p.Results.FigureWidth, p.Results.FigureHeight]);
    
    % Prepare data for plotting
    means = zeros(1, num_groups);
    sems = zeros(1, num_groups);
    all_data = cell(1, num_groups);
    
    for g = 1:num_groups
        group_name = group_order{g};
        means(g) = diff_data.(group_name).mean;
        sems(g) = diff_data.(group_name).sem;
        all_data{g} = diff_data.(group_name).differences;
    end
    
    % Create bar plot
    b = bar(1:num_groups, means, 0.6, 'FaceColor', 'flat');
    
    % Set colors for each bar
    for g = 1:num_groups
        b.CData(g,:) = colors(g,:);
    end
    
    hold on;
    
    % Add error bars (SEM) with styling matching your preferences
    for g = 1:num_groups
        errorbar(g, means(g), sems(g), symbols{g}, ...
            'LineWidth', 1.5, 'MarkerSize', 10, ...
            'Color', colors(g,:), ...
            'MarkerFaceColor', colors(g,:), ...
            'MarkerEdgeColor', 'none', ...
            'CapSize', 10);
    end
    
    % Add individual data points with jitter using group-specific symbols and colors
    for g = 1:num_groups
        x_pos = g + 0.15 * (rand(size(all_data{g})) - 0.5); % Add jitter
        
        % Determine color for this group's scatter points
        if use_group_colors
            current_scatter_color = colors(g,:);
        else
            current_scatter_color = scatter_color;
        end
        
        % Extract just the marker symbol (e.g., 'o' from '-o')
        marker_symbol = extract_marker(symbols{g});
        
        % Plot scatter points with group-specific symbol
        scatter(x_pos, all_data{g}, 50, current_scatter_color, marker_symbol, ...
                'filled', 'MarkerFaceAlpha', p.Results.ScatterAlpha);
    end
    
    % Add horizontal line at y=0
    plot([0.5, num_groups+0.5], [0, 0], 'k--', 'LineWidth', 1);
    
    hold off;
    
    % Formatting
    xlim([0.5, num_groups+0.5]);
    xticks(1:num_groups);
    xticklabels(group_order);
    ylabel(p.Results.YLabel);
    
    % Set Y-axis limits if provided
    if ~isempty(p.Results.YLim)
        ylim(p.Results.YLim);
    end
    
    if ~isempty(p.Results.Title)
        title(p.Results.Title);
    end
    
    % Add sample sizes to x-axis labels
    xticklabels_with_n = cell(1, num_groups);
    for g = 1:num_groups
        group_name = group_order{g};
        xticklabels_with_n{g} = sprintf('%s\n(n=%d)', group_name, diff_data.(group_name).n);
    end
    xticklabels(xticklabels_with_n);
    
    % Display summary statistics
    fprintf('\n=== WITHIN-SUBJECT LASER DIFFERENCES ===\n');
    fprintf('Measure: %s\n\n', measure_name);
    for g = 1:num_groups
        group_name = group_order{g};
        fprintf('%s (n=%d):\n', group_name, diff_data.(group_name).n);
        fprintf('  Mean difference: %.3f ± %.3f (SEM)\n', ...
                diff_data.(group_name).mean, diff_data.(group_name).sem);
        fprintf('  Range: [%.3f, %.3f]\n', ...
                min(diff_data.(group_name).differences), ...
                max(diff_data.(group_name).differences));
        fprintf('\n');
    end
end