function create_auc_barplot_new(data_first, data_second, plot_indiv, title_str, ylabel_str, shapes, axis_downsample_factor, ylimits, num_pairs, plot_connecting_lines)
    
    % Set default number of pairs if not specified
    if nargin < 9
        num_pairs = size(data_first, 1); % Use all available rows by default
    end
    
    % Set default for connecting lines if not specified
    if nargin < 10
        plot_connecting_lines = 1; % Plot connecting lines by default
    end
    
    % Ensure num_pairs doesn't exceed available data
    num_pairs = min(num_pairs, min(size(data_first, 1), size(data_second, 1)));
    
    % Calculate means and SEMs for each row
    data_first = data_first/axis_downsample_factor;
    data_second = data_second/axis_downsample_factor;

    means1 = mean(data_first(1:num_pairs, :), 2);
    sems1 = std(data_first(1:num_pairs, :), 0, 2) / sqrt(size(data_first, 2));
    
    means2 = mean(data_second(1:num_pairs, :), 2);
    sems2 = std(data_second(1:num_pairs, :), 0, 2) / sqrt(size(data_second, 2));
    
    % Create figure - adjust width based on number of pairs
    fig_width = 150 + num_pairs * 75; % Dynamic width
    figure('Position', [100, 100, fig_width, 400]);
    
    % X positions for bars - dynamically spaced
    x_spacing = 1.5;
    x_pos = 1:x_spacing:(num_pairs * x_spacing);
    bar_width = 0.35;
    division_factor = 2;

    % Create grouped bar plot positions
    x1 = x_pos - bar_width/division_factor;
    x2 = x_pos + bar_width/division_factor;
    
    % Plot bars
    b1 = bar(x1, means1, bar_width, 'EdgeColor', 'g', 'LineWidth', 1);
    hold on;
    b2 = bar(x2, means2, bar_width, 'EdgeColor', 'r', 'LineWidth', 1);
    
    % Add error bars
    errorbar(x1, means1, sems1, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 5);
    errorbar(x2, means2, sems2, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 5);
    
    % Add individual data points as scatter
    if plot_indiv == 1
        for i = 1:num_pairs
            % File 1 data points
            y_data1 = data_first(i, :);
            x_scatter1 = x1(i) + (rand(size(y_data1)) - 0.5) * bar_width * 0.6;
            first_shape = shapes{1};
            scatter(x_scatter1, y_data1, 20, 'g', 'filled', first_shape);

            % File 2 data points
            y_data2 = data_second(i, :);
            x_scatter2 = x2(i) + (rand(size(y_data2)) - 0.5) * bar_width * 0.6;
            second_shape = shapes{2};
            scatter(x_scatter2, y_data2, 20, 'r', 'filled', second_shape);
        end
    end

    % Connect paired points within each dataset
    if plot_connecting_lines == 1
        if num_pairs == 1
            % Special case: connect the two bars within the single pair
            for i = 1:size(data_first, 2)
                point1 = data_first(1, i);
                point2 = data_second(1, i);
                if ~isnan(point1) && ~isnan(point2)
                    plot([x1(1), x2(1)], [point1, point2], 'Color', [0.5 0.5 0.5], ...
                        'MarkerSize', 4, 'LineWidth', 0.8);
                end
            end
        else
            % Multiple pairs: connect adjacent pairs within each dataset
            for i = 1:size(data_first, 2)
                for j = 1:num_pairs-1
                    mouse_current = data_first(j, i);
                    mouse_next = data_first(j+1, i);
                    if ~isnan(mouse_current) && ~isnan(mouse_next)
                        plot([x1(j), x1(j+1)], [mouse_current, mouse_next], 'Color', 'g', ...
                            'MarkerSize', 4, 'LineWidth', 0.8, 'MarkerFaceColor', 'g');
                    end
                end
            end

            for i = 1:size(data_second, 2)
                for j = 1:num_pairs-1
                    mouse_current = data_second(j, i);
                    mouse_next = data_second(j+1, i);
                    if ~isnan(mouse_current) && ~isnan(mouse_next)
                        plot([x2(j), x2(j+1)], [mouse_current, mouse_next], 'Color', 'r', ...
                            'MarkerSize', 4, 'LineWidth', 0.8, 'MarkerFaceColor', 'r');
                    end
                end
            end
        end
    end
    
    % Formatting
    ylabel(ylabel_str);
    title(title_str);
    
    % Generate appropriate x-tick labels
    x_labels = cellstr(string(1:num_pairs));
    set(gca, 'XTick', x_pos, 'XTickLabel', x_labels);
    
    % Set axis properties
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    box on;
    ylim([ylimits])
    ytickformat('%.2f');
    hold off;
end