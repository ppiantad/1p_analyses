%% Function to create bar plot with error bars and scatter points
function create_auc_barplot(data_first, data_second, plot_indiv, title_str, ylabel_str, shapes, axis_downsample_factor, ylimits)
    
    % Calculate means and SEMs for each row
    data_first = data_first/axis_downsample_factor;
    data_second = data_second/axis_downsample_factor;

    means1 = mean(data_first, 2, 'omitnan');
    sems1 = std(data_first, 0, 2, 'omitnan') / sqrt(size(data_first, 2) - sum(any(isnan(data_first))));
    
    means2 = mean(data_second, 2, 'omitnan');
    sems2 = std(data_second, 0, 2, 'omitnan') / sqrt(size(data_second, 2)- sum(any(isnan(data_second))));
    
    % Create figure
    figure('Position', [100, 100, 200, 400]);
    
    % X positions for bars
    % x_pos = [1, 2];
    % bar_width = 0.35;
    % division_factor = 2;

    x_pos = [1, 2.5];
    bar_width = 0.65;
    division_factor = 2;

    % Create grouped bar plot
    x1 = [x_pos(1) - bar_width(1)/division_factor,  x_pos(1) + bar_width(1)/division_factor];
    x2 = [x_pos(2) - bar_width(1)/division_factor,  x_pos(2) + bar_width(1)/division_factor];

    % x1 = [x_pos - bar_width/division_factor];
    % x2 = [x_pos + bar_width/division_factor];

    % x_points = [x_pos - bar_width/2, x_pos + bar_width/2;]
    
    % Plot bars
    b1 = bar([x1], means1, bar_width, 'EdgeColor', 'g', 'LineWidth', 1);
    hold on;
    b2 = bar([x2], means2, bar_width, 'EdgeColor', 'r', 'LineWidth', 1);
    
    % Add error bars
    errorbar([x1], means1, sems1, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 5);
    errorbar([x2], means2, sems2, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 5);
    
    % Add individual data points as scatter
    if plot_indiv == 1
        for i = 1:2
            % File 1 data points
            y_data1 = data_first(i, :);
            x_scatter1 = x1(i) + (rand(size(y_data1)) - 0.5) * bar_width * 0.3; % Add jitter
            first_shape = shapes{1};
            scatter(x_scatter1, y_data1, 20, 'g', 'filled', first_shape);

            % File 2 data points
            y_data2 = data_second(i, :);
            x_scatter2 = x2(i) + (rand(size(y_data2)) - 0.5) * bar_width * 0.3; % Add jitter
            second_shape = shapes{2};
            scatter(x_scatter2, y_data2, 20, 'r', 'filled', second_shape);
        end
    else
    end


    for i = 1:size(data_first, 2)
        mouse_large = mean(data_first(1, i), 'omitnan');
        mouse_small = mean(data_first(2, i), 'omitnan');
        if ~isnan(mouse_large) && ~isnan(mouse_small)
            plot([x1(1), x1(2)], [mouse_large, mouse_small], 'Color', 'g', ...
                'MarkerSize', 4, 'LineWidth', 0.8, 'MarkerFaceColor', 'g');
        end
    end


    for i = 1:size(data_second, 2)
        mouse_large = mean(data_second(1, i), 'omitnan');
        mouse_small = mean(data_second(2, i), 'omitnan');
        if ~isnan(mouse_large) && ~isnan(mouse_small)
            plot([x2(1), x2(2)], [mouse_large, mouse_small], 'Color', 'r', ...
                'MarkerSize', 4, 'LineWidth', 0.8, 'MarkerFaceColor', 'r');
        end
    end
    % Formatting

    ylabel(ylabel_str);
    title(title_str);
    set(gca, 'XTick', x_pos, 'XTickLabel', {'1', '2', '3'});
    % legend([b1, b2], {'File 1', 'File 2'}, 'Location', 'best');
    % grid on;
    % grid alpha 0.3;
    
    % Set axis properties
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    box on;
    ylim([ylimits])
    ytickformat('%.2f');
    hold off;
end
