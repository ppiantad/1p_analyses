function plot_line(filtered_motion, normalized_velocity_bounded, plot_index, title_str)
    % Extract X-Y coordinates
    x = filtered_motion{plot_index}(1, :);
    y = filtered_motion{plot_index}(2, :);

    % Extract velocity data
    velocity_values = normalized_velocity_bounded{plot_index};

    % Normalize velocity_values to the range [0, 1] for colormap mapping
    normalized_velocity = (velocity_values - min(velocity_values)) / (max(velocity_values) - min(velocity_values));

    % Set marker colors based on normalized_velocity
    marker_colors = colormap('jet');
    marker_colors_mapped = interp1(linspace(0, 1, size(marker_colors, 1)), marker_colors, normalized_velocity);

    % Plot the line using scatter with varying marker colors
    scatter(x, y, 50, marker_colors_mapped, 'filled');
    plot(x, y, 'Color', 'k');

    % Add a black square at the start of the line
    scatter(x(1), y(1), 200, 'k', 's', 'filled');

    % Add a black triangle at the end of the line
    scatter(x(end), y(end), 200, 'k', '^', 'filled');

    xlabel('X');
    ylabel('Y');
    title(['Plot ' num2str(plot_index) ' - ' title_str]);
end