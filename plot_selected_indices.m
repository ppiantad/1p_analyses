function plot_selected_indices(filtered_motion, normalized_velocity_bounded, selected_indices, title_str)
    figure;

    % Set the colormap for the heatmap
    colormap('jet');

    % Loop through each line and plot
    for j = 1:length(selected_indices)
        subplot(1, length(selected_indices), j);

        % Extract X-Y coordinates
        x = filtered_motion{selected_indices(j)}(1, :);
        y = filtered_motion{selected_indices(j)}(2, :);

        % Extract velocity data
        velocity_values = normalized_velocity_bounded{selected_indices(j)};

        % Normalize velocity_values to the range [0, 1] for colormap mapping
        normalized_velocity = (velocity_values - min(velocity_values)) / (max(velocity_values) - min(velocity_values));

        % Set marker colors based on normalized_velocity
        marker_colors = colormap('jet');
        marker_colors_mapped = interp1(linspace(0, 1, size(marker_colors, 1)), marker_colors, normalized_velocity);

        % Plot the line using scatter with varying marker colors
        scatter(x, y, 50, marker_colors_mapped, 'filled');
        hold on;
        plot(x, y, 'Color', 'k');

        % Add a black square at the start of the line
        scatter(x(1), y(1), 200, 'k', 's', 'filled');

        % Add a black triangle at the end of the line
        scatter(x(end), y(end), 200, 'k', '^', 'filled');

        title(['Selected Indices ' title_str ' - Plot ' num2str(j)]);

        % Set consistent axis limits across subplots
        if j == 1
            xlim_auto = xlim;
            ylim_auto = ylim;
        else
            xlim(xlim_auto);
            ylim(ylim_auto);
        end
    end
end
