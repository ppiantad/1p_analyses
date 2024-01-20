function set_axes_limits(filtered_motion, selected_indices)
    % Initialize variables to store min and max values
    x_min = inf;
    x_max = -inf;
    y_min = inf;
    y_max = -inf;

    % Loop through selected indices to find min and max values
    for j = 1:length(selected_indices)
        x = filtered_motion{selected_indices(j)}(1, :);
        y = filtered_motion{selected_indices(j)}(2, :);

        x_min = min(x_min, min(x));
        x_max = max(x_max, max(x));
        y_min = min(y_min, min(y));
        y_max = max(y_max, max(y));
    end

    % Set axis limits based on min and max values
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
end