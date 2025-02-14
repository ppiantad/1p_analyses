function [x_fig, y_fig] = ds2nfu(x_data, y_data)
    ax = gca; % Get current axis
    x_limits = ax.XLim;
    y_limits = ax.YLim;
    
    % Convert data coordinates to normalized figure units (0 to 1)
    x_fig = (x_data - x_limits(1)) / (x_limits(2) - x_limits(1));
    y_fig = (y_data - y_limits(1)) / (y_limits(2) - y_limits(1));
end
