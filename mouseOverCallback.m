function mouseOverCallback(~, eventdata)
    % Check if the mouse pointer is over the figure
    if ishandle(gcf)
        % Get the current mouse position
        currentPoint = get(gca, 'CurrentPoint');
        xCoord = currentPoint(1, 1);
        yCoord = currentPoint(1, 2);

        % Display the X and Y coordinates
        disp(['X: ' num2str(xCoord) ', Y: ' num2str(yCoord)]);
    end
end