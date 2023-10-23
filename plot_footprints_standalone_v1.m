% Load "_final.mat" for a given session from a single mouse
% Can load that session, then run the eventRelatedActivity etc script to
% get cell activity indices, then plot the FOV and color code accordingly. 
figure;
hold on;

for i = 1:numel(Coor)
    % Get the coordinates of the current circle
    circleCoords = Coor{i};
    
    % Extract x and y coordinates
    y = circleCoords(1, :);
    x = circleCoords(2, :);
    
    % currently need to load a mouse's data, and select these filters etc.
    % manually. can probably integrate with the
    % eventRelatedActivitAndClassification script to make it more automatic
    if evalin( 'base', 'exist(''respClass_mouse'', ''var'') == 1')
        if respClass_mouse.BLA_Insc_30.RDT_D1.choiceTime.Outcome_0to2.SHK_1.activated(i) == 1
            plot_color = "red";
        elseif respClass_mouse.BLA_Insc_30.RDT_D1.choiceTime.Outcome_0to2.SHK_1.inhibited(i) == 1
            plot_color = "blue";
        elseif respClass_mouse.BLA_Insc_30.RDT_D1.choiceTime.Outcome_0to2.SHK_1.neutral(i) == 1
            plot_color = "black";
        end
    else
    end
    % Plot the circle
    plot(x, y, 'Color', plot_color);
end

hold off;

% Add labels and title
xlabel('X');
ylabel('Y');
% title('Plot of Circles');

% Adjust the aspect ratio if needed
axis equal;