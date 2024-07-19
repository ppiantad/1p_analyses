% Create the 3D plot
figure;
plot3(PCScore{1, 1}(1, :), PCScore{1, 1}(2, :), PCScore{1, 1}(3, :), '-o'); % '-o' adds circle markers to the line
title('3D Plot of PCScore');
xlabel('X values (Row 1)');
ylabel('Y values (Row 2)');
zlabel('Z values (Row 3)');
grid on; % Optional: adds a grid to the plot

hold on; plot3(PCScore{1, 2}(1, :), PCScore{1, 2}(2, :), PCScore{1, 2}(3, :), '-o'); % '-o' adds circle markers to the line
% hold on; plot3(PCScore{1, 3}(1, :), PCScore{1, 3}(2, :), PCScore{1, 3}(3, :), '-o'); % '-o' adds circle markers to the line

%%