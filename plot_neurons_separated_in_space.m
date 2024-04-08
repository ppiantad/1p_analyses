% Define spacing and line thickness
row_spacing = 5; % Adjust as needed
line_thickness = 1.0; % Adjust as needed

% Create a figure
figure;

% Loop through each row and plot the neural activity for the combined neurons
for i = 1:size(caTraceTrials, 1)
    % Extract neural activity for the current neuron
    neuron_activity = caTraceTrials(i, :);
    
    % Calculate the y-coordinate for the current neuron's plot
    y_coordinate = (i - 1) * row_spacing;
    
    % Plot the neural activity with the specified x and y coordinates
    plot(ts1, neuron_activity + y_coordinate, 'k', 'LineWidth', line_thickness);
    
    hold on; % To overlay all plots on the same figure
end

figure;
null_data = nullDistTrace(1:30, :);

for i = 1:size(null_data, 1)
    % Extract neural activity for the current neuron
    neuron_activity = null_data(i, :);
    
    % Calculate the y-coordinate for the current neuron's plot
    y_coordinate = (i - 1) * row_spacing;
    
    % Plot the neural activity with the specified x and y coordinates
    plot(ts1, neuron_activity + y_coordinate, 'k', 'LineWidth', line_thickness);
    
    hold on; % To overlay all plots on the same figure
end