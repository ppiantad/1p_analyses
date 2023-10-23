% Define the number of neurons responsive to each event and their overlap
responsive_to_event1 = 527;
responsive_to_event2 = 274;
overlap = 203; % Responsive to both events

% Calculate the number of neurons exclusively responsive to each event
exclusive_to_event1 = responsive_to_event1 - overlap;
exclusive_to_event2 = responsive_to_event2 - overlap;
K = [exclusive_to_event1 exclusive_to_event2 overlap]
% Create a Venn diagram
figure;
[H, S] = venn([exclusive_to_event1, exclusive_to_event2, overlap], 'FaceColor', {'r', 'g'});
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end
% Set labels and title
labels = {'Event 1', 'Event 2'};
title('Neurons Responsive to Events');
legend(labels, 'Location', 'Best');

% Optionally, you can save the figure as an image file
% saveas(gcf, 'venn_diagram.png');

% Display the figure
axis off;


responsive_to_event1 = total_modulated(1);
responsive_to_event2 = total_modulated(2);
overlap = total_co_modulated;

% Calculate the number of neurons exclusively responsive to each event
exclusive_to_event1 = responsive_to_event1 - overlap;
exclusive_to_event2 = responsive_to_event2 - overlap;
K = [exclusive_to_event1 exclusive_to_event2 overlap]
% Create a Venn diagram
figure;
[H, S] = venn([exclusive_to_event1, exclusive_to_event2, overlap], 'FaceColor', {'r', 'g'});
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end
% Set labels and title
labels = {'Event 1', 'Event 2'};
title('Neurons Responsive to Events');
legend(labels, 'Location', 'Best');

% Optionally, you can save the figure as an image file
% saveas(gcf, 'venn_diagram.png');

% Display the figure
axis off;



responsive_to_event1 = (total_modulated(1)/neuron_num)*100;
responsive_to_event2 = (total_modulated(2)/neuron_num)*100;
overlap = (total_co_modulated/neuron_num)*100;

% Calculate the number of neurons exclusively responsive to each event
exclusive_to_event1 = responsive_to_event1 - overlap;
exclusive_to_event2 = responsive_to_event2 - overlap;
K = [exclusive_to_event1 exclusive_to_event2 overlap]
% Create a Venn diagram
figure;
[H, S] = venn([exclusive_to_event1, exclusive_to_event2, overlap], 'FaceColor', {'r', 'g'});
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end
% Set labels and title
labels = {'Event 1', 'Event 2'};
title('Neurons Responsive to Events');
legend(labels, 'Location', 'Best');

% Optionally, you can save the figure as an image file
% saveas(gcf, 'venn_diagram.png');

% Display the figure
axis off;
