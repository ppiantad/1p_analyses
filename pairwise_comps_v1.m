test = [neuron_mean_array{1, 1}(respClass_all_array{1, 1}==1, :)];
test = [test; neuron_mean_array{1, 1}(respClass_all_array{1, 2}==1, :)];

%%
data = test;

% Initialize a matrix to store correlation coefficients
correlation_matrix = zeros(size(data, 1));

% Calculate correlation coefficients between rows
for i = 1:size(data, 1)
    for j = 1:size(data, 1)
        correlation_matrix(i, j) = corr(data(i, :)', data(j, :)');
    end
end

% Plot the correlation matrix
figure;
imagesc(correlation_matrix);
colorbar; % Add a colorbar to the plot
axis square; % Make the plot square for better visualization
title('Neuronal Correlation Matrix');
xlabel('Neuron Index');
ylabel('Neuron Index');

% Show row and column indices on the plot
xticks(1:size(data, 1));
yticks(1:size(data, 1));

% If you want to customize the color map, you can use colormap function
% For example, using a blue-white-red colormap:
colormap(bluewhitered);

% If you want to limit the color scale to the range [0, 1]
caxis([-1 1]); % Assu