% Load aggregated dataset, which should contain cell coordinates (assuming
% .mat file is from later in 2023) for each session
% Run the eventRelatedActivity etc script to
% get cell activity indices, then plot the FOV and color code accordingly. 




% String to compare
targetAnimal = 'BLA_Insc_25';

% Perform element-wise comparison
animal_index_to_plot = find(strcmp(animalIDs, targetAnimal));



Coor = final.(targetAnimal).(session_to_analyze).CNMFe_data.Coor;  

% these variables need to be customized based on what you want to plot
peri_choice_activated = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) == 1;
consumption_activated = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) == 1;
neutral = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) ~= 1 & respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) ~= 1;

% Create a combined array
Identity = zeros(size(Coor));

% Assign values based on conditions
Identity(peri_choice_activated) = 1;
Identity(consumption_activated) = 2;
Identity(neutral) = 3;



%%
figure;
hold on;



for i = 1:numel(Coor)
    % Get the coordinates of the current circle
    circleCoords = Coor{i};
    
    % Extract x and y coordinates
    y = circleCoords(1, :);
    x = circleCoords(2, :);

    % Calculate centroid
    centroid = [mean(x), mean(y)];

    % Store centroid in the array
    centroids(i, :) = round(centroid);
    

    if evalin( 'base', 'exist(''respClass_mouse'', ''var'') == 1')
        if peri_choice_activated(i) == 1
            plot_color = "red";
        elseif consumption_activated(i) == 1
            plot_color = "blue";
        elseif neutral(i) == 1
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

%% from ChatGPT, attempting to get it to model what Luthi does
% Initialize an array to store pairwise distances
pairwiseDistances = zeros(numel(Coor));

% Calculate pairwise Euclidean distances using a loop
for i = 1:numel(Coor)
    for j = i+1:numel(Coor)
        % Euclidean distance formula
        pairwiseDistances(i, j) = sqrt(sum((centroids(i, :) - centroids(j, :)).^2));
        % Since distances are symmetric, we can fill both sides of the matrix
        pairwiseDistances(j, i) = pairwiseDistances(i, j);
    end
end


% Plot the resulting distance metric
figure; imagesc(pairwiseDistances, [0 max(pairwiseDistances(:))]);

%% from Luthi lab code
% now quantify distance between Identities
Distances = nan(numel(unique(Identity)), numel(unique(Identity)));% preallocate space
for k = 1:numel(unique(Identity)) % run through all identities
    IDSubselect = pairwiseDistances(Identity == k, :); %select all distances of that ID
    for kk = 1:numel(unique(Identity)) 
        % now of that subset calculate the mean distance of all neurons to
        % all other neurons. This step is variable and people do different
        % things here, so maybe read the literature a bit. If you want to
        % quantify, repeat this step for every animal that leaves you with
        % a 3d instead of a 2d matrix (ID, ID, animal) and test the
        % differences against each other
        Distances(k, kk) = nanmean(nanmean(IDSubselect(:, Identity == kk)));
    end
end


% Now plot the mean distance metric
figure; imagesc(Distances, [0 max(pairwiseDistances(:))]);
colormap(bluewhitered); set(gca, 'YTick', 1:3, 'XTick', 1:3, ...
    'XTickLabel', {'ID1', 'ID2', 'ID3'}, 'YTickLabel', {'ID1', 'ID2', 'ID3'});
box off;colorbar


%%
% Assuming you already have centroids calculated and stored in 'centroids'
numCategories = max(Identity);

% Initialize arrays to store distances
withinIdentityDistances = zeros(size(centroids, 1), size(centroids, 1), numCategories); % 3rd dimension is for Identity categories
acrossIdentityDistances = zeros(size(centroids, 1), size(centroids, 1), numCategories);
neuron = 0;
% Loop through each identity category for within-identity distances
for identityCategory = 1:numCategories
    % Find indices of points in the current identity category
    withinIdentityIndices = find(Identity == identityCategory);
    
    % Extract coordinates of points within the current identity category
    withinIdentityCoords = centroids(withinIdentityIndices, :);
    withinIdentityIndices_array(:, identityCategory) = {withinIdentityIndices};
    for tt = 1:size(withinIdentityCoords, 1)
        neuron = neuron + 1;
        for qq = 1:size(withinIdentityCoords, 1)
            % Calculate pairwise distances within the current identity category
            withinIdentityDistances(tt, qq) = pdist([centroids(tt, :); centroids(qq,:)], 'euclidean');
            withinIdentityDistances_array(:, identityCategory) = {withinIdentityDistances}; 
            withinIdentityDistances_mean(neuron) = mean(withinIdentityDistances(tt, qq));
            % DistanceMatrix(k, kk) = pdist([Centroids(k, :); Centroids(kk, :)], 'euclidean'); 
        end
    end
    clear withinIdentityDistances withinIdentityIndices withinIdentityCoords
end



%%
% Assuming you already have centroids calculated and stored in 'centroids'
numNeurons = size(centroids, 1);


% Loop through each neuron
for neuron = 1:numNeurons
    % Get the Identity value for the current neuron
    currentIdentity = Identity(neuron);
    
    % Find indices of neurons with different Identity values
    otherIdentityIndices = find(Identity ~= currentIdentity);
    

    % add the index for the neuron so that there is a 0 in the array (same
    % as withinIdentityDistances - can always filter out the non-zero
    % values later when calculating the mean! 
    otherIdentityIndices = [neuron; otherIdentityIndices];

    % Extract coordinates of the current neuron
    currentNeuronCoords = centroids(neuron, :);
    
    % Extract coordinates of neurons with different Identity values
    otherIdentityCoords = centroids(otherIdentityIndices, :);
    
    % Calculate pairwise distances between the current neuron and neurons with different Identity values
    % acrossIdentityDistancesPerNeuron(neuron, :) = pdist2(currentNeuronCoords, otherIdentityCoords);
    acrossIdentityDistances_array(:, neuron) = {pdist2(currentNeuronCoords, otherIdentityCoords,'euclidean')}; 
    acrossIdentityDistances_mean(neuron) = mean(acrossIdentityDistances_array{1, neuron});
end


%%
% Assuming you have two arrays: withinIdentityDistances_mean and acrossIdentityDistances_mean
% Also, assuming you have an array called Identity
% Scatter plot

% Scatter plot
figure;
scatter(acrossIdentityDistances_mean, withinIdentityDistances_mean, [], Identity, 'filled');
hold on;

% Plotting the diagonal line
plot([min(acrossIdentityDistances_mean), max(acrossIdentityDistances_mean)], [min(withinIdentityDistances_mean), max(withinIdentityDistances_mean)], '--k');  % Dotted line

colormap('jet');  % You can use a different colormap if desired
title('Scatter Plot of Distances');
xlabel('Across-Identity Distances Mean');
ylabel('Within-Identity Distances Mean');
colorbar;

hold off;

%% gathering ALL cells
animalIDs = (fieldnames(final));
Coor_all_mice = {};
for qq = 1:size(fieldnames(final),1)
    animal = animalIDs{qq};
    Coor_indiv = final.(animal).(session_to_analyze).CNMFe_data.Coor;
    Coor_all_mice = [Coor_all_mice; Coor_indiv ];
end

%% plot ALL cells
figure;
hold on;



for i = 1:numel(Coor_all_mice)
    % Get the coordinates of the current circle
    circleCoords = Coor_all_mice{i};
    
    % Extract x and y coordinates
    y = circleCoords(1, :);
    x = circleCoords(2, :);

    % Calculate centroid
    centroid = [mean(x), mean(y)];

    % Store centroid in the array
    centroids(i, :) = round(centroid);
    
    % currently need to load a mouse's data, and select these filters etc.
    % manually. can probably integrate with the
    % eventRelatedActivitAndClassification script to make it more automatic
    if evalin( 'base', 'exist(''respClass_mouse'', ''var'') == 1')
        if respClass_all_array{1, 1}(1, i) == 1
            plot_color = "red";
        elseif respClass_all_array{1, 2}(1, i) == 2
            plot_color = "blue";
        elseif respClass_all_array{1, 1}(1, i) ~= 1 & respClass_all_array{1, 2}(1, i) ~= 1
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


%% withinIdentity for ALL cells
% Assuming you already have centroids calculated and stored in 'centroids'
numCategories = max(respClass_all_array{1, 1});

% Initialize arrays to store distances
% withinIdentityDistances = zeros(size(centroids, 1), size(centroids, 1), numCategories); % 3rd dimension is for Identity categories
% acrossIdentityDistances = zeros(size(centroids, 1), size(centroids, 1), numCategories);
neuron = 0;
% Loop through each identity category for within-identity distances
for identityCategory = 1:numCategories
    % Find indices of points in the current identity category
    withinIdentityIndices = find(respClass_all_array{1, 1} == identityCategory);
    
    % Extract coordinates of points within the current identity category
    withinIdentityCoords = centroids(withinIdentityIndices, :);
    withinIdentityIndices_array(:, identityCategory) = {withinIdentityIndices};
    for tt = 1:size(withinIdentityCoords, 1)
        neuron = neuron + 1;
        for qq = 1:size(withinIdentityCoords, 1)
            % Calculate pairwise distances within the current identity category
            withinIdentityDistances(tt, qq) = pdist([centroids(tt, :); centroids(qq,:)], 'euclidean');
            withinIdentityDistances_array(:, identityCategory) = {withinIdentityDistances}; 
            withinIdentityDistances_mean(neuron) = mean(withinIdentityDistances(tt, qq));
            % DistanceMatrix(k, kk) = pdist([Centroids(k, :); Centroids(kk, :)], 'euclidean'); 
        end
    end
    clear withinIdentityDistances withinIdentityIndices withinIdentityCoords
end

%% acrossIdentity for ALL cells
% Assuming you already have centroids calculated and stored in 'centroids'
numNeurons = size(centroids, 1);


% Loop through each neuron
for neuron = 1:numNeurons
    % Get the Identity value for the current neuron
    currentIdentity = respClass_all_array{1, 1}(neuron);
    
    % Find indices of neurons with different Identity values
    otherIdentityIndices = find(respClass_all_array{1, 1} ~= currentIdentity)';
    

    % add the index for the neuron so that there is a 0 in the array (same
    % as withinIdentityDistances - can always filter out the non-zero
    % values later when calculating the mean! 
    otherIdentityIndices = [neuron; otherIdentityIndices];

    % Extract coordinates of the current neuron
    currentNeuronCoords = centroids(neuron, :);
    
    % Extract coordinates of neurons with different Identity values
    otherIdentityCoords = centroids(otherIdentityIndices, :);
    
    % Calculate pairwise distances between the current neuron and neurons with different Identity values
    % acrossIdentityDistancesPerNeuron(neuron, :) = pdist2(currentNeuronCoords, otherIdentityCoords);
    acrossIdentityDistances_array(:, neuron) = {pdist2(currentNeuronCoords, otherIdentityCoords,'euclidean')}; 
    acrossIdentityDistances_mean(neuron) = mean(acrossIdentityDistances_array{1, neuron});
end

%%
% Assuming you have two arrays: withinIdentityDistances_mean and acrossIdentityDistances_mean
% Also, assuming you have an array called Identity
% Scatter plot

% Scatter plot
figure;
scatter(acrossIdentityDistances_mean, withinIdentityDistances_mean, [], respClass_all_array{1, 1}, 'filled');
hold on;

% Plotting the diagonal line
plot([min(acrossIdentityDistances_mean), max(acrossIdentityDistances_mean)], [min(withinIdentityDistances_mean), max(withinIdentityDistances_mean)], '--k');  % Dotted line

colormap('jet');  % You can use a different colormap if desired
title('Scatter Plot of Distances');
xlabel('Across-Identity Distances Mean');
ylabel('Within-Identity Distances Mean');
colorbar;

hold off;

%%
% Assuming you have two arrays: withinIdentityDistances_mean and acrossIdentityDistances_mean

% Calculate cumulative distribution functions (CDF)
[unique_within, ~, ic_within] = unique(withinIdentityDistances_mean);
cdf_within = cumsum(hist(ic_within, 1:max(ic_within)) / numel(withinIdentityDistances_mean));

[unique_across, ~, ic_across] = unique(acrossIdentityDistances_mean);
cdf_across = cumsum(hist(ic_across, 1:max(ic_across)) / numel(acrossIdentityDistances_mean));

% Combine the distances for the entire population
allDistances = [withinIdentityDistances_mean, acrossIdentityDistances_mean];
[unique_all, ~, ic_all] = unique(allDistances);
cdf_all = cumsum(hist(ic_all, 1:max(ic_all)) / numel(allDistances));

% Plot cumulative probability distributions
figure;
plot(unique_within, cdf_within, 'LineWidth', 2, 'DisplayName', 'Within Identity');
hold on;
plot(unique_across, cdf_across, 'LineWidth', 2, 'DisplayName', 'Across Identity');
plot(unique_all, cdf_all, 'LineWidth', 2, 'DisplayName', 'Entire Population');

title('Cumulative Probability Distribution of Mean Pairwise Distances');
xlabel('Mean Pairwise Distance');
ylabel('Cumulative Probability');
legend('Location', 'best');
grid on;