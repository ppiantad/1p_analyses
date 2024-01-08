%% Run eventRelatedActivityAndClassification first. This code will loop through the respClass_mouse array to find all the "activated" substructures. These can be used to calculate the % activated neurons per mouse, which is a necessary supplemental figure. 

[activatedStructs] = findAndSaveActivatedFields(respClass_mouse);


duplicatedAnimalIDs = cellfun(@(x) {x; x}, animalIDs, 'UniformOutput', false);

% Convert the cell array to a single column cell array
duplicatedAnimalIDs = vertcat(duplicatedAnimalIDs{:});

cells_per_mouse = arrayfun(@(row) size(neuron_mean_mouse{row, 2}, 1), 1:size(neuron_mean_mouse, 1), 'UniformOutput', false);

duplicate_cels_per_mouse = cellfun(@(x) {x; x}, cells_per_mouse, 'UniformOutput', false);
duplicate_cels_per_mouse = vertcat(duplicate_cels_per_mouse {:});

for bb = 1:size(activatedStructs, 2)
    activated_temp = activatedStructs{bb};
    sum_activated_mouse(bb) = sum(activated_temp == 1);
    
end


% Initialize a matrix for the stacked bar plot
stacked_data = zeros(14, iter+1); % Assuming 3 segments (adjust as needed)

% Create stacked bar plots for each mouse
figure;
for i = 1:14
    % Calculate the third segment (grey) for each mouse
    stacked_data(i, 3) = cells_per_mouse{i} - sum(mouse_data(i, 1:2));
    
    % Combine the mouse_data and stacked_data
    stacked_data(i, 1:2) = mouse_data(i, :);
    X = 1;
    % Plot the stacked bar
    subplot(7, 2, i);  % Adjust subplot layout as needed
    bar(X, stacked_data(i, :), 'stacked');
    title(['Mouse ', num2str(i)]);
end

function activatedStructs = findAndSaveActivatedFields(structure)
    % Get field names of the current structure
    fields = fieldnames(structure);
    
    % Initialize the activatedStructs cell array to store multiple structures
    activatedStructs = cell(0);

    % Check if the current structure has a field named 'activated'
    if isfield(structure, 'activated')
        disp('Found ''activated'' field:');
        disp(structure.activated);
        % Save the 'activated' substructure
        activatedStructs{end+1} = structure.activated;
        
        % Your custom action goes here
        % Replace the above lines with your desired action
    end
    
    % Recursively call the function for each field that is a structure
    for i = 1:length(fields)
        if isstruct(structure.(fields{i}))
            % Recursively call the function and concatenate activatedStructs
            activatedStructs = [activatedStructs, findAndSaveActivatedFields(structure.(fields{i}))];
        end
    end
end

