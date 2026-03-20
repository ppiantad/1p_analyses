% Load the two .mat files
data1 = load('Re_DREADD_XY_and_freeze_data_09052025.mat');
data2 = load('PL_DREADD_XY_and_freeze_data_05232025.mat');

% Extract the final_DLC structs from each file
final_DLC_1 = data1.final_DLC;
final_DLC_2 = data2.final_DLC;

% Initialize the combined struct
final_DLC_combined = struct();

% Get field names from both structs
fields1 = fieldnames(final_DLC_1);
fields2 = fieldnames(final_DLC_2);

% Combine fields from the first struct
for i = 1:length(fields1)
    final_DLC_combined.(fields1{i}) = final_DLC_1.(fields1{i});
end

% Add fields from the second struct
for i = 1:length(fields2)
    if isfield(final_DLC_combined, fields2{i})
        % If field already exists, you might want to rename it or handle the collision
        warning('Field %s already exists in combined struct. Adding with suffix _PL', fields2{i});
        final_DLC_combined.([fields2{i} '_PL']) = final_DLC_2.(fields2{i});
    else
        final_DLC_combined.(fields2{i}) = final_DLC_2.(fields2{i});
    end
end

% Save the combined struct with the original variable name
final_DLC = final_DLC_combined;

% Display the structure to verify
disp('Combined final_DLC structure:');
disp(fieldnames(final_DLC));

% Optionally save to a new file
% save('Combined_DREADD_XY_and_freeze_data.mat', 'final_DLC');