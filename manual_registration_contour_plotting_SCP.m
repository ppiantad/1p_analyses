%Drag in final file from the first session then run this section
CoorD1 = Coor;
CnD1 = neuron.Cn;
figure;
imagesc(CnD1) %only one session will look good here typically

%% %Drag in final file from the first session then run this section

CoorD2 = Coor;
CnD2 = neuron.Cn;
figure;
imagesc(CnD2) %only one session will look good here typically

%%
cell_registered_struct_temp = struct;

cell_registered_struct_temp.cell_to_index_map = (1:length(CoorD1(:,1)))';

cell_registered_struct_temp.cell_to_index_map = [cell_registered_struct_temp.cell_to_index_map; zeros(abs(length(CoorD1(:,1)) - length(CoorD2(:,1))), 1)];

cell_registered_struct_temp.cell_to_index_map(:, 2) = zeros(size(cell_registered_struct_temp.cell_to_index_map, 1), 1);

%%



% Assuming you have already plotted the neurons as described in your code

% Create a figure
figure;
hold on;
imagesc(neuron.Cn) %only one session will look good here typically

% Initialize variables to store clicked cell indices
selectedCellSession1 = [];
selectedCellSession2 = [];
% Plot Session 1 cells
for m = 1:length(CoorD1)
    h_fill1(m) = fill(CoorD1{m, 1}(1,:), CoorD1{m, 1}(end,:), 'blue');
    text(CoorD1{m,1}(1,1), CoorD1{m, 1}(end,1), ['Sess1' ' ' num2str(m)], 'FontSize', 16);
    set(h_fill1(m), 'ButtonDownFcn', @(src, event) cellClickCallback(src, event, m, CoorD1, 'Session1'));
end

% Plot Session 2 cells
for m = 1:length(CoorD2)
    h_fill2(m) = fill(CoorD2{m, 1}(1,:), CoorD2{m, 1}(end,:), 'red');
    text(CoorD2{m,1}(1,1), CoorD2{m, 1}(end,1), ['Sess2' ' ' num2str(m)], 'FontSize', 16);
    set(h_fill2(m), 'ButtonDownFcn', @(src, event) cellClickCallback(src, event, m, CoorD2, 'Session2'));
end



%%
% Now need to manually assign values that MATCH based on the image above
% Change the value in (this spot, 2) according to map and update = this spot
cell_registered_struct_temp.cell_to_index_map(6, 2) = 4;

%%
% use the below line if you want to try clicking on the cells using the UI,
% still a work in progress
cell_registered_struct_temp.cell_to_index_map(selectedCellSession1, 2) = selectedCellSession2;

%% set the rest of the cells that aren't matched to the 0 spots that were created in Column 1 (Session 1)

session_1_cell_num = (length(CoorD1(:,1)))';

session_2_cell_IDs_total = (1:length(CoorD2(:,1)))';

% Create a logical index for values greater than 0
logical_index = cell_registered_struct_temp.cell_to_index_map(:, 2) > 0;

% Use the logical index to extract values
values_greater_than_zero = cell_registered_struct_temp.cell_to_index_map(logical_index, 2);

session_2_cell_IDs_missing = setdiff(session_2_cell_IDs_total, values_greater_than_zero);

% Check the condition for each row and update the array if needed
for row = 1:size(cell_registered_struct_temp.cell_to_index_map, 1)
    if cell_registered_struct_temp.cell_to_index_map(row, 1) ~= 0 && cell_registered_struct_temp.cell_to_index_map(row, 2) == 0
        % Move the value from column 1 to the final row + 1
        cell_registered_struct_temp.cell_to_index_map(end+1, 1) = cell_registered_struct_temp.cell_to_index_map(row, 1);
        
        % Add a 0 in column 2 of the final row + 1
        cell_registered_struct_temp.cell_to_index_map(end, 2) = 0;
        cell_registered_struct_temp.cell_to_index_map(row, :) = 0;
        % You can break here if you want to handle one such case at a time
        % break;
    end
end


% Iterate through session_2_cell_IDs_missing and fill them into column 2
missing_index = 1;
for row = 1:size(cell_registered_struct_temp.cell_to_index_map, 1)
    if cell_registered_struct_temp.cell_to_index_map(row, 1) == 0 && cell_registered_struct_temp.cell_to_index_map(row, 2) == 0
        % If both column 1 and column 2 are 0, add the value from session_2_cell_IDs_missing
        if missing_index <= numel(session_2_cell_IDs_missing)
            cell_registered_struct_temp.cell_to_index_map(row, 2) = session_2_cell_IDs_missing(missing_index);
            missing_index = missing_index + 1;
        else
            % Break if there are no more missing values
            break;
        end
    end
end

% Create a logical index to identify rows where both column 1 and column 2 are not 0
logical_index = ~((cell_registered_struct_temp.cell_to_index_map(:, 1) == 0) & (cell_registered_struct_temp.cell_to_index_map(:, 2) == 0));

% Use the logical index to keep only the rows where the condition is not met
cell_registered_struct_temp.cell_to_index_map = cell_registered_struct_temp.cell_to_index_map(logical_index, :);



cell_registered_struct = cell_registered_struct_temp;

%%
% Create the desired filename
filename = sprintf('cellRegistered_%s_manualReg.mat', datestr(now, 'yyyymmdd'));

% Open a save dialog box to choose the location and filename
[filename, path] = uiputfile(filename, 'Save cell_registered_struct_temp as .mat file');

% Check if the user canceled the operation
if isequal(filename,0) || isequal(path,0)
    disp('User canceled the save operation');
else
    % Save the array with the chosen filename
    full_filename = fullfile(path, filename);
    save(full_filename, 'cell_registered_struct');
    disp(['Saved as ' full_filename]);

    % Save the figure as a .png with the desired filename
    png_filename = sprintf('footprints_overlaid_%s.png', datestr(now, 'yyyymmdd'));
    full_png_filename = fullfile(path, png_filename);
    saveas(gcf, full_png_filename);
    fig_filename = sprintf('footprints_overlaid_%s.fig', datestr(now, 'yyyymmdd'));
    full_fig_filename = fullfile(path, fig_filename);
    savefig(gcf, full_fig_filename);
    

end

