% This code currently will utilize data loaded in to the script
% load_data_4_cellreg.m and will loop through the sessions & CellReg files
% to create or add to the cellreg_struct data structure. 

load('BLA_panneuronal_Risk_2023_07_06.mat') %Load completely unmatched dataset - BLA panneuronal

load('BLA_panneuronal_Risk_MATCHED_04232023.mat') %Load matched dataset (to be added to)

%%

% cellreg_struct = struct;

ts1 = (-10:.1:10);

for qq = 1:size(paired_sessions, 1)
    load(cell_reg_data{qq});
%     idx = sort(cell_registered_struct.cell_to_index_map);
    % Extract the two columns of data
    column1 = cell_registered_struct.cell_to_index_map(:, 1);
    column2 = cell_registered_struct.cell_to_index_map(:, 2);
    % Sort the first column and save the indices
    [sorted_column1, sorted_indices] = sort(column1);
    % Rearrange the second column using the saved indices
    sorted_column2 = column2(sorted_indices);
    % Combine the sorted columns into a new array
    sorted_map = [sorted_column1, sorted_column2];
    % Filter out rows where either column has a value of zero
    nonzero_rows = (sorted_map(:,1) ~= 0) & (sorted_map(:,2) ~= 0);
    % Apply the filter to the sorted map
    filtered_map = sorted_map(nonzero_rows,:);
    paired_sessions_struct_name = [paired_sessions{qq, 1}, '_vs_', paired_sessions{qq, 2}];

    for zz = 1:size(paired_sessions(qq,:), 2)
        paired_sessions_current = paired_sessions{qq, zz};
        
        for i = 1:length(uv.behav)
            alignment_event = char(uv.behav(i));
            for ii = 1:size(fieldnames(final),1)
%                 current_animal = char(animalIDs(ii));

                if isfield(final.(current_animal), paired_sessions_current)
%                     current_session = paired_sessions{1};

                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(paired_sessions_current).(alignment_event).time = final.(current_animal).(paired_sessions_current).(alignment_event).time; %final(i).time = caTime;
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(paired_sessions_current).(alignment_event).unitAVG.caTraces = final.(current_animal).(paired_sessions_current).(alignment_event).unitAVG.caTraces(filtered_map(:,zz),:);
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(paired_sessions_current).(alignment_event).unitXTrials = final.(current_animal).(paired_sessions_current).(alignment_event).unitXTrials(1,filtered_map(:,zz),:);
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(paired_sessions_current).(alignment_event).uv = final.(current_animal).(paired_sessions_current).(alignment_event).uv;
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(paired_sessions_current).(alignment_event).unitSEM = final.(current_animal).(paired_sessions_current).(alignment_event).unitSEM.caTraces(filtered_map(:,zz),:);
                    %             new_struct.(current_animal).(current_session).(alignment_event).uv = final.(current_animal).(current_session).(alignment_event).uv;

               

                elseif ~isfield(final.(current_animal), session_to_analyze)


                end
            end
        end
    end
    clear cell_registered_struct column1 column2 sorted_column1 sort_indices sorted_column2 sorted_map nonzero_rows filtered_map
end