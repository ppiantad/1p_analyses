
uv.behav = {'stTime','choiceTime','collectionTime'}; %which behavior/timestamp to look at
animalIDs = (fieldnames(final));
neuron_num = 0;

paired_sessions = {'RM_D1', 'Pre_RDT_RM'};
paired_sessions_struct_name = [paired_sessions{1}, '_vs_', paired_sessions{2}];

load('cellRegistered_20230422_152526.mat');
idx = sort(cell_registered_struct.cell_to_index_map);

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



cellreg_struct = struct;

ts1 = (-10:.1:10);

for qq = 1:size(paired_sessions, 1)
    for zz = 1:size(paired_sessions(qq,:))
        paired_sessions_current = paired_sessions{qq, 1};
        for i = 1:length(uv.behav)
            alignment_event = char(uv.behav(i));
            for ii = 1:size(fieldnames(final),1)
                current_animal = char(animalIDs(ii));

                if isfield(final.(current_animal), paired_sessions{1})
                    current_session = paired_sessions{1};

                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).time = final.(current_animal).(current_session).(alignment_event).time; %final(i).time = caTime;
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).unitAVG.caTraces = final.(current_animal).(current_session).(alignment_event).unitAVG.caTraces(filtered_map(:,1),:);
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).unitXTrials = final.(current_animal).(current_session).(alignment_event).unitXTrials(1,filtered_map(:,1),:);
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).uv = final.(current_animal).(current_session).(alignment_event).uv;
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).unitSEM = final.(current_animal).(current_session).(alignment_event).unitSEM.caTraces(filtered_map(:,1),:);
                    %             new_struct.(current_animal).(current_session).(alignment_event).uv = final.(current_animal).(current_session).(alignment_event).uv;

                elseif isfield(final.(currentanimal), paired_sessions{2})
                    current_session = paired_sessions{2};
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).time = final.(current_animal).(current_session).(alignment_event).time; %final(i).time = caTime;
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).unitAVG.caTraces = final.(current_animal).(current_session).(alignment_event).unitAVG.caTraces(filtered_map(:,2),:);
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).unitXTrials = final.(current_animal).(current_session).(alignment_event).unitXTrials(1,filtered_map(:,2),:);
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).uv = final.(current_animal).(current_session).(alignment_event).uv;
                    cellreg_struct.(paired_sessions_struct_name).(current_animal).(current_session).(alignment_event).unitSEM = final.(current_animal).(current_session).(alignment_event).unitSEM.caTraces(filtered_map(:,2),:);



                    clear ;




                elseif ~isfield(final.(current_animal), session_to_analyze)


                end
            end
        end
    end
end