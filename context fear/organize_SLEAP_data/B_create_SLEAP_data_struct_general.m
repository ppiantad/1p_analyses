% Define the top-level directory
top_level_directory = 'F:\Maddy Pilot\full_pilot'; % Update this to your top-level directory

SLEAP_and_freeze_combined_filepaths = dir(fullfile(top_level_directory, '**', '*SLEAP_and_freezing_combined*.csv'));

for k = 1:length(SLEAP_and_freeze_combined_filepaths)
    file_path = fullfile(SLEAP_and_freeze_combined_filepaths(k).folder, SLEAP_and_freeze_combined_filepaths(k).name);
    file_path_string = strsplit(file_path, '\');

    current_animal = file_path_string{4}; % 
    current_session_raw = file_path_string{5};

    current_session = strsplit(current_session_raw, '_');

    if strcmp(current_session{end}, 'EPM')
        current_session = current_session{end};

    elseif ~strcmp(current_session{end}, 'EPM')
        current_session = append(current_session{2}, '_', current_session{end});

    end

    if contains(lower(current_session), 'aft') || contains(lower(current_session), 'noon')
        if contains(lower(current_session), '2')
            current_session_renamed = 'D1_Afternoon';

        elseif contains(lower(current_session), '3')
            current_session_renamed = 'D2_Afternoon';
        end
    elseif contains(lower(current_session), 'morn')
        if contains(lower(current_session), '2')
            current_session_renamed = 'D1_Morning';

        elseif contains(lower(current_session), '3')
            current_session_renamed = 'D2_Morning';
        end
    elseif contains(lower(current_session), '4')
        current_session_renamed = 'D3';

    elseif contains(lower(current_session), '5')
        current_session_renamed = 'D4';
    else
        current_session_renamed = current_session;
    end

    
    if ~isempty(file_path)
        % Read the CSV file
        data = readtable(file_path);
        final_DLC.(current_animal).(current_session_renamed).movement_data = data;
        final_DLC.(current_animal).(current_session_renamed).filename = SLEAP_and_freeze_combined_filepaths(k).name;

    end
end


%%
experimental_grps = readtable('E:\MATLAB\my_repo\context fear\organize_DLC_data\PFC mice.xlsx');
animalIDs = fieldnames(final_DLC);

for dd = 1:size(experimental_grps, 1)
    current_mouse = experimental_grps.mouse{dd};
    current_mouse_condition = experimental_grps.group{dd};
    for hh = 1:size(animalIDs, 1)
        if strcmp(current_mouse, animalIDs(hh))
            final_DLC.(current_mouse).experimental_grp = current_mouse_condition;
            if any("sex" == string(experimental_grps.Properties.VariableNames))
                final_DLC.(current_mouse).sex = experimental_grps.sex{hh};
            end
        end

    end

end