

path_to_save_data = 'D:\MATLAB\my_repo\context fear\data_for_maddy_decoding';


animalIDs = fieldnames(final);

animals_data_to_grab = {'B51618'};
sessions_to_grab = {'D1_Afternoon'};
for i = 1:numel(animals_data_to_grab)
    animalID = animals_data_to_grab{i};

    animalFolder = fullfile(path_to_save_data, animalID);

    if ~exist(animalFolder, 'dir')
        mkdir(animalFolder);
    end

    sessions = fieldnames(final.(animalID));
    for j = 1:numel(sessions_to_grab)
        session = sessions_to_grab{j};
        sessionFolder = fullfile(animalFolder, session);

        if ~exist(sessionFolder, 'dir')
            mkdir(sessionFolder);
        end


        % if isfield(final.(animalID).(session), 'uv')
        %     % Save BehavData to CSV
        %     behavData = final_behavior.(animalID).(session).uv.BehavData;
        %     behavDataFile = fullfile(sessionFolder, ['BehavData_', animalID, '_', session, '.csv']);
        %     writetable(behavData, behavDataFile);
        % end

        % if isfield(final.(animalID).(session), 'CNMFe_data')
        % time_array = final.(animalID).(session).time';

        C_raw = final.(animalID).(session).CNMFe_data.C_raw;
        C = final.(animalID).(session).CNMFe_data.C;
        % if size(time_array, 2) < size(C_raw, 2)
        %
        %     time_array = [0 time_array];
        % elseif size(time_array, 2) > size(C_raw, 2)
        %     time_array = time_array(2, :);
        %
        % end

        C_raw_File = fullfile(sessionFolder, ['C_raw_', animalID, '_', session, '.csv']);
        writematrix(C_raw, C_raw_File);


        C_File = fullfile(sessionFolder, ['C_', animalID, '_', session, '.csv']);
        writematrix(C, C_File);
    end


end


%%