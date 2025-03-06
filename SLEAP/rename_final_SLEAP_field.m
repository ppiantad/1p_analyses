% Get top-level field names (animal IDs)
animalIDs = fieldnames(final_SLEAP);

for i = 1:numel(animalIDs)
    animalID = animalIDs{i};
    
    % Get second-level field names (behavioral sessions)
    sessions = fieldnames(final_SLEAP.(animalID));
    
    for j = 1:numel(sessions)
        session = sessions{j};
        
        % Check if the session is named "RDT_D3_SAL"
        if strcmp(session, 'RDT_D3_SAL')
            % Rename the session by copying data and deleting the old field
            final_SLEAP.(animalID).('RDT_D3_SALINE') = final_SLEAP.(animalID).(session);
            final_SLEAP.(animalID) = rmfield(final_SLEAP.(animalID), session);
        end
    end
end
