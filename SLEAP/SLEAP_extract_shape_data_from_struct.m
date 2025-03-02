% Initialize the new struct
final_SLEAP_shapeData = struct();

% Get top-level field names (animal IDs)
animalIDs = fieldnames(final_SLEAP);

for i = 1:numel(animalIDs)
    animalID = animalIDs{i};
    
    % Get second-level field names (behavioral sessions)
    sessions = fieldnames(final_SLEAP.(animalID));
    
    for j = 1:numel(sessions)
        session = sessions{j};
        
        % Check if 'shapeData' exists directly at this level
        if isfield(final_SLEAP.(animalID).(session), 'shapeData')
            final_SLEAP_shapeData.(animalID).(session).shapeData = ...
                final_SLEAP.(animalID).(session).shapeData;
        end
        
        % Get third-level field names (data fields within the session)
        dataFields = fieldnames(final_SLEAP.(animalID).(session));
        
        for k = 1:numel(dataFields)
            field = dataFields{k};
            
            % Check if this subfield is a struct
            if isstruct(final_SLEAP.(animalID).(session).(field)) 
                % Check for 'shapeData' inside this struct
                if isfield(final_SLEAP.(animalID).(session).(field), 'shapeData')
                    final_SLEAP_shapeData.(animalID).(session).(field) = ...
                        final_SLEAP.(animalID).(session).(field).shapeData;
                end
            end
        end
    end
end

%% add shapeData back to UPDATED STRUCT (clear old final_SLEAP, load new one)

% Get top-level field names (animal IDs)
animalIDs = fieldnames(final_SLEAP_shapeData);

for i = 1:numel(animalIDs)
    animalID = animalIDs{i};
    
    % Get second-level field names (behavioral sessions)
    sessions = fieldnames(final_SLEAP_shapeData.(animalID));
    
    for j = 1:numel(sessions)
        session = sessions{j};
        
        % Check if shapeData exists at this level
        if isfield(final_SLEAP_shapeData.(animalID).(session), 'shapeData')
            final_SLEAP.(animalID).(session).shapeData = ...
                final_SLEAP_shapeData.(animalID).(session).shapeData;
        end
        

    end
end