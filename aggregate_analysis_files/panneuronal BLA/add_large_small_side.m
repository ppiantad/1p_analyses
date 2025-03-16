
imaging_group = {...,
'BLA_Insc_24';
'BLA_Insc_25';
'BLA_Insc_26';
'BLA_Insc_27';
'BLA_Insc_30';
'BLA_Insc_32';
'BLA_Insc_34';
'BLA_Insc_35';
'BLA_Insc_37';
'BLA_Insc_40';
};






large_screen_side = {...,
'left';
'left';
'right';
'right';
'left';
'left';
'right';
'left';
'left';
'left';
};

%%
% Get all first-level field names (animal IDs) from final_SLEAP
animalIDs = fieldnames(final_SLEAP);

% Loop through each animal ID
for i = 1:numel(animalIDs)
    animalID = animalIDs{i};
    
    % Find index of the animal ID in imaging_group
    match_idx = find(strcmp(imaging_group, animalID));
    
    % If a match is found, assign the corresponding large_screen_side value
    if ~isempty(match_idx)
        final_SLEAP.(animalID).large_screen_side = large_screen_side{match_idx};
    end
end

