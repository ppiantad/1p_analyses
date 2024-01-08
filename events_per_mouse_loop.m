for qq = 1:size(animalIDs, 1)
    current_mouse = animalIDs{qq};
    mouse_fieldnames = fieldnames(respClass_mouse)
    if strcmp(current_mouse, mouse_fieldnames)
        for mm = 1:size(fieldnames(respClass_mouse.(current_mouse)))
            session_field = cell2mat(fieldnames(respClass_mouse.(current_mouse)(mm)))
            for zz = 1:size(fieldnames(respClass_mouse.(current_mouse).(session_field)), 1)
                epoch_fields = fieldnames(respClass_mouse.(current_mouse).(session_field)(zz))
                epoch_field_current = epoch_fields{zz};
                for pp = 1:size(fieldnames(respClass_mouse.(current_mouse).(session_field).(epoch_field_current)), 1)
                    
                end
            end
        end

    end
end



