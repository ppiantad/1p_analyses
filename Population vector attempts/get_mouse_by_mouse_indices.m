
if size(respClass_all_array, 2) == 10
    comparison_arrays = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays = [1 2 3; 4 5 6]
end



arrays_to_examine = [1 8];

inhib_or_excite = 1;

event_for_figures = 1; 


for hh = 1:size(respClass_all_array_mouse, 1)
    % respClass_all_array_mouse_prechoice_block_1_data = respClass_all_array_mouse{hh, 1};
    % respClass_all_array_mouse_prechoice_block_2_3_data = respClass_all_array_mouse{hh, 8};
    % 
    % respClass_all_array_mouse_postchoice_block_1_data = respClass_all_array_mouse{hh, 2};
    % respClass_all_array_mouse_postchoice_block_2_3_data = respClass_all_array_mouse{hh, 9};
    % 
    % respClass_all_array_mouse_consumption_block_1_data = respClass_all_array_mouse{hh, 3};
    % respClass_all_array_mouse_consumption_block_2_3_data = respClass_all_array_mouse{hh, 10};

    prechoice_block_1_not_block_2_3{hh} = respClass_all_array_mouse{hh, 1} == inhib_or_excite & respClass_all_array_mouse{hh, 8} ~= inhib_or_excite;
    prechoice_block_2_and_3_not_block_1{hh} = respClass_all_array_mouse{hh, 8} == inhib_or_excite & respClass_all_array_mouse{hh, 1} ~= inhib_or_excite;

    postchoice_block_1_not_block_2_3{hh} = respClass_all_array_mouse{hh, 2} == inhib_or_excite & respClass_all_array_mouse{hh, 9} ~= inhib_or_excite;
    postchoice_block_2_and_3_not_block_1{hh} = respClass_all_array_mouse{hh, 9} == inhib_or_excite & respClass_all_array_mouse{hh, 2} ~= inhib_or_excite;

    consumption_block_1_not_block_2_3{hh} = respClass_all_array_mouse{hh, 3} == inhib_or_excite & respClass_all_array_mouse{hh, 10} ~= inhib_or_excite;
    consumption_block_2_and_3_not_block_1{hh} = respClass_all_array_mouse{hh, 10} == inhib_or_excite & respClass_all_array_mouse{hh, 3} ~= inhib_or_excite;
end

prechoice_indices_for_PV = [prechoice_block_1_not_block_2_3; prechoice_block_2_and_3_not_block_1];

postchoice_indices_for_PV = [postchoice_block_1_not_block_2_3; postchoice_block_2_and_3_not_block_1];

consumption_indices_for_PV = [consumption_block_1_not_block_2_3; consumption_block_2_and_3_not_block_1];