
data_to_load = {'BLA_C_raw_no_additional_filtering_Pre_RDT_RM_only_completed_sessions_zall_window_base_workspace_10_categories.mat', 'BLA_C_raw_no_additional_filtering_RDT_D1_only_completed_sessions_zall_window_base_workspace_10_categories.mat'}

for qq = 1:size(data_to_load, 2)
    
    load(data_to_load{qq})
    if size(respClass_all_array, 2) == 10 | size(respClass_all_array, 2) == 11
        comparison_arrays = [1 2 3; 8 9 10]
    elseif size(respClass_all_array, 2) == 6
        comparison_arrays = [1 2 3; 4 5 6]
    elseif size(respClass_all_array, 2) == 7
        comparison_arrays = [1 2 3; 5 6 7]
    end



    arrays_to_examine = [1 8];

    inhib_or_excite = 1;

    event_for_figures = 1;
    %%
    for kk = 1:size(animalIDs, 1)
        % prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
        prechoice_block_1_mouse{kk, :} = respClass_all_array_mouse{kk, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 3)} ~= inhib_or_excite;
        prechoice_blocks_2_and_3_mouse{kk, :} = respClass_all_array_mouse{kk, comparison_arrays(2, 1)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 2)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 3)} ~= inhib_or_excite;

        % postchoice_reward_block_1 = respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
        postchoice_reward_block_1_mouse{kk, :}  = respClass_all_array_mouse{kk, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 3)} ~= inhib_or_excite;
        postchoice_reward_blocks_2_and_3_mouse{kk, :}  = respClass_all_array_mouse{kk, comparison_arrays(2, 2)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 3)} ~= inhib_or_excite;

        % collect_block_1 = respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
        collect_block_1_mouse{kk, :}  = respClass_all_array_mouse{kk, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(1, 2)} ~= inhib_or_excite;
        collect_blocks_2_and_3_mouse{kk, :}  = respClass_all_array_mouse{kk, comparison_arrays(2, 3)} == inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 1)} ~= inhib_or_excite & respClass_all_array_mouse{kk, comparison_arrays(2, 2)} ~= inhib_or_excite;

    end




    %%
    for kk = 1:size(animalIDs, 1)
        % prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
        prechoice_conserved_mouse{kk, :} = prechoice_block_1_mouse{kk, :} == event_for_figures & prechoice_blocks_2_and_3_mouse{kk, :} == event_for_figures;
        prechoice_lost_mouse{kk, :} = prechoice_block_1_mouse{kk, :} == event_for_figures & prechoice_blocks_2_and_3_mouse{kk, :} ~= event_for_figures;
        prechoice_remapped_mouse{kk, :} = prechoice_block_1_mouse{kk, :} ~= event_for_figures & prechoice_blocks_2_and_3_mouse{kk, :} == event_for_figures;

        postchoice_conserved_mouse{kk, :} = postchoice_reward_block_1_mouse{kk, :} == event_for_figures & postchoice_reward_blocks_2_and_3_mouse{kk, :} == event_for_figures;
        postchoice_lost_mouse{kk, :} = postchoice_reward_block_1_mouse{kk, :} == event_for_figures & postchoice_reward_blocks_2_and_3_mouse{kk, :} ~= event_for_figures;
        postchoice_remapped_mouse{kk, :} = postchoice_reward_block_1_mouse{kk, :} ~= event_for_figures & postchoice_reward_blocks_2_and_3_mouse{kk, :} == event_for_figures;

        collection_conserved_mouse{kk, :} = collect_block_1_mouse{kk, :} == event_for_figures & collect_blocks_2_and_3_mouse{kk, :} == event_for_figures;
        collection_lost_mouse{kk, :} = collect_block_1_mouse{kk, :} == event_for_figures & collect_blocks_2_and_3_mouse{kk, :} ~= event_for_figures;
        collection_remapped_mouse{kk, :} = collect_block_1_mouse{kk, :} ~= event_for_figures & collect_blocks_2_and_3_mouse{kk, :} == event_for_figures;

    end

    %%
    for kk = 1:size(animalIDs, 1)


        conserved_sum(kk) = sum([sum(prechoice_conserved_mouse{kk, :}), sum(postchoice_conserved_mouse{kk, :}), sum(collection_conserved_mouse{kk, :})]);
        lost_sum(kk) = sum([sum(prechoice_lost_mouse{kk, :}), sum(postchoice_lost_mouse{kk, :}), sum(collection_lost_mouse{kk, :})]);
        remapped_sum(kk) = sum([sum(prechoice_remapped_mouse{kk, :}), sum(postchoice_remapped_mouse{kk, :}), sum(collection_remapped_mouse{kk, :})]);

        conserved_ratio(kk) = conserved_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);
        lost_ratio(kk) = lost_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);
        remapped_ratio(kk) = remapped_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);

        conserved_prechoice_sum(kk) = sum(prechoice_conserved_mouse{kk, :});
        remapped_prechoice_sum(kk) = sum(prechoice_remapped_mouse{kk, :});
        lost_prechoice_sum(kk) = sum(prechoice_lost_mouse{kk, :});

        conserved_prechoice_ratio(kk) = conserved_prechoice_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);
        remapped_prechoice_ratio(kk) = remapped_prechoice_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);
        lost_prechoice_ratio(kk) = lost_prechoice_sum(kk)/size(prechoice_conserved_mouse{kk, :}, 2);


        conserved_postchoice_sum(kk) = sum(postchoice_conserved_mouse{kk, :});
        remapped_postchoice_sum(kk) = sum(postchoice_remapped_mouse{kk, :});
        lost_postchoice_sum(kk) = sum(postchoice_lost_mouse{kk, :});

        conserved_postchoice_ratio(kk) = conserved_postchoice_sum(kk)/size(postchoice_conserved_mouse{kk, :}, 2);
        remapped_postchoice_ratio(kk) = remapped_postchoice_sum(kk)/size(postchoice_conserved_mouse{kk, :}, 2);
        lost_postchoice_ratio(kk) = lost_postchoice_sum(kk)/size(postchoice_conserved_mouse{kk, :}, 2);

        conserved_collection_sum(kk) = sum(collection_conserved_mouse{kk, :});
        remapped_collection_sum(kk) = sum(collection_remapped_mouse{kk, :});
        lost_collection_sum(kk) = sum(collection_lost_mouse{kk, :});

        conserved_collection_ratio(kk) = conserved_collection_sum(kk)/size(collection_conserved_mouse{kk, :}, 2);
        remapped_collection_ratio(kk) = remapped_collection_sum(kk)/size(collection_conserved_mouse{kk, :}, 2);
        lost_collection_ratio(kk) = lost_collection_sum(kk)/size(collection_conserved_mouse{kk, :}, 2);
        conserved_ratio_all(kk, qq) = conserved_ratio(kk);
        lost_ratio_all(kk, qq) = lost_ratio(kk);
        remapped_ratio_all(kk, qq) = remapped_ratio(kk);


        conserved_sum_all(kk, qq) = conserved_sum(kk);
        lost_sum_all(kk, qq) = lost_sum(kk);
        remapped_sum_all(kk, qq) = remapped_sum(kk);
    end


end