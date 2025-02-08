% Pre_RDT_RM & RDT D1 matched
% RDT D1
% collectionTime
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 1
% 'OMITALL', 0, 'BLANK_TOUCH', 0, 'SHK', 0, 'BLOCK', 2, 'BLOCK', 3

% actually, load overall workspace (e.g. workspace with 5 variables
% identified, then identify the one that you want to see change, (e.g.,
% collectionTime Blocks 2 & 3), then change the variable below to do that
% comparison with the appropriate array #s

% use this to specify which arrays should be compared below. this is
% critical - you must know which arrays correspond to variables to be
% compared, i..e, prechoice B1 vs prechoice B2/3, which may correspond to
% respClass_all_arrays 1 and 4, in the case of arrays of size 6. this code
% will break and could make invalid comparisons if you do not adjust these values appropriately 
if size(respClass_all_array, 2) == 10 | size(respClass_all_array, 2) == 11 | size(respClass_all_array, 2) == 12
    comparison_arrays = [1 2 3; 8 9 10]
elseif size(respClass_all_array, 2) == 6
    comparison_arrays = [1 2 3; 4 5 6]
elseif size(respClass_all_array, 2) == 7
    comparison_arrays = [1 2 3; 5 6 7]
elseif size(respClass_all_array, 2) == 3
    comparison_arrays = [1, 2, 3; 1, 2, 3]
end



arrays_to_examine = [1 8];

inhib_or_excite = 1;

event_for_figures = 1; 




%%
% Initialize output cell arrays
trial_peak_sum_all = cell(size(trial_peaks_mouse_array));
trial_peak_num_all = cell(size(trial_peaks_mouse_array));

% Loop over rows and columns of trial_peaks_mouse_array
for row = 1:size(trial_peaks_mouse_array, 1)
    for col = 1:size(trial_peaks_mouse_array, 2)
        
        % Extract the current trial peaks cell array (1-row, many columns)
        current_trial_peaks_array = trial_peaks_mouse_array{row, col};
        
        % Ensure it's not empty before processing
        if ~isempty(current_trial_peaks_array)
            % Initialize arrays to store results for each column
            trial_peak_sum_temp = zeros(1, size(current_trial_peaks_array, 2));
            trial_peak_num_temp = zeros(1, size(current_trial_peaks_array, 2));

            % Loop over each column in the current cell array
            for ff = 1:size(current_trial_peaks_array, 2)
                
                % Extract the current column data (double array)
                current_data = current_trial_peaks_array{1, ff};
                
                % Identify peaks (values >= 1)
                current_data_ind = current_data >= 1;

                % Compute sum and peak count, store as double
                trial_peak_sum_temp(1, ff) = sum(current_data);
                trial_peak_num_temp(1, ff) = sum(current_data_ind);
            end
            
            % Store results in the main cell array (as doubles)
            trial_peak_sum_all{row, col} = trial_peak_sum_temp;
            trial_peak_num_all{row, col} = trial_peak_num_temp;
        else
            % Store empty arrays if input is empty
            trial_peak_sum_all{row, col} = [];
            trial_peak_num_all{row, col} = [];
        end
    end
end




%%
for kk = 1:size(animalIDs, 1)
    % prechoice_block_1 = respClass_all_array{1, comparison_arrays(1, 1)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
    prechoice_block_1_mouse{kk, :} = trial_peak_num_all{kk, comparison_arrays(1, 1)} == inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(1, 2)} ~= inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(1, 3)} ~= inhib_or_excite;
    prechoice_blocks_2_and_3_mouse{kk, :} = trial_peak_num_all{kk, comparison_arrays(2, 1)} == inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(2, 2)} ~= inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(2, 3)} ~= inhib_or_excite;

    % postchoice_reward_block_1 = respClass_all_array{1, comparison_arrays(1, 2)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 3)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
    postchoice_reward_block_1_mouse{kk, :}  = trial_peak_num_all{kk, comparison_arrays(1, 2)} == inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(1, 1)} ~= inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(1, 3)} ~= inhib_or_excite;
    postchoice_reward_blocks_2_and_3_mouse{kk, :}  = trial_peak_num_all{kk, comparison_arrays(2, 2)} == inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(2, 1)} ~= inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(2, 3)} ~= inhib_or_excite;

    % collect_block_1 = respClass_all_array{1, comparison_arrays(1, 3)} == inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 1)} ~= inhib_or_excite & respClass_all_array{1, comparison_arrays(1, 2)} ~= inhib_or_excite & respClass_all_array{1, 4} ~= inhib_or_excite;
    collect_block_1_mouse{kk, :}  = trial_peak_num_all{kk, comparison_arrays(1, 3)} == inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(1, 1)} ~= inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(1, 2)} ~= inhib_or_excite;
    collect_blocks_2_and_3_mouse{kk, :}  = trial_peak_num_all{kk, comparison_arrays(2, 3)} == inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(2, 1)} ~= inhib_or_excite & trial_peak_num_all{kk, comparison_arrays(2, 2)} ~= inhib_or_excite;

    num_cells_mouse(kk) = size(trial_peak_num_all{kk, 1}, 2);
    true_neutral_block_1_mouse{kk, :} = trial_peak_num_all{kk, comparison_arrays(1, 1)} == 3 & trial_peak_num_all{kk, comparison_arrays(1, 2)} == 3 & trial_peak_num_all{kk, comparison_arrays(1, 3)} == 3;
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
end