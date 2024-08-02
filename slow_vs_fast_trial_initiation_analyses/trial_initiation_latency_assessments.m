
load('BLA_C_raw_no_additional_filtering_Pre_RDT_RM_only_completed_sessions_zall_window_base_workspace_3_categories.mat')

% arrays_to_examine = [1 8];

inhib_or_excite = 1;

event_for_figures = 1; 

%%

prechoice = respClass_all_array{1, 1} == inhib_or_excite & respClass_all_array{1, 2} ~= inhib_or_excite & respClass_all_array{1, 3} ~= inhib_or_excite;

postchoice = respClass_all_array{1, 2} == inhib_or_excite & respClass_all_array{1, 1} ~= inhib_or_excite & respClass_all_array{1, 3} ~= inhib_or_excite;

collect = respClass_all_array{1, 3} == inhib_or_excite & respClass_all_array{1, 1} ~= inhib_or_excite & respClass_all_array{1, 2} ~= inhib_or_excite;






%%

relevant_behav_array = behav_tbl_iter{1, 1};


for zz = 1:size(relevant_behav_array, 1)
    BehavData = relevant_behav_array{zz};
    
    % Create a logi0.10cal index array based on your conditions
    logical_index{zz} = BehavData.stTime - BehavData.TrialPossible >= 10 & BehavData.stTime - BehavData.TrialPossible <= 50;

    % Use the logical index array to subset BehavData
    % BehavData = BehavData(logical_index,: );
    % trials = trials(logical_index);

end

%%
only_prechoice_neurons = zall_array(1, prechoice == 1);

only_choice_aligned_array =  zall_mouse(:, 1);

mean_data_from_cells = []; 

count = 1; 
for pp = 1:size(only_choice_aligned_array, 1)
    mouse_data = only_choice_aligned_array{pp, 1};
    mouse_data = mouse_data(respClass_all_array_mouse{pp, 3} == 1);
    for gg = 1:size(mouse_data, 2)
        cell_data = mouse_data{gg};
        cell_data = cell_data(logical_index{1, pp}, :);
        mean_data_from_cells(count, :) = mean(cell_data);
        count = count+1;
    end
    
end


figure; plot(ts1, mean(mean_data_from_cells));