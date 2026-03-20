% run generage_behav_figs_RDT.m first 
%%

% Define block structure
trials_per_block = 30;
forced_trials = 8;
free_trials = 22;
free_choice_to_analyze = 10;

% Initialize results array (10 mice x 3 blocks)
percent_large_choice = zeros(10, 3);

% Calculate for each mouse and block
for mouse = 1:size(large_sequences_mouse, 1)
    for block = 1:3
        % Find the start of this block
        block_start = (block - 1) * trials_per_block + 1;
        
        % Get the first 10 free choice trials (skip first 8 forced trials)
        free_choice_start = block_start + forced_trials;
        free_choice_end = free_choice_start + free_choice_to_analyze - 1;
        
        % Extract the relevant trials
        trials = large_sequences_mouse(mouse, free_choice_start:free_choice_end);
        
        % Calculate percentage (1s are large choices)
        percent_large_choice(mouse, block) = (sum(trials) / free_choice_to_analyze) * 100;
    end
end