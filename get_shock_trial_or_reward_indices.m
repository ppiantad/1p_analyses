load('BLA_RDT_11_variables_incl_AA_use_for_AA_PV.mat')

%%

for ii = 1:size(animalIDs, 1)
    currentanimal = char(animalIDs(ii));
    if isfield(final_behavior.(currentanimal), session_to_analyze)
        % Extract the table for easy reference
        data = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData  ;
        [BehavData,trials,varargin_identity_class]=TrialFilter_test(data, 'OMITALL', 0, 'BLANK_TOUCH', 0);
        shock_trials{ii} = find(BehavData.shock == 1)+1; %+1 because trials are indexed in ABET starting from 2
        num_shocks(ii) = size(shock_trials{ii}, 1);
    end
end

shock_matrix = padcat(shock_trials{:}); % from File Exchange

% Transpose so that each *column* corresponds to one animal
shock_matrix = shock_matrix; 

% Convert numeric matrix to a table
shock_table = array2table(shock_matrix);

% Add headers based on animal IDs
shock_table.Properties.VariableNames = matlab.lang.makeValidName(cellstr(animalIDs));

% Save as CSV
% output_filename = fullfile(pwd, 'shock_trials.csv');
% writetable(shock_table, 'shock_trials.csv');

%%

% MATLAB R2017a
LB = 2;   % starting @ 2 because trials are indexed in ABET starting at 2 (because trial counter increments @ choice)
UB = 91;  % upper bound of range (integer),  UB > LB
total_trials = 2:91;

for ii = 1:size(animalIDs, 1)
    currentanimal = char(animalIDs(ii));
    if isfield(final_behavior.(currentanimal), session_to_analyze)
        % Extract the table for easy reference
        data = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData  ;
        [BehavData,trials,varargin_identity_class]=TrialFilter_test(data, 'OMITALL', 0, 'BLANK_TOUCH', 0);
        large_rew_trials(ii) = sum(BehavData.bigSmall == 1.2);
        small_rew_trials(ii) = sum(BehavData.bigSmall == 0.3);
        trials_for_large_delivery{ii} = sort(randsample(LB:UB,large_rew_trials(ii),'false'));
        trials_for_small_delivery{ii} = setdiff(total_trials, trials_for_large_delivery{ii});
    end
end

large_trials_matrix = padcat(trials_for_large_delivery{:}); % from File Exchange

% Transpose so that each *column* corresponds to one animal
large_trials_matrix = large_trials_matrix'; 

% Convert numeric matrix to a table
large_trials_table = array2table(large_trials_matrix);

% Add headers based on animal IDs
large_trials_table.Properties.VariableNames = matlab.lang.makeValidName(cellstr(animalIDs));

% Save as CSV
% output_filename = fullfile(pwd, 'shock_trials.csv');
writetable(large_trials_table, 'large_trials.csv');


small_trials_matrix = padcat(trials_for_small_delivery{:}); % from File Exchange

% Transpose so that each *column* corresponds to one animal
small_trials_matrix = small_trials_matrix'; 

% Convert numeric matrix to a table
small_trials_table = array2table(small_trials_matrix);

% Add headers based on animal IDs
small_trials_table.Properties.VariableNames = matlab.lang.makeValidName(cellstr(animalIDs));

% Save as CSV
% output_filename = fullfile(pwd, 'shock_trials.csv');
writetable(small_trials_table, 'small_trials.csv');
