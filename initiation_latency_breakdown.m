initiation_to_choice_short = 0;
initiation_to_choice_long = 0;
total_trials = 0;
% Initialize arrays to store per-animal data
initiation_to_choice_short_per_animal = zeros(size(valid_animalIDs, 1), 1);
initiation_to_choice_long_per_animal = zeros(size(valid_animalIDs, 1), 1);
animal_names = cell(size(valid_animalIDs, 1), 1);

for ii = 1:size(valid_animalIDs,1) % 1:size(fieldnames(final),1)
    currentanimal = char(valid_animalIDs(ii));
    animal_names{ii} = currentanimal;  % Store animal name for reference
    
    % Initialize counters for this animal
    animal_initiation_to_choice_short = 0;
    animal_initiation_to_choice_long = 0;
    
    if isfield(final_behavior.(currentanimal), session_to_analyze)
        % if contains(session_to_analyze, 'CNO')
        %     BehavData = final_behavior.(currentanimal).(session_to_analyze).BehavData;
        % else
        %     BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        % end
        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        BehavData=TrialFilter(BehavData,'OMITALL', 0, 'BLANK_TOUCH', 0);
        for hh = 1:size(BehavData, 1)
            % Check if current row has type_binary of 1 or 2
            if (BehavData.choiceTime(hh) - BehavData.stTime(hh)) <= 4
                animal_initiation_to_choice_short = animal_initiation_to_choice_short + 1;
                initiation_to_choice_short = initiation_to_choice_short + 1;
                total_trials = total_trials + 1;
            elseif (BehavData.choiceTime(hh) - BehavData.stTime(hh)) >= 4
                animal_initiation_to_choice_long = animal_initiation_to_choice_long + 1;
                initiation_to_choice_long = initiation_to_choice_long + 1;
                total_trials = total_trials + 1;
            end
        end
    end
    
    % Store this animal's data
    initiation_to_choice_short_per_animal(ii) = animal_initiation_to_choice_short;
    initiation_to_choice_long_per_animal(ii) = animal_initiation_to_choice_long;
    
    % Print individual animal results (optional)
    fprintf('Animal %s - Short: %d, Long: %d\n', currentanimal, animal_initiation_to_choice_short, animal_initiation_to_choice_long);
end
percent_short_trials = (initiation_to_choice_short/total_trials)*100;
percent_long_trials = (initiation_to_choice_long/total_trials)*100;
% Print summary across all animals
fprintf('\n--- Summary Across All Animals ---\n');
fprintf('Total trials where initiation to choice < 4 sec: %d\n', initiation_to_choice_short);
fprintf('Total trials where initiation to choice > 4 sec: %d\n', initiation_to_choice_long);
fprintf('Percent of trials where initiation to choice < 4 sec: %g\n', percent_short_trials);
fprintf('Percent of trials where initiation to choice > 4 sec: %g\n', percent_long_trials);


