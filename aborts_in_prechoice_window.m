aborts_outside_prechoice = 0;
aborts_inside_prechoice = 0;

% Initialize arrays to store per-animal data
aborts_outside_per_animal = zeros(size(valid_animalIDs, 1), 1);
aborts_inside_per_animal = zeros(size(valid_animalIDs, 1), 1);
animal_names = cell(size(valid_animalIDs, 1), 1);

for ii = 1:size(valid_animalIDs,1) % 1:size(fieldnames(final),1)
    currentanimal = char(valid_animalIDs(ii));
    animal_names{ii} = currentanimal;  % Store animal name for reference
    
    % Initialize counters for this animal
    animal_aborts_outside = 0;
    animal_aborts_inside = 0;
    
    if isfield(final_behavior.(currentanimal), session_to_analyze)
        % if contains(session_to_analyze, 'CNO')
        %     BehavData = final_behavior.(currentanimal).(session_to_analyze).BehavData;
        % else
        %     BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        % end
        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        
        for hh = 1:size(BehavData, 1)
            % Check if current row has type_binary of 1 or 2
            if BehavData.type_binary(hh) == 1 || BehavData.type_binary(hh) == 2
                % Look at all subsequent rows
                for jj = (hh+1):size(BehavData, 1)
                    % Check if this subsequent row has bigSmall of 1.2 or 0.3
                    if BehavData.bigSmall(jj) == 1.2 || BehavData.bigSmall(jj) == 0.3
                        % Calculate time difference
                        abort_to_choice_latency = BehavData.choiceTime(jj) - BehavData.choiceTime(hh);
                        
                        % Categorize based on latency
                        if abort_to_choice_latency > 4
                            aborts_outside_prechoice = aborts_outside_prechoice + 1;
                            animal_aborts_outside = animal_aborts_outside + 1;
                        else
                            aborts_inside_prechoice = aborts_inside_prechoice + 1;
                            animal_aborts_inside = animal_aborts_inside + 1;
                        end
                        
                        % Found a match, so break out of inner loop to avoid double-counting
                        break;
                    end
                end
            end
        end
    end
    
    % Store this animal's data
    aborts_outside_per_animal(ii) = animal_aborts_outside;
    aborts_inside_per_animal(ii) = animal_aborts_inside;
    
    % Print individual animal results (optional)
    fprintf('Animal %s - Outside: %d, Inside: %d\n', currentanimal, animal_aborts_outside, animal_aborts_inside);
end

% Print summary across all animals
fprintf('\n--- Summary Across All Animals ---\n');
fprintf('Total aborts outside prechoice window (>4s): %d\n', aborts_outside_prechoice);
fprintf('Total aborts inside prechoice window (≤4s): %d\n', aborts_inside_prechoice);

% Optional: Create a summary table for easier viewing
summary_table = table(animal_names, aborts_outside_per_animal, aborts_inside_per_animal, ...
    'VariableNames', {'AnimalID', 'AbortsOutside', 'AbortsInside'});
disp(summary_table);

% Optional: Calculate and display percentages
total_per_animal = aborts_outside_per_animal + aborts_inside_per_animal;
percent_inside = (aborts_inside_per_animal ./ total_per_animal) * 100;
percent_inside(isnan(percent_inside)) = 0;  % Handle division by zero

summary_table_with_percent = table(animal_names, aborts_outside_per_animal, aborts_inside_per_animal, ...
    total_per_animal, percent_inside, ...
    'VariableNames', {'AnimalID', 'AbortsOutside', 'AbortsInside', 'Total', 'PercentInside'});
disp(summary_table_with_percent);