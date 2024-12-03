% Load data
data_path = 'final_DLC.B46837.D3.DLC_data_raw'; % Replace with the actual variable assignment as a string
test_data = eval(data_path);

% Define the columns that need to be updated
columns_to_update = {'left_ear_x', 'left_ear_y', 'right_ear_x', ...
                     'right_ear_y', 'tail_x', 'tail_y', 'mean_x_pix', 'mean_y_pix'};

% Loop through each column for jump identification and correction
for col_idx = 1:numel(columns_to_update)
    col_name = columns_to_update{col_idx};
    disp(['Processing column: ', col_name]);

    % Extract the data for the current column
    fixedX = test_data.(col_name);

    % Identify jumps in the current column
    jumps = abs([0; diff(fixedX)]); % Include a leading zero for alignment
    jump_ind = find(jumps > 30); % Adjust threshold as needed

    % Initialize jump correction loop for the current column
    idx = 1; % Start with the first jump index
    while idx <= length(jump_ind)
        current_idx = jump_ind(idx);

        % Store the value prior to the first jump for reference
        if idx == 1
            value_prior_to_first_jump = fixedX(max(1, current_idx - 1));
        end

        % Inner loop to handle decisions for the current jump
        while true
            % Extract values just prior to, at, and just after the current index
            value_prior = fixedX(max(1, current_idx - 1));
            value_jump = fixedX(current_idx);
            value_after = fixedX(min(length(fixedX), current_idx + 1));

            % Display information and prompt for action
            disp(['Column: ', col_name, ', Index: ', num2str(current_idx)]);
            disp(['Value prior to first jump: ', num2str(value_prior_to_first_jump)]);
            disp(['Value prior to current index: ', num2str(value_prior)]);
            disp(['Value at current index: ', num2str(value_jump)]);
            disp(['Value after current index: ', num2str(value_after)]);

            disp('Options:');
            disp('1: Replace with value prior to the first jump');
            disp('2: Replace with value prior to the current index');
            disp('3: Replace with the current value');
            disp('4: Replace with value after the current index');
            disp('5: Move to the next index of this jump');
            disp('6: Skip to the next jump');
            choice = input('Select an option: ');

            % Take action based on user choice
            if choice == 1
                fixedX(current_idx) = value_prior_to_first_jump;
                disp('Replaced with value prior to the first jump.');
                current_idx = current_idx + 1;
            elseif choice == 2
                fixedX(current_idx) = value_prior;
                disp('Replaced with value prior to the current index.');
                current_idx = current_idx + 1;
            elseif choice == 3
                fixedX(current_idx) = value_jump;
                disp('Replaced with the current value.');
                current_idx = current_idx + 1;
            elseif choice == 4
                fixedX(current_idx) = value_after;
                disp('Replaced with value after the current index.');
                current_idx = current_idx + 1;
            elseif choice == 5
                disp('Moving to the next index of this jump...');
                current_idx = current_idx + 1;
            elseif choice == 6
                disp('Skipping to the next jump...');
                idx = idx + 1;
                if idx > length(jump_ind)
                    disp('No more jumps to process in this column.');
                    break; % Exit the inner loop
                else
                    current_idx = jump_ind(idx); % Update to the next jump index
                end
            else
                disp('Invalid choice. Please try again.');
            end

            % Check if we have reached the end of the data for this jump
            if current_idx > length(fixedX)
                disp('Reached the end of the data for this jump.');
                break; % Exit the inner loop
            end
        end

        if idx > length(jump_ind)
            break; % Exit the outer loop for the column
        end
    end

    % Save the corrected data back to the table for this column
    test_data.(col_name) = fixedX;

    % Display progress
    disp(['Finished processing column: ', col_name]);
end

% Convert the variable name string dynamically
corrected_path = strrep(data_path, 'DLC_data_raw', 'DLC_raw_data_corrected');

% Use `eval` to save the corrected data
eval([corrected_path, ' = test_data;']);

% Optional: Confirm that the data has been saved
disp(['Corrected data saved to ', corrected_path]);
