%%


test_data = final_DLC.C68604.D3.movement_data  ;

fixedX = test_data.mean_x_pix; 
plot(fixedX)
pks = findpeaks(fixedX)

for qq = 1:size(fixedX, 1)
    if qq > 1
        jumps(qq) = abs(fixedX(qq-1) - fixedX(qq));
    end
end


jump_ind = find(jumps > 30);

% Loop through each jump index
idx = 1; % Start with the first jump index
while idx <= length(jump_ind)
    current_idx = jump_ind(idx);
    
    % Store the value prior to the first jump for reference
    if idx == 1
        value_prior_to_first_jump = fixedX(max(1, current_idx - 1));
    end
    
    % Inner loop to handle decisions for the current jump and subsequent indices
    while true
        % Extract values just prior to, at, and just after the current index
        value_prior = fixedX(max(1, current_idx - 1));
        value_jump = fixedX(current_idx);
        value_after = fixedX(min(size(fixedX, 1), current_idx + 1));
        
        % Display information and prompt for action
        disp(['Index: ', num2str(current_idx)]);
        disp(['Value just prior to the first jump: ', num2str(value_prior_to_first_jump)]);
        disp(['Value just prior to the current index: ', num2str(value_prior)]);
        disp(['Value at the current index: ', num2str(value_jump)]);
        disp(['Value just after the current index: ', num2str(value_after)]);
        
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
            current_idx = current_idx + 1; % Move to the next index automatically
        elseif choice == 2
            fixedX(current_idx) = value_prior;
            disp('Replaced with value prior to the current index.');
            current_idx = current_idx + 1; % Move to the next index automatically
        elseif choice == 3
            fixedX(current_idx) = value_jump;
            disp('Replaced with the current value.');
            current_idx = current_idx + 1; % Move to the next index automatically
        elseif choice == 4
            fixedX(current_idx) = value_after;
            disp('Replaced with value after the current index.');
            current_idx = current_idx + 1; % Move to the next index automatically
        elseif choice == 5
            disp('Moving to the next index of this jump...');
            current_idx = current_idx + 1;
        elseif choice == 6
            disp('Skipping to the next jump...');
            idx = idx + 1; % Move to the next jump index
            if idx > length(jump_ind)
                disp('No more jumps to process.');
                break; % Exit the outer loop
            else
                current_idx = jump_ind(idx); % Update to the new jump index
            end
        else
            disp('Invalid choice. Please try again.');
        end
        
        % Check if we have reached the end of the current jump sequence
        if current_idx > size(fixedX, 1)
            disp('Reached the end of the data for this jump.');
            break; % Exit the inner loop
        end
    end
    
    if idx > length(jump_ind)
        break; % Exit the outer loop if all jumps are processed
    end
end

% Display or use fixedX as needed
disp('Final corrected values:');
disp(fixedX);



%%

% this code somewhat works, though can easily run for infinite time if it
% detects a large jump in X/Y coordinates that that isn't "fixed" by
% averaging across the prior 5 frames

% Initialize variables
jump_threshold = 10; % Set your threshold for defining a jump
% fixedX = correctedX; % Initialize a copy of correctedX for modifications
fixedX = SLEAP_data.x_pix; 



% Flag to track if jumps are still being fixed
jumps_fixed = true;

while jumps_fixed
    jumps_fixed = false; % Reset jumps_fixed flag for each iteration
    
    % Initialize start index for the loop
    start_index = 2;
    
    % Identify jumps and fix them
    for qq = start_index:size(fixedX, 1)
        jump = abs(fixedX(qq-1) - fixedX(qq));
        if jump > jump_threshold
            % Fix the jump by replacing it with the average of the prior 5 samples
            fixedX(qq) = mean(fixedX(max(1, qq-5):qq-1));
            jumps_fixed = true; % Set jumps_fixed flag to indicate a jump was fixed
            start_index = qq; % Update start index for the next iteration
        end
    end
end

% Display or use fixedX as needed
disp(fixedX);


SLEAP_data.x_pix = fixedX;


%%

% Initialize variables
jump_threshold = 10; % Set your threshold for defining a jump
% fixedX = correctedX; % Initialize a copy of correctedX for modifications
fixedY = SLEAP_data.y_pix; 



% Flag to track if jumps are still being fixed
jumps_fixed = true;

while jumps_fixed
    jumps_fixed = false; % Reset jumps_fixed flag for each iteration
    
    % Initialize start index for the loop
    start_index = 2;
    
    % Identify jumps and fix them
    for qq = start_index:size(fixedX, 1)
        jump = abs(fixedY(qq-1) - fixedY(qq));
        if jump > jump_threshold
            % Fix the jump by replacing it with the average of the prior 5 samples
            fixedY(qq) = mean(fixedY(max(1, qq-5):qq-1));
            jumps_fixed = true; % Set jumps_fixed flag to indicate a jump was fixed
            start_index = qq; % Update start index for the next iteration
        end
    end
end

% Display or use fixedX as needed
disp(fixedY);

SLEAP_data.y_pix = fixedY;


%%
for bb = 1:size(fixedX, 1)
    if bb > 1
        jumps(bb) = abs(fixedX(bb-1) - fixedX(bb));
    end
end

ind = find(jumps > 10);