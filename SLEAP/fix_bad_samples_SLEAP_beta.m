%%




fixedX = SLEAP_data.x_pix'; 
plot(fixedX)
pks = findpeaks(fixedX)

for qq = 1:size(fixedX, 1)
    if qq > 1
        jumps(qq) = abs(fixedX(qq-1) - fixedX(qq));
    end
end


jump_ind = find(jumps > 70);

% Loop through each jump index
for idx = 1:length(jump_ind)
    current_idx = jump_ind(idx);
    
    % Extract values just prior to, at, and just after the jump index
    value_prior = fixedX(max(1, current_idx - 1));
    value_jump = fixedX(current_idx);
    value_after = fixedX(min(size(fixedX, 1), current_idx + 1));
    
    % Display options to the user and prompt for input
    disp(['Index: ', num2str(current_idx)]);
    disp(['Value just prior to the jump: ', num2str(value_prior)]);
    disp(['Value at the jump index: ', num2str(value_jump)]);
    disp(['Value just after the jump: ', num2str(value_after)]);
    choice = input('Which value to replace the jump with? (1: prior, 2: jump, 3: after): ');
    
    % Replace the jump with the chosen value
    if choice == 1
        fixedX(current_idx) = value_prior;
    elseif choice == 2
        fixedX(current_idx) = value_jump;
    elseif choice == 3
        fixedX(current_idx) = value_after;
    else
        disp('Invalid choice. Skipping...');
    end
end

% Display or use fixedX as needed
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