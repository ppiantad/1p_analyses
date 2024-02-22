% this code somewhat works, though can easily run for infinite time if it
% detects a large jump in X/Y coordinates that that isn't "fixed" by
% averaging across the prior 5 frames

% Initialize variables
jump_threshold = 20; % Set your threshold for defining a jump
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
jump_threshold = 20; % Set your threshold for defining a jump
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