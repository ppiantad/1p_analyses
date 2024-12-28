valid_session_names = {'D1_Morning', 'D1_Afternoon', 'D2_Morning', 'D2_Afternoon', 'D3', 'D4'};

% Example inputs
struct_data = final_DLC.B46857  ; % Replace with your actual struct
% final_DLC.B57422: experimental, ok to plot
% final_DLC.B46852: No Shock, good to plot
% final_DLC.B46857: One Context, good to plot

% Get field names of the top-level struct
fields = fieldnames(struct_data);

% Filter the fields based on valid_session_names
matching_fields = fields(ismember(fields, valid_session_names));

% Initialize variables to store global min and max
global_x_min = inf;
global_x_max = -inf;
global_y_min = inf;
global_y_max = -inf;

% First pass: Calculate global min and max
for i = 1:numel(valid_session_names)
    session_name = valid_session_names{i};
    
    if isfield(struct_data, session_name) && isfield(struct_data.(session_name), 'movement_data')
        movement_data = struct_data.(session_name).movement_data;
        
        % Extract X and Y positions, omitting the first 100 rows
        x_positions = movement_data.mean_x_pix(101:end);
        y_positions = movement_data.mean_y_pix(101:end);
        
        % Update global min and max
        global_x_min = min(global_x_min, min(x_positions));
        global_x_max = max(global_x_max, max(x_positions));
        global_y_min = min(global_y_min, min(y_positions));
        global_y_max = max(global_y_max, max(y_positions));
    end
end

% Prepare figure
num_sessions = numel(matching_fields);
figure;
tiledlayout(1, num_sessions); % Use tiledlayout for 1-row subplots

% Second pass: Plot data with consistent axes
for i = 1:num_sessions
    session_name = valid_session_names{i};
    
    % Check if the session contains the expected field
    if isfield(struct_data, session_name) && isfield(struct_data.(session_name), 'movement_data')
        movement_data = struct_data.(session_name).movement_data;
        
        % Extract X and Y positions, omitting the first 100 rows
        x_positions = movement_data.body_x_pix(101:end);
        y_positions = movement_data.body_y_pix(101:end);
        
        % Create a subplot
        nexttile;
        plot(x_positions, y_positions, 'LineWidth', 1.5);
        title(sprintf('Session: %s', session_name), 'Interpreter', 'none');
        xlabel('X Position (pixels)');
        ylabel('Y Position (pixels)');
        grid on;
        
        % Set consistent axis limits
        xlim([global_x_min, global_x_max]);
        ylim([global_y_min, global_y_max]);
    else
        warning('No "movement_data" field in session: %s', session_name);
    end
end

% Adjust the layout
sgtitle('Mouse Movement Across Sessions');

