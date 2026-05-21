%% Export D3 behavioral and calcium data to Python-readable .mat files
% Converts MATLAB table objects in context_fear_PFC_data.mat to plain arrays.

load('I:\MATLAB\my_repo\context fear\aggregate_data\context_fear_PFC_data.mat');
load('I:\MATLAB\my_repo\context fear\aggregate_data\PFC alone imaging\calcium\PFC_alone_imaging_data_07182025.mat');

imaging_animals = fieldnames(final);
behav_animals   = fieldnames(final_DLC);

% Animals that have both imaging and behavior
common_animals = intersect(imaging_animals, behav_animals);

export = struct();

for ii = 1:length(common_animals)
    aname = common_animals{ii};

    % Skip if no D3 imaging
    if ~isfield(final.(aname), 'D3'), continue; end
    if ~isfield(final_DLC.(aname), 'D3'), continue; end

    % --- Calcium data ---
    C      = final.(aname).D3.CNMFe_data.C;       % neurons x frames
    C_raw  = final.(aname).D3.CNMFe_data.C_raw;
    S      = full(final.(aname).D3.CNMFe_data.S);
    sr     = final.(aname).D3.uv.sampling_rate;

    % --- Behavioral data (MATLAB table) ---
    T = final_DLC.(aname).D3.DLC_data_raw;

    % Print column names for first animal
    if ii == 1
        disp('DLC table variables:');
        disp(T.Properties.VariableNames);
    end

    % Position
    mean_x = T.mean_x_pix;
    mean_y = T.mean_y_pix;

    % Freezing — look for a freezing or immobility column
    var_names_lower = lower(T.Properties.VariableNames);
    freeze_col = find(contains(var_names_lower, 'freez') | contains(var_names_lower, 'immob') | contains(var_names_lower, 'freeze'));
    if ~isempty(freeze_col)
        freezing = T{:, freeze_col(1)};
    else
        freezing = nan(height(T), 1);
        disp(['WARNING: no freezing column found for ' aname]);
        disp('Available columns:');
        disp(T.Properties.VariableNames);
    end

    % Velocity / motion column
    vel_col = find(contains(var_names_lower, 'speed') | contains(var_names_lower, 'veloc') | contains(var_names_lower, 'locom') | contains(var_names_lower, 'motion'));
    if ~isempty(vel_col)
        velocity = T{:, vel_col(1)};
    else
        velocity = nan(height(T), 1);
    end

    % Store in export struct
    export.(aname).C          = C;
    export.(aname).C_raw      = C_raw;
    export.(aname).S          = S;
    export.(aname).sampling_rate = sr;
    export.(aname).mean_x     = mean_x;
    export.(aname).mean_y     = mean_y;
    export.(aname).freezing   = freezing;
    export.(aname).velocity   = velocity;
    export.(aname).n_frames_behav = height(T);
    export.(aname).experimental_grp = final.(aname).experimental_grp;

    % All column names and full table as numeric where possible
    export.(aname).dlc_varnames = strjoin(T.Properties.VariableNames, ',');

    fprintf('Exported %s: %d neurons x %d imaging frames, %d behav frames\n', ...
        aname, size(C,1), size(C,2), height(T));
end

save_path = 'I:\MATLAB\my_repo\context fear\aggregate_data\D3_export_for_python.mat';
save(save_path, 'export', '-v7');
fprintf('\nSaved to %s\n', save_path);
