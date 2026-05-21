% locomotor_trajectories_session_normalized_GoPro.m
%
% Plots XY trajectories with velocity colormapped dots for randomly selected
% large and small reward trials (Block 1, 2, 3 — one subplot each).
%
% Key difference from v2: velocity is normalized to [0,1] across the ENTIRE
% SESSION before trial extraction, so the colormap is consistent across all
% trials and all subplots. A single shared colorbar applies to all 3 subplots.
%
% Requires in workspace: final_SLEAP, final_behavior

%% ---- Settings ----
select_mouse       = 'RDT_M_1';
session_to_analyze = 'RDT_D1';
lines_to_plot      = 1;       % random trials per block (capped by available)
cmap_name          = 'jet';   % colormap applied to all subplots
vel_clim           = [0, 0.6]; % colormap range in session-normalized velocity units

%% ---- Load data ----
shapeData  = final_SLEAP.(select_mouse).(session_to_analyze).shapeData;
SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;
BehavData  = final_behavior.(select_mouse).(session_to_analyze).uv.BehavData;

onset_trials  = BehavData.stTime';
choice_trials = BehavData.choiceTime';
offset_trials = BehavData.collectionTime';
fs_cam = 30;
time_ranges_trials = [onset_trials; choice_trials; offset_trials];

% Smoothed XY
X_data = sgolayfilt(SLEAP_data.x_pix, 9, 33);
Y_data = sgolayfilt(SLEAP_data.y_pix, 9, 33);

%% ---- Session-level velocity normalization (done ONCE, before trial extraction) ----
vel_raw = SLEAP_data.vel_cm_s;
vel_min = min(vel_raw);
vel_max = max(vel_raw);
vel_norm = (vel_raw - vel_min) / max(vel_max - vel_min, eps);  % [0, 1] across session

%% ---- Extract per-trial XY, velocity, and choice point ----
filtered_motion   = {};
filtered_velocity = {};   % values already in [0,1] — no further scaling needed
choice_times      = {};

max_ind    = SLEAP_data.idx_frame(end);
good_index = 1;

for j = 1:size(time_ranges_trials, 2)
    onset  = round(time_ranges_trials(1,j) * fs_cam) + 1;
    offset = round(time_ranges_trials(3,j) * fs_cam) + 1;

    if isinf(offset)
        if onset <= max_ind && onset > 0
            filtered_motion{good_index} = SLEAP_data(onset:end); %#ok<SAGROW>
            break
        end
    else
        mask = SLEAP_data.idx_time > time_ranges_trials(1,j) & ...
               SLEAP_data.idx_time < time_ranges_trials(3,j);

        filtered_motion{j}   = [X_data(mask)'; Y_data(mask)'];   % 2 × k
        filtered_velocity{j} = vel_norm(mask)';                   % 1 × k, [0,1]

        % Choice-point XY (nearest frame to choiceTime)
        ci = interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), ...
                     time_ranges_trials(2,j), 'nearest');
        if ~isnan(ci)
            choice_times{j} = [X_data(ci); Y_data(ci)];
        else
            choice_times{j} = [nan; nan];
        end

        good_index = good_index + 1;
    end
end

%% ---- Trial indices by reward × block ----
large_b1 = find(BehavData.bigSmall == 1.2 & BehavData.Block == 1);
large_b2 = find(BehavData.bigSmall == 1.2 & BehavData.Block == 2);
large_b3 = find(BehavData.bigSmall == 1.2 & BehavData.Block == 3);
small_b1 = find(BehavData.bigSmall == 0.3 & BehavData.Block == 1);
small_b2 = find(BehavData.bigSmall == 0.3 & BehavData.Block == 2);
small_b3 = find(BehavData.bigSmall == 0.3 & BehavData.Block == 3);

%% ---- Helper: draw arena shapes ----
% (circles and left/right squares from shapeData)

%% ---- Figure: LARGE reward trials ----
n_large = min(lines_to_plot, min([numel(large_b1), numel(large_b2), numel(large_b3)]));
sel_L = { large_b1(randperm(numel(large_b1), n_large)), ...
          large_b2(randperm(numel(large_b2), n_large)), ...
          large_b3(randperm(numel(large_b3), n_large)) };
titles_L = {'Large Reward — Block 1', 'Large Reward — Block 2', 'Large Reward — Block 3'};

plot_trajectory_figure(sel_L, titles_L, filtered_motion, filtered_velocity, ...
    choice_times, shapeData, cmap_name, vel_clim, ...
    sprintf('%s  %s  — Large Reward', strrep(select_mouse,'_','\_'), ...
            strrep(session_to_analyze,'_','\_')));

%% ---- Figure: SMALL reward trials ----
n_small = min(lines_to_plot, min([numel(small_b1), numel(small_b2), numel(small_b3)]));
sel_S = { small_b1(randperm(numel(small_b1), n_small)), ...
          small_b2(randperm(numel(small_b2), n_small)), ...
          small_b3(randperm(numel(small_b3), n_small)) };
titles_S = {'Small Reward — Block 1', 'Small Reward — Block 2', 'Small Reward — Block 3'};

plot_trajectory_figure(sel_S, titles_S, filtered_motion, filtered_velocity, ...
    choice_times, shapeData, cmap_name, vel_clim, ...
    sprintf('%s  %s  — Small Reward', strrep(select_mouse,'_','\_'), ...
            strrep(session_to_analyze,'_','\_')));


%% ================================================================
%  Local functions
%  ================================================================

function plot_trajectory_figure(sel_by_block, subplot_titles, ...
        filtered_motion, filtered_velocity, choice_times, shapeData, ...
        cmap_name, vel_clim, fig_title)
% Creates a 3-subplot figure (one per block) with velocity-colormapped
% trajectories. clim is fixed at [0 1] on every subplot so all three share
% the same colormap scale. A single colorbar is added at the figure level.

    figure;
    set(gcf, 'Position', [100, 100, 620, 900]);
    colormap(gcf, cmap_name);      % set once for the whole figure
    cmap = colormap(gcf);          % retrieve for manual color indexing

    ax = gobjects(3, 1);

    for sp = 1:3
        ax(sp) = subplot(3, 1, sp);
        hold on;

        sel = sel_by_block{sp};
        for j = 1:numel(sel)
            idx = sel(j);
            if idx > numel(filtered_motion) || isempty(filtered_motion{idx}); continue; end

            x = filtered_motion{idx}(1, :);
            y = filtered_motion{idx}(2, :);
            v = filtered_velocity{idx};       % session-normalised [0,1]

            % Clamp to vel_clim, then rescale that range to [0,1] for
            % colormap indexing — values above vel_clim(2) saturate at top colour
            v_clamped = v;
            v_clamped(v_clamped < vel_clim(1)) = vel_clim(1);
            v_clamped(v_clamped > vel_clim(2)) = vel_clim(2);
            v_scaled  = (v_clamped - vel_clim(1)) / (vel_clim(2) - vel_clim(1));
            cmap_rows  = max(1, round(v_scaled * (size(cmap,1) - 1)) + 1);
            dot_colors = cmap(cmap_rows, :);

            % Thin black line underneath for continuity
            plot(x, y, 'k-', 'LineWidth', 0.5);

            % Velocity-coloured dots
            scatter(x, y, 25, dot_colors, 'filled', 'MarkerEdgeColor', 'none');

            % Start marker: filled black square
            scatter(x(1), y(1), 150, 'k', 's', 'filled');

            % Choice marker: filled black circle
            cx = choice_times{idx}(1);
            cy = choice_times{idx}(2);
            if ~isnan(cx)
                scatter(cx, cy, 150, 'k', 'o', 'filled');
            end

            % End marker: filled black triangle
            scatter(x(end), y(end), 150, 'k', '^', 'filled');
        end

        % Arena shapes
        draw_arena(shapeData);

        % Fix colormap range — same on every subplot
        set(gca, 'CLim', vel_clim);
        title(subplot_titles{sp}, 'FontSize', 11);
        axis equal;
        box off;
        hold off;
    end

    % Compress subplot widths to make room for a shared colorbar
    for sp = 1:3
        pos = ax(sp).Position;
        ax(sp).Position = [pos(1), pos(2), pos(3) * 0.87, pos(4)];
    end

    % Single colorbar for the whole figure
    cb = colorbar(ax(2), 'Position', [0.89, 0.11, 0.025, 0.77]);
    cb.Limits      = vel_clim;
    cb.Label.String = 'Velocity (session min–max normalised)';
    cb.Label.FontSize = 10;

    sgtitle(fig_title, 'FontSize', 12);
end


function draw_arena(shapeData)
% Draws the arena circles and reward-port rectangles from shapeData.
    for k = 1:numel(shapeData)
        if strcmp(shapeData{k}.Type, 'Circle')
            viscircles(shapeData{k}.Center, shapeData{k}.Radius, ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 1);

        elseif strcmp(shapeData{k}.Type, 'Square')
            sc = shapeData{k}.Center;
            ss = shapeData{k}.Size;
            if strcmp(shapeData{k}.Location, 'left screen')
                ec = 'b';
            else
                ec = 'r';
            end
            % Position = [x, y, width, height] — matches original code exactly
            rectangle('Position', [sc - ss/2, ss], 'EdgeColor', ec, 'LineWidth', 1.5);
        end
    end
end
