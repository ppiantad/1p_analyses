% locomotor_velocity_heatmap_GoPro.m
%
% Heatmap of session-normalised velocity — one row per trial.
%   X-axis : time from trial onset (s)
%   Y-axis : trial number, reversed (trial 1 at top)
%   Colour : velocity mapped through cmap_name, clamped to vel_clim
%   Overlaid scatter points:
%       filled black diamond  — choice time
%       filled white diamond  — reward collection time  (black edge)
%
% Two figures are produced: one for Large reward trials, one for Small,
% each with 3 subplots (Block 1 / 2 / 3).  A single shared colorbar
% appears at the right of every figure.
%
% Requires in workspace: final_SLEAP, final_behavior

%% ---- Settings ----
select_mouse       = 'RDT_M_4';
session_to_analyze = 'RDT_D1';
cmap_name          = 'jet';
vel_clim           = [0, 0.6];   % colormap range in session-normalised units
fs_cam             = 30;         % camera frame rate (Hz)

%% ---- Load data ----
SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;
BehavData  = final_behavior.(select_mouse).(session_to_analyze).uv.BehavData;

onset_trials       = BehavData.stTime';
choice_trials      = BehavData.choiceTime';
offset_trials      = BehavData.collectionTime';
time_ranges_trials = [onset_trials; choice_trials; offset_trials];

%% ---- Session-level velocity normalisation (once, before trial extraction) ----
vel_raw  = SLEAP_data.vel_cm_s;
vel_min  = min(vel_raw);
vel_max  = max(vel_raw);
vel_norm = (vel_raw - vel_min) / max(vel_max - vel_min, eps);   % [0,1]

%% ---- Extract per-trial velocity and event times (relative to onset) ----
n_trials_total = size(time_ranges_trials, 2);
filtered_velocity = cell(1, n_trials_total);
t_choice_rel      = nan(1, n_trials_total);   % seconds from onset
t_collect_rel     = nan(1, n_trials_total);   % seconds from onset

for j = 1:n_trials_total
    t_onset  = time_ranges_trials(1, j);
    t_choice = time_ranges_trials(2, j);
    t_offset = time_ranges_trials(3, j);

    if isinf(t_offset) || isnan(t_offset); continue; end

    mask = SLEAP_data.idx_time > t_onset & SLEAP_data.idx_time < t_offset;
    if ~any(mask); continue; end

    filtered_velocity{j} = vel_norm(mask)';           % row vector, [0,1]
    t_choice_rel(j)      = t_choice - t_onset;        % s from onset
    t_collect_rel(j)     = t_offset - t_onset;        % s from onset
end

%% ---- Trial indices by reward × block ----
large_b1 = find(BehavData.bigSmall == 1.2 & BehavData.Block == 1);
large_b2 = find(BehavData.bigSmall == 1.2 & BehavData.Block == 2);
large_b3 = find(BehavData.bigSmall == 1.2 & BehavData.Block == 3);
small_b1 = find(BehavData.bigSmall == 0.3 & BehavData.Block == 1);
small_b2 = find(BehavData.bigSmall == 0.3 & BehavData.Block == 2);
small_b3 = find(BehavData.bigSmall == 0.3 & BehavData.Block == 3);

%% ---- Figure: LARGE reward ----
plot_velocity_heatmap( ...
    {large_b1, large_b2, large_b3}, ...
    {'Large Reward — Block 1', 'Large Reward — Block 2', 'Large Reward — Block 3'}, ...
    filtered_velocity, t_choice_rel, t_collect_rel, fs_cam, cmap_name, vel_clim, ...
    sprintf('%s  %s  — Large Reward Velocity Heatmap', ...
        strrep(select_mouse,'_','\_'), strrep(session_to_analyze,'_','\_')));

%% ---- Figure: SMALL reward ----
plot_velocity_heatmap( ...
    {small_b1, small_b2, small_b3}, ...
    {'Small Reward — Block 1', 'Small Reward — Block 2', 'Small Reward — Block 3'}, ...
    filtered_velocity, t_choice_rel, t_collect_rel, fs_cam, cmap_name, vel_clim, ...
    sprintf('%s  %s  — Small Reward Velocity Heatmap', ...
        strrep(select_mouse,'_','\_'), strrep(session_to_analyze,'_','\_')));


%% ================================================================
%  Local functions
%  ================================================================

function plot_velocity_heatmap(sel_by_block, subplot_titles, ...
        filtered_velocity, t_choice_rel, t_collect_rel, ...
        fs_cam, cmap_name, vel_clim, fig_title)
% Creates a 3-subplot figure (one per block).
% Each row of the heatmap = one trial's session-normalised velocity.
% Overlaid markers: filled black diamond = choice time;
%                   filled white diamond (black edge) = collection time.

    col_choice  = [0.05 0.05 0.05];   % near-black fill for choice marker
    col_collect = [1.00 1.00 1.00];   % white fill for collection marker

    figure;
    set(gcf, 'Position', [100, 100, 720, 900]);
    colormap(gcf, cmap_name);

    ax = gobjects(3, 1);

    for sp = 1:3
        ax(sp) = subplot(3, 1, sp);
        hold on;

        trial_idxs = sel_by_block{sp};
        % Keep only trials that have extracted velocity data
        has_data = ~cellfun(@isempty, filtered_velocity(trial_idxs));
        valid    = trial_idxs(has_data);
        n_t      = numel(valid);

        if n_t == 0
            title(subplot_titles{sp}, 'FontSize', 11);
            axis off;
            hold off;
            continue;
        end

        % ---- Build heatmap matrix (n_t × max_frames), NaN-padded ----
        lengths = cellfun(@numel, filtered_velocity(valid));
        max_len = max(lengths);
        hmap    = nan(n_t, max_len);

        for ti = 1:n_t
            v = filtered_velocity{valid(ti)};
            v_clamped = v;
            v_clamped(v_clamped < vel_clim(1)) = vel_clim(1);
            v_clamped(v_clamped > vel_clim(2)) = vel_clim(2);
            hmap(ti, 1:numel(v_clamped)) = v_clamped;
        end

        % Time axis in seconds (frame 1 → t = 0)
        t_axis = (0 : max_len - 1) / fs_cam;

        % ---- Draw heatmap ----
        % NaN padding regions appear in the axes background colour (grey)
        set(ax(sp), 'Color', [0.82 0.82 0.82]);
        imagesc(t_axis, 1:n_t, hmap);
        set(ax(sp), 'CLim', vel_clim);
        set(ax(sp), 'YDir', 'reverse');    % trial 1 at top

        % ---- Overlay event markers ----
        for ti = 1:n_t
            idx = valid(ti);

            tc = t_choice_rel(idx);
            tl = t_collect_rel(idx);

            % Choice time: filled black diamond
            if ~isnan(tc) && ~isinf(tc) && tc >= 0
                scatter(tc, ti, 55, col_choice, 'd', 'filled', ...
                    'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
            end

            % Collection time: filled white diamond with black edge
            if ~isnan(tl) && ~isinf(tl) && tl >= 0
                scatter(tl, ti, 55, col_collect, 'd', 'filled', ...
                    'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
            end
        end

        % ---- Axis formatting ----
        xlim([0, t_axis(end) + 0.5]);
        ylim([0.5, n_t + 0.5]);
        xlabel('Time from trial onset (s)', 'FontSize', 10);
        ylabel('Trial',                     'FontSize', 10);
        title(subplot_titles{sp}, 'FontSize', 11);
        box off;
        hold off;
    end

    % ---- Compress subplot widths to make room for shared colorbar ----
    for sp = 1:3
        pos = ax(sp).Position;
        ax(sp).Position = [pos(1), pos(2), pos(3) * 0.87, pos(4)];
    end

    % ---- Single colorbar for the whole figure ----
    cb = colorbar(ax(2), 'Position', [0.89, 0.11, 0.025, 0.77]);
    cb.Limits            = vel_clim;
    cb.Label.String      = 'Velocity (session min–max normalised)';
    cb.Label.FontSize    = 10;

    % ---- Legend (drawn on middle subplot, hold re-enabled to avoid wiping heatmap) ----
    hold(ax(2), 'on');
    hl(1) = scatter(ax(2), nan, nan, 55, col_choice,  'd', 'filled', ...
        'MarkerEdgeColor', 'w');
    hl(2) = scatter(ax(2), nan, nan, 55, col_collect, 'd', 'filled', ...
        'MarkerEdgeColor', [0.1 0.1 0.1]);
    legend(ax(2), hl, {'Choice time', 'Collection time'}, ...
        'Location', 'northeast', 'FontSize', 8);
    hold(ax(2), 'off');

    sgtitle(fig_title, 'FontSize', 12);
end
