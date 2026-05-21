% generate_abort_raster_RDT_sex_comparison.m
% Abort raster — one figure per animal × session (RDT D1 and D2)
%
% Y-axis (reversed): proper trial 1 at top, trial N at bottom
%   Proper trials = rows where bigSmall == 1.2 or 0.3 (Large or Small reward choice)
% X-axis: seconds from proper trial onset (stTime = 0 for every row)
%
% Aborts are identified by BehavData.type_binary:
%   type_binary == 1  →  Large reward abort  (red square)
%   type_binary == 2  →  Small reward abort  (blue square)
%
% An abort is assigned to proper trial n if its stTime falls within:
%   [stTime(trial n),  stTime(trial n+1))
% For the final proper trial the window extends to the end of the session.
%
% Requires in workspace (run females_vs_males_groups_correct.m first):
%   ChrimsonR_IDs, ChrimsonR_treatment_groups, valid_sessions
% Requires: final_behavior (load new_females_and_males_behavior_05012026.mat)

%% Settings
col_large_abort = [0.15, 0.35, 0.80];   % blue
col_small_abort =  [0.80, 0.15, 0.15];   % red

fig_px_w     = 760;
px_per_trial = 8;
fig_px_h_min = 450;
fig_px_h_max = 1000;

%% Loop: one figure per animal × session
n_animals = numel(ChrimsonR_IDs);
raster_sessions = {'RDT_D1', 'RDT_D2'};

for rs = 1:numel(raster_sessions)
    session_to_analyze = raster_sessions{rs};

    for ii = 1:n_animals
        currentanimal = ChrimsonR_IDs{ii};

        % Skip if session not valid or data missing
        if ~any(strcmp(valid_sessions{ii}, session_to_analyze)); continue; end
        if ~isfield(final_behavior, currentanimal) || ...
           ~isfield(final_behavior.(currentanimal), session_to_analyze)
            continue;
        end

        BD = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        % Check that abort column exists
        if ~ismember('type_binary', BD.Properties.VariableNames)
            fprintf('[%s %s] type_binary column missing — skipping abort raster\n', ...
                currentanimal, session_to_analyze);
            continue;
        end

        % ---- Define proper trials (Y-axis rows) ----
        proper_mask = BD.bigSmall == 1.2 | BD.bigSmall == 0.3;
        proper_rows = find(proper_mask);
        BD_proper   = BD(proper_rows, :);
        n_trials    = height(BD_proper);
        if n_trials == 0; continue; end

        % ---- Define abort rows ----
        is_large_abort = BD.type_binary == 1;
        is_small_abort = BD.type_binary == 2;
        is_any_abort   = is_large_abort | is_small_abort;

        % ---- Build time-window boundaries for each proper trial ----
        % Window for trial n: [stTime(n), stTime(n+1))
        % Last trial: [stTime(end), Inf)
        t_starts = BD_proper.stTime;
        t_ends   = [BD_proper.stTime(2:end); Inf];

        % ---- Figure ----
        fig_px_h = min(max(n_trials * px_per_trial + 100, fig_px_h_min), fig_px_h_max);
        figure;
        set(gcf, 'Position', [100, 100, fig_px_w, fig_px_h]);
        hold on;

        % Grey track spanning each trial's window (x = 0 to window length)
        for t = 1:n_trials
            win_len = t_ends(t) - t_starts(t);
            if ~isinf(win_len)
                plot([0, win_len], [t, t], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
            end
        end

        % Trial-start reference dots at x = 0
        scatter(zeros(n_trials, 1), (1:n_trials)', 12, [0.65 0.65 0.65], 'filled', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);

        % ---- Plot aborts within each proper trial's window ----
        ax_large = [];  ay_large = [];
        ax_small = [];  ay_small = [];

        for t = 1:n_trials
            t0  = t_starts(t);
            t1  = t_ends(t);

            % Abort rows in this window
            in_window = is_any_abort & BD.stTime >= t0 & BD.stTime < t1;

            % Large aborts
            la_times = BD.stTime(in_window & is_large_abort) - t0;
            if ~isempty(la_times)
                ax_large = [ax_large; la_times(:)];           %#ok<AGROW>
                ay_large = [ay_large; repmat(t, numel(la_times), 1)]; %#ok<AGROW>
            end

            % Small aborts
            sa_times = BD.stTime(in_window & is_small_abort) - t0;
            if ~isempty(sa_times)
                ax_small = [ax_small; sa_times(:)];           %#ok<AGROW>
                ay_small = [ay_small; repmat(t, numel(sa_times), 1)]; %#ok<AGROW>
            end
        end

        % Draw abort markers — filled squares
        if ~isempty(ax_large)
            scatter(ax_large, ay_large, 38, col_large_abort, 's', 'filled', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.88);
        end
        if ~isempty(ax_small)
            scatter(ax_small, ay_small, 38, col_small_abort, 's', 'filled', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.88);
        end

        % Block divider lines
        if ismember('Block', BD_proper.Properties.VariableNames)
            blk_edges = find(diff(BD_proper.Block) ~= 0);
            for bb = 1:numel(blk_edges)
                yline(blk_edges(bb) + 0.5, '--', ...
                    'Color', [0.40, 0.40, 0.40], 'LineWidth', 0.9);
                text(xlim * [0.98; 0.02], blk_edges(bb) + 0.5, ...
                    sprintf('  B%d', BD_proper.Block(blk_edges(bb) + 1)), ...
                    'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
                    'VerticalAlignment', 'bottom');
            end
        end

        % Axes formatting
        set(gca, 'YDir', 'reverse', 'FontSize', 11);
        ylim([0.5, n_trials + 0.5]);
        xlabel('Time from trial onset (s)', 'FontSize', 12);
        ylabel('Trial',                     'FontSize', 12);

        n_la = numel(ax_large);
        n_sa = numel(ax_small);
        title(sprintf('%s  —  %s  (%d large aborts, %d small aborts)', ...
            strrep(currentanimal,      '_', '\_'), ...
            strrep(session_to_analyze, '_', '\_'), n_la, n_sa), 'FontSize', 11);
        box off;

        % Legend
        hl(1) = scatter(nan, nan, 38, col_large_abort, 's', 'filled');
        hl(2) = scatter(nan, nan, 38, col_small_abort, 's', 'filled');
        legend(hl, {'Large abort', 'Small abort'}, ...
            'Location', 'eastoutside', 'FontSize', 9);

        hold off;
        clear hl ax_large ay_large ax_small ay_small
    end
end
