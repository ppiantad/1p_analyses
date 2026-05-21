% generate_raster_figs_RDT_sex_comparison.m
% Trial timing raster — one figure per animal × session (RDT D1 and D2)
%
% Y-axis (reversed): trial 1 at top, trial N at bottom
% X-axis: seconds from trial onset (stTime = 0 for every row)
% Diamonds   — choice time:    Large = blue,       Small = red
% Pentagrams — collect time:   Large = light-blue, Small = light-red
% Yellow '*' — shock trial (closest MATLAB marker to a lightning bolt;
%              placed at collection time, or choice time if no collection)
% Thin grey line connects onset → collection for each trial as a track.
% Dashed lines separate blocks.
%
% Requires in workspace (run females_vs_males_groups_correct.m first):
%   ChrimsonR_IDs, ChrimsonR_treatment_groups, valid_sessions
% Requires: final_behavior (load new_females_and_males_behavior_05012026.mat)

%% Settings — adjust colours and figure width here
raster_col_large_choice  = [0.18, 0.36, 0.82];   % blue
raster_col_small_choice  = [0.82, 0.18, 0.18];   % red
raster_col_large_collect = [0.62, 0.78, 1.00];   % light blue
raster_col_small_collect = [1.00, 0.62, 0.62];   % light red
raster_col_shock         = [1.00, 0.84, 0.00];   % yellow

fig_px_w    = 760;    % figure width (pixels)
px_per_trial = 8;     % figure height scales with trial count
fig_px_h_min = 450;
fig_px_h_max = 1000;

%% Loop: one figure per animal × session
n_animals = numel(ChrimsonR_IDs);

raster_sessions = {'RDT_D1', 'RDT_D2'};

for rs = 1:numel(raster_sessions)
    session_to_analyze = raster_sessions{rs};

    for ii = 1:n_animals
        currentanimal = ChrimsonR_IDs{ii};

        % Skip if session not valid for this animal
        if ~any(strcmp(valid_sessions{ii}, session_to_analyze)); continue; end

        % Skip if behavioral data missing
        if ~isfield(final_behavior, currentanimal) || ...
           ~isfield(final_behavior.(currentanimal), session_to_analyze)
            continue;
        end

        BD_raw = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        % Proper trials only: Large (1.2) or Small (0.3) reward choice
        % Excludes omissions (bigSmall = NaN) and any non-reward rows
        proper = BD_raw.bigSmall == 1.2 | BD_raw.bigSmall == 0.3;
        BD     = BD_raw(proper, :);
        if height(BD) == 0; continue; end

        n_trials = height(BD);
        trial_y  = (1:n_trials)';

        % Times relative to each trial's own onset
        t_choice  = BD.choiceTime     - BD.stTime;
        t_collect = BD.collectionTime - BD.stTime;

        is_large = BD.bigSmall == 1.2;
        is_small = BD.bigSmall == 0.3;
        is_shock = BD.shock    == 1;
        ok_col   = ~isinf(t_collect) & ~isnan(t_collect);   % valid collection time

        % ---- Figure ----
        fig_px_h = min(max(n_trials * px_per_trial + 100, fig_px_h_min), fig_px_h_max);
        figure;
        set(gcf, 'Position', [100, 100, fig_px_w, fig_px_h]);
        hold on;

        % Grey tracks: onset (x=0) → collection time (or choice if no collection)
        lx = [];  ly = [];
        for t = 1:n_trials
            if ok_col(t)
                x_end = t_collect(t);
            else
                x_end = t_choice(t);
            end
            lx = [lx,  0, x_end,      NaN]; %#ok<AGROW>
            ly = [ly, trial_y(t), trial_y(t), NaN]; %#ok<AGROW>
        end
        plot(lx, ly, 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.5);

        % Trial-start reference dots (grey, x = 0)
        scatter(zeros(n_trials, 1), trial_y, 12, [0.65 0.65 0.65], 'filled', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);

        % Choice markers — filled diamonds
        if any(is_large)
            scatter(t_choice(is_large), trial_y(is_large), 38, ...
                raster_col_large_choice, 'd', 'filled', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.88);
        end
        if any(is_small)
            scatter(t_choice(is_small), trial_y(is_small), 38, ...
                raster_col_small_choice, 'd', 'filled', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.88);
        end

        % Shock markers — yellow '*' at choice time (shock occurs at moment of choice)
        if any(is_shock)
            scatter(t_choice(is_shock), trial_y(is_shock), 120, ...
                raster_col_shock, '*', 'LineWidth', 1.8);
        end

        % Collection markers — filled pentagrams (stars)
        % Shock collection times are unobscured because shock is plotted at choice time
        if any(is_large & ok_col)
            scatter(t_collect(is_large & ok_col), trial_y(is_large & ok_col), 38, ...
                raster_col_large_collect, 'p', 'filled', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.88);
        end
        if any(is_small & ok_col)
            scatter(t_collect(is_small & ok_col), trial_y(is_small & ok_col), 38, ...
                raster_col_small_collect, 'p', 'filled', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.88);
        end

        % Block divider lines
        if ismember('Block', BD.Properties.VariableNames)
            blk_edges = find(diff(BD.Block) ~= 0);
            for bb = 1:numel(blk_edges)
                yline(blk_edges(bb) + 0.5, '--', ...
                    'Color', [0.40, 0.40, 0.40], 'LineWidth', 0.9);
                text(xlim * [0.98; 0.02], blk_edges(bb) + 0.5, ...
                    sprintf('  B%d', BD.Block(blk_edges(bb) + 1)), ...
                    'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
                    'VerticalAlignment', 'bottom');
            end
        end

        % Axes formatting
        set(gca, 'YDir', 'reverse', 'FontSize', 11);
        ylim([0.5, n_trials + 0.5]);
        xlabel('Time from trial onset (s)', 'FontSize', 12);
        ylabel('Trial',                     'FontSize', 12);
        title(sprintf('%s  —  %s  (n = %d trials)', ...
            strrep(currentanimal,      '_', '\_'), ...
            strrep(session_to_analyze, '_', '\_'), n_trials), 'FontSize', 12);
        box off;

        % Legend
        hl(1) = scatter(nan, nan, 38,  raster_col_large_choice,  'd', 'filled');
        hl(2) = scatter(nan, nan, 38,  raster_col_small_choice,  'd', 'filled');
        hl(3) = scatter(nan, nan, 38,  raster_col_large_collect, 'p', 'filled');
        hl(4) = scatter(nan, nan, 38,  raster_col_small_collect, 'p', 'filled');
        hl(5) = scatter(nan, nan, 100, raster_col_shock, '*', 'LineWidth', 1.8);
        legend(hl, {'Large choice', 'Small choice', 'Large collection', ...
            'Small collection', 'Shock'}, ...
            'Location', 'eastoutside', 'FontSize', 9);

        hold off;
        clear hl
    end
end
