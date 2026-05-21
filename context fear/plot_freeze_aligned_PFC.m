%% PFC Pseudopopulation Activity Aligned to Freeze Bout Onset and Offset
%
% Data sources:
%   Imaging  : PFC_alone_imaging_data_07182025.mat  → 'final'
%              final.(mouse).(session).CNMFe_data.C_raw  [neurons × frames, 10 Hz]
%              final.(mouse).experimental_grp             'Experimental' | 'No Shock' | 'One Context'
%
%   Behavior : context_fear_PFC_data_12162024.mat   → 'final_DLC'
%              final_DLC.(mouse).(session).movement_data  table with freeze column [30 fps]
%
% Session context-block definitions (frame numbers, 1-indexed, 10 Hz):
%   D3  Neutral  : [1,1200 ; 2401,3600 ; 4801,6000]
%   D3  Aversive : [1201,2400 ; 3601,4800 ; 6001,7200]
%   D4  Aversive : [1,1200 ; 2401,3600 ; 4801,6000]
%   D4  Neutral  : [1201,2400 ; 3601,4800 ; 6001,7200]
%
% respclass encoding:  0 = unmodulated  |  1 = activated  |  2 = inhibited
%
% Outputs:
%   pseudopop.(grp_field).(evt).(ctx)   [n_neurons × n_win]  trial-averaged z-scores
%   respclass.(grp_field).(evt).(ctx)   [n_neurons × 1]      Wilcoxon classification
%
%   Figures:
%     - Heatmap figure per event type  (rows = contexts, cols = groups)
%     - Mean-trace figure per event type (one panel per context, groups overlaid)
%     - Classification summary (% activated / inhibited bar charts)

%% ── 1. Load data ─────────────────────────────────────────────────────────
% clear; clc;

% imaging_file = 'I:\MATLAB\my_repo\context fear\aggregate_data\PFC alone imaging\calcium\PFC_alone_imaging_data_07182025.mat';
% behav_file   = 'I:\MATLAB\my_repo\context fear\aggregate_data\context_fear_PFC_data_12162024.mat';
% fprintf('Loading imaging data...\n');
% load(imaging_file);   % → 'final'
% fprintf('Loading behavioral data...\n');
% load(behav_file);     % → 'final_DLC'

%% ── 2. Parameters ────────────────────────────────────────────────────────
sessions_to_analyze = {'D3', 'D4'};
imaging_fps         = 10;    % Hz
behav_fps           = 30;    % Hz
ca_data_type        = 'C_raw';

pre_s         = 5;   % seconds before event (window + baseline)
post_s        = 10;  % seconds after event
bl_end_offset = 1;   % baseline ends this many seconds before event

min_bout_s = 2;      % minimum freeze bout duration (s) to include

% Wilcoxon test windows (seconds, relative to event)
wilcox_pre_win  = [-pre_s, -bl_end_offset];  % e.g. [-5, -1]
wilcox_post_win = [0, 3];                    % e.g. [0,  3]

% Context block definitions (1-indexed imaging frames, 10 Hz)
session_ctx.D3.neutral  = [1,1200 ; 2401,3600 ; 4801,6000];
session_ctx.D3.aversive = [1201,2400 ; 3601,4800 ; 6001,7200];
session_ctx.D4.aversive = [1,1200 ; 2401,3600 ; 4801,6000];
session_ctx.D4.neutral  = [1201,2400 ; 3601,4800 ; 6001,7200];

% Group labels (as stored in experimental_grp) and safe struct fieldnames
grp_labels   = {'Experimental', 'No Shock',  'One Context'};
grp_fields   = {'Experimental', 'No_Shock',  'One_Context'};
group_colors = {[0.85 0.15 0.15], [0.15 0.45 0.85], [0.20 0.70 0.35]};

ctx_names  = {'aversive', 'neutral', 'combined'};
ctx_labels = {'Aversive context', 'Neutral context', 'Both contexts'};
evt_names  = {'onset', 'offset'};
evt_labels = {'Freeze onset', 'Freeze offset'};

pre_frames  = pre_s  * imaging_fps;
post_frames = post_s * imaging_fps;
n_win       = pre_frames + post_frames;
ts          = (-pre_s : 1/imaging_fps : post_s - 1/imaging_fps);

bl_idx          = 1 : (pre_frames - bl_end_offset * imaging_fps);
min_bout_frames = round(min_bout_s * imaging_fps);

%% ── 3. Initialise storage ────────────────────────────────────────────────
for gg = 1:numel(grp_fields)
    for ev = 1:numel(evt_names)
        for cx = 1:numel(ctx_names)
            pseudopop.(grp_fields{gg}).(evt_names{ev}).(ctx_names{cx}) = [];
            respclass.(grp_fields{gg}).(evt_names{ev}).(ctx_names{cx}) = [];
        end
    end
end

%% ── 4. Main extraction loop ───────────────────────────────────────────────
imaging_animals = fieldnames(final);
behav_animals   = fieldnames(final_DLC);
common_animals  = intersect(imaging_animals, behav_animals);

fprintf('\nProcessing %d animals...\n\n', numel(common_animals));

for aa = 1:numel(common_animals)
    aname = common_animals{aa};

    if ~isfield(final.(aname), 'experimental_grp'), continue; end
    grp_label = final.(aname).experimental_grp;
    grp_idx   = find(strcmp(grp_label, grp_labels), 1);
    if isempty(grp_idx)
        fprintf('  Skipping %s — unknown group "%s"\n', aname, grp_label);
        continue;
    end
    grp_field = grp_fields{grp_idx};

    for ss = 1:numel(sessions_to_analyze)
        sess = sessions_to_analyze{ss};
        if ~isfield(final.(aname),     sess), continue; end
        if ~isfield(final_DLC.(aname), sess), continue; end

        % ── Imaging traces ─────────────────────────────────────────────────
        ca        = final.(aname).(sess).CNMFe_data.(ca_data_type);
        n_neurons = size(ca, 1);
        n_frames  = size(ca, 2);

        % ── Behavioral freeze vector ───────────────────────────────────────
        T = final_DLC.(aname).(sess).movement_data;
        col_lower  = lower(T.Properties.VariableNames);
        freeze_col = find(contains(col_lower, 'freez') | contains(col_lower, 'immob'));
        if isempty(freeze_col)
            warning('No freeze column for %s %s — skipping.', aname, sess);
            continue;
        end
        freeze_behav   = double(T{:, freeze_col(1)});
        freeze_imaging = downsample_freeze(freeze_behav, behav_fps, imaging_fps, n_frames);

        % ── Find all freeze bouts ──────────────────────────────────────────
        d         = diff([0; freeze_imaging(:); 0]);
        starts    = find(d ==  1);
        ends      = find(d == -1) - 1;
        durations = ends - starts + 1;
        keep      = durations >= min_bout_frames;
        all_starts = starts(keep)';
        all_ends   = ends(keep)';
        if isempty(all_starts), continue; end

        % ── Context block sets for this session ────────────────────────────
        aversive_blocks = session_ctx.(sess).aversive;
        neutral_blocks  = session_ctx.(sess).neutral;
        ctx_block_sets  = {aversive_blocks, neutral_blocks, [aversive_blocks; neutral_blocks]};

        % ── Per-context, per-event extraction + classification ─────────────
        for cx = 1:numel(ctx_names)
            blks   = ctx_block_sets{cx};
            in_ctx = false(1, numel(all_starts));
            for bb = 1:size(blks, 1)
                in_ctx = in_ctx | ...
                    (all_starts >= blks(bb,1) & all_starts <= blks(bb,2));
            end
            cx_starts = all_starts(in_ctx);
            cx_ends   = all_ends(in_ctx);
            if isempty(cx_starts), continue; end

            for ev = 1:numel(evt_names)
                ev_frames = cx_starts;
                if ev == 2, ev_frames = cx_ends; end

                % Remove events too close to session boundaries
                valid = ev_frames > pre_frames & ...
                        (ev_frames + post_frames - 1) <= n_frames;
                ev_frames = ev_frames(valid);
                if isempty(ev_frames), continue; end

                % Build per-trial z-scored traces for every neuron
                zall_array   = cell(n_neurons, 1);
                neuron_means = NaN(n_neurons, n_win);
                for nn = 1:n_neurons
                    raw     = extract_trials(ca(nn,:), ev_frames, pre_frames, post_frames);
                    ztrials = zscore_to_baseline(raw, bl_idx);
                    zall_array{nn}      = ztrials;
                    neuron_means(nn,:)  = nanmean(ztrials, 1);
                end

                % Wilcoxon classification (requires ≥ 2 trials per neuron)
                modulation = wilcoxon_analyzeNeuronalResponses_fn( ...
                    zall_array, ts, wilcox_pre_win, wilcox_post_win);

                % Append to pseudopop and respclass
                pseudopop.(grp_field).(evt_names{ev}).(ctx_names{cx}) = ...
                    [pseudopop.(grp_field).(evt_names{ev}).(ctx_names{cx}); neuron_means];
                respclass.(grp_field).(evt_names{ev}).(ctx_names{cx}) = ...
                    [respclass.(grp_field).(evt_names{ev}).(ctx_names{cx}); modulation];
            end
        end

        fprintf('  %s | %-12s | %s | %d neurons | %d bouts\n', ...
            aname, grp_label, sess, n_neurons, numel(all_starts));
    end
end

fprintf('\nExtraction complete.\n');

%% ── 5. Print classification summary to command window ────────────────────
fprintf('\n%s\n', repmat('─',1,72));
fprintf('%-14s  %-8s  %-10s  %6s  %7s  %7s  %7s\n', ...
    'Group','Event','Context','Total','Act (%)','Inh (%)','None (%)');
fprintf('%s\n', repmat('─',1,72));
for gg = 1:numel(grp_fields)
    for ev = 1:numel(evt_names)
        for cx = 1:numel(ctx_names)
            rc  = respclass.(grp_fields{gg}).(evt_names{ev}).(ctx_names{cx});
            if isempty(rc), continue; end
            n   = numel(rc);
            na  = sum(rc == 1);
            ni  = sum(rc == 2);
            nn0 = sum(rc == 0);
            fprintf('%-14s  %-8s  %-10s  %6d  %6.1f%%  %6.1f%%  %6.1f%%\n', ...
                grp_labels{gg}, evt_names{ev}, ctx_names{cx}, ...
                n, 100*na/n, 100*ni/n, 100*nn0/n);
        end
    end
end
fprintf('%s\n\n', repmat('─',1,72));

%% ── 6. Build colormaps ───────────────────────────────────────────────────
half     = 128;
bwr_cmap = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1); ...
            ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];

%% ── 7. Heatmap figures: one per event type ───────────────────────────────
% Layout: rows = contexts, cols = groups
% Neurons sorted: activated (top) → unmodulated (middle) → inhibited (bottom)

for ev = 1:numel(evt_names)
    figure('Name', sprintf('Heatmaps — %s', evt_labels{ev}), ...
           'Position', [30 30 1400 900]);
    tiledlayout(numel(ctx_names), numel(grp_fields), ...
                'TileSpacing', 'compact', 'Padding', 'compact');

    for cx = 1:numel(ctx_names)
        for gg = 1:numel(grp_fields)
            ax  = nexttile;
            dat = pseudopop.(grp_fields{gg}).(evt_names{ev}).(ctx_names{cx});
            rc  = respclass.(grp_fields{gg}).(evt_names{ev}).(ctx_names{cx});

            if isempty(dat)
                title(ax, sprintf('%s\n%s\n(no data)', grp_labels{gg}, ctx_labels{cx}));
                axis(ax, 'off');
                continue;
            end

            % Sort by class then by post-event mean within each class
            sort_win = ts >= 0 & ts < wilcox_post_win(2);
            post_mean = nanmean(dat(:, sort_win), 2);

            act_idx  = find(rc == 1);  [~,o] = sort(post_mean(act_idx),  'descend'); act_idx  = act_idx(o);
            none_idx = find(rc == 0);  [~,o] = sort(post_mean(none_idx), 'descend'); none_idx = none_idx(o);
            inh_idx  = find(rc == 2);  [~,o] = sort(post_mean(inh_idx),  'ascend');  inh_idx  = inh_idx(o);
            sort_idx = [act_idx; none_idx; inh_idx];

            imagesc(ax, ts, 1:size(dat,1), dat(sort_idx, :));
            colormap(ax, bwr_cmap);
            clim_val = max(prctile(abs(dat(:)), 98), 0.1);
            clim(ax, [-clim_val, clim_val]);
            cb = colorbar(ax, 'eastoutside');
            cb.Label.String = 'Z-score';

            % Draw lines separating classes
            n_act  = numel(act_idx);
            n_none = numel(none_idx);
            if n_act > 0 && (n_none + numel(inh_idx)) > 0
                yline(ax, n_act + 0.5, 'w-', 'LineWidth', 1);
            end
            if n_act + n_none > 0 && numel(inh_idx) > 0
                yline(ax, n_act + n_none + 0.5, 'w-', 'LineWidth', 1);
            end

            xline(ax, 0, 'w--', 'LineWidth', 1.5);
            xlabel(ax, 'Time (s)');
            ylabel(ax, 'Neuron');
            title(ax, sprintf('%s | %s\nn=%d  act=%d  inh=%d', ...
                grp_labels{gg}, ctx_labels{cx}, ...
                numel(rc), sum(rc==1), sum(rc==2)));
            set(ax, 'FontSize', 8);
        end
    end

    sgtitle(sprintf('PFC Pseudopopulation — %s', evt_labels{ev}), ...
            'FontSize', 13, 'FontWeight', 'bold');
end

%% ── 8. Mean-trace figures: one per event type ────────────────────────────
% Panels: all groups overlaid per context
% Separate lines for activated / all neurons within each group

for ev = 1:numel(evt_names)
    figure('Name', sprintf('Mean traces — %s', evt_labels{ev}), ...
           'Position', [60 60 1300 450]);
    tiledlayout(1, numel(ctx_names), 'TileSpacing', 'compact', 'Padding', 'compact');

    for cx = 1:numel(ctx_names)
        ax = nexttile;
        hold(ax, 'on');
        leg_handles = [];
        leg_strs    = {};

        for gg = 1:numel(grp_fields)
            dat = pseudopop.(grp_fields{gg}).(evt_names{ev}).(ctx_names{cx});
            if isempty(dat), continue; end

            mu  = nanmean(dat, 1);
            sem = nanstd(dat, 0, 1) / sqrt(size(dat, 1));
            col = group_colors{gg};

            fill(ax, [ts, fliplr(ts)], [mu+sem, fliplr(mu-sem)], ...
                col, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
            h = plot(ax, ts, mu, 'Color', col, 'LineWidth', 2);

            leg_handles(end+1) = h; %#ok<AGROW>
            leg_strs{end+1}    = sprintf('%s (n=%d)', grp_labels{gg}, size(dat,1)); %#ok<AGROW>
        end

        xline(ax, 0, 'k--', 'LineWidth', 1.2);
        yline(ax, 0, 'k:',  'LineWidth', 0.8);
        xlabel(ax, 'Time from event (s)');
        ylabel(ax, 'Z-score');
        title(ax, ctx_labels{cx});
        if ~isempty(leg_handles)
            legend(ax, leg_handles, leg_strs, 'Location', 'best', 'FontSize', 8);
        end
        set(ax, 'FontSize', 10);
        hold(ax, 'off');
    end

    sgtitle(sprintf('PFC Population Mean ± SEM — %s', evt_labels{ev}), ...
            'FontSize', 13, 'FontWeight', 'bold');
end

%% ── 9. Classification summary bar charts ─────────────────────────────────
% One figure per context, rows = events, cols = groups
% Bars: % activated (filled) and % inhibited (hatched / lighter)

for cx = 1:numel(ctx_names)
    figure('Name', sprintf('Classification — %s', ctx_labels{cx}), ...
           'Position', [90 90 900 500]);
    tiledlayout(1, numel(evt_names), 'TileSpacing', 'compact', 'Padding', 'compact');

    for ev = 1:numel(evt_names)
        ax = nexttile;
        hold(ax, 'on');

        pct_act = zeros(1, numel(grp_fields));
        pct_inh = zeros(1, numel(grp_fields));
        ns      = zeros(1, numel(grp_fields));

        for gg = 1:numel(grp_fields)
            rc = respclass.(grp_fields{gg}).(evt_names{ev}).(ctx_names{cx});
            if isempty(rc), continue; end
            ns(gg)      = numel(rc);
            pct_act(gg) = 100 * sum(rc == 1) / numel(rc);
            pct_inh(gg) = 100 * sum(rc == 2) / numel(rc);
        end

        x = 1:numel(grp_fields);
        b1 = bar(ax, x - 0.2, pct_act, 0.35, 'FaceColor', 'flat');
        b2 = bar(ax, x + 0.2, pct_inh, 0.35, 'FaceColor', 'flat');

        for gg = 1:numel(grp_fields)
            b1.CData(gg,:) = group_colors{gg};
            b2.CData(gg,:) = group_colors{gg} * 0.45 + 0.55;  % lighter shade for inhibited
        end

        % Neuron count labels above each pair
        for gg = 1:numel(grp_fields)
            if ns(gg) > 0
                text(ax, gg, max(pct_act(gg), pct_inh(gg)) + 1, ...
                    sprintf('n=%d', ns(gg)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 7);
            end
        end

        set(ax, 'XTick', x, 'XTickLabel', grp_labels, 'FontSize', 9);
        xtickangle(ax, 20);
        ylabel(ax, '% neurons');
        ylim(ax, [0, max([pct_act, pct_inh]) * 1.25 + 5]);
        title(ax, evt_labels{ev});
        legend(ax, [b1, b2], {'Activated', 'Inhibited'}, 'Location', 'northeast', 'FontSize', 8);
        hold(ax, 'off');
    end

    sgtitle(sprintf('Freeze-responsive neurons — %s', ctx_labels{cx}), ...
            'FontSize', 12, 'FontWeight', 'bold');
end

fprintf('Plotting complete.\n');

%% ── Local functions ──────────────────────────────────────────────────────

function freeze_ds = downsample_freeze(freeze_behav, behav_fps, imaging_fps, n_imaging_frames)
    ds_factor = behav_fps / imaging_fps;
    freeze_ds = zeros(n_imaging_frames, 1);
    n_behav   = numel(freeze_behav);
    for fi = 1:n_imaging_frames
        b_start = round((fi - 1) * ds_factor) + 1;
        b_end   = min(round(fi * ds_factor), n_behav);
        if b_start > n_behav, break; end
        freeze_ds(fi) = any(freeze_behav(b_start:b_end) > 0);
    end
end

function trials = extract_trials(trace, event_frames, pre_frames, post_frames)
    n_win  = pre_frames + post_frames;
    trials = NaN(numel(event_frames), n_win);
    for tt = 1:numel(event_frames)
        s = event_frames(tt) - pre_frames;
        e = event_frames(tt) + post_frames - 1;
        if s >= 1 && e <= numel(trace)
            trials(tt, :) = trace(s:e);
        end
    end
end

function ztrials = zscore_to_baseline(trials, bl_idx)
    bl_data = trials(:, bl_idx);
    zb      = nanmean(bl_data(:));
    zsd     = nanstd(bl_data(:));
    if isnan(zsd) || zsd < 1e-10, zsd = 1; end
    ztrials = (trials - zb) ./ zsd;
end
