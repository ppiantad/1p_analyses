% generate_velocity_aligned_figs_RDT_sex_comparison.m
%
% Locomotor velocity (cm/s) aligned to choice time, split by:
%   Session  : RDT_D1, RDT_D2
%   Block    : 1, 2, 3
%   Reward   : Large (1.2 pellets), Small (0.3 pellets)
%   Sex      : Female (green), Male (grey)
%
% Two figures produced (Large reward / Small reward), each with
%   2 rows (D1 top, D2 bottom) × 3 columns (Block 1, 2, 3).
% Each subplot: group mean ± SEM via shadedErrorBar.
% Vertical dotted line at t = 0 (choice); dashed lines at median
%   collection latency per sex.
% Shared Y-axis within each figure, derived from the data.
%
% Requires workspace from females_vs_males_groups_correct.m:
%   ChrimsonR_IDs, ChrimsonR_treatment_groups, large_screen_side,
%   valid_sessions, session_weight_raw_grams, freefeed_weight
% Requires: final_SLEAP, final_behavior

clc

%% ---- Settings ----
evtWin = [-8,  8];        % s around choice time
fs_cam = 30;              % camera frame rate (Hz)
dt     = 1 / fs_cam;
ts     = evtWin(1) : dt : evtWin(2) - dt;   % common time axis

sg_ord = 9;               % Savitzky-Golay polynomial order
sg_len = 21;              % Savitzky-Golay window length (odd)

reward_types = [1.2,     0.3];
reward_names = {'Large', 'Small'};
sessions     = {'RDT_D1', 'RDT_D2'};

%% ---- Figure / colour settings ----
fig_width  = 900;
fig_height = 520;

col_f_d1 = [0.55, 0.88, 0.55];
col_f_d2 = [0.05, 0.42, 0.12];
col_m_d1 = [0.72, 0.72, 0.72];
col_m_d2 = [0.28, 0.28, 0.28];
% day_colors{sex_row, day_col}:  row 1 = Female, row 2 = Male
day_colors = {col_f_d1, col_f_d2; ...
              col_m_d1, col_m_d2};

%% ---- Body-weight filter (matches other scripts) ----
bw_threshold      = 75;   % % of free-feed weight; 0 = include all
session_weight_mat = cell2mat(session_weight_raw_grams');   % n × 2
% weight_pct_D1     = (session_weight_mat(:,1) ./ freefeed_weight') * 100;
% weight_pct_D2     = (session_weight_mat(:,2) ./ freefeed_weight') * 100;

weight_pct_D1 = [(session_weight_mat(:,1) ./ freefeed_weight') * 100]+5;
weight_pct_D2 = [(session_weight_mat(:,2) ./ freefeed_weight') * 100]+5;

bw_include        = (isnan(weight_pct_D1) | weight_pct_D1 >= bw_threshold) & ...
                    (isnan(weight_pct_D2) | weight_pct_D2 >= bw_threshold);

%% SLEAP-specific exclusions (video quality — separate from valid_sessions)
% Format: each row is {session, animal_ID}
vel_exclusions = { ...
    'RDT_D2', 'RDT_M_4'; ...
    'RDT_D1', 'RDT_M_5'; ...
};

%% ---- Extraction: per-animal mean velocity aligned to choice ----
n_animals = numel(ChrimsonR_IDs);
is_female = strcmp(ChrimsonR_treatment_groups, 'Female');
n_ts      = numel(ts);

% Preallocate storage: vel_mean{session, reward, block}(animal, time)
vel_mean = cell(2, 2, 3);
coll_med = cell(2, 2, 3);   % median collection latency per animal (s from choice)
for ss = 1:2; for rr = 1:2; for bb = 1:3
    vel_mean{ss,rr,bb} = nan(n_animals, n_ts);
    coll_med{ss,rr,bb} = nan(n_animals, 1);
end; end; end

for s_idx = 1:2
    session_to_analyze = sessions{s_idx};

    for ii = 1:n_animals
        currentanimal = ChrimsonR_IDs{ii};

        % ---- Skip checks ----
        if ~any(strcmp(valid_sessions{ii}, session_to_analyze)); continue; end
        if ~bw_include(ii); continue; end

        % SLEAP-specific exclusion
        if any(strcmp(vel_exclusions(:,1), session_to_analyze) & ...
               strcmp(vel_exclusions(:,2), currentanimal))
            fprintf('  [Excluding %s %s — SLEAP exclusion]\n', currentanimal, session_to_analyze);
            continue;
        end

        if ~isfield(final_SLEAP, currentanimal) || ...
           ~isfield(final_SLEAP.(currentanimal), session_to_analyze); continue; end
        if ~isfield(final_behavior, currentanimal) || ...
           ~isfield(final_behavior.(currentanimal), session_to_analyze); continue; end

        SLEAP_data = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data_raw;
        BehavData  = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        % Smooth velocity and grab time vector
        vel_s    = sgolayfilt(SLEAP_data.vel_cm_s(:), sg_ord, sg_len);
        idx_time = SLEAP_data.idx_time(:);

        % Valid trial mask (exclude omissions, blank touches, and shock trials)
        valid_trial = true(height(BehavData), 1);
        if ismember('omissionALL', BehavData.Properties.VariableNames)
            valid_trial = valid_trial & BehavData.omissionALL ~= 1;
        end
        if ismember('Blank_Touch', BehavData.Properties.VariableNames)
            valid_trial = valid_trial & BehavData.Blank_Touch ~= 1;
        end
        if ismember('shock', BehavData.Properties.VariableNames)
            valid_trial = valid_trial & BehavData.shock ~= 1;
        end

        for ri = 1:2
            rew_val = reward_types(ri);
            for blk = 1:3
                t_mask = valid_trial & ...
                         BehavData.bigSmall == rew_val & ...
                         BehavData.Block    == blk;
                t_idxs = find(t_mask);
                if isempty(t_idxs); continue; end

                trial_mat = nan(numel(t_idxs), n_ts);
                coll_lats = nan(numel(t_idxs), 1);

                for ti = 1:numel(t_idxs)
                    jj       = t_idxs(ti);
                    choice_t = BehavData.choiceTime(jj);
                    coll_t   = BehavData.collectionTime(jj);
                    if isnan(choice_t) || isinf(choice_t); continue; end

                    seg_mask = idx_time >= choice_t + evtWin(1) & ...
                               idx_time <= choice_t + evtWin(2);
                    if sum(seg_mask) < 5; continue; end

                    t_rel = idx_time(seg_mask) - choice_t;
                    v_seg = vel_s(seg_mask);
                    trial_mat(ti, :) = interp1(t_rel, v_seg, ts, 'linear', NaN);

                    if ~isnan(coll_t) && ~isinf(coll_t)
                        coll_lats(ti) = coll_t - choice_t;
                    end
                end

                vel_mean{s_idx, ri, blk}(ii, :) = nanmean(trial_mat, 1);
                coll_med{s_idx, ri, blk}(ii)    = nanmedian(coll_lats);
            end
        end
    end
    fprintf('Finished extraction: %s\n', session_to_analyze);
end

%% ---- Sex / BW masks ----
f_mask = is_female  & bw_include;
m_mask = ~is_female & bw_include;
fprintf('Females: %d  |  Males: %d  (after BW filter)\n', sum(f_mask), sum(m_mask));

%% ---- Line plots: one figure per day, rows = reward type, cols = block ----
% Row 1 = Large reward, Row 2 = Small reward
for s_idx = 1:2
    figure;
    set(gcf, 'Position', [100, 100, fig_width, fig_height]);

    ax_handles   = gobjects(2, 3);
    all_vals_fig = [];

    col_f = day_colors{1, s_idx};
    col_m = day_colors{2, s_idx};

    for ri = 1:2
        for blk = 1:3
            sp = (ri - 1) * 3 + blk;
            ax_handles(ri, blk) = subplot(2, 3, sp);
            hold on;

            vm = vel_mean{s_idx, ri, blk};

            f_data = vm(f_mask, :);
            m_data = vm(m_mask, :);

            nf = max(sum(~isnan(f_data(:, 1))), 1);
            nm = max(sum(~isnan(m_data(:, 1))), 1);

            f_mean = nanmean(f_data, 1);
            m_mean = nanmean(m_data, 1);
            f_sem  = nanstd(f_data, 0, 1) ./ sqrt(nf);
            m_sem  = nanstd(m_data, 0, 1) ./ sqrt(nm);

            if any(~isnan(f_mean))
                shadedErrorBar(ts, f_mean, f_sem, ...
                    'lineProps', {'Color', col_f, 'LineWidth', 1.8}, 'transparent', 1);
            end
            if any(~isnan(m_mean))
                shadedErrorBar(ts, m_mean, m_sem, ...
                    'lineProps', {'Color', col_m, 'LineWidth', 1.8}, 'transparent', 1);
            end

            xline(0, 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            f_coll = nanmedian(coll_med{s_idx, ri, blk}(f_mask));
            m_coll = nanmedian(coll_med{s_idx, ri, blk}(m_mask));
            if ~isnan(f_coll) && f_coll >= evtWin(1) && f_coll <= evtWin(2)
                xline(f_coll, '--', 'Color', col_f, 'LineWidth', 1.0, 'HandleVisibility', 'off');
            end
            if ~isnan(m_coll) && m_coll >= evtWin(1) && m_coll <= evtWin(2)
                xline(m_coll, '--', 'Color', col_m, 'LineWidth', 1.0, 'HandleVisibility', 'off');
            end

            all_vals_fig = [all_vals_fig, ...
                f_mean + f_sem, f_mean - f_sem, ...
                m_mean + m_sem, m_mean - m_sem]; %#ok<AGROW>

            xlabel('Time from choice (s)', 'FontSize', 10);
            ylabel('Velocity (cm/s)',       'FontSize', 10);
            title(sprintf('%s — Block %d', reward_names{ri}, blk), 'FontSize', 11);
            xlim(evtWin);
            set(gca, 'XTick', evtWin(1):4:evtWin(2));
            box off;
            hold off;
        end
    end

    % Shared Y-axis
    finite_vals = all_vals_fig(isfinite(all_vals_fig));
    if ~isempty(finite_vals)
        ypad = range(finite_vals) * 0.06;
        y_lo = max(0, min(finite_vals) - ypad);
        y_hi = max(finite_vals) + ypad;
        for ri = 1:2
            for blk = 1:3
                ylim(ax_handles(ri, blk), [y_lo, y_hi]);
            end
        end
    end

    % Legend on top-right subplot
    hold(ax_handles(1, 3), 'on');
    lh(1) = plot(ax_handles(1,3), nan, nan, '-', 'Color', col_f, 'LineWidth', 2);
    lh(2) = plot(ax_handles(1,3), nan, nan, '-', 'Color', col_m, 'LineWidth', 2);
    legend(ax_handles(1,3), lh, {'Female', 'Male'}, 'Location', 'best', 'FontSize', 9);
    hold(ax_handles(1,3), 'off');

    sgtitle(sprintf('RDT D%d — velocity aligned to choice time', s_idx), 'FontSize', 13);
end

%% ---- Publication figure: Large reward D1, lines only, scale bars ----
% 1 row × 3 columns (one per block). No axes, no labels.
% Female (D1 green) and male (D1 grey) mean ± SEM lines.
% Vertical dotted line at t = 0 in each subplot.
% L-shaped scale bars (time + velocity) on the last subplot.

pub_fig_w   = 680;
pub_fig_h   = 190;
pub_x_scale = 2;     % time scale bar (s)
pub_y_scale = NaN;   % velocity scale bar (cm/s); NaN = auto ~25% of y range

figure;
set(gcf, 'Position', [100, 550, pub_fig_w, pub_fig_h]);

ax_pub    = gobjects(1, 3);
all_v_pub = [];

for blk = 1:3
    ax_pub(blk) = subplot(1, 3, blk);
    hold on;

    vm     = vel_mean{1, 1, blk};   % D1, Large
    f_data = vm(f_mask, :);
    m_data = vm(m_mask, :);

    nf = max(sum(~isnan(f_data(:,1))), 1);
    nm = max(sum(~isnan(m_data(:,1))), 1);

    f_mn = nanmean(f_data, 1);
    m_mn = nanmean(m_data, 1);
    f_se = nanstd(f_data, 0, 1) ./ sqrt(nf);
    m_se = nanstd(m_data, 0, 1) ./ sqrt(nm);

    if any(~isnan(f_mn))
        shadedErrorBar(ts, f_mn, f_se, ...
            'lineProps', {'Color', col_f_d1, 'LineWidth', 1.8}, 'transparent', 1);
    end
    if any(~isnan(m_mn))
        shadedErrorBar(ts, m_mn, m_se, ...
            'lineProps', {'Color', col_m_d1, 'LineWidth', 1.8}, 'transparent', 1);
    end

    % Dotted reference lines at t = 0 (choice) and y = 0
    xline(0, 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    yline(0, 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');

    xlim(evtWin);
    all_v_pub = [all_v_pub, f_mn+f_se, f_mn-f_se, m_mn+m_se, m_mn-m_se]; %#ok<AGROW>
    hold off;
end

% Shared y limits
finite_v_pub = all_v_pub(isfinite(all_v_pub));
if ~isempty(finite_v_pub)
    ypad_pub = range(finite_v_pub) * 0.08;
    y_lo_pub = max(0, min(finite_v_pub) - ypad_pub);
    y_hi_pub = max(finite_v_pub) + ypad_pub;
else
    y_lo_pub = 0;  y_hi_pub = 1;
end
for blk = 1:3
    ylim(ax_pub(blk), [y_lo_pub, y_hi_pub]);
end

% Auto y scale bar: round number covering ~20–40% of y range
if isnan(pub_y_scale)
    y_rng_pub   = y_hi_pub - y_lo_pub;
    y_sb_opts   = [0.5, 1, 2, 5, 10, 15, 20, 25, 50];
    pub_y_scale = y_sb_opts(end);
    for ys = y_sb_opts
        if ys / y_rng_pub >= 0.18 && ys / y_rng_pub <= 0.45
            pub_y_scale = ys;
            break;
        end
    end
end

% Strip all axis decorations
for blk = 1:3
    axis(ax_pub(blk), 'off');
end

% Draw L-shaped scale bars on last subplot (data coordinates)
axes(ax_pub(3));
hold(ax_pub(3), 'on');

x_span   = evtWin(2) - evtWin(1);
y_span   = y_hi_pub  - y_lo_pub;
sb_xpad  = x_span * 0.04;
sb_ypad  = y_span * 0.05;

x_sb_r = evtWin(2) - sb_xpad;           % right end of horizontal bar
x_sb_l = x_sb_r    - pub_x_scale;       % left  end
y_sb_b = y_lo_pub  + sb_ypad;           % bottom of vertical bar
y_sb_t = y_sb_b    + pub_y_scale;       % top

line([x_sb_l, x_sb_r], [y_sb_b, y_sb_b], 'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');
line([x_sb_r, x_sb_r], [y_sb_b, y_sb_t], 'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');

text(x_sb_l + pub_x_scale/2, y_sb_b - sb_ypad*1.2, ...
    sprintf('%g s', pub_x_scale), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 9, 'Clipping', 'off');
text(x_sb_r + sb_xpad*0.7, y_sb_b + pub_y_scale/2, ...
    sprintf('%g cm/s', pub_y_scale), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 9, 'Clipping', 'off');

hold(ax_pub(3), 'off');

sgtitle('Large reward — D1 — velocity aligned to choice  [pub]', 'FontSize', 11);

%% ---- Publication figure: Large reward D1, rows = sex, cols = block ----
% 2 rows × 3 columns. No axes. Scale bars. X-axis clipped to ±4 s.
% Row 1 (top)    : Male   — D1 per block
% Row 2 (bottom) : Female — D1 per block

pub2_fig_w  = 250;
pub2_fig_h  = 310;
pub2_evtWin = [-4, 4];   % display window (s around choice)

% Indices into ts that fall within the display window
pub2_ts_mask = ts >= pub2_evtWin(1) & ts <= pub2_evtWin(2);
pub2_ts      = ts(pub2_ts_mask);

figure;
set(gcf, 'Position', [100, 300, pub2_fig_w, pub2_fig_h]);

ax_pub2    = gobjects(2, 3);
all_v_pub2 = [];

% Layout: row 1 = Male, row 2 = Female
sex_masks_pub2 = {m_mask,   f_mask  };
sex_cols_pub2  = {col_m_d1, col_f_d1};   % D1 colour per row

for row = 1:2
    smask = sex_masks_pub2{row};

    for blk = 1:3
        ax_pub2(row, blk) = subplot(2, 3, (row-1)*3 + blk);

        vm = vel_mean{1, 1, blk};            % D1, Large reward
        d  = vm(smask, pub2_ts_mask);        % crop to display window
        n  = max(sum(~isnan(d(:,1))), 1);
        mn = nanmean(d, 1);
        se = nanstd(d, 0, 1) ./ sqrt(n);

        axes(ax_pub2(row, blk)); hold on; %#ok<LAXES>
        if any(~isnan(mn))
            shadedErrorBar(pub2_ts, mn, se, ...
                'lineProps', {'Color', sex_cols_pub2{row}, 'LineWidth', 1.8}, 'transparent', 1);
        end
        all_v_pub2 = [all_v_pub2, mn+se, mn-se]; %#ok<AGROW>

        xline(0, 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        yline(0, 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        xlim(pub2_evtWin);
        hold off;
    end
end

% Shared y limits
finite_v_pub2 = all_v_pub2(isfinite(all_v_pub2));
if ~isempty(finite_v_pub2)
    ypad_pub2 = range(finite_v_pub2) * 0.08;
    y_lo_pub2 = min(finite_v_pub2) - ypad_pub2;
    y_hi_pub2 = max(finite_v_pub2) + ypad_pub2;
else
    y_lo_pub2 = 0;  y_hi_pub2 = 1;
end
for row = 1:2
    for blk = 1:3
        ylim(ax_pub2(row, blk), [y_lo_pub2, y_hi_pub2]);
    end
end

% Auto y scale bar (~20–40% of y range)
y_rng_pub2   = y_hi_pub2 - y_lo_pub2;
pub2_y_scale = 50;
for ys = [0.5, 1, 2, 5, 10, 15, 20, 25, 50]
    if ys / y_rng_pub2 >= 0.18 && ys / y_rng_pub2 <= 0.45
        pub2_y_scale = ys;
        break;
    end
end

% Strip all axis decorations
for row = 1:2
    for blk = 1:3
        axis(ax_pub2(row, blk), 'off');
    end
end

% L-shaped scale bars on bottom-right subplot
axes(ax_pub2(2, 3));
hold(ax_pub2(2, 3), 'on');

x_span2  = pub2_evtWin(2) - pub2_evtWin(1);
y_span2  = y_hi_pub2 - y_lo_pub2;
sb_xpad2 = x_span2 * 0.04;
sb_ypad2 = y_span2 * 0.05;

x_sb_r2 = pub2_evtWin(2) - sb_xpad2;
x_sb_l2 = x_sb_r2        - pub_x_scale;
y_sb_b2 = y_lo_pub2      + sb_ypad2;
y_sb_t2 = y_sb_b2        + pub2_y_scale;

line([x_sb_l2, x_sb_r2], [y_sb_b2, y_sb_b2], 'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');
line([x_sb_r2, x_sb_r2], [y_sb_b2, y_sb_t2], 'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');

text(x_sb_l2 + pub_x_scale/2, y_sb_b2 - sb_ypad2*1.2, ...
    sprintf('%g s', pub_x_scale), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 9, 'Clipping', 'off');
text(x_sb_r2 + sb_xpad2*0.7, y_sb_b2 + pub2_y_scale/2, ...
    sprintf('%g cm/s', pub2_y_scale), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 9, 'Clipping', 'off');

hold(ax_pub2(2, 3), 'off');

sgtitle('Large reward — D1 by sex  [pub]', 'FontSize', 11);

%% ---- AUC: mean velocity in pre- and post-choice windows ----
pre_win  = [-1,  0];   % s relative to choice
post_win = [ 0,  1];

pre_idx  = ts >= pre_win(1)  & ts <= pre_win(2);
post_idx = ts >= post_win(1) & ts <= post_win(2);

% auc_pre/post{s_idx, ri, blk}(n_animals × 1) = mean velocity in window
auc_pre  = cell(2, 2, 3);
auc_post = cell(2, 2, 3);
for ss = 1:2; for rr = 1:2; for bb = 1:3
    auc_pre{ss,rr,bb}  = nanmean(vel_mean{ss,rr,bb}(:, pre_idx),  2);
    auc_post{ss,rr,bb} = nanmean(vel_mean{ss,rr,bb}(:, post_idx), 2);
end; end; end

% Bar plot colours (F-D1, M-D1, F-D2, M-D2 — matches behav figs order)
bar_colors   = {col_f_d1, col_m_d1, col_f_d2, col_m_d2};
sc_colors    = cellfun(@(c) min(c * 0.75, 1), bar_colors, 'UniformOutput', false);
sc_markers   = {'s', 'o', 's', 'o'};
bar_labels   = {'Female D1','Male D1','Female D2','Male D2'};
nbars        = 4;
ngroups      = 3;
groupwidth   = min(0.8, nbars / (nbars + 1.5));
bar_x_offsets = zeros(ngroups, nbars);
for b = 1:nbars
    bar_x_offsets(:, b) = (1:ngroups)' - groupwidth/2 + (2*b-1)*groupwidth/(2*nbars);
end

% Global y limits — derived from every individual animal value across
% both reward types, both periods, all sessions and blocks.
% This ensures all 4 AUC subplots (2 figures × 2 subplots) share identical axes.
global_auc_pts = [];
for ss = 1:2; for rr = 1:2; for bb = 1:3
    vp = auc_pre{ss,rr,bb};   global_auc_pts = [global_auc_pts; vp(isfinite(vp))]; %#ok<AGROW>
    vq = auc_post{ss,rr,bb};  global_auc_pts = [global_auc_pts; vq(isfinite(vq))]; %#ok<AGROW>
end; end; end

if ~isempty(global_auc_pts)
    ypad_g          = range(global_auc_pts) * 0.06;
    global_ylim_auc = [min(global_auc_pts) - ypad_g, ...
                        max(global_auc_pts) + ypad_g];
else
    global_ylim_auc = [0, 1];
end

% Compute y-ticks anchored to 0 so it is always present.
% Pick the smallest step that gives 3–7 ticks across the range.
auc_step_opts = [0.5, 1, 2, 5, 10, 15, 20, 25, 50];
auc_yt_step   = auc_step_opts(end);
for ys = auc_step_opts
    if diff(global_ylim_auc) / ys >= 3 && diff(global_ylim_auc) / ys <= 7
        auc_yt_step = ys;
        break;
    end
end
% Always start axis at 0 so it appears as a tick reference
global_ylim_auc(1) = 0;
auc_yt_lo         = 0;
auc_yt_hi         = floor(global_ylim_auc(2) / auc_yt_step) * auc_yt_step;
global_auc_yticks  = auc_yt_lo : auc_yt_step : auc_yt_hi;

% One figure per reward type; 2 subplots (pre | post)
auc_fig_w = 560;   auc_fig_h = 380;

for ri = 1:2
    figure;
    set(gcf, 'Position', [100, 100, auc_fig_w, auc_fig_h]);

    auc_data_sets = {auc_pre, auc_post};
    period_labels = {sprintf('Pre-choice AUC  [%.0f to 0 s]', pre_win(1)), ...
                     sprintf('Post-choice AUC  [0 to %.0f s]', post_win(2))};

    all_auc_vals = [];
    ax_auc = gobjects(1, 2);

    for pp = 1:2
        ax_auc(pp) = subplot(1, 2, pp);
        hold on;

        % Build (n_animals × 3) matrices for each of the 4 bar groups
        % Order: F-D1, M-D1, F-D2, M-D2
        grp_data = { ...
            cell2mat(cellfun(@(c) c(f_mask), {auc_data_sets{pp}{1,ri,1}, auc_data_sets{pp}{1,ri,2}, auc_data_sets{pp}{1,ri,3}}, 'UniformOutput', false)), ...
            cell2mat(cellfun(@(c) c(m_mask), {auc_data_sets{pp}{1,ri,1}, auc_data_sets{pp}{1,ri,2}, auc_data_sets{pp}{1,ri,3}}, 'UniformOutput', false)), ...
            cell2mat(cellfun(@(c) c(f_mask), {auc_data_sets{pp}{2,ri,1}, auc_data_sets{pp}{2,ri,2}, auc_data_sets{pp}{2,ri,3}}, 'UniformOutput', false)), ...
            cell2mat(cellfun(@(c) c(m_mask), {auc_data_sets{pp}{2,ri,1}, auc_data_sets{pp}{2,ri,2}, auc_data_sets{pp}{2,ri,3}}, 'UniformOutput', false))};

        means_mat = zeros(ngroups, nbars);
        sems_mat  = zeros(ngroups, nbars);
        for g = 1:nbars
            d = grp_data{g};   % animals × blocks
            n = max(sum(~isnan(d), 1), 1);
            means_mat(:, g) = nanmean(d, 1)';
            sems_mat(:, g)  = nanstd(d, 0, 1)' ./ sqrt(n)';
        end

        bar_h = bar(1:ngroups, means_mat, 'grouped', 'BarWidth', 0.85);
        for g = 1:nbars
            bar_h(g).FaceColor = bar_colors{g};
            bar_h(g).EdgeColor = 'none';
        end

        for g = 1:nbars
            d = grp_data{g};
            for blk = 1:ngroups
                xc   = bar_x_offsets(blk, g);
                errorbar(xc, means_mat(blk, g), sems_mat(blk, g), ...
                    'k', 'LineWidth', 1.2, 'CapSize', 5, 'LineStyle', 'none');
                vals = d(~isnan(d(:, blk)), blk);
                if ~isempty(vals)
                    jit = (rand(numel(vals), 1) - 0.5) * groupwidth / (nbars * 3);
                    scatter(xc + jit, vals, 26, sc_colors{g}, sc_markers{g}, ...
                        'filled', 'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', 'none');
                end
            end
            all_auc_vals = [all_auc_vals; means_mat(:,g)+sems_mat(:,g); means_mat(:,g)-sems_mat(:,g)]; %#ok<AGROW>
        end

        set(gca, 'XTick', 1:3, 'XTickLabel', {'Block 1','Block 2','Block 3'}, 'FontSize', 11);
        ylabel('Mean velocity (cm/s)', 'FontSize', 11);
        title(period_labels{pp}, 'FontSize', 11);
        if pp == 2
            legend(bar_h, bar_labels, 'Location', 'best', 'FontSize', 8);
        end
        box off;
        hold off;
    end

    % Apply global ylim and explicit ticks (same across all AUC figures)
    set(ax_auc(1), 'YLim', global_ylim_auc, 'YTick', global_auc_yticks);
    set(ax_auc(2), 'YLim', global_ylim_auc, 'YTick', global_auc_yticks);

    sgtitle(sprintf('%s reward — pre/post-choice mean velocity', reward_names{ri}), 'FontSize', 12);
end

%% ---- AUC statistics: 3-way mixed ANOVA (Sex × Day × Block) ----
auc_sources = {auc_pre, auc_post};
period_strs = {sprintf('Pre-choice AUC [%g to 0 s]', pre_win(1)), ...
               sprintf('Post-choice AUC [0 to %g s]', post_win(2))};

for ri = 1:2
    for pp = 1:2
        src = auc_sources{pp};
        f_d1 = [src{1,ri,1}(f_mask), src{1,ri,2}(f_mask), src{1,ri,3}(f_mask)];
        f_d2 = [src{2,ri,1}(f_mask), src{2,ri,2}(f_mask), src{2,ri,3}(f_mask)];
        m_d1 = [src{1,ri,1}(m_mask), src{1,ri,2}(m_mask), src{1,ri,3}(m_mask)];
        m_d2 = [src{2,ri,1}(m_mask), src{2,ri,2}(m_mask), src{2,ri,3}(m_mask)];
        run_3way_anova_vel(f_d1, f_d2, m_d1, m_d2, ...
            sprintf('%s — %s reward', period_strs{pp}, reward_names{ri}));
    end
end

%% ---- Shock-trial velocity: one figure per RDT day, all blocks pooled ----
% Extract per-animal mean velocity on shock trials only (all blocks combined)
vel_shk = cell(1, 2);   % {session}(animal, time)
for s_idx = 1:2
    vel_shk{s_idx} = nan(n_animals, n_ts);
end

for s_idx = 1:2
    session_to_analyze = sessions{s_idx};

    for ii = 1:n_animals
        currentanimal = ChrimsonR_IDs{ii};

        if ~any(strcmp(valid_sessions{ii}, session_to_analyze)); continue; end
        if ~bw_include(ii); continue; end
        if any(strcmp(vel_exclusions(:,1), session_to_analyze) & ...
               strcmp(vel_exclusions(:,2), currentanimal))
            continue;
        end
        if ~isfield(final_SLEAP, currentanimal) || ...
           ~isfield(final_SLEAP.(currentanimal), session_to_analyze); continue; end
        if ~isfield(final_behavior, currentanimal) || ...
           ~isfield(final_behavior.(currentanimal), session_to_analyze); continue; end

        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        % Must have shock column and at least one shock trial
        if ~ismember('shock', BehavData.Properties.VariableNames); continue; end
        shk_mask = BehavData.shock == 1;
        if ~any(shk_mask); continue; end

        SLEAP_data = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data_raw;
        vel_s      = sgolayfilt(SLEAP_data.vel_cm_s(:), sg_ord, sg_len);
        idx_time   = SLEAP_data.idx_time(:);

        t_idxs    = find(shk_mask);
        trial_mat = nan(numel(t_idxs), n_ts);

        for ti = 1:numel(t_idxs)
            jj       = t_idxs(ti);
            choice_t = BehavData.choiceTime(jj);
            if isnan(choice_t) || isinf(choice_t); continue; end

            seg_mask = idx_time >= choice_t + evtWin(1) & ...
                       idx_time <= choice_t + evtWin(2);
            if sum(seg_mask) < 5; continue; end

            t_rel = idx_time(seg_mask) - choice_t;
            v_seg = vel_s(seg_mask);
            trial_mat(ti, :) = interp1(t_rel, v_seg, ts, 'linear', NaN);
        end

        vel_shk{s_idx}(ii, :) = nanmean(trial_mat, 1);
    end
end

% ---- Plot: one figure per day ----
day_labels = {'RDT D1', 'RDT D2'};
shk_fig_colors = {col_f_d1, col_f_d2; col_m_d1, col_m_d2};  % (sex, day)

for s_idx = 1:2
    figure;
    set(gcf, 'Position', [100, 100, 380, 380]);
    hold on;

    f_data = vel_shk{s_idx}(f_mask, :);
    m_data = vel_shk{s_idx}(m_mask, :);

    nf = max(sum(~isnan(f_data(:, 1))), 1);
    nm = max(sum(~isnan(m_data(:, 1))), 1);

    f_mean = nanmean(f_data, 1);
    m_mean = nanmean(m_data, 1);
    f_sem  = nanstd(f_data, 0, 1) ./ sqrt(nf);
    m_sem  = nanstd(m_data, 0, 1) ./ sqrt(nm);

    col_f = shk_fig_colors{1, s_idx};
    col_m = shk_fig_colors{2, s_idx};

    if any(~isnan(f_mean))
        shadedErrorBar(ts, f_mean, f_sem, ...
            'lineProps', {'Color', col_f, 'LineWidth', 1.8}, ...
            'transparent', 1);
    end
    if any(~isnan(m_mean))
        shadedErrorBar(ts, m_mean, m_sem, ...
            'lineProps', {'Color', col_m, 'LineWidth', 1.8}, ...
            'transparent', 1);
    end

    xline(0, 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    % Legend
    lh_s(1) = plot(nan, nan, '-', 'Color', col_f, 'LineWidth', 2);
    lh_s(2) = plot(nan, nan, '-', 'Color', col_m, 'LineWidth', 2);
    legend(lh_s, {sprintf('Female (n=%d)', nf), sprintf('Male (n=%d)', nm)}, ...
        'Location', 'best', 'FontSize', 9);

    % Shared ylim from data
    all_shk_vals = [f_mean + f_sem, f_mean - f_sem, m_mean + m_sem, m_mean - m_sem];
    finite_shk   = all_shk_vals(isfinite(all_shk_vals));
    if ~isempty(finite_shk)
        ypad = range(finite_shk) * 0.06;
        ylim([max(0, min(finite_shk) - ypad),  max(finite_shk) + ypad]);
    end

    xlim(evtWin);
    set(gca, 'XTick', evtWin(1):4:evtWin(2));
    xlabel('Time from choice (s)', 'FontSize', 11);
    ylabel('Velocity (cm/s)',      'FontSize', 11);
    title(sprintf('Shock trials — %s', day_labels{s_idx}), 'FontSize', 12);
    box off;
    hold off;
end

%% ---- Shock-trial velocity: D1 + D2 combined per animal ----
% Average the per-session means (nanmean ignores days with no shock data)
vel_shk_combined = nanmean(cat(3, vel_shk{1}, vel_shk{2}), 3);   % n_animals × n_ts

figure;
set(gcf, 'Position', [100, 100, 380, 380]);
hold on;

f_data_c = vel_shk_combined(f_mask, :);
m_data_c = vel_shk_combined(m_mask, :);

nf_c = max(sum(~isnan(f_data_c(:, 1))), 1);
nm_c = max(sum(~isnan(m_data_c(:, 1))), 1);

f_mean_c = nanmean(f_data_c, 1);
m_mean_c = nanmean(m_data_c, 1);
f_sem_c  = nanstd(f_data_c, 0, 1) ./ sqrt(nf_c);
m_sem_c  = nanstd(m_data_c, 0, 1) ./ sqrt(nm_c);

% Use mid-range green/grey (between D1 and D2 shades)
col_f_c = (col_f_d1 + col_f_d2) / 2;
col_m_c = (col_m_d1 + col_m_d2) / 2;

if any(~isnan(f_mean_c))
    shadedErrorBar(ts, f_mean_c, f_sem_c, ...
        'lineProps', {'Color', col_f_c, 'LineWidth', 1.8}, ...
        'transparent', 1);
end
if any(~isnan(m_mean_c))
    shadedErrorBar(ts, m_mean_c, m_sem_c, ...
        'lineProps', {'Color', col_m_c, 'LineWidth', 1.8}, ...
        'transparent', 1);
end

xline(0, 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off');

lh_c(1) = plot(nan, nan, '-', 'Color', col_f_c, 'LineWidth', 2);
lh_c(2) = plot(nan, nan, '-', 'Color', col_m_c, 'LineWidth', 2);
legend(lh_c, {sprintf('Female (n=%d)', nf_c), sprintf('Male (n=%d)', nm_c)}, ...
    'Location', 'best', 'FontSize', 9);

all_c_vals = [f_mean_c + f_sem_c, f_mean_c - f_sem_c, ...
              m_mean_c + m_sem_c, m_mean_c - m_sem_c];
finite_c   = all_c_vals(isfinite(all_c_vals));
if ~isempty(finite_c)
    ypad_c = range(finite_c) * 0.06;
    ylim([max(0, min(finite_c) - ypad_c),  max(finite_c) + ypad_c]);
end

xlim(evtWin);
set(gca, 'XTick', evtWin(1):4:evtWin(2));
xlabel('Time from choice (s)', 'FontSize', 11);
ylabel('Velocity (cm/s)',      'FontSize', 11);
title('Shock trials — D1 + D2 combined', 'FontSize', 12);
box off;
hold off;


%% ================================================================
%  Local functions
%  ================================================================

function run_3way_anova_vel(f_d1, f_d2, m_d1, m_d2, measure_name)
% 3-way mixed ANOVA: Sex (between) × Day × Block (within).
% Inputs: each matrix is (n_animals × 3 blocks), one row per animal.

nF = size(f_d1, 1);
nM = size(m_d1, 1);

all_wide   = [f_d1, f_d2; m_d1, m_d2];
sex_labels = categorical([repmat({'Female'}, nF, 1); repmat({'Male'}, nM, 1)]);

% Apply listwise deletion ourselves so fitrm sees both Sex levels explicitly.
% If MATLAB does it internally the Sex column can lose one level, silently
% dropping the between-subjects term from the model.
complete     = all(~isnan(all_wide), 2);
n_complete   = sum(complete);
n_excluded   = (nF + nM) - n_complete;
nF_c         = sum(complete(1:nF));
nM_c         = sum(complete(nF+1:end));

all_wide_c   = all_wide(complete, :);
sex_labels_c = sex_labels(complete);

tbl     = array2table(all_wide_c, 'VariableNames', {'D1B1','D1B2','D1B3','D2B1','D2B2','D2B3'});
tbl.Sex = sex_labels_c;
within  = table( ...
    categorical({'D1';'D1';'D1';'D2';'D2';'D2'}), ...
    categorical({'B1';'B2';'B3';'B1';'B2';'B3'}), ...
    'VariableNames', {'Day','Block'});

fprintf('\n========== 3-Way Mixed ANOVA: %s ==========\n', measure_name);
fprintf('N = %d subjects (%d female, %d male)', n_complete, nF_c, nM_c);
if n_excluded > 0
    fprintf('  [%d excluded — missing session data]', n_excluded);
end
fprintf('\n');

if nF_c == 0 || nM_c == 0
    fprintf('  [Only one sex has complete data — between-subjects test unavailable]\n');
    return;
end

try
    rm = fitrm(tbl, 'D1B1,D1B2,D1B3,D2B1,D2B2,D2B3 ~ Sex', 'WithinDesign', within);
catch ME
    fprintf('  [fitrm failed: %s]\n', ME.message);
    return;
end

try
    bs_tbl = rm.anova();
catch ME
    fprintf('  [rm.anova() failed: %s]\n', ME.message);
    bs_tbl = table();
end

try
    ws_tbl = ranova(rm, 'WithinModel', 'Day*Block');
catch ME
    fprintf('  [ranova failed: %s]\n', ME.message);
    return;
end

try
    ws_rows = ws_tbl.Properties.RowNames;
    bs_rows = bs_tbl.Properties.RowNames;
    bs_vars = bs_tbl.Properties.VariableNames;

    fprintf('\n  Between-subjects:\n');

    % rm.anova() has two layouts depending on MATLAB version:
    %   Format A (older): source names are row names; no Between/Within columns
    %   Format B (newer): no row names; source stored in a 'Between' column
    err_labels_bs = {'Error(Sex)', 'Error', 'Residuals', 'Residual'};
    sex_idx = [];  err_idx = [];
    if ~isempty(bs_rows)
        % Format A — row-name based lookup
        sex_idx = find(strcmp(bs_rows, 'Sex'));
        for lbl = err_labels_bs
            idx = find(strcmp(bs_rows, lbl{1}));
            if ~isempty(idx); err_idx = idx; break; end
        end
    elseif ismember('Between', bs_vars)
        % Format B — 'Between' column identifies the source
        btw_vals = cellstr(bs_tbl.Between);
        sex_idx  = find(strcmp(btw_vals, 'Sex'));
        for lbl = err_labels_bs
            idx = find(strcmp(btw_vals, lbl{1}));
            if ~isempty(idx); err_idx = idx; break; end
        end
    end

    p_sex = NaN;
    if isempty(sex_idx)
        fprintf('    [Sex row not found — bs_tbl rows: {%s}  cols: {%s}]\n', ...
            strjoin(bs_rows, ', '), strjoin(bs_vars, ', '));
    else
        F_val = NaN;  p_val = NaN;  df1 = NaN;  df2 = NaN;
        if ismember('F',      bs_vars); F_val = bs_tbl.F(sex_idx);      end
        if ismember('pValue', bs_vars); p_val = bs_tbl.pValue(sex_idx); end
        if ismember('DF',     bs_vars); df1   = bs_tbl.DF(sex_idx);     end
        if ~isempty(err_idx) && ismember('DF', bs_vars)
            df2 = bs_tbl.DF(err_idx);
        end
        if ~isnan(F_val) && ~isnan(p_val)
            print_effect_v('Sex', F_val, df1, df2, p_val);
            p_sex = p_val;
        else
            fprintf('    [Sex F/p are NaN — F=%s  p=%s]\n', mat2str(F_val), mat2str(p_val));
        end
    end

    fprintf('  Within-subjects and interactions:\n');
    effects_ws = { ...
        '(Intercept):Day',       'Day',             'Error(Day)'; ...
        'Sex:Day',               'Sex x Day',        'Error(Day)'; ...
        '(Intercept):Block',     'Block',            'Error(Block)'; ...
        'Sex:Block',             'Sex x Block',      'Error(Block)'; ...
        '(Intercept):Day:Block', 'Day x Block',      'Error(Day:Block)'; ...
        'Sex:Day:Block',         'Sex x Day x Block','Error(Day:Block)'};

    for e = 1:size(effects_ws, 1)
        eff_idx = find(strcmp(ws_rows, effects_ws{e,1}));
        err_idx_w = find(strcmp(ws_rows, effects_ws{e,3}));
        if ~isempty(eff_idx) && ~isempty(err_idx_w)
            print_effect_v(effects_ws{e,2}, ws_tbl.F(eff_idx), ...
                ws_tbl.DF(eff_idx), ws_tbl.DF(err_idx_w), ws_tbl.pValue(eff_idx));
        end
    end

    p_3way      = get_p_v(ws_tbl, 'Sex:Day:Block');
    p_sex_block = get_p_v(ws_tbl, 'Sex:Block');
    p_sex_day   = get_p_v(ws_tbl, 'Sex:Day');
    p_day_block = get_p_v(ws_tbl, '(Intercept):Day:Block');
    p_block     = get_p_v(ws_tbl, '(Intercept):Block');
    p_day       = get_p_v(ws_tbl, '(Intercept):Day');

    ran_posthoc = false;
    if p_3way < 0.05
        fprintf('\n  Significant 3-way — post-hocs (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp_v(mc);
        mc = multcompare(rm, 'Day', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp_v(mc);
        ran_posthoc = true;
    end
    if p_sex_block < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex x Block (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp_v(mc);
        ran_posthoc = true;
    end
    if p_sex_day < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex x Day (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
        print_multcomp_v(mc);
        ran_posthoc = true;
    end
    if p_day_block < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Day x Block (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Day', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp_v(mc);
        ran_posthoc = true;
    end
    if p_block < 0.05 && p_sex_block >= 0.05 && p_day_block >= 0.05 && p_3way >= 0.05
        fprintf('\n  Significant Block main effect (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp_v(mc);
    end
    if p_day < 0.05 && p_sex_day >= 0.05 && p_day_block >= 0.05 && p_3way >= 0.05
        fprintf('\n  Significant Day main effect (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Day', 'ComparisonType', 'tukey-kramer');
        print_multcomp_v(mc);
    end
    if ~isnan(p_sex) && p_sex < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex main effect (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'ComparisonType', 'tukey-kramer');
        print_multcomp_v(mc);
    end

catch ME
    fprintf('  [ANOVA failed: %s]\n', ME.message);
end
end


function p = get_p_v(tbl, row_name)
    try
        p = tbl{row_name, 'pValue'};
    catch
        idx = find(strcmp(tbl.Properties.RowNames, row_name));
        if ~isempty(idx); p = tbl.pValue(idx(1)); else; p = NaN; end
    end
end

function print_effect_v(label, F, df1, df2, p)
    fprintf('    %-26s F(%d,%d) = %7.3f,  p = %.4f%s\n', label, df1, df2, F, p, sig_stars_v(p));
end

function print_multcomp_v(mc)
    if ~istable(mc)
        for r = 1:size(mc, 1)
            p_val = mc(r, 6);
            fprintf('      Group %d vs %d: diff=%7.3f, 95%%CI[%6.3f,%6.3f], p=%.4f%s\n', ...
                mc(r,1), mc(r,2), mc(r,4), mc(r,3), mc(r,5), p_val, sig_stars_v(p_val));
        end
        return;
    end
    if ~ismember('pValue', mc.Properties.VariableNames); disp(mc); return; end

    % Count leading categorical / cell columns (before the first numeric column).
    % multcompare(..., 'By', X) puts X in col 1, compared-factor levels in col 2 & 3.
    % multcompare without 'By' puts compared-factor levels in col 1 & 2.
    n_cat = 0;
    for c = 1:width(mc)
        if iscategorical(mc{:,c}) || iscell(mc{:,c})
            n_cat = n_cat + 1;
        else
            break;
        end
    end

    seen_pairs = {};
    for r = 1:height(mc)
        p_val = mc.pValue(r);
        if n_cat >= 3
            % 'By' form: col1 = By-level, col2/col3 = compared factor levels
            by_lv  = char(mc{r, 1});
            lv1    = char(mc{r, 2});
            lv2    = char(mc{r, 3});
            by_str = [' [' by_lv ']'];
        elseif n_cat >= 2
            % No 'By': col1/col2 = compared factor levels
            lv1    = char(mc{r, 1});
            lv2    = char(mc{r, 2});
            by_lv  = '';
            by_str = '';
        else
            lv1 = char(mc{r,1}); lv2 = '?'; by_lv = ''; by_str = '';
        end

        % Suppress duplicate rows — MATLAB returns both A-vs-B and B-vs-A
        sorted_pair = sort({lv1, lv2});
        pair_key    = [sorted_pair{1} '|' sorted_pair{2} '|' by_lv];
        if any(strcmp(seen_pairs, pair_key)); continue; end
        seen_pairs{end+1} = pair_key; %#ok<AGROW>

        fprintf('      %-42s p = %.4f%s\n', [lv1 ' vs ' lv2 by_str], p_val, sig_stars_v(p_val));
    end
end

function s = sig_stars_v(p)
    if p < 0.001;    s = ' ***';
    elseif p < 0.01; s = ' **';
    elseif p < 0.05; s = ' *';
    else;            s = '';
    end
end
