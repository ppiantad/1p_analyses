% locomotor_velocity_heatmap_choice_aligned.m
%
% Heatmap of session-normalised velocity — one row per valid trial,
% aligned to the animal's choice time.
%
%   X-axis      : time relative to choice (s); t=0 = choice
%   Y-axis      : trial number, chronological (trial 1 at top)
%   Colour      : Savitzky-Golay-smoothed velocity, session min–max
%                 normalised, rendered with a Crameri Scientific Colour Map
%   Left strip  : reward type per trial (Large = dark / Small = light)
%   Overlaid markers:
%       white  triangle (▲) — trial onset  (always left of t=0)
%       black  diamond  (◆) — reward collection
%       white dotted line   — choice time (t=0)
%
% "Proper" trials = Large or Small reward choices only.
%   Excluded: omissions, blank touches, shock trials.
%
% Requires workspace: final_SLEAP, final_behavior

%% ---- Settings ----
select_mouse       = 'RDT_M_4';
session_to_analyze = 'RDT_D1';

% Crameri Scientific Colour Map
%   Any subfolder name under scm_root: batlow, hawaii, lajolla, tokyo,
%   devon, oslo, nuuk, davos, bilbao, acton, cork, …
scm_root  = 'I:\MATLAB\Packages\ScientificColourMaps7';
cmap_name = 'batlow';

% Velocity colour-axis limits (session-normalised 0–1 units).
% Set to [NaN, NaN] to derive automatically from the 2nd–98th percentile.
vel_clim = [0, 0.65];

% Time window around choice (s)
evtWin = [-8, 8];

% Camera frame rate (Hz)
fs_cam = 30;

% Savitzky-Golay smoothing parameters
sg_ord = 9;
sg_len = 21;

% Event marker colours
col_onset_fill    = [1.00, 1.00, 1.00];   % white fill  — trial onset triangle
col_onset_edge    = [0.25, 0.25, 0.25];   % dark edge
col_collect_large = [0.15, 0.35, 0.80];   % blue  — Large reward collection diamond
col_collect_small = [0.85, 0.15, 0.15];   % red   — Small reward collection diamond
col_collect_edge  = [0.92, 0.92, 0.92];   % light edge on collection markers

% Figure size (px)
fig_width = 720;

%% ---- Load Crameri colourmap ----
cmap_file = fullfile(scm_root, cmap_name, [cmap_name '.mat']);
if ~exist(cmap_file, 'file')
    error('Colourmap file not found:\n  %s\nCheck scm_root and cmap_name.', cmap_file);
end
cmap_loaded = load(cmap_file);
cmap_rgb    = cmap_loaded.(cmap_name);   % 256×3 RGB matrix

%% ---- Load session data ----
if ~isfield(final_SLEAP, select_mouse) || ...
   ~isfield(final_SLEAP.(select_mouse), session_to_analyze)
    error('No SLEAP data for %s / %s.', select_mouse, session_to_analyze);
end
if ~isfield(final_behavior, select_mouse) || ...
   ~isfield(final_behavior.(select_mouse), session_to_analyze)
    error('No behaviour data for %s / %s.', select_mouse, session_to_analyze);
end

SLEAP_data = final_SLEAP.(select_mouse).(session_to_analyze).SLEAP_data_raw;
BehavData  = final_behavior.(select_mouse).(session_to_analyze).uv.BehavData;

%% ---- Smooth and session-normalise velocity ----
vel_s    = sgolayfilt(SLEAP_data.vel_cm_s(:), sg_ord, sg_len);
idx_time = SLEAP_data.idx_time(:);

vel_min  = min(vel_s);
vel_max  = max(vel_s);
vel_norm = (vel_s - vel_min) / max(vel_max - vel_min, eps);

%% ---- Valid trial mask ----
valid_trial = ismember(BehavData.bigSmall, [1.2, 0.3]);
if ismember('omissionALL', BehavData.Properties.VariableNames)
    valid_trial = valid_trial & BehavData.omissionALL ~= 1;
end
if ismember('Blank_Touch', BehavData.Properties.VariableNames)
    valid_trial = valid_trial & BehavData.Blank_Touch ~= 1;
end
if ismember('shock', BehavData.Properties.VariableNames)
    valid_trial = valid_trial & BehavData.shock ~= 1;
end

t_idxs = find(valid_trial);
n_t    = numel(t_idxs);
fprintf('%d valid trials found (%d total).\n', n_t, height(BehavData));

%% ---- Build heatmap matrix and collect event times ----
dt   = 1 / fs_cam;
ts   = evtWin(1) : dt : evtWin(2) - dt;
n_ts = numel(ts);

hmap          = nan(n_t, n_ts);
t_onset_rel   = nan(n_t, 1);    % trial onset relative to choice (s)
t_collect_rel = nan(n_t, 1);    % reward collection relative to choice (s)
is_large      = false(n_t, 1);  % true = Large reward trial
trial_block   = nan(n_t, 1);   % block number (1/2/3) for each trial row

for ti = 1:n_t
    jj       = t_idxs(ti);
    choice_t = BehavData.choiceTime(jj);
    onset_t  = BehavData.stTime(jj);
    coll_t   = BehavData.collectionTime(jj);

    if isnan(choice_t) || isinf(choice_t); continue; end

    seg_mask = idx_time >= choice_t + evtWin(1) & ...
               idx_time <= choice_t + evtWin(2);
    if sum(seg_mask) < 5; continue; end

    t_rel = idx_time(seg_mask) - choice_t;
    v_seg = vel_norm(seg_mask);
    hmap(ti, :) = interp1(t_rel, v_seg, ts, 'linear', NaN);

    if ~isnan(onset_t) && ~isinf(onset_t)
        t_onset_rel(ti) = onset_t - choice_t;
    end
    if ~isnan(coll_t) && ~isinf(coll_t)
        t_collect_rel(ti) = coll_t - choice_t;
    end
    is_large(ti)    = BehavData.bigSmall(jj) == 1.2;
    trial_block(ti) = BehavData.Block(jj);
end

%% ---- Auto colour limits if requested ----
if any(isnan(vel_clim))
    finite_v = hmap(isfinite(hmap));
    if ~isempty(finite_v)
        vel_clim = [prctile(finite_v, 2), prctile(finite_v, 98)];
    end
end

%% ---- Figure layout ----
fig_h = max(420, n_t * 7 + 150);

hfig = figure;
set(hfig, 'Position', [80, 80, fig_width, fig_h]);
colormap(hfig, cmap_rgb);   % set figure-level default (also covers the colorbar)

margin_l = 0.10;   % left margin
hmap_w   = 0.73;   % heatmap width
cb_gap   = 0.015;  % gap between heatmap and colorbar
cb_w     = 0.022;  % colorbar width
bot      = 0.10;
ht       = 0.82;

pos_main = [margin_l,                  bot, hmap_w, ht];
pos_cb   = [margin_l + hmap_w + cb_gap, bot, cb_w,  ht];

ax_main = axes('Position', pos_main);

%% ---- Draw main heatmap ----
axes(ax_main);
colormap(ax_main, cmap_rgb);
imagesc(ax_main, ts, 1:n_t, hmap);
set(ax_main, 'CLim', vel_clim, 'YDir', 'reverse', ...
    'Color', [0.78, 0.78, 0.78], ...   % background for NaN regions
    'YLim', [0.5, n_t + 0.5], 'YTick', get_yticks_hm(n_t), ...
    'XLim', evtWin, 'XTick', evtWin(1):2:evtWin(2));
xlabel(ax_main, 'Time from choice (s)', 'FontSize', 11);
ylabel(ax_main, 'Trial number', 'FontSize', 11);
hold(ax_main, 'on');

% Vertical reference line at choice (t = 0)
xline(ax_main, 0, 'w:', 'LineWidth', 1.4, 'HandleVisibility', 'off');

% Horizontal dashed lines at block boundaries (B1→B2 and B2→B3)
for blk = 1:2
    last_row = find(trial_block == blk, 1, 'last');
    if ~isempty(last_row)
        yline(ax_main, last_row + 0.5, 'w--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end

% Overlay event markers
for ti = 1:n_t
    to = t_onset_rel(ti);
    tc = t_collect_rel(ti);

    if ~isnan(to) && to >= evtWin(1) && to <= evtWin(2)
        scatter(ax_main, to, ti, 36, col_onset_fill, '^', 'filled', ...
            'MarkerEdgeColor', col_onset_edge, 'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end

    if ~isnan(tc) && tc >= evtWin(1) && tc <= evtWin(2)
        c_fill = col_collect_large * is_large(ti) + col_collect_small * ~is_large(ti);
        scatter(ax_main, tc, ti, 36, c_fill, 'd', 'filled', ...
            'MarkerEdgeColor', col_collect_edge, 'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end
end

% Legend (dummy handles)
hl(1) = scatter(ax_main, nan, nan, 36, col_onset_fill, '^', 'filled', ...
    'MarkerEdgeColor', col_onset_edge);
hl(2) = scatter(ax_main, nan, nan, 36, col_collect_large, 'd', 'filled', ...
    'MarkerEdgeColor', col_collect_edge);
hl(3) = scatter(ax_main, nan, nan, 36, col_collect_small, 'd', 'filled', ...
    'MarkerEdgeColor', col_collect_edge);
legend(ax_main, hl, {'Trial onset', 'Large reward collection', 'Small reward collection'}, ...
    'Location', 'northeast', 'FontSize', 8, 'Box', 'on');

box(ax_main, 'off');
hold(ax_main, 'off');

%% ---- Colourbar ----
% Create linked to ax_main first, then reposition — passing 'Position'
% inline can silently detach the bar from the axes-level colormap.
cb = colorbar(ax_main);
cb.Position          = pos_cb;
cb.Limits            = vel_clim;
cb.Label.String      = sprintf('Velocity — %s  (session min–max normalised)', cmap_name);
cb.Label.FontSize    = 9;
cb.FontSize          = 9;

%% ---- Title ----
n_large = sum(is_large);
n_small = sum(~is_large);
sgtitle(sprintf('%s  |  %s  —  velocity aligned to choice  (n=%d: %d large, %d small)', ...
    strrep(select_mouse, '_', '\_'), ...
    strrep(session_to_analyze, '_', '\_'), ...
    n_t, n_large, n_small), 'FontSize', 11);


%% ================================================================
%  Local functions
%  ================================================================

function tks = get_yticks_hm(n)
% Return a reasonable set of Y-tick values for n trials.
    steps = [1, 2, 5, 10, 20, 25, 50];
    step  = steps(end);
    for s = steps
        if ceil(n / s) <= 8
            step = s;
            break;
        end
    end
    tks = step : step : n;
end
