% generate_loco_trajectory_overlay_GoPro.m
%
% Overlays one representative locomotor trajectory per mouse on a shared figure.
%
% Trajectories are ARENA-NORMALISED:
%   - centred on the largest circle in shapeData
%   - scaled by that circle's radius  →  units of arena radii
%   - mice whose large reward port is on the RIGHT have X negated so that
%     all mice appear with large reward on the LEFT
%
% Overlaid per trial:
%   Filled square  — trial start
%   Open circle    — choice time position
%   Filled triangle — trial end (collection)
%
% Reward port rectangles (from shapeData, normalised):
%   Blue outline  — large reward screen (always plotted on the left)
%   Red  outline  — small reward screen (always plotted on the right)
%   Rectangles are taken from the first mouse that has valid shapeData.
%
% One randomly-selected trial per mouse per block is extracted.
% Two figures produced:
%   Figure 1 — Large reward  : Block 1 (left subplot) | Block 3 (right subplot)
%   Figure 2 — Small reward  : Block 1 (left subplot) | Block 3 (right subplot)
%
% Requires in workspace (run females_vs_males_groups_correct.m first):
%   ChrimsonR_IDs, large_screen_side, valid_sessions
% Requires: final_SLEAP, final_behavior

%% ---- Settings ----
session_to_analyze = 'RDT_D1';
line_alpha         = 0.75;    % trajectory line transparency [0, 1]
line_width         = 1.4;     % trajectory line width (pts)
marker_size        = 55;      % start / end marker size
choice_marker_size = 60;      % open circle at choice time

col_large_rect = [0.18, 0.36, 0.82];   % blue  — large reward port outline
col_small_rect = [0.82, 0.18, 0.18];   % red   — small reward port outline

%% ---- Per-animal colour map: females = greens, males = greys ----
n_animals = numel(ChrimsonR_IDs);
is_female = strcmp(ChrimsonR_treatment_groups, 'Female');
n_female  = sum(is_female);
n_male    = sum(~is_female);

% Green ramp: light green → dark green
green_light = [0.55, 0.88, 0.55];
green_dark  = [0.05, 0.42, 0.12];
green_ramp  = [linspace(green_light(1), green_dark(1), max(n_female,1))', ...
               linspace(green_light(2), green_dark(2), max(n_female,1))', ...
               linspace(green_light(3), green_dark(3), max(n_female,1))'];

% Grey ramp: light grey → dark grey
grey_light  = [0.78, 0.78, 0.78];
grey_dark   = [0.18, 0.18, 0.18];
grey_ramp   = [linspace(grey_light(1), grey_dark(1), max(n_male,1))', ...
               linspace(grey_light(2), grey_dark(2), max(n_male,1))', ...
               linspace(grey_light(3), grey_dark(3), max(n_male,1))'];

% Assemble cmap_all: one colour per animal in ChrimsonR_IDs order
cmap_all   = zeros(n_animals, 3);
f_count    = 0;
m_count    = 0;
for ii = 1:n_animals
    if is_female(ii)
        f_count = f_count + 1;
        cmap_all(ii, :) = green_ramp(f_count, :);
    else
        m_count = m_count + 1;
        cmap_all(ii, :) = grey_ramp(m_count, :);
    end
end

%% ---- Extract representative trials for every mouse ----
% Trajectories: 2 × n_frames (normalised XY)
% Choice positions: 2 × 1 (normalised XY at choice time)
traj_large_b1   = cell(n_animals, 1);
traj_large_b3   = cell(n_animals, 1);
traj_small_b1   = cell(n_animals, 1);
traj_small_b3   = cell(n_animals, 1);

choice_large_b1 = cell(n_animals, 1);
choice_large_b3 = cell(n_animals, 1);
choice_small_b1 = cell(n_animals, 1);
choice_small_b3 = cell(n_animals, 1);

% Normalised reward port rectangles — taken from first valid mouse.
% Each struct: .pos_norm [x y w h], .is_large (logical)
norm_squares = [];   % filled on first successful mouse

for ii = 1:n_animals
    currentanimal = ChrimsonR_IDs{ii};

    % --- Skip if session or data not available ---
    if ~any(strcmp(valid_sessions{ii}, session_to_analyze))
        fprintf('[%s] session %s not valid — skipping\n', currentanimal, session_to_analyze);
        continue;
    end
    if ~isfield(final_SLEAP, currentanimal) || ...
       ~isfield(final_SLEAP.(currentanimal), session_to_analyze)
        fprintf('[%s] SLEAP data missing — skipping\n', currentanimal);
        continue;
    end
    if ~isfield(final_behavior, currentanimal) || ...
       ~isfield(final_behavior.(currentanimal), session_to_analyze)
        fprintf('[%s] behavior data missing — skipping\n', currentanimal);
        continue;
    end

    SLEAP_data = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data_raw;
    shapeData  = final_SLEAP.(currentanimal).(session_to_analyze).shapeData;
    BehavData  = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

    % --- Find the main (largest) arena circle for normalisation ---
    arena_r  = 0;
    arena_cx = 0;
    arena_cy = 0;
    for k = 1:numel(shapeData)
        if strcmp(shapeData{k}.Type, 'Circle') && shapeData{k}.Radius > arena_r
            arena_r  = shapeData{k}.Radius;
            arena_cx = shapeData{k}.Center(1);
            arena_cy = shapeData{k}.Center(2);
        end
    end
    if arena_r == 0
        fprintf('[%s] No circle found in shapeData — skipping\n', currentanimal);
        continue;
    end

    % --- Mirror flag: negate X if large reward is on the right ---
    do_mirror = strcmp(large_screen_side{ii}, 'right');

    % --- Smooth XY and normalise to arena coordinates ---
    X_sm   = sgolayfilt(SLEAP_data.x_pix, 9, 33);
    Y_sm   = sgolayfilt(SLEAP_data.y_pix, 9, 33);
    X_norm = (X_sm - arena_cx) / arena_r;
    Y_norm = (Y_sm - arena_cy) / arena_r;
    if do_mirror
        X_norm = -X_norm;   % flip so large reward is always on the left
    end

    % --- Collect normalised reward-port rectangles from the first valid mouse ---
    if isempty(norm_squares)
        sq_list = struct('pos_norm', {}, 'is_large', {});
        for k = 1:numel(shapeData)
            if strcmp(shapeData{k}.Type, 'Square')
                sc = shapeData{k}.Center;  % [cx, cy] pixels
                ss = shapeData{k}.Size;    % [w, h]   pixels
                % Normalise centre
                sc_n = (sc - [arena_cx, arena_cy]) / arena_r;
                ss_n = ss / arena_r;
                if do_mirror
                    sc_n(1) = -sc_n(1);   % mirror x
                end
                % After mirroring, large reward is always on the left (x < 0)
                is_lg = sc_n(1) < 0;
                entry.pos_norm = [sc_n(1) - ss_n(1)/2, sc_n(2) - ss_n(2)/2, ss_n(1), ss_n(2)];
                entry.is_large = is_lg;
                sq_list(end+1) = entry; %#ok<AGROW>
            end
        end
        if ~isempty(sq_list)
            norm_squares = sq_list;
        end
    end

    % --- Extract one random trial per reward-type × block combination ---
    combos    = [1.2, 1;  1.2, 3;  0.3, 1;  0.3, 3];
    traj_res  = cell(4, 1);
    choice_res = cell(4, 1);
    found_str  = '';

    for c = 1:4
        rv  = combos(c, 1);
        blk = combos(c, 2);

        rows = find(BehavData.bigSmall == rv & BehavData.Block == blk);
        if isempty(rows); continue; end

        pick     = rows(randi(numel(rows)));
        t_on     = BehavData.stTime(pick);
        t_off    = BehavData.collectionTime(pick);
        t_choice = BehavData.choiceTime(pick);
        if isinf(t_off) || isnan(t_off); continue; end

        time_mask = SLEAP_data.idx_time > t_on & SLEAP_data.idx_time < t_off;
        if ~any(time_mask); continue; end

        traj_res{c} = [X_norm(time_mask)'; Y_norm(time_mask)'];   % 2 × n_frames

        % Choice position: nearest frame to choiceTime
        ci = interp1(SLEAP_data.idx_time, 1:numel(SLEAP_data.idx_time), ...
                     t_choice, 'nearest', 'extrap');
        ci = max(1, min(ci, numel(X_norm)));
        choice_res{c} = [X_norm(ci); Y_norm(ci)];   % 2 × 1

        label_map = {'L-B1','L-B3','S-B1','S-B3'};
        found_str = [found_str, label_map{c}, ' ']; %#ok<AGROW>
    end

    traj_large_b1{ii}   = traj_res{1};
    traj_large_b3{ii}   = traj_res{2};
    traj_small_b1{ii}   = traj_res{3};
    traj_small_b3{ii}   = traj_res{4};

    choice_large_b1{ii} = choice_res{1};
    choice_large_b3{ii} = choice_res{2};
    choice_small_b1{ii} = choice_res{3};
    choice_small_b3{ii} = choice_res{4};

    mirror_str = '';
    if do_mirror; mirror_str = ' [mirrored]'; end
    fprintf('[%s] extracted: %s%s\n', currentanimal, found_str, mirror_str);
end

%% ---- Compute global axis limits across ALL trajectories and ALL figures ----
all_traj = [traj_large_b1; traj_large_b3; traj_small_b1; traj_small_b3];
all_x = [];
all_y = [];
for tt = 1:numel(all_traj)
    tr = all_traj{tt};
    if ~isempty(tr)
        all_x = [all_x, tr(1, :)]; %#ok<AGROW>
        all_y = [all_y, tr(2, :)]; %#ok<AGROW>
    end
end
% Also include rectangle corners so the reward port squares are never clipped
if ~isempty(norm_squares)
    for k = 1:numel(norm_squares)
        p = norm_squares(k).pos_norm;  % [x y w h]
        all_x = [all_x, p(1), p(1)+p(3)]; %#ok<AGROW>
        all_y = [all_y, p(2), p(2)+p(4)]; %#ok<AGROW>
    end
end

axis_pad = 0.08;   % fractional padding added to each side of the data range
if ~isempty(all_x)
    x_range = max(all_x) - min(all_x);
    y_range = max(all_y) - min(all_y);
    global_xlim = [min(all_x) - axis_pad*x_range,  max(all_x) + axis_pad*x_range];
    global_ylim = [min(all_y) - axis_pad*y_range,  max(all_y) + axis_pad*y_range];
    % Square and symmetric — use the larger span for both axes
    half_span   = max(range(global_xlim), range(global_ylim)) / 2;
    global_xlim = [mean(global_xlim) - half_span, mean(global_xlim) + half_span];
    global_ylim = [mean(global_ylim) - half_span, mean(global_ylim) + half_span];
else
    global_xlim = [-1.35, 1.35];
    global_ylim = [-1.35, 1.35];
end

%% ---- Figures ----
plot_overlay_figure( ...
    {traj_large_b1, traj_large_b3}, ...
    {choice_large_b1, choice_large_b3}, ...
    {'Large Reward — Block 1', 'Large Reward — Block 3'}, ...
    ChrimsonR_IDs, is_female, cmap_all, line_alpha, line_width, marker_size, choice_marker_size, ...
    norm_squares, col_large_rect, col_small_rect, ...
    global_xlim, global_ylim, ...
    sprintf('Large Reward Trajectories — %s  (LargeRew normalised to left)', ...
        strrep(session_to_analyze, '_', '\_')));

plot_overlay_figure( ...
    {traj_small_b1, traj_small_b3}, ...
    {choice_small_b1, choice_small_b3}, ...
    {'Small Reward — Block 1', 'Small Reward — Block 3'}, ...
    ChrimsonR_IDs, is_female, cmap_all, line_alpha, line_width, marker_size, choice_marker_size, ...
    norm_squares, col_large_rect, col_small_rect, ...
    global_xlim, global_ylim, ...
    sprintf('Small Reward Trajectories — %s  (LargeRew normalised to left)', ...
        strrep(session_to_analyze, '_', '\_')));


%% ================================================================
%  Local functions
%  ================================================================

function plot_overlay_figure(traj_by_block, choice_by_block, subplot_titles, ...
        animal_ids, is_female, cmap_all, line_alpha, line_width, marker_size, choice_marker_size, ...
        norm_squares, col_large_rect, col_small_rect, ...
        shared_xlim, shared_ylim, fig_title)
% Draws 2 subplots side-by-side with a heatmap-style legend to the right.
% Axes positions are hard-coded in pixels so layout is identical across figures.
%
% Legend: two stacked columns of coloured squares, one column per sex.
%   Left column  = Females (greens),  Right column = Males (greys).
%   Number beside each square = trailing digits extracted from animal ID.

% ---- Fixed layout (pixels) ----
fig_w            = 1080;
fig_h            =  540;
sp_w_px          = 340;
left_margin_px   =  70;
bottom_margin_px =  60;
gap_px           =  65;
leg_gap_px       =  20;   % gap between subplot 2 right edge and legend
legend_w_px      = 200;   % width of legend axes area

w_n     = sp_w_px / fig_w;
h_n     = sp_w_px / fig_h;                            % square axes (w = h in px)
b_n     = bottom_margin_px / fig_h;
l1_n    = left_margin_px / fig_w;
l2_n    = (left_margin_px + sp_w_px + gap_px) / fig_w;
leg_l_n = (left_margin_px + 2*sp_w_px + gap_px + leg_gap_px) / fig_w;
leg_w_n = legend_w_px / fig_w;

n_animals = numel(animal_ids);

figure;
set(gcf, 'Position', [80, 80, fig_w, fig_h]);

ax(1) = axes('Position', [l1_n, b_n, w_n, h_n]);
ax(2) = axes('Position', [l2_n, b_n, w_n, h_n]);

for sp = 1:2
    axes(ax(sp)); %#ok<LAXES>
    hold on;

    % Arena outline (normalised unit circle)
    viscircles([0, 0], 1, 'Color', [0.70 0.70 0.70], 'LineWidth', 1.5);
    plot(0, 0, '+', 'Color', [0.75 0.75 0.75], 'MarkerSize', 6, 'LineWidth', 0.8);

    % Reward port rectangles
    for k = 1:numel(norm_squares)
        ec = col_small_rect;
        if norm_squares(k).is_large; ec = col_large_rect; end
        rectangle('Position', norm_squares(k).pos_norm, ...
            'EdgeColor', ec, 'LineWidth', 2.0, 'FaceColor', 'none');
    end

    % Trajectories
    for ii = 1:n_animals
        traj = traj_by_block{sp}{ii};
        if isempty(traj); continue; end
        x   = traj(1, :);
        y   = traj(2, :);
        col = cmap_all(ii, :);

        plot(x, y, '-', 'Color', [col, line_alpha], 'LineWidth', line_width);
        scatter(x(1),   y(1),   marker_size, col, 's', 'filled', 'MarkerEdgeColor', 'none');
        scatter(x(end), y(end), marker_size, col, '^', 'filled', 'MarkerEdgeColor', 'none');

        cp = choice_by_block{sp}{ii};
        if ~isempty(cp) && ~any(isnan(cp))
            scatter(cp(1), cp(2), choice_marker_size, col, 'o', 'filled', ...
                'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', col, 'LineWidth', 1.5);
        end
    end

    xlim(shared_xlim);
    ylim(shared_ylim);
    xlabel('Arena position (radii)', 'FontSize', 10);
    if sp == 1; ylabel('Arena position (radii)', 'FontSize', 10); end
    title(subplot_titles{sp}, 'FontSize', 11);
    box off;
    hold off;
end

% ---- Heatmap-style legend ----
f_colors = cmap_all(is_female,  :);   % greens in animal order
m_colors = cmap_all(~is_female, :);   % greys  in animal order
f_ids    = animal_ids(is_female);
m_ids    = animal_ids(~is_female);
nf       = size(f_colors, 1);
nm       = size(m_colors, 1);
n_max    = max(nf, nm);

leg_ax = axes('Position', [leg_l_n, b_n, leg_w_n, h_n]);
set(leg_ax, 'Visible', 'off', 'XLim', [0, 5], 'YLim', [-0.3, n_max + 1.2]);
hold(leg_ax, 'on');

sq_w = 0.85;   % square width  in legend data units
sq_h = 0.80;   % square height in legend data units
f_x  = 0.0;    % female column left edge
m_x  = 2.5;    % male column left edge

% Column headers
text(leg_ax, f_x + sq_w/2, n_max + 0.85, 'F', ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', ...
    'Color', f_colors(max(end,1), :));
text(leg_ax, m_x + sq_w/2, n_max + 0.85, 'M', ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', ...
    'Color', m_colors(max(end,1), :));

% Female squares (animal 1 at top)
for fi = 1:nf
    y_bot = n_max - fi;
    rectangle(leg_ax, 'Position', [f_x, y_bot, sq_w, sq_h], ...
        'FaceColor', f_colors(fi, :), 'EdgeColor', 'none');
    tok = regexp(f_ids{fi}, '\d+$', 'match');
    num_str = tok{1};
    text(leg_ax, f_x + sq_w + 0.12, y_bot + sq_h/2, num_str, ...
        'FontSize', 7, 'VerticalAlignment', 'middle', 'Color', [0.25 0.25 0.25]);
end

% Male squares
for mi = 1:nm
    y_bot = n_max - mi;
    rectangle(leg_ax, 'Position', [m_x, y_bot, sq_w, sq_h], ...
        'FaceColor', m_colors(mi, :), 'EdgeColor', 'none');
    tok = regexp(m_ids{mi}, '\d+$', 'match');
    num_str = tok{1};
    text(leg_ax, m_x + sq_w + 0.12, y_bot + sq_h/2, num_str, ...
        'FontSize', 7, 'VerticalAlignment', 'middle', 'Color', [0.25 0.25 0.25]);
end

hold(leg_ax, 'off');
sgtitle(fig_title, 'FontSize', 12);
end
