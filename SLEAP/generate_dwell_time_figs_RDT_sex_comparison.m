% generate_dwell_time_figs_RDT_sex_comparison.m
%
% Extracts per-trial zone dwell times from SLEAP tracking data and plots
% Female vs Male bar graphs (same style as generate_behav_figs_RDT_sex_comparison.m).
%
% Zones (detected by type from shapeData):
%   Circle     — reward cup area
%   Square 1   — one reward screen (left or right)
%   Square 2   — other reward screen
%   Zone 0     — outside all defined shapes (corridor / other)
%
% Large/small reward screen assignment uses large_screen_side per animal.
% Dwell times expressed as % of trial duration.
%
% Four figures produced (one per zone metric), each with bars:
%   X-axis  : Block 1 | Block 2 | Block 3
%   Per block: Female-D1 | Male-D1 | Female-D2 | Male-D2
%
% Requires in workspace (run females_vs_males_groups_correct.m first):
%   ChrimsonR_IDs, ChrimsonR_treatment_groups, large_screen_side, valid_sessions
% Requires: final_SLEAP, final_behavior
clc
%% ---- Settings ----
sessions   = {'RDT_D1', 'RDT_D2'};
fs_cam     = 30;

%% ---- Figure / colour settings (matches generate_behav_figs) ----
fig_width  = 320;
fig_height = 400;

color_f_d1 = [0.60, 0.88, 0.60];
color_f_d2 = [0.13, 0.52, 0.13];
color_m_d1 = [0.72, 0.72, 0.72];
color_m_d2 = [0.28, 0.28, 0.28];
all_colors     = {color_f_d1, color_f_d2, color_m_d1, color_m_d2};
scatter_colors = cellfun(@(c) min(c * 0.75, 1), all_colors, 'UniformOutput', false);
scatter_markers = {'s', 's', 'o', 'o'};

%% ---- Body-weight inclusion threshold (matches generate_behav_figs) ----
bw_threshold = 75;   % % of free-feed weight; set to 0 to include all

%% SLEAP-specific exclusions (video quality issues — separate from valid_sessions)
% Format: each row is {session, animal_ID}
dwell_exclusions = { ...
    'RDT_D2', 'RDT_M_4'; ...
    'RDT_D1', 'RDT_M_5'; ...
};

%% ---- Extract dwell times for RDT_D1 and RDT_D2 ----
n_animals    = numel(ChrimsonR_IDs);
is_female    = strcmp(ChrimsonR_treatment_groups, 'Female');

% Pre-compute BW filter (same logic as generate_behav_figs)
session_weight_mat = cell2mat(session_weight_raw_grams');   % n_animals × 2
% weight_pct_D1 = (session_weight_mat(:,1) ./ freefeed_weight') * 100;
% weight_pct_D2 = (session_weight_mat(:,2) ./ freefeed_weight') * 100;


weight_pct_D1 = [(session_weight_mat(:,1) ./ freefeed_weight') * 100]+5;
weight_pct_D2 = [(session_weight_mat(:,2) ./ freefeed_weight') * 100]+5;

bw_include = (isnan(weight_pct_D1) | weight_pct_D1 >= bw_threshold) & ...
             (isnan(weight_pct_D2) | weight_pct_D2 >= bw_threshold);

dwell_tables = cell(1, 2);

for s_idx = 1:2
    session_to_analyze = sessions{s_idx};

    % Per-animal output arrays (NaN = missing)
    large_screen_B1 = nan(n_animals, 1);
    large_screen_B2 = nan(n_animals, 1);
    large_screen_B3 = nan(n_animals, 1);
    small_screen_B1 = nan(n_animals, 1);
    small_screen_B2 = nan(n_animals, 1);
    small_screen_B3 = nan(n_animals, 1);
    reward_cup_B1   = nan(n_animals, 1);
    reward_cup_B2   = nan(n_animals, 1);
    reward_cup_B3   = nan(n_animals, 1);
    other_zone_B1   = nan(n_animals, 1);
    other_zone_B2   = nan(n_animals, 1);
    other_zone_B3   = nan(n_animals, 1);

    for ii = 1:n_animals
        currentanimal = ChrimsonR_IDs{ii};

        % Skip if session not valid for this animal
        if ~any(strcmp(valid_sessions{ii}, session_to_analyze)); continue; end

        % Skip SLEAP-specific exclusions (e.g. bad video)
        if any(strcmp(dwell_exclusions(:,1), session_to_analyze) & ...
               strcmp(dwell_exclusions(:,2), currentanimal))
            fprintf('  [Excluding %s %s — SLEAP exclusion list]\n', currentanimal, session_to_analyze);
            continue;
        end

        if ~isfield(final_SLEAP, currentanimal) || ...
           ~isfield(final_SLEAP.(currentanimal), session_to_analyze); continue; end
        if ~isfield(final_behavior, currentanimal) || ...
           ~isfield(final_behavior.(currentanimal), session_to_analyze); continue; end

        SLEAP_data = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data_raw;
        shapeData  = final_SLEAP.(currentanimal).(session_to_analyze).shapeData;
        BehavData  = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        % ---- Identify zone shapes by type ----
        circle_idx = [];
        square_idx = [];
        for k = 1:numel(shapeData)
            if strcmp(shapeData{k}.Type, 'Circle')
                % Use the largest circle as the reward cup zone
                if isempty(circle_idx) || shapeData{k}.Radius > shapeData{circle_idx}.Radius
                    circle_idx = k;
                end
            elseif strcmp(shapeData{k}.Type, 'Square')
                square_idx(end+1) = k; %#ok<AGROW>
            end
        end
        if isempty(circle_idx) || numel(square_idx) < 2
            fprintf('[%s %s] shapeData incomplete — skipping\n', currentanimal, session_to_analyze);
            continue;
        end

        % ---- Identify which square is large / small reward screen ----
        % Squares are labelled by Location ('left screen' / 'right screen')
        sq1 = shapeData{square_idx(1)};
        sq2 = shapeData{square_idx(2)};
        large_rew_side = large_screen_side{ii};
        if strcmp(sq1.Location, [large_rew_side ' screen'])
            large_sq = sq1;  small_sq = sq2;
        else
            large_sq = sq2;  small_sq = sq1;
        end

        % ---- Smooth XY ----
        X_data = sgolayfilt(SLEAP_data.x_pix, 9, 33);
        Y_data = sgolayfilt(SLEAP_data.y_pix, 9, 33);

        % ---- Trial timing ----
        onset_t  = BehavData.stTime';
        offset_t = BehavData.collectionTime';

        % ---- Valid trial mask (exclude omissions and blank touches) ----
        has_omit  = ismember('omissionALL', BehavData.Properties.VariableNames);
        has_blank = ismember('Blank_Touch', BehavData.Properties.VariableNames);
        valid_trial = true(height(BehavData), 1);
        if has_omit;  valid_trial = valid_trial & BehavData.omissionALL ~= 1; end
        if has_blank; valid_trial = valid_trial & BehavData.Blank_Touch  ~= 1; end

        % ---- Per-trial zone classification ----
        circ  = shapeData{circle_idx};
        n_trials = height(BehavData);
        frac_cup   = nan(n_trials, 1);
        frac_large = nan(n_trials, 1);
        frac_small = nan(n_trials, 1);
        frac_other = nan(n_trials, 1);

        for jj = 1:n_trials
            if ~valid_trial(jj); continue; end
            t_on  = onset_t(jj);
            t_off = offset_t(jj);
            if isinf(t_off) || isnan(t_off); continue; end

            mask = SLEAP_data.idx_time > t_on & SLEAP_data.idx_time < t_off;
            if ~any(mask); continue; end

            x = X_data(mask);
            y = Y_data(mask);
            n_frames = numel(x);

            zone = zeros(n_frames, 1);
            for hh = 1:n_frames
                % Reward cup circle
                if sqrt((x(hh)-circ.Center(1))^2 + (y(hh)-circ.Center(2))^2) <= circ.Radius
                    zone(hh) = 1;
                % Large reward screen
                elseif x(hh) >= large_sq.Center(1) - large_sq.Size(1)/2 && ...
                       x(hh) <= large_sq.Center(1) + large_sq.Size(1)/2 && ...
                       y(hh) >= large_sq.Center(2) - large_sq.Size(2)/2 && ...
                       y(hh) <= large_sq.Center(2) + large_sq.Size(2)/2
                    zone(hh) = 2;
                % Small reward screen
                elseif x(hh) >= small_sq.Center(1) - small_sq.Size(1)/2 && ...
                       x(hh) <= small_sq.Center(1) + small_sq.Size(1)/2 && ...
                       y(hh) >= small_sq.Center(2) - small_sq.Size(2)/2 && ...
                       y(hh) <= small_sq.Center(2) + small_sq.Size(2)/2
                    zone(hh) = 3;
                % else zone = 0 (other)
                end
            end

            frac_cup(jj)   = sum(zone == 1) / n_frames;
            frac_large(jj) = sum(zone == 2) / n_frames;
            frac_small(jj) = sum(zone == 3) / n_frames;
            frac_other(jj) = sum(zone == 0) / n_frames;
        end

        % ---- Average per block ----
        for blk = 1:3
            blk_mask = valid_trial & BehavData.Block == blk;
            switch blk
                case 1
                    large_screen_B1(ii) = nanmean(frac_large(blk_mask)) * 100;
                    small_screen_B1(ii) = nanmean(frac_small(blk_mask)) * 100;
                    reward_cup_B1(ii)   = nanmean(frac_cup(blk_mask))   * 100;
                    other_zone_B1(ii)   = nanmean(frac_other(blk_mask)) * 100;
                case 2
                    large_screen_B2(ii) = nanmean(frac_large(blk_mask)) * 100;
                    small_screen_B2(ii) = nanmean(frac_small(blk_mask)) * 100;
                    reward_cup_B2(ii)   = nanmean(frac_cup(blk_mask))   * 100;
                    other_zone_B2(ii)   = nanmean(frac_other(blk_mask)) * 100;
                case 3
                    large_screen_B3(ii) = nanmean(frac_large(blk_mask)) * 100;
                    small_screen_B3(ii) = nanmean(frac_small(blk_mask)) * 100;
                    reward_cup_B3(ii)   = nanmean(frac_cup(blk_mask))   * 100;
                    other_zone_B3(ii)   = nanmean(frac_other(blk_mask)) * 100;
            end
        end
    end % animal loop

    % Build table
    dt = table();
    dt.animalID = ChrimsonR_IDs;
    dt.sex      = ChrimsonR_treatment_groups;
    dt.large_screen_B1 = large_screen_B1;
    dt.large_screen_B2 = large_screen_B2;
    dt.large_screen_B3 = large_screen_B3;
    dt.small_screen_B1 = small_screen_B1;
    dt.small_screen_B2 = small_screen_B2;
    dt.small_screen_B3 = small_screen_B3;
    dt.reward_cup_B1   = reward_cup_B1;
    dt.reward_cup_B2   = reward_cup_B2;
    dt.reward_cup_B3   = reward_cup_B3;
    dt.other_zone_B1   = other_zone_B1;
    dt.other_zone_B2   = other_zone_B2;
    dt.other_zone_B3   = other_zone_B3;

    dwell_tables{s_idx} = dt;
    fprintf('Finished extraction: %s\n', session_to_analyze);
end

dt_D1 = dwell_tables{1};
dt_D2 = dwell_tables{2};

%% ---- Apply BW filter and separate by sex ----
f_mask = is_female  & bw_include;
m_mask = ~is_female & bw_include;

female_D1 = dt_D1(f_mask, :);
female_D2 = dt_D2(f_mask, :);
male_D1   = dt_D1(m_mask, :);
male_D2   = dt_D2(m_mask, :);

fprintf('Females: %d  |  Males: %d  (after BW filter)\n', sum(f_mask), sum(m_mask));

%% ---- Shared Y-axis limits (derived from data across all zones and sessions) ----
all_zone_vars = { ...
    'large_screen_B1','large_screen_B2','large_screen_B3', ...
    'small_screen_B1','small_screen_B2','small_screen_B3', ...
    'reward_cup_B1',  'reward_cup_B2',  'reward_cup_B3', ...
    'other_zone_B1',  'other_zone_B2',  'other_zone_B3'};

all_dwell_vals = [ ...
    table2array(female_D1(:, all_zone_vars)); ...
    table2array(female_D2(:, all_zone_vars)); ...
    table2array(male_D1(:,   all_zone_vars)); ...
    table2array(male_D2(:,   all_zone_vars))];

global_max = max(all_dwell_vals(:), [], 'omitnan');

% Pick the smallest 'nice' tick step that gives no more than 6 ticks
nice_steps = [1, 2, 5, 10, 20, 25, 50];
shared_step = nice_steps(end);
for ss = nice_steps
    ylim_top = ceil(global_max / ss) * ss;
    if ylim_top / ss <= 6
        shared_step = ss;
        break;
    end
end
ylim_top      = ceil(global_max / shared_step) * shared_step;
shared_ylim   = [0, ylim_top];
shared_yticks = 0 : shared_step : ylim_top;

%% ---- Figures ----
raw = 1;   % values already in %

% Large reward screen dwell time
vars = {'large_screen_B1','large_screen_B2','large_screen_B3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Large screen dwell time (%)', ...
    shared_ylim, shared_yticks, fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Large screen dwell time');

% Small reward screen dwell time
vars = {'small_screen_B1','small_screen_B2','small_screen_B3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Small screen dwell time (%)', ...
    shared_ylim, shared_yticks, fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Small screen dwell time');

% Reward cup dwell time
vars = {'reward_cup_B1','reward_cup_B2','reward_cup_B3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Reward cup dwell time (%)', ...
    shared_ylim, shared_yticks, fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Reward cup dwell time');

% Other zone dwell time
vars = {'other_zone_B1','other_zone_B2','other_zone_B3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Other zone dwell time (%)', ...
    shared_ylim, shared_yticks, fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Other zone dwell time');


%% ================================================================
%  Local functions  (identical to generate_behav_figs_RDT_sex_comparison.m)
%  ================================================================

function plot_sex_session_bars(f_d1, f_d2, m_d1, m_d2, var_names, scale, ...
    all_colors, scatter_colors, scatter_markers, ylabel_str, ylim_vals, ytick_vals, fig_w, fig_h)

day_order = [1, 3, 2, 4];
all_colors      = all_colors(day_order);
scatter_colors  = scatter_colors(day_order);
scatter_markers = scatter_markers(day_order);

all_data_sex = { ...
    table2array(f_d1(:, var_names)) * scale, ...
    table2array(f_d2(:, var_names)) * scale, ...
    table2array(m_d1(:, var_names)) * scale, ...
    table2array(m_d2(:, var_names)) * scale};
all_data = all_data_sex(day_order);

ngroups    = 3;
nbars      = 4;
groupwidth = min(0.8, nbars / (nbars + 1.5));

bar_x = zeros(ngroups, nbars);
for b = 1:nbars
    bar_x(:, b) = (1:ngroups)' - groupwidth/2 + (2*b-1) * groupwidth / (2*nbars);
end

means_mat = zeros(ngroups, nbars);
sems_mat  = zeros(ngroups, nbars);
for g = 1:nbars
    d = all_data{g};
    n = max(sum(~isnan(d), 1), 1);
    means_mat(:, g) = nanmean(d, 1)';
    sems_mat(:, g)  = nanstd(d, 0, 1)' ./ sqrt(n)';
end

figure;
hold on;
set(gcf, 'Position', [100, 100, fig_w, fig_h]);

bar_h = bar(1:ngroups, means_mat, 'grouped', 'BarWidth', 0.85);
for g = 1:nbars
    bar_h(g).FaceColor = all_colors{g};
    bar_h(g).EdgeColor = 'none';
end

for g = 1:nbars
    d = all_data{g};
    for blk = 1:ngroups
        xc = bar_x(blk, g);
        errorbar(xc, means_mat(blk, g), sems_mat(blk, g), ...
            'k', 'LineWidth', 1.2, 'CapSize', 5, 'LineStyle', 'none');
        vals = d(~isnan(d(:, blk)), blk);
        if ~isempty(vals)
            jitter = (rand(numel(vals), 1) - 0.5) * (groupwidth / (nbars * 3));
            scatter(xc + jitter, vals, 28, scatter_colors{g}, ...
                scatter_markers{g}, 'filled', ...
                'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', 'none');
        end
    end
end

set(gca, 'XTick', 1:3, 'XTickLabel', {'Block 1','Block 2','Block 3'}, 'FontSize', 12);
ylabel(ylabel_str, 'FontSize', 12);
if ~isempty(ylim_vals);  ylim(ylim_vals);   end
if ~isempty(ytick_vals); yticks(ytick_vals); end
legend(bar_h, {'Female D1','Male D1','Female D2','Male D2'}, ...
    'Location', 'best', 'FontSize', 9);
box off;

margin_left = 65;  margin_bottom = 55;
margin_right = 15; margin_top   = 25;
ax = gca;
ax.Units = 'pixels';
ax.Position = [margin_left, margin_bottom, ...
               fig_w - margin_left - margin_right, ...
               fig_h - margin_bottom - margin_top];
hold off;
end


function run_3way_anova(f_d1, f_d2, m_d1, m_d2, var_names, scale, measure_name)

data_fd1 = table2array(f_d1(:, var_names)) * scale;
data_fd2 = table2array(f_d2(:, var_names)) * scale;
data_md1 = table2array(m_d1(:, var_names)) * scale;
data_md2 = table2array(m_d2(:, var_names)) * scale;

nF = size(data_fd1, 1);
nM = size(data_md1, 1);

all_wide   = [data_fd1, data_fd2; data_md1, data_md2];
sex_labels = categorical([repmat({'Female'}, nF, 1); repmat({'Male'}, nM, 1)]);

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
    rm     = fitrm(tbl, 'D1B1,D1B2,D1B3,D2B1,D2B2,D2B3 ~ Sex', 'WithinDesign', within);
    bs_tbl = rm.anova();
    ws_tbl = ranova(rm, 'WithinModel', 'Day*Block');

    ws_rows = ws_tbl.Properties.RowNames;
    bs_rows = bs_tbl.Properties.RowNames;

    fprintf('\n  Between-subjects:\n');
    bs_vars = bs_tbl.Properties.VariableNames;

    % rm.anova() layout differs by MATLAB version:
    %   Format A: source names as row names
    %   Format B: no row names; source in a 'Between' column
    err_labels_bs = {'Error(Sex)', 'Error', 'Residuals', 'Residual'};
    sex_idx = [];  err_idx = [];
    if ~isempty(bs_rows)
        sex_idx = find(strcmp(bs_rows, 'Sex'));
        for lbl = err_labels_bs
            idx = find(strcmp(bs_rows, lbl{1}));
            if ~isempty(idx); err_idx = idx; break; end
        end
    elseif ismember('Between', bs_vars)
        btw_vals = cellstr(bs_tbl.Between);
        sex_idx  = find(strcmp(btw_vals, 'Sex'));
        for lbl = err_labels_bs
            idx = find(strcmp(btw_vals, lbl{1}));
            if ~isempty(idx); err_idx = idx; break; end
        end
    end

    p_sex = NaN;
    if ~isempty(sex_idx)
        F_val = NaN;  p_val = NaN;  df1 = NaN;  df2 = NaN;
        if ismember('F',      bs_vars); F_val = bs_tbl.F(sex_idx);      end
        if ismember('pValue', bs_vars); p_val = bs_tbl.pValue(sex_idx); end
        if ismember('DF',     bs_vars); df1   = bs_tbl.DF(sex_idx);     end
        if ~isempty(err_idx) && ismember('DF', bs_vars)
            df2 = bs_tbl.DF(err_idx);
        end
        if ~isnan(F_val) && ~isnan(p_val)
            print_effect('Sex', F_val, df1, df2, p_val);
            p_sex = p_val;
        end
    else
        fprintf('    [Sex row not found — bs_tbl rows: {%s}  cols: {%s}]\n', ...
            strjoin(bs_rows, ', '), strjoin(bs_vars, ', '));
    end

    fprintf('  Within-subjects and interactions:\n');
    effects_ws = { ...
        '(Intercept):Day',       'Day',               'Error(Day)'; ...
        'Sex:Day',               'Sex x Day',          'Error(Day)'; ...
        '(Intercept):Block',     'Block',              'Error(Block)'; ...
        'Sex:Block',             'Sex x Block',         'Error(Block)'; ...
        '(Intercept):Day:Block', 'Day x Block',         'Error(Day:Block)'; ...
        'Sex:Day:Block',         'Sex x Day x Block',   'Error(Day:Block)'};

    for e = 1:size(effects_ws, 1)
        eff_idx = find(strcmp(ws_rows, effects_ws{e,1}));
        err_idx_w = find(strcmp(ws_rows, effects_ws{e,3}));
        if ~isempty(eff_idx) && ~isempty(err_idx_w)
            print_effect(effects_ws{e,2}, ws_tbl.F(eff_idx), ...
                         ws_tbl.DF(eff_idx), ws_tbl.DF(err_idx_w), ws_tbl.pValue(eff_idx));
        end
    end

    p_3way      = get_p(ws_tbl, 'Sex:Day:Block');
    p_sex_block = get_p(ws_tbl, 'Sex:Block');
    p_sex_day   = get_p(ws_tbl, 'Sex:Day');
    p_day_block = get_p(ws_tbl, '(Intercept):Day:Block');
    p_block     = get_p(ws_tbl, '(Intercept):Block');
    p_day       = get_p(ws_tbl, '(Intercept):Day');

    ran_posthoc = false;
    if p_3way < 0.05
        fprintf('\n  Significant 3-way — Sex x Day within each Block (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        mc = multcompare(rm, 'Day', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if p_sex_block < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex x Block (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if p_sex_day < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex x Day (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if p_day_block < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Day x Block (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Day', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if p_block < 0.05 && p_sex_block >= 0.05 && p_day_block >= 0.05 && p_3way >= 0.05
        fprintf('\n  Significant Block main effect (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
    end
    if p_day < 0.05 && p_sex_day >= 0.05 && p_day_block >= 0.05 && p_3way >= 0.05
        fprintf('\n  Significant Day main effect (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Day', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
    end
    if ~isnan(p_sex) && p_sex < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex main effect (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
    end

catch ME
    fprintf('  [ANOVA failed: %s]\n', ME.message);
end
end


function p = get_p(tbl, row_name)
    try
        p = tbl{row_name, 'pValue'};
    catch
        idx = find(strcmp(tbl.Properties.RowNames, row_name));
        if ~isempty(idx); p = tbl.pValue(idx(1)); else; p = NaN; end
    end
end

function print_effect(label, F, df1, df2, p)
    fprintf('    %-26s F(%d,%d) = %7.3f,  p = %.4f%s\n', label, df1, df2, F, p, sig_stars(p));
end

function print_multcomp(mc)
    if ~istable(mc)
        for r = 1:size(mc, 1)
            p_val = mc(r, 6);
            fprintf('      Group %d vs %d: diff=%7.3f, 95%%CI[%6.3f,%6.3f], p=%.4f%s\n', ...
                mc(r,1), mc(r,2), mc(r,4), mc(r,3), mc(r,5), p_val, sig_stars(p_val));
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

        fprintf('      %-42s p = %.4f%s\n', [lv1 ' vs ' lv2 by_str], p_val, sig_stars(p_val));
    end
end

function s = sig_stars(p)
    if p < 0.001;    s = ' ***';
    elseif p < 0.01; s = ' **';
    elseif p < 0.05; s = ' *';
    else;            s = '';
    end
end
