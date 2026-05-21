% generate_behav_figs_PR_sex_comparison.m
%
% Bar graphs and cumulative press curves comparing male vs female behaviour
% across two progressive ratio sessions (PR_D1 and PR_D2).
%
% PR sessions are trimmed to a maximum of 1800 s (30 min) to remove
% accidental session overruns.
%
% Variables extracted per animal per session:
%   breakpoint           — final ratio successfully completed
%   total_presses        — total lever presses
%   rewards_obtained     — number of rewards collected
%   press_rate           — presses per minute (total / session duration)
%   mean_collect_latency — mean time from reward delivery to collection (s)
%
% Figures produced:
%   1 per variable  — 4-bar plot (F-D1, M-D1, F-D2, M-D2)
%   1 cumulative press curve — group mean ± SEM over session time
%
% Requires in workspace (run females_vs_males_groups_correct.m first):
%   ChrimsonR_IDs, ChrimsonR_treatment_groups, valid_sessions,
%   session_weight_raw_grams, freefeed_weight
% Requires: final_behavior
clc

%% Animal list and sex labels
n_animals = numel(ChrimsonR_IDs);
is_female = strcmp(ChrimsonR_treatment_groups, 'Female');
fprintf('Including %d females and %d males\n', sum(is_female), sum(~is_female));

%% Settings
max_session_time = 1800;        % PR session cap (s) — trim any data beyond this
cum_dt           = 5;           % time-grid resolution for cumulative press curve (s)
cum_t_grid       = 0 : cum_dt : max_session_time;   % shared time axis (361 points)

%% Figure dimensions
fig_width  = 280;
fig_height = 400;

%% Colour scheme (females = greens, males = greys; D1 light / D2 dark)
color_f_d1 = [0.60, 0.88, 0.60];
color_f_d2 = [0.13, 0.52, 0.13];
color_m_d1 = [0.72, 0.72, 0.72];
color_m_d2 = [0.28, 0.28, 0.28];
all_colors     = {color_f_d1, color_f_d2, color_m_d1, color_m_d2};
scatter_colors = cellfun(@(c) min(c * 0.75, 1), all_colors, 'UniformOutput', false);
scatter_markers = {'s', 's', 'o', 'o'};   % females = squares, males = circles

%% Behavioural variable names
variable_names = { ...
    'total_presses', ...
    'breakpoint', ...
    'rewards_obtained', ...
    'press_rate', ...
    'mean_collect_latency'};
n_vars = numel(variable_names);

%% Extract data for PR_D1 and PR_D2
sessions    = {'PR_D1', 'PR_D2'};
pr_tables   = cell(1, 2);
cum_presses = cell(2, n_animals);   % cumulative press curves: {session, animal}

for s_idx = 1:2
    session_to_analyze = sessions{s_idx};
    valid_for_session  = cellfun(@(s) any(strcmp(s, session_to_analyze)), valid_sessions);

    data_matrix = nan(n_animals, n_vars);

    for ii = 1:n_animals
        if ~valid_for_session(ii); continue; end

        currentanimal = ChrimsonR_IDs{ii};
        if ~isfield(final_behavior, currentanimal) || ...
           ~isfield(final_behavior.(currentanimal), session_to_analyze)
            continue;
        end

        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        % --- Trim to max_session_time (30 min) ---
        BehavData = BehavData(BehavData.PressTime <= max_session_time, :);
        if isempty(BehavData); continue; end

        % --- Core metrics ---
        total_presses    = BehavData.Trial(end);
        breakpoint       = max(BehavData.Ratio) - 1;    % last completed ratio
        rewards_obtained = sum(BehavData.collectTrial);
        session_length   = BehavData.PressTime(end);
        press_rate       = total_presses / (session_length / 60);  % presses / min

        % --- Collection latency (reward delivery → collection) ---
        if ismember('rewDeliveryTime', BehavData.Properties.VariableNames)
            collected = BehavData.collectTrial == 1;
            if any(collected)
                cl_vals = BehavData.collectionTime(collected) - ...
                          BehavData.rewDeliveryTime(collected);
                mean_collect_latency = nanmean(cl_vals);
            else
                mean_collect_latency = NaN;
            end
        else
            mean_collect_latency = NaN;
        end

        data_matrix(ii, :) = [total_presses, breakpoint, rewards_obtained, ...
                               press_rate, mean_collect_latency];

        % --- Cumulative press curve on shared time grid ---
        press_times = BehavData.PressTime;
        cum_vec = arrayfun(@(t) sum(press_times <= t), cum_t_grid);
        cum_presses{s_idx, ii} = cum_vec;
    end

    pr_table             = array2table(data_matrix, 'VariableNames', variable_names);
    pr_table.animalID    = ChrimsonR_IDs;
    pr_table.sex         = ChrimsonR_treatment_groups;
    pr_table.valid_PR_D1 = cellfun(@(s) any(strcmp(s, 'PR_D1')), valid_sessions);
    pr_table.valid_PR_D2 = cellfun(@(s) any(strcmp(s, 'PR_D2')), valid_sessions);
    pr_tables{s_idx}     = pr_table;
end

pr_table_D1 = pr_tables{1};
pr_table_D2 = pr_tables{2};

%% Body-weight filter
% Re-use bw_include from the RDT workspace if available; otherwise recompute.
if ~exist('bw_include', 'var')
    session_weight_mat = cell2mat(session_weight_raw_grams');
    weight_pct_D1_bw   = (session_weight_mat(:,1) ./ freefeed_weight') * 100 + 5;
    weight_pct_D2_bw   = (session_weight_mat(:,2) ./ freefeed_weight') * 100 + 5;
    bw_threshold       = 70;
    bw_include = (isnan(weight_pct_D1_bw) | weight_pct_D1_bw >= bw_threshold) & ...
                 (isnan(weight_pct_D2_bw) | weight_pct_D2_bw >= bw_threshold);
end

%% Separate by sex (BW filter applied)
female_D1 = pr_table_D1(is_female  & bw_include, :);
female_D2 = pr_table_D2(is_female  & bw_include, :);
male_D1   = pr_table_D1(~is_female & bw_include, :);
male_D2   = pr_table_D2(~is_female & bw_include, :);

f_idx = find(is_female  & bw_include);
m_idx = find(~is_female & bw_include);

fprintf('PR D1 — females: %d valid | males: %d valid\n', ...
    sum(~isnan(female_D1.total_presses)), sum(~isnan(male_D1.total_presses)));
fprintf('PR D2 — females: %d valid | males: %d valid\n', ...
    sum(~isnan(female_D2.total_presses)), sum(~isnan(male_D2.total_presses)));

raw = 1;

%% ---- Bar figures ----

% Breakpoint (final ratio completed)
plot_pr_bars(female_D1, female_D2, male_D1, male_D2, 'breakpoint', raw, ...
    all_colors, scatter_colors, scatter_markers, 'Breakpoint (final ratio)', ...
    [], [], fig_width, fig_height);
run_2way_anova_pr(female_D1, female_D2, male_D1, male_D2, 'breakpoint', 'Breakpoint');

% Total presses
plot_pr_bars(female_D1, female_D2, male_D1, male_D2, 'total_presses', raw, ...
    all_colors, scatter_colors, scatter_markers, 'Total lever presses', ...
    [], [], fig_width, fig_height);
run_2way_anova_pr(female_D1, female_D2, male_D1, male_D2, 'total_presses', 'Total presses');

% Rewards obtained
plot_pr_bars(female_D1, female_D2, male_D1, male_D2, 'rewards_obtained', raw, ...
    all_colors, scatter_colors, scatter_markers, 'Rewards obtained', ...
    [], [], fig_width, fig_height);
run_2way_anova_pr(female_D1, female_D2, male_D1, male_D2, 'rewards_obtained', 'Rewards obtained');

% Press rate
plot_pr_bars(female_D1, female_D2, male_D1, male_D2, 'press_rate', raw, ...
    all_colors, scatter_colors, scatter_markers, 'Press rate (presses/min)', ...
    [], [], fig_width, fig_height);
run_2way_anova_pr(female_D1, female_D2, male_D1, male_D2, 'press_rate', 'Press rate');

% Collection latency
plot_pr_bars(female_D1, female_D2, male_D1, male_D2, 'mean_collect_latency', raw, ...
    all_colors, scatter_colors, scatter_markers, 'Reward collection latency (s)', ...
    [], [], fig_width, fig_height);
run_2way_anova_pr(female_D1, female_D2, male_D1, male_D2, 'mean_collect_latency', 'Collection latency');

%% ---- Cumulative press curve ----
% Mean ± SEM cumulative presses over session time, one line per sex × session.

% Aggregate into matrices [animals × time points], BW filter applied
f_cum_d1 = []; f_cum_d2 = []; m_cum_d1 = []; m_cum_d2 = [];
for k = 1:numel(f_idx)
    ii = f_idx(k);
    if ~isempty(cum_presses{1,ii}); f_cum_d1 = [f_cum_d1; cum_presses{1,ii}]; end %#ok<AGROW>
    if ~isempty(cum_presses{2,ii}); f_cum_d2 = [f_cum_d2; cum_presses{2,ii}]; end %#ok<AGROW>
end
for k = 1:numel(m_idx)
    ii = m_idx(k);
    if ~isempty(cum_presses{1,ii}); m_cum_d1 = [m_cum_d1; cum_presses{1,ii}]; end %#ok<AGROW>
    if ~isempty(cum_presses{2,ii}); m_cum_d2 = [m_cum_d2; cum_presses{2,ii}]; end %#ok<AGROW>
end

% Per-column SEM (handles animals that ran out of time — trailing NaNs)
safe_sem = @(mat) nanstd(mat, 0, 1) ./ sqrt(max(sum(~isnan(mat), 1), 1));

figure('Position', [100, 100, 500, 360]);

% Explicit hold on before each shadedErrorBar call — shadedErrorBar restores
% hold state on exit, so omitting this causes subsequent lines to overwrite.
% 'HandleVisibility','off' keeps the patch and line out of the auto-legend.
hold on;
shadedErrorBar(cum_t_grid, nanmean(f_cum_d1, 1), safe_sem(f_cum_d1), ...
    'lineProps', {'Color', color_f_d1, 'LineWidth', 1.8, 'HandleVisibility', 'off'});
hold on;
shadedErrorBar(cum_t_grid, nanmean(f_cum_d2, 1), safe_sem(f_cum_d2), ...
    'lineProps', {'Color', color_f_d2, 'LineWidth', 1.8, 'HandleVisibility', 'off'});
hold on;
shadedErrorBar(cum_t_grid, nanmean(m_cum_d1, 1), safe_sem(m_cum_d1), ...
    'lineProps', {'Color', color_m_d1, 'LineWidth', 1.8, 'HandleVisibility', 'off'});
hold on;
shadedErrorBar(cum_t_grid, nanmean(m_cum_d2, 1), safe_sem(m_cum_d2), ...
    'lineProps', {'Color', color_m_d2, 'LineWidth', 1.8, 'HandleVisibility', 'off'});

% Manual legend — one dummy line per group with the correct colour
hl(1) = plot(NaN, NaN, '-', 'Color', color_f_d1, 'LineWidth', 2);
hl(2) = plot(NaN, NaN, '-', 'Color', color_f_d2, 'LineWidth', 2);
hl(3) = plot(NaN, NaN, '-', 'Color', color_m_d1, 'LineWidth', 2);
hl(4) = plot(NaN, NaN, '-', 'Color', color_m_d2, 'LineWidth', 2);
legend(hl, {'Female D1','Female D2','Male D1','Male D2'}, ...
    'Location', 'northwest', 'FontSize', 9);

xlabel('Session time (min)', 'FontSize', 12);
ylabel('Cumulative presses', 'FontSize', 12);
title('Progressive Ratio — cumulative press curves', 'FontSize', 12);
xlim([0, max_session_time]);
xticks(0 : 300 : max_session_time);
xticklabels({'0','5','10','15','20','25','30'});
box off;
hold off;


%% ================================================================
%  Local functions
%  ================================================================

% ----- Bar plot: 4 bars (F-D1, M-D1, F-D2, M-D2), x-axis = PR D1 | PR D2 -----
function plot_pr_bars(f_d1, f_d2, m_d1, m_d2, var_name, scale, ...
    all_colors, scatter_colors, scatter_markers, ylabel_str, ylim_vals, ytick_vals, fig_w, fig_h)

% Day-grouped order: F-D1, M-D1, F-D2, M-D2 (indices into input cell arrays)
day_order    = [1, 3, 2, 4];
plot_colors  = all_colors(day_order);
plot_scatter = scatter_colors(day_order);
plot_markers = scatter_markers(day_order);

all_vals = { ...
    f_d1.(var_name) * scale, ...
    m_d1.(var_name) * scale, ...
    f_d2.(var_name) * scale, ...
    m_d2.(var_name) * scale};

% Manual x-positions: D1 pair centred on 1, D2 pair centred on 2
bar_w = 0.35;
x_pos = [0.82, 1.18, 1.82, 2.18];

means_vec = cellfun(@(v) nanmean(v),                               all_vals);
sems_vec  = cellfun(@(v) nanstd(v) / max(sqrt(sum(~isnan(v))),1), all_vals);

figure;
hold on;
set(gcf, 'Position', [100, 100, fig_w, fig_h]);

bar_h = gobjects(4, 1);
for b = 1:4
    bar_h(b) = bar(x_pos(b), means_vec(b), bar_w, ...
        'FaceColor', plot_colors{b}, 'EdgeColor', 'none');
    errorbar(x_pos(b), means_vec(b), sems_vec(b), ...
        'k', 'LineWidth', 1.2, 'CapSize', 5, 'LineStyle', 'none');
    vals = all_vals{b};
    vals = vals(~isnan(vals));
    if ~isempty(vals)
        jitter = (rand(numel(vals), 1) - 0.5) * bar_w * 0.4;
        scatter(x_pos(b) + jitter, vals, 28, plot_scatter{b}, ...
            plot_markers{b}, 'filled', ...
            'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', 'none');
    end
end

set(gca, 'XTick', [1, 2], 'XTickLabel', {'PR D1', 'PR D2'}, 'FontSize', 12);
xlim([0.4, 2.6]);
ylabel(ylabel_str, 'FontSize', 12);
if ~isempty(ylim_vals);  ylim(ylim_vals);    end
if ~isempty(ytick_vals); yticks(ytick_vals); end
legend(bar_h, {'Female D1','Male D1','Female D2','Male D2'}, ...
    'Location', 'best', 'FontSize', 9);
box off;

margin_left   = 65;
margin_bottom = 55;
margin_right  = 15;
margin_top    = 25;
ax = gca;
ax.Units = 'pixels';
ax.Position = [margin_left, margin_bottom, ...
               fig_w - margin_left - margin_right, ...
               fig_h - margin_bottom - margin_top];
hold off;
end


% ----- 2-way mixed ANOVA: Sex (between) × Session (within) -----
function run_2way_anova_pr(f_d1, f_d2, m_d1, m_d2, var_name, measure_name)

f_d1_vals = f_d1.(var_name)(:);
f_d2_vals = f_d2.(var_name)(:);
m_d1_vals = m_d1.(var_name)(:);
m_d2_vals = m_d2.(var_name)(:);

nF = numel(f_d1_vals);
nM = numel(m_d1_vals);

all_wide   = [f_d1_vals, f_d2_vals; m_d1_vals, m_d2_vals];
sex_labels = categorical([repmat({'Female'}, nF, 1); repmat({'Male'}, nM, 1)]);

% Listwise deletion
complete   = all(~isnan(all_wide), 2);
n_complete = sum(complete);
n_excluded = (nF + nM) - n_complete;
nF_c       = sum(complete(1:nF));
nM_c       = sum(complete(nF+1:end));

all_wide_c   = all_wide(complete, :);
sex_labels_c = sex_labels(complete);

fprintf('\n========== 2-Way Mixed ANOVA: %s ==========\n', measure_name);
fprintf('N = %d subjects (%d female, %d male)', n_complete, nF_c, nM_c);
if n_excluded > 0
    fprintf('  [%d excluded — missing session data]', n_excluded);
end
fprintf('\n');

if nF_c < 2 || nM_c < 2
    fprintf('  [Insufficient data for ANOVA]\n');
    return;
end

try
    tbl     = array2table(all_wide_c, 'VariableNames', {'D1','D2'});
    tbl.Sex = sex_labels_c;
    within  = table(categorical({'D1';'D2'}), 'VariableNames', {'Session'});
    rm      = fitrm(tbl, 'D1,D2 ~ Sex', 'WithinDesign', within);

    % --- Between-subjects: Sex ---
    fprintf('\n  Between-subjects:\n');
    bs_tbl  = rm.anova();
    bs_rows = bs_tbl.Properties.RowNames;
    bs_vars = bs_tbl.Properties.VariableNames;
    err_labels_bs = {'Error(Sex)','Error','Residuals','Residual'};

    sex_idx = [];  err_idx_bs = [];
    if ~isempty(bs_rows)
        sex_idx = find(strcmp(bs_rows, 'Sex'));
        for lbl = err_labels_bs
            idx = find(strcmp(bs_rows, lbl{1}));
            if ~isempty(idx); err_idx_bs = idx; break; end
        end
    elseif ismember('Between', bs_vars)
        btw_vals = cellstr(bs_tbl.Between);
        sex_idx  = find(strcmp(btw_vals, 'Sex'));
        for lbl = err_labels_bs
            idx = find(strcmp(btw_vals, lbl{1}));
            if ~isempty(idx); err_idx_bs = idx; break; end
        end
    end

    p_sex = NaN;
    if ~isempty(sex_idx)
        F_val = NaN; p_val = NaN; df1 = NaN; df2 = NaN;
        if ismember('F',      bs_vars); F_val = bs_tbl.F(sex_idx);      end
        if ismember('pValue', bs_vars); p_val = bs_tbl.pValue(sex_idx); end
        if ismember('DF',     bs_vars); df1   = bs_tbl.DF(sex_idx);     end
        if ~isempty(err_idx_bs) && ismember('DF', bs_vars)
            df2 = bs_tbl.DF(err_idx_bs);
        end
        if ~isnan(F_val) && ~isnan(p_val)
            print_effect_pr('Sex', F_val, df1, df2, p_val);
            p_sex = p_val;
        end
    else
        fprintf('    [Sex row not found — bs_tbl rows: {%s}  cols: {%s}]\n', ...
            strjoin(bs_rows, ', '), strjoin(bs_vars, ', '));
    end

    % --- Within-subjects: Session, Sex × Session ---
    fprintf('\n  Within-subjects:\n');
    ws_tbl  = ranova(rm, 'WithinModel', 'Session');
    ws_rows = ws_tbl.Properties.RowNames;

    effects_ws = { ...
        '(Intercept):Session', 'Session',        'Error(Session)'; ...
        'Sex:Session',         'Sex x Session',  'Error(Session)'};

    for e = 1:size(effects_ws, 1)
        eff_idx = find(strcmp(ws_rows, effects_ws{e,1}));
        err_idx = find(strcmp(ws_rows, effects_ws{e,3}));
        if ~isempty(eff_idx) && ~isempty(err_idx)
            print_effect_pr(effects_ws{e,2}, ws_tbl.F(eff_idx), ...
                ws_tbl.DF(eff_idx), ws_tbl.DF(err_idx), ws_tbl.pValue(eff_idx));
        end
    end

    % --- Post-hoc ---
    p_session     = get_p_pr(ws_tbl, '(Intercept):Session');
    p_sex_session = get_p_pr(ws_tbl, 'Sex:Session');

    if p_sex_session < 0.05
        fprintf('\n  Significant Sex x Session — Sex within each Session (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Session', 'ComparisonType', 'tukey-kramer');
        print_multcomp_pr(mc);
        fprintf('\n  Session within each Sex (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Session', 'By', 'Sex', 'ComparisonType', 'tukey-kramer');
        print_multcomp_pr(mc);
    elseif p_session < 0.05
        fprintf('\n  Significant Session — D1 vs D2 (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Session', 'ComparisonType', 'tukey-kramer');
        print_multcomp_pr(mc);
    end

    if ~isnan(p_sex) && p_sex < 0.05 && p_sex_session >= 0.05
        fprintf('\n  Significant Sex main effect — Female vs Male (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'ComparisonType', 'tukey-kramer');
        print_multcomp_pr(mc);
    end

catch ME
    fprintf('  [ANOVA failed: %s]\n', ME.message);
end
end


% ----- Helpers -----
function p = get_p_pr(tbl, row_name)
    try
        p = tbl{row_name, 'pValue'};
    catch
        idx = find(strcmp(tbl.Properties.RowNames, row_name));
        if ~isempty(idx); p = tbl.pValue(idx(1)); else; p = NaN; end
    end
end

function print_effect_pr(label, F, df1, df2, p)
    fprintf('    %-26s F(%d,%d) = %7.3f,  p = %.4f%s\n', ...
        label, df1, df2, F, p, sig_stars_pr(p));
end

function print_multcomp_pr(mc)
    if ~istable(mc)
        for r = 1:size(mc, 1)
            p_val = mc(r, 6);
            fprintf('      Group %d vs %d: diff=%7.3f, 95%%CI[%6.3f,%6.3f], p=%.4f%s\n', ...
                mc(r,1), mc(r,2), mc(r,4), mc(r,3), mc(r,5), p_val, sig_stars_pr(p_val));
        end
        return;
    end
    if ~ismember('pValue', mc.Properties.VariableNames); disp(mc); return; end

    % Count leading categorical/cell columns to determine 'By' form
    n_cat = 0;
    for c = 1:width(mc)
        if iscategorical(mc{:,c}) || iscell(mc{:,c}); n_cat = n_cat + 1; else; break; end
    end

    seen_pairs = {};
    for r = 1:height(mc)
        p_val = mc.pValue(r);
        if n_cat >= 3
            by_lv  = char(mc{r,1}); lv1 = char(mc{r,2}); lv2 = char(mc{r,3});
            by_str = [' [' by_lv ']'];
        elseif n_cat >= 2
            lv1 = char(mc{r,1}); lv2 = char(mc{r,2}); by_lv = ''; by_str = '';
        else
            lv1 = char(mc{r,1}); lv2 = '?'; by_lv = ''; by_str = '';
        end
        sorted_pair = sort({lv1, lv2});
        pair_key    = [sorted_pair{1} '|' sorted_pair{2} '|' by_lv];
        if any(strcmp(seen_pairs, pair_key)); continue; end
        seen_pairs{end+1} = pair_key; %#ok<AGROW>
        fprintf('      %-42s p = %.4f%s\n', ...
            [lv1 ' vs ' lv2 by_str], p_val, sig_stars_pr(p_val));
    end
end

function s = sig_stars_pr(p)
    if     p < 0.001; s = ' ***';
    elseif p < 0.01;  s = ' **';
    elseif p < 0.05;  s = ' *';
    else;             s = '';
    end
end
