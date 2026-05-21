% generate_loco_figs_RDT_sex_comparison.m
% Bar graphs comparing male vs female locomotor path lengths on RDT D1 and RDT D2
% Variables: mean XY path length (pixels) during large / small reward trials, per block
% Figures and 3-way ANOVAs match the style of generate_behav_figs_RDT_sex_comparison.m
%
% Requires in workspace (run females_vs_males_groups_correct.m first):
%   ChrimsonR_IDs, ChrimsonR_treatment_groups, valid_sessions
% Requires: final_SLEAP   (locomotor data)
%           final_behavior (trial timing — stTime, collectionTime, bigSmall, Block)

%% Sex labels and animal list
n_animals = numel(ChrimsonR_IDs);
is_female = strcmp(ChrimsonR_treatment_groups, 'Female');
fprintf('Including %d females and %d males\n', sum(is_female), sum(~is_female));

%% Path-length variable names
loco_var_names = { ...
    'b1_large_path', 'b2_large_path', 'b3_large_path', ...
    'b1_small_path', 'b2_small_path', 'b3_small_path'};
n_loco_vars = numel(loco_var_names);

%% SLEAP-specific exclusions (video quality issues — separate from valid_sessions)
% Format: each row is {session, animal_ID}
loco_exclusions = { ...
    'RDT_D2', 'RDT_M_4'; ...
    'RDT_D1', 'RDT_M_5'; ...
};

%% Extract path lengths for RDT_D1 and RDT_D2
sessions    = {'RDT_D1', 'RDT_D2'};
loco_tables = cell(1, 2);

for s_idx = 1:2
    session_to_analyze = sessions{s_idx};
    valid_for_session  = cellfun(@(s) any(strcmp(s, session_to_analyze)), valid_sessions);

    data_matrix = nan(n_animals, n_loco_vars);

    for ii = 1:n_animals
        if ~valid_for_session(ii); continue; end

        currentanimal = ChrimsonR_IDs{ii};

        % Skip SLEAP-specific exclusions (e.g. bad video)
        if any(strcmp(loco_exclusions(:,1), session_to_analyze) & ...
               strcmp(loco_exclusions(:,2), currentanimal))
            fprintf('  [Excluding %s %s — SLEAP exclusion list]\n', currentanimal, session_to_analyze);
            continue;
        end

        % Both SLEAP and behavioral data must exist for this animal/session
        if ~isfield(final_SLEAP, currentanimal) || ...
           ~isfield(final_SLEAP.(currentanimal), session_to_analyze) || ...
           ~isfield(final_behavior, currentanimal) || ...
           ~isfield(final_behavior.(currentanimal), session_to_analyze)
            fprintf('  [Skipping %s %s — SLEAP or BehavData missing]\n', ...
                currentanimal, session_to_analyze);
            continue;
        end

        % ----- SLEAP: load and smooth XY coordinates (identical to get_path_lengths_loop_GoPro) -----
        SLEAP_data = final_SLEAP.(currentanimal).(session_to_analyze).SLEAP_data_raw;
        X_data = sgolayfilt(SLEAP_data.x_pix, 9, 33);
        Y_data = sgolayfilt(SLEAP_data.y_pix, 9, 33);

        % ----- Trial timing from behavioral data -----
        BehavData     = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
        onset_trials  = BehavData.stTime';
        choice_trials = BehavData.choiceTime';
        offset_trials = BehavData.collectionTime';
        fs_cam = 30;
        time_ranges_trials = [onset_trials; choice_trials; offset_trials];

        % ----- Pass 1: build filtered_motion cell (exact copy of original) -----
        filtered_motion = [];
        max_ind    = SLEAP_data.idx_frame(end);
        good_index = 1;
        for j = 1:size(time_ranges_trials, 2)
            onset  = round(time_ranges_trials(1,j) * fs_cam) + 1;
            offset = round(time_ranges_trials(3,j) * fs_cam) + 1;

            if isinf(offset)
                if onset <= max_ind && onset > 0
                    filtered_motion{good_index} = SLEAP_data(onset:end); %#ok<AGROW>
                    break
                end
            else
                mask = SLEAP_data.idx_time > time_ranges_trials(1,j) & ...
                       SLEAP_data.idx_time < time_ranges_trials(3,j);
                % Store as 2×k matrix [x-row; y-row] — identical to original
                filtered_motion{j} = [X_data(mask)'; Y_data(mask)']; %#ok<AGROW>
                good_index = good_index + 1;
            end
        end

        % ----- Pass 2: compute path length per trial (exact copy of original) -----
        path_length_array = [];
        for qq = 1:size(filtered_motion, 2)
            try
                % Transpose 2×k → k×2, identical to original
                coordinates = filtered_motion{1, qq}';
                if isempty(coordinates) || ~isnumeric(coordinates)
                    path_length_array(qq) = 0; %#ok<AGROW>
                    continue;
                end
                distances = pdist2(coordinates, coordinates);
                path_length_array(qq) = sum(diag(distances, 1)); %#ok<AGROW>
            catch
                path_length_array(qq) = 0; %#ok<AGROW>
            end
        end

        % ----- Average by block × reward type (exact copy of original) -----
        if numel(path_length_array) >= height(BehavData)
            b1L = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 1));
            b2L = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 2));
            b3L = mean(path_length_array(1, BehavData.bigSmall == 1.2 & BehavData.Block == 3));
            b1S = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 1));
            b2S = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 2));
            b3S = mean(path_length_array(1, BehavData.bigSmall == 0.3 & BehavData.Block == 3));
            data_matrix(ii, :) = [b1L, b2L, b3L, b1S, b2S, b3S];
        else
            fprintf('  [%s %s: path_length_array (%d) shorter than BehavData (%d) — skipping]\n', ...
                currentanimal, session_to_analyze, numel(path_length_array), height(BehavData));
        end
    end

    loco_tbl          = array2table(data_matrix, 'VariableNames', loco_var_names);
    loco_tbl.animalID = ChrimsonR_IDs;
    loco_tbl.sex      = ChrimsonR_treatment_groups;
    loco_tables{s_idx} = loco_tbl;
end

loco_D1 = loco_tables{1};
loco_D2 = loco_tables{2};

%% Separate by sex
female_loco_D1 = loco_D1(is_female,  :);
female_loco_D2 = loco_D2(is_female,  :);
male_loco_D1   = loco_D1(~is_female, :);
male_loco_D2   = loco_D2(~is_female, :);

fprintf('Loco D1 — females: %d valid | males: %d valid\n', ...
    sum(~isnan(female_loco_D1.b1_large_path)), sum(~isnan(male_loco_D1.b1_large_path)));
fprintf('Loco D2 — females: %d valid | males: %d valid\n', ...
    sum(~isnan(female_loco_D2.b1_large_path)), sum(~isnan(male_loco_D2.b1_large_path)));

%% Figure dimensions and colour scheme  (match behavioral script)
fig_width  = 320;
fig_height = 400;

color_f_d1 = [0.60, 0.88, 0.60];   % light green
color_f_d2 = [0.13, 0.52, 0.13];   % dark green
color_m_d1 = [0.72, 0.72, 0.72];   % light grey
color_m_d2 = [0.28, 0.28, 0.28];   % dark grey
all_colors     = {color_f_d1, color_f_d2, color_m_d1, color_m_d2};
scatter_colors = cellfun(@(c) min(c * 0.75, 1), all_colors, 'UniformOutput', false);
scatter_markers = {'s', 's', 'o', 'o'};  % females = squares, males = circles

%% Figures and 3-way ANOVAs
% Large vs small path length plotted with shared Y-axis so they are comparable.
% Adjust ylim_vals / ytick_vals manually after seeing the data.

vars_large = {'b1_large_path', 'b2_large_path', 'b3_large_path'};
vars_small = {'b1_small_path', 'b2_small_path', 'b3_small_path'};

plot_sex_session_bars(female_loco_D1, female_loco_D2, male_loco_D1, male_loco_D2, ...
    vars_large, 1, all_colors, scatter_colors, scatter_markers, ...
    'Large reward path length (px)', [0 2200], [0:500:2000], fig_width, fig_height);
run_3way_anova(female_loco_D1, female_loco_D2, male_loco_D1, male_loco_D2, ...
    vars_large, 1, 'Large reward path length');

plot_sex_session_bars(female_loco_D1, female_loco_D2, male_loco_D1, male_loco_D2, ...
    vars_small, 1, all_colors, scatter_colors, scatter_markers, ...
    'Small reward path length (px)', [0 2200], [0:500:2000], fig_width, fig_height);
run_3way_anova(female_loco_D1, female_loco_D2, male_loco_D1, male_loco_D2, ...
    vars_small, 1, 'Small reward path length');

%% ================================================================
%  Local functions  (identical interface to generate_behav_figs_RDT_sex_comparison)
%  ================================================================

% ----- Grouped bar plot with scatter overlay -----
function plot_sex_session_bars(f_d1, f_d2, m_d1, m_d2, var_names, scale, ...
    all_colors, scatter_colors, scatter_markers, ylabel_str, ylim_vals, ytick_vals, fig_w, fig_h)
% Bar order per block: Female-D1 | Male-D1 | Female-D2 | Male-D2  (day-grouped)
% Input colour/marker arrays are in sex order {F-D1,F-D2,M-D1,M-D2}; reindexed here.
day_order       = [1, 3, 2, 4];
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
    bar_x(:, b) = (1:ngroups)' - groupwidth/2 + (2*b - 1) * groupwidth / (2*nbars);
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

% ----- 3-way mixed ANOVA: Sex (between) x Day (within) x Block (within) -----
function run_3way_anova(f_d1, f_d2, m_d1, m_d2, var_names, scale, measure_name)

data_fd1 = table2array(f_d1(:, var_names)) * scale;
data_fd2 = table2array(f_d2(:, var_names)) * scale;
data_md1 = table2array(m_d1(:, var_names)) * scale;
data_md2 = table2array(m_d2(:, var_names)) * scale;

nF = size(data_fd1, 1);
nM = size(data_md1, 1);

all_wide   = [data_fd1, data_fd2; data_md1, data_md2];
sex_labels = categorical([repmat({'Female'}, nF, 1); repmat({'Male'}, nM, 1)]);

complete   = all(~isnan(all_wide), 2);
n_complete = sum(complete);
n_excluded = (nF + nM) - n_complete;

tbl     = array2table(all_wide, 'VariableNames', {'D1B1','D1B2','D1B3','D2B1','D2B2','D2B3'});
tbl.Sex = sex_labels;
within  = table( ...
    categorical({'D1';'D1';'D1';'D2';'D2';'D2'}), ...
    categorical({'B1';'B2';'B3';'B1';'B2';'B3'}), ...
    'VariableNames', {'Day','Block'});

fprintf('\n========== 3-Way Mixed ANOVA: %s ==========\n', measure_name);
fprintf('N = %d subjects (%d female, %d male)', n_complete, sum(complete(1:nF)), sum(complete(nF+1:end)));
if n_excluded > 0
    fprintf('  [%d excluded — missing session data]', n_excluded);
end
fprintf('\n');

try
    rm     = fitrm(tbl, 'D1B1,D1B2,D1B3,D2B1,D2B2,D2B3 ~ Sex', 'WithinDesign', within);
    bs_tbl = rm.anova();
    ws_tbl = ranova(rm, 'WithinModel', 'Day*Block');

    ws_rows = ws_tbl.Properties.RowNames;
    bs_rows = bs_tbl.Properties.RowNames;

    fprintf('\n  Between-subjects:\n');
    sex_idx      = find(strcmp(bs_rows, 'Sex'));
    err_candidates = {'Error(Sex)', 'Error'};
    err_idx = [];
    for ci = 1:numel(err_candidates)
        idx = find(strcmp(bs_rows, err_candidates{ci}));
        if ~isempty(idx); err_idx = idx; break; end
    end
    p_sex = NaN;
    if ~isempty(sex_idx) && ~isempty(err_idx)
        print_effect('Sex', bs_tbl.F(sex_idx), bs_tbl.DF(sex_idx), ...
                     bs_tbl.DF(err_idx), bs_tbl.pValue(sex_idx));
        p_sex = bs_tbl.pValue(sex_idx);
    else
        fprintf('    [Sex row not found — bs_tbl rows: %s]\n', strjoin(bs_rows, ', '));
        disp(bs_tbl);
    end

    fprintf('  Within-subjects and interactions:\n');
    effects_ws = { ...
        '(Intercept):Day',       'Day',              'Error(Day)'; ...
        'Sex:Day',               'Sex x Day',         'Error(Day)'; ...
        '(Intercept):Block',     'Block',             'Error(Block)'; ...
        'Sex:Block',             'Sex x Block',        'Error(Block)'; ...
        '(Intercept):Day:Block', 'Day x Block',        'Error(Day:Block)'; ...
        'Sex:Day:Block',         'Sex x Day x Block',  'Error(Day:Block)'};

    for e = 1:size(effects_ws, 1)
        eff_idx = find(strcmp(ws_rows, effects_ws{e,1}));
        err_idx = find(strcmp(ws_rows, effects_ws{e,3}));
        if ~isempty(eff_idx) && ~isempty(err_idx)
            print_effect(effects_ws{e,2}, ws_tbl.F(eff_idx), ...
                         ws_tbl.DF(eff_idx), ws_tbl.DF(err_idx), ws_tbl.pValue(eff_idx));
        end
    end

    % Post-hoc (cascading)
    p_3way      = get_p(ws_tbl, 'Sex:Day:Block');
    p_sex_block = get_p(ws_tbl, 'Sex:Block');
    p_sex_day   = get_p(ws_tbl, 'Sex:Day');
    p_day_block = get_p(ws_tbl, '(Intercept):Day:Block');
    p_block     = get_p(ws_tbl, '(Intercept):Block');
    p_day       = get_p(ws_tbl, '(Intercept):Day');

    ran_posthoc = false;
    if p_3way < 0.05
        fprintf('\n  Significant 3-way interaction — Sex x Day within each Block:\n');
        mc = multcompare(rm, 'Sex', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        fprintf('\n  Day within each Block:\n');
        mc = multcompare(rm, 'Day', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if p_sex_block < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex x Block — Sex within each Block:\n');
        mc = multcompare(rm, 'Sex', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if p_sex_day < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex x Day — Sex within each Day:\n');
        mc = multcompare(rm, 'Sex', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if p_day_block < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Day x Block — Day within each Block:\n');
        mc = multcompare(rm, 'Day', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if p_block < 0.05 && p_sex_block >= 0.05 && p_day_block >= 0.05 && p_3way >= 0.05
        fprintf('\n  Significant Block main effect — pairwise block comparisons:\n');
        mc = multcompare(rm, 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
    end
    if p_day < 0.05 && p_sex_day >= 0.05 && p_day_block >= 0.05 && p_3way >= 0.05
        fprintf('\n  Significant Day main effect — D1 vs D2:\n');
        mc = multcompare(rm, 'Day', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
    end

    % Between-subjects Sex main effect — post-hoc if no interaction covered it
    if ~isnan(p_sex) && p_sex < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex main effect — Female vs Male (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
    end

catch ME
    fprintf('  [ANOVA failed: %s]\n', ME.message);
end
end

% ----- Helper: extract p-value from ranova table by row name -----
function p = get_p(tbl, row_name)
    try
        p = tbl{row_name, 'pValue'};
    catch
        rows = tbl.Properties.RowNames;
        idx  = find(strcmp(rows, row_name));
        if ~isempty(idx); p = tbl.pValue(idx(1)); else; p = NaN; end
    end
end

% ----- Helper: print one ANOVA effect line -----
function print_effect(label, F, df1, df2, p)
    fprintf('    %-26s F(%d,%d) = %7.3f,  p = %.4f%s\n', label, df1, df2, F, p, sig_stars(p));
end

% ----- Helper: print multcompare output -----
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
    for r = 1:height(mc)
        p_val = mc.pValue(r);
        lv1   = char(mc{r,1});
        lv2   = char(mc{r,2});
        by_str = '';
        if width(mc) >= 3
            v3 = mc{r,3};
            if iscategorical(v3) || ischar(v3) || iscell(v3)
                by_str = [' [' char(v3) ']'];
            end
        end
        fprintf('      %-42s p = %.4f%s\n', [lv1 ' vs ' lv2 by_str], p_val, sig_stars(p_val));
    end
end

% ----- Helper: significance stars -----
function s = sig_stars(p)
    if p < 0.001;    s = ' ***';
    elseif p < 0.01; s = ' **';
    elseif p < 0.05; s = ' *';
    else;            s = '';
    end
end
