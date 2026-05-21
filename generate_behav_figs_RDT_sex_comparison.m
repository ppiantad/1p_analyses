% generate_behav_figs_RDT_sex_comparison.m
% Bar graphs comparing male vs female behavior on RDT D1 and RDT D2
% X-axis: Block 1, 2, 3 | Within each block: Female-D1, Male-D1, Female-D2, Male-D2
% Individual scatter points overlaid on each bar
%
% Requires in workspace (run females_vs_males_groups_correct.m first):
%   ChrimsonR_IDs, ChrimsonR_treatment_groups, valid_sessions,
%   session_weight_raw_grams, freefeed_weight, shock_sens,
%   pretraining_sessions, sessions_to_rm_criterion
% Requires: final_behavior (load new_females_and_males_behavior_05012026.mat)
clc

%% Use animal list and sex labels from groups file
n_animals  = numel(ChrimsonR_IDs);
is_female  = strcmp(ChrimsonR_treatment_groups, 'Female');
fprintf('Including %d females and %d males\n', sum(is_female), sum(~is_female));

%% Build per-animal metadata table (for scatter plot covariates)
pretraining_sessions_sum = cellfun(@sum, pretraining_sessions);

% Sum all 3 staircase values per animal (NaN if missing)
shock_level_sum = nan(n_animals, 1);
for ii = 1:n_animals
    s = shock_sens{ii};
    if isnumeric(s) && ~all(isnan(s(:)))
        shock_level_sum(ii) = sum(s(:));
    end
end

% Session weights as % of free-feed body weight
session_weight_mat = cell2mat(session_weight_raw_grams');  % n_animals x 2
% weight_pct_D1 = (session_weight_mat(:,1) ./ freefeed_weight') * 100;
% weight_pct_D2 = (session_weight_mat(:,2) ./ freefeed_weight') * 100;

weight_pct_D1 = [(session_weight_mat(:,1) ./ freefeed_weight') * 100]+5;
weight_pct_D2 = [(session_weight_mat(:,2) ./ freefeed_weight') * 100]+5;

metadata_table = table(ChrimsonR_IDs, ChrimsonR_treatment_groups, ...
    weight_pct_D1, weight_pct_D2, shock_level_sum, ...
    pretraining_sessions_sum', sessions_to_rm_criterion', ...
    'VariableNames', {'animalID','sex','weight_pct_D1','weight_pct_D2', ...
                      'shock_level_sum','pretraining_sessions','sessions_to_criterion'});

%% Define behavioral variable names (46 total)
variable_names = { ...
    'block_1_large','block_2_large','block_3_large', ...
    'block_1_small','block_2_small','block_3_small', ...
    'block_1_shocks','block_2_shocks','block_3_shocks', ...
    'block_1_omissions','block_2_omissions','block_3_omissions', ...
    'win_stay_b1','win_stay_b2','win_stay_b3', ...
    'lose_shift_b1','lose_shift_b2','lose_shift_b3', ...
    'choice_lat_all_b1','choice_lat_all_b2','choice_lat_all_b3', ...
    'choice_lat_large_b1','choice_lat_large_b2','choice_lat_large_b3', ...
    'choice_lat_small_b1','choice_lat_small_b2','choice_lat_small_b3', ...
    'collect_lat_large_b1','collect_lat_large_b2','collect_lat_large_b3', ...
    'collect_lat_large_shock_b1','collect_lat_large_shock_b2','collect_lat_large_shock_b3', ...
    'collect_lat_large_noshock_b1','collect_lat_large_noshock_b2','collect_lat_large_noshock_b3', ...
    'collect_lat_small_b1','collect_lat_small_b2','collect_lat_small_b3', ...
    'consum_dur_large_b1','consum_dur_large_b2','consum_dur_large_b3', ...
    'consum_dur_small_b1','consum_dur_small_b2','consum_dur_small_b3', ...
    'large_aborts_b1','large_aborts_b2','large_aborts_b3', ...
    'small_aborts_b1','small_aborts_b2','small_aborts_b3', ...
    'trials_completed', ...
    'blank_touch_b1','blank_touch_b2','blank_touch_b3', ...
    'lose_omit_b1','lose_omit_b2','lose_omit_b3', ...
    'lose_stay_b1','lose_stay_b2','lose_stay_b3', ...
    'win_shift_b1','win_shift_b2','win_shift_b3', ...
    'win_omit_b1','win_omit_b2','win_omit_b3'};

n_vars = numel(variable_names);  % 61

%% Extract behavioral data for RDT_D1 and RDT_D2
sessions    = {'RDT_D1', 'RDT_D2'};
risk_tables = cell(1, 2);

for s_idx = 1:2
    session_to_analyze = sessions{s_idx};

    % Which animals have this session listed as valid
    valid_for_session = cellfun(@(s) any(strcmp(s, session_to_analyze)), valid_sessions);

    data_matrix = nan(n_animals, n_vars);

    for ii = 1:n_animals
        % Skip if this session is not valid for this animal
        if ~valid_for_session(ii)
            continue;
        end

        currentanimal = ChrimsonR_IDs{ii};
        if ~isfield(final_behavior, currentanimal) || ...
           ~isfield(final_behavior.(currentanimal), session_to_analyze)
            continue;
        end

        BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;

        % --- trial_after_shk ---
        if ~ismember('trial_after_shk', BehavData.Properties.VariableNames)
            BehavData.trial_after_shk = zeros(height(BehavData), 1);
        end
        for row = 1:height(BehavData)
            if BehavData.shock(row) == 1
                kk = 1;
                while true
                    if (row + kk) > height(BehavData); break; end
                    if ~isnan(BehavData.bigSmall(row+kk)) && BehavData.ForceFree(row+kk) ~= 999
                        BehavData.trial_after_shk(row+kk) = 1;
                        break;
                    end
                    kk = kk + 1;
                end
            end
        end

        [BehavData, ~, ~] = TrialFilter_test(BehavData, 'ALL', 1);

        % --- % large / small choice (free-choice trials only) ---
        fc = BehavData.ForceFree == 0;
        ch = (BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3);
        bl1L = sum(BehavData.bigSmall==1.2 & BehavData.Block==1 & fc) / sum(ch & fc & BehavData.Block==1);
        bl2L = sum(BehavData.bigSmall==1.2 & BehavData.Block==2 & fc) / sum(ch & fc & BehavData.Block==2);
        bl3L = sum(BehavData.bigSmall==1.2 & BehavData.Block==3 & fc) / sum(ch & fc & BehavData.Block==3);
        bl1S = sum(BehavData.bigSmall==0.3 & BehavData.Block==1 & fc) / sum(ch & fc & BehavData.Block==1);
        bl2S = sum(BehavData.bigSmall==0.3 & BehavData.Block==2 & fc) / sum(ch & fc & BehavData.Block==2);
        bl3S = sum(BehavData.bigSmall==0.3 & BehavData.Block==3 & fc) / sum(ch & fc & BehavData.Block==3);

        % --- shocks and omissions ---
        shk1  = sum(BehavData.shock==1       & BehavData.Block==1);
        shk2  = sum(BehavData.shock==1       & BehavData.Block==2);
        shk3  = sum(BehavData.shock==1       & BehavData.Block==3);
        omit1 = sum(BehavData.omissionALL==1 & BehavData.Block==1);
        omit2 = sum(BehavData.omissionALL==1 & BehavData.Block==2);
        omit3 = sum(BehavData.omissionALL==1 & BehavData.Block==3);

        % --- win-stay / lose-shift by block ---
        ws1 = sum(BehavData.win_stay==1   & BehavData.Block==1) / sum(BehavData.WL==1 & fc & BehavData.Block==1);
        ws2 = sum(BehavData.win_stay==1   & BehavData.Block==2) / sum(BehavData.WL==1 & fc & BehavData.Block==2);
        ws3 = sum(BehavData.win_stay==1   & BehavData.Block==3) / sum(BehavData.WL==1 & fc & BehavData.Block==3);
        ls1 = 0;  % no losses in Block 1 by design
        ls2 = sum(BehavData.lose_shift==1 & BehavData.Block==2) / sum(BehavData.WL==3 & fc & BehavData.Block==2);
        ls3 = sum(BehavData.lose_shift==1 & BehavData.Block==3) / sum(BehavData.WL==3 & fc & BehavData.Block==3);

        % --- latencies ---
        BehavData.choice_latency  = BehavData.choiceTime    - BehavData.stTime;
        BehavData.collect_latency = BehavData.collectionTime - BehavData.choiceTime;

        cl_all1 = nanmean(BehavData.choice_latency(BehavData.Block==1 & ch));
        cl_all2 = nanmean(BehavData.choice_latency(BehavData.Block==2 & ch));
        cl_all3 = nanmean(BehavData.choice_latency(BehavData.Block==3 & ch));
        cl_L1   = nanmean(BehavData.choice_latency(BehavData.bigSmall==1.2 & BehavData.Block==1));
        cl_L2   = nanmean(BehavData.choice_latency(BehavData.bigSmall==1.2 & BehavData.Block==2));
        cl_L3   = nanmean(BehavData.choice_latency(BehavData.bigSmall==1.2 & BehavData.Block==3));
        cl_S1   = nanmean(BehavData.choice_latency(BehavData.bigSmall==0.3 & BehavData.Block==1));
        cl_S2   = nanmean(BehavData.choice_latency(BehavData.bigSmall==0.3 & BehavData.Block==2));
        cl_S3   = nanmean(BehavData.choice_latency(BehavData.bigSmall==0.3 & BehavData.Block==3));

        col_L1 = nanmean(BehavData.collect_latency(BehavData.bigSmall==1.2 & BehavData.Block==1));
        col_L2 = nanmean(BehavData.collect_latency(BehavData.bigSmall==1.2 & BehavData.Block==2));
        col_L3 = nanmean(BehavData.collect_latency(BehavData.bigSmall==1.2 & BehavData.Block==3));

        % --- collect latency: large reward, shock trials only ---
        col_L_shk1 = nanmean(BehavData.collect_latency(BehavData.bigSmall==1.2 & BehavData.shock==1 & BehavData.Block==1));
        col_L_shk2 = nanmean(BehavData.collect_latency(BehavData.bigSmall==1.2 & BehavData.shock==1 & BehavData.Block==2));
        col_L_shk3 = nanmean(BehavData.collect_latency(BehavData.bigSmall==1.2 & BehavData.shock==1 & BehavData.Block==3));

        % --- collect latency: large reward, no-shock trials only ---
        col_L_noshk1 = nanmean(BehavData.collect_latency(BehavData.bigSmall==1.2 & BehavData.shock==0 & BehavData.Block==1));
        col_L_noshk2 = nanmean(BehavData.collect_latency(BehavData.bigSmall==1.2 & BehavData.shock==0 & BehavData.Block==2));
        col_L_noshk3 = nanmean(BehavData.collect_latency(BehavData.bigSmall==1.2 & BehavData.shock==0 & BehavData.Block==3));

        col_S1 = nanmean(BehavData.collect_latency(BehavData.bigSmall==0.3 & BehavData.Block==1));
        col_S2 = nanmean(BehavData.collect_latency(BehavData.bigSmall==0.3 & BehavData.Block==2));
        col_S3 = nanmean(BehavData.collect_latency(BehavData.bigSmall==0.3 & BehavData.Block==3));

        % --- consumption duration ---
        if ismember('collectionTime_end', BehavData.Properties.VariableNames)
            BehavData.consum_dur = BehavData.collectionTime_end - BehavData.collectionTime;
            cd_L1 = nanmean(BehavData.consum_dur(BehavData.bigSmall==1.2 & BehavData.Block==1));
            cd_L2 = nanmean(BehavData.consum_dur(BehavData.bigSmall==1.2 & BehavData.Block==2));
            cd_L3 = nanmean(BehavData.consum_dur(BehavData.bigSmall==1.2 & BehavData.Block==3));
            cd_S1 = nanmean(BehavData.consum_dur(BehavData.bigSmall==0.3 & BehavData.Block==1));
            cd_S2 = nanmean(BehavData.consum_dur(BehavData.bigSmall==0.3 & BehavData.Block==2));
            cd_S3 = nanmean(BehavData.consum_dur(BehavData.bigSmall==0.3 & BehavData.Block==3));
        else
            [cd_L1, cd_L2, cd_L3, cd_S1, cd_S2, cd_S3] = deal(NaN);
        end

        % --- aborts ---
        if ismember('type_binary', BehavData.Properties.VariableNames)
            ab_L1 = sum(BehavData.type_binary==1 & BehavData.Block==1);
            ab_L2 = sum(BehavData.type_binary==1 & BehavData.Block==2);
            ab_L3 = sum(BehavData.type_binary==1 & BehavData.Block==3);
            ab_S1 = sum(BehavData.type_binary==2 & BehavData.Block==1);
            ab_S2 = sum(BehavData.type_binary==2 & BehavData.Block==2);
            ab_S3 = sum(BehavData.type_binary==2 & BehavData.Block==3);
        else
            [ab_L1, ab_L2, ab_L3, ab_S1, ab_S2, ab_S3] = deal(0);
        end

        % --- blank touches ---
        has_blank_touch = ismember('Blank_Touch', BehavData.Properties.VariableNames);
        if has_blank_touch
            bt1 = sum(BehavData.Blank_Touch == 1 & BehavData.Block == 1);
            bt2 = sum(BehavData.Blank_Touch == 1 & BehavData.Block == 2);
            bt3 = sum(BehavData.Blank_Touch == 1 & BehavData.Block == 3);
        else
            [bt1, bt2, bt3] = deal(NaN);
        end

        % --- lose-omit / lose-stay (proportion of loss trials) ---
        % Denominator = free-choice loss trials per block (WL==3 & fc).
        % Block 1 has no losses by design → hardcode to 0 (not NaN) so stats work.
        n_loss_b2 = sum(BehavData.WL==3 & fc & BehavData.Block==2);
        n_loss_b3 = sum(BehavData.WL==3 & fc & BehavData.Block==3);

        lo1  = 0;
        lo2  = sum(BehavData.lose_omit==1 & BehavData.Block==2) / n_loss_b2;
        lo3  = sum(BehavData.lose_omit==1 & BehavData.Block==3) / n_loss_b3;

        lsy1 = 0;
        lsy2 = sum(BehavData.lose_stay==1 & BehavData.Block==2) / n_loss_b2;
        lsy3 = sum(BehavData.lose_stay==1 & BehavData.Block==3) / n_loss_b3;

        % --- win-shift / win-omit ---
        % For each win trial, advance past blank touches and aborts to find the next
        % real trial, then check whether the animal shifted to small (win-shift)
        % or omitted (win-omit).  Denominator = total wins per block (same as win-stay).
        has_type_binary = ismember('type_binary', BehavData.Properties.VariableNames);
        n_rows_bd = height(BehavData);
        wsh_vals = nan(1, 3);
        wo_vals  = nan(1, 3);
        for blk_w = 1:3
            win_rows = find(BehavData.WL==1 & fc & BehavData.Block==blk_w);
            n_wins = numel(win_rows);
            if n_wins == 0; continue; end
            n_wsh = 0;  n_wo = 0;
            for wr = 1:n_wins
                nr = win_rows(wr) + 1;
                % advance past blank touches and aborts
                while nr <= n_rows_bd
                    is_blank = has_blank_touch && BehavData.Blank_Touch(nr) == 1;
                    is_abort = has_type_binary  && (BehavData.type_binary(nr)==1 || BehavData.type_binary(nr)==2);
                    if ~is_blank && ~is_abort; break; end
                    nr = nr + 1;
                end
                if nr > n_rows_bd; continue; end
                if BehavData.bigSmall(nr) == 0.3
                    n_wsh = n_wsh + 1;
                elseif BehavData.omissionALL(nr) == 1
                    n_wo = n_wo + 1;
                end
            end
            wsh_vals(blk_w) = n_wsh / n_wins;
            wo_vals(blk_w)  = n_wo  / n_wins;
        end
        wsh1 = wsh_vals(1); wsh2 = wsh_vals(2); wsh3 = wsh_vals(3);
        wo1  = wo_vals(1);  wo2  = wo_vals(2);  wo3  = wo_vals(3);

        trials_completed = sum(ch);

        data_matrix(ii, :) = [ ...
            bl1L, bl2L, bl3L, bl1S, bl2S, bl3S, ...
            shk1, shk2, shk3, omit1, omit2, omit3, ...
            ws1, ws2, ws3, ls1, ls2, ls3, ...
            cl_all1, cl_all2, cl_all3, ...
            cl_L1, cl_L2, cl_L3, cl_S1, cl_S2, cl_S3, ...
            col_L1, col_L2, col_L3, ...
            col_L_shk1, col_L_shk2, col_L_shk3, ...
            col_L_noshk1, col_L_noshk2, col_L_noshk3, ...
            col_S1, col_S2, col_S3, ...
            cd_L1, cd_L2, cd_L3, cd_S1, cd_S2, cd_S3, ...
            ab_L1, ab_L2, ab_L3, ab_S1, ab_S2, ab_S3, ...
            trials_completed, bt1, bt2, bt3, lo1, lo2, lo3, lsy1, lsy2, lsy3, ...
            wsh1, wsh2, wsh3, wo1, wo2, wo3];
    end

    risk_table = array2table(data_matrix, 'VariableNames', variable_names);
    risk_table.animalID    = ChrimsonR_IDs;
    risk_table.sex         = ChrimsonR_treatment_groups;
    risk_table.valid_D1    = cellfun(@(s) any(strcmp(s,'RDT_D1')), valid_sessions);
    risk_table.valid_D2    = cellfun(@(s) any(strcmp(s,'RDT_D2')), valid_sessions);
    % Attach metadata for scatter plot covariates
    % weight_pct_D1/D2 already computed once at the top — just index by session
    weight_pcts_all = {weight_pct_D1, weight_pct_D2};
    risk_table.weight_pct       = weight_pcts_all{s_idx};
    risk_table.session_weight_g = session_weight_mat(:, s_idx);
    risk_table.freefeed_weight_g = freefeed_weight';
    risk_table.shock_level      = shock_level_sum;
    risk_table.sessions_crit    = sessions_to_rm_criterion';
    risk_table.pretrain_total   = pretraining_sessions_sum';

    risk_tables{s_idx} = risk_table;
end

risk_table_D1 = risk_tables{1};
risk_table_D2 = risk_tables{2};

%% Body-weight inclusion threshold
% Animals below this % of free-feed weight on EITHER RDT session are excluded.
% Set to 0 to include all animals regardless of weight.
bw_threshold = 70;   % % of free-feed weight
%% Body-weight filter

% Exclude animals whose session weight is below bw_threshold % of free-feed
% weight on EITHER RDT session. NaN weight (session not run) does not trigger
% exclusion for that session.
bw_ok_D1   = isnan(weight_pct_D1) | weight_pct_D1 >= bw_threshold;
bw_ok_D2   = isnan(weight_pct_D2) | weight_pct_D2 >= bw_threshold;
bw_include = bw_ok_D1 & bw_ok_D2;   % n_animals × 1 logical

n_excl_bw = sum(~bw_include);
if n_excl_bw > 0
    excl_ids = ChrimsonR_IDs(~bw_include);
    fprintf('BW filter (>= %.0f%% free-feed): excluding %d animal(s): %s\n', ...
        bw_threshold, n_excl_bw, strjoin(excl_ids, ', '));
else
    fprintf('BW filter (>= %.0f%% free-feed): no animals excluded\n', bw_threshold);
end

%% Separate by sex (BW filter applied)
female_D1 = risk_table_D1(is_female  & bw_include, :);
female_D2 = risk_table_D2(is_female  & bw_include, :);
male_D1   = risk_table_D1(~is_female & bw_include, :);
male_D2   = risk_table_D2(~is_female & bw_include, :);

fprintf('D1 females: %d valid | D2 females: %d valid\n', ...
    sum(~isnan(female_D1.block_1_large)), sum(~isnan(female_D2.block_1_large)));
fprintf('D1 males:   %d valid | D2 males:   %d valid\n', ...
    sum(~isnan(male_D1.block_1_large)),   sum(~isnan(male_D2.block_1_large)));

%% Figure dimensions  — adjust here to resize all figures at once
fig_width  = 320;   % pixels  (original 560; set to 280 for half-width)
fig_height = 400;   % pixels
scatter_fig_size = 350;   % pixels — scatter plots are square



%% Color scheme  (females = greens, D1 light / D2 dark; males = greys, D1 light / D2 dark)
color_f_d1 = [0.60, 0.88, 0.60];   % light green
color_f_d2 = [0.13, 0.52, 0.13];   % dark green
color_m_d1 = [0.72, 0.72, 0.72];   % light grey
color_m_d2 = [0.28, 0.28, 0.28];   % dark grey
all_colors     = {color_f_d1, color_f_d2, color_m_d1, color_m_d2};
scatter_colors = cellfun(@(c) min(c * 0.75, 1), all_colors, 'UniformOutput', false);
% Females = squares ('s'), Males = circles ('o')
scatter_markers = {'s', 's', 'o', 'o'};

%% Generate figures and run 3-way ANOVAs
% For each call: plot_sex_session_bars(..., ylim_vals, ytick_vals, fig_w, fig_h)
%   ylim_vals  — [lo hi] or [] for auto
%   ytick_vals — e.g. [0 25 50 75 100] or [] for auto
pct = 100;  % scale proportions to %
raw = 1;

% --- % choice ---
vars = {'block_1_large','block_2_large','block_3_large'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, pct, ...
    all_colors, scatter_colors, scatter_markers, '% Large choice', ...
    [0 100], [0 25 50 75 100], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, pct, '% Large choice');

vars = {'block_1_small','block_2_small','block_3_small'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, pct, ...
    all_colors, scatter_colors, scatter_markers, '% Small choice', ...
    [0 100], [0 25 50 75 100], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, pct, '% Small choice');

% --- Shocks and omissions ---
vars = {'block_1_shocks','block_2_shocks','block_3_shocks'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Shocks per block', ...
    [], [], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Shocks per block');

vars = {'block_1_omissions','block_2_omissions','block_3_omissions'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Omissions per block', ...
    [], [], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Omissions per block');

% --- Choice latency: all ---
vars = {'choice_lat_all_b1','choice_lat_all_b2','choice_lat_all_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Choice latency - all (s)', ...
    [0 20], [0:5:20], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Choice latency - all');

% --- Choice latency: large vs small ---
vars = {'choice_lat_large_b1','choice_lat_large_b2','choice_lat_large_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Large choice latency (s)', ...
    [0 20], [0:5:20], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Large choice latency');

vars = {'choice_lat_small_b1','choice_lat_small_b2','choice_lat_small_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Small choice latency (s)', ...
    [0 20], [0:5:20], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Small choice latency');

% --- Collect latency: large vs small ---
vars = {'collect_lat_large_b1','collect_lat_large_b2','collect_lat_large_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Large collect latency — all trials (s)', ...
    [0 6], [0:2:6], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Large collect latency — all trials');

vars = {'collect_lat_large_shock_b1','collect_lat_large_shock_b2','collect_lat_large_shock_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Large collect latency — shock trials (s)', ...
    [0 6], [0:2:6], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Large collect latency — shock trials');

vars = {'collect_lat_large_noshock_b1','collect_lat_large_noshock_b2','collect_lat_large_noshock_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Large collect latency — no-shock trials (s)', ...
    [0 6], [0:2:6], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Large collect latency — no-shock trials');

vars = {'collect_lat_small_b1','collect_lat_small_b2','collect_lat_small_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Small collect latency (s)', ...
    [0 6], [0:2:6], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Small collect latency');

% --- Consumption duration: large vs small ---
vars = {'consum_dur_large_b1','consum_dur_large_b2','consum_dur_large_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Large consumption duration (s)', ...
    [0 6], [0:2:6], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Large consumption duration');

vars = {'consum_dur_small_b1','consum_dur_small_b2','consum_dur_small_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Small consumption duration (s)', ...
    [0 6], [0:2:6], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Small consumption duration');

% --- Strategy analyses (all as %; lose metrics = 0 in Block 1 by design) ---
vars = {'win_stay_b1','win_stay_b2','win_stay_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, pct, ...
    all_colors, scatter_colors, scatter_markers, '% win-stay', ...
    [0 100], [0 25 50 75 100], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, pct, '% win-stay');

vars = {'win_shift_b1','win_shift_b2','win_shift_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, pct, ...
    all_colors, scatter_colors, scatter_markers, '% win-shift', ...
    [0 100], [0 25 50 75 100], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, pct, '% win-shift');

vars = {'win_omit_b1','win_omit_b2','win_omit_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, pct, ...
    all_colors, scatter_colors, scatter_markers, '% win-omit', ...
    [0 100], [0 25 50 75 100], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, pct, '% win-omit');

vars = {'lose_shift_b1','lose_shift_b2','lose_shift_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, pct, ...
    all_colors, scatter_colors, scatter_markers, '% lose-shift', ...
    [0 100], [0 25 50 75 100], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, pct, '% lose-shift');

vars = {'lose_omit_b1','lose_omit_b2','lose_omit_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, pct, ...
    all_colors, scatter_colors, scatter_markers, '% lose-omit', ...
    [0 100], [0 25 50 75 100], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, pct, '% lose-omit');

vars = {'lose_stay_b1','lose_stay_b2','lose_stay_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, pct, ...
    all_colors, scatter_colors, scatter_markers, '% lose-stay', ...
    [0 100], [0 25 50 75 100], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, pct, '% lose-stay');

% --- Aborts: large vs small ---
vars = {'large_aborts_b1','large_aborts_b2','large_aborts_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Large aborts per block', ...
    [0 160], [0:40:160], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Large aborts');

vars = {'small_aborts_b1','small_aborts_b2','small_aborts_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Small aborts per block', ...
    [0 160], [0:40:160], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Small aborts');

% --- Blank touches per block ---
vars = {'blank_touch_b1','blank_touch_b2','blank_touch_b3'};
plot_sex_session_bars(female_D1, female_D2, male_D1, male_D2, vars, raw, ...
    all_colors, scatter_colors, scatter_markers, 'Blank touches per block', ...
    [], [], fig_width, fig_height);
run_3way_anova(female_D1, female_D2, male_D1, male_D2, vars, raw, 'Blank touches per block');

% --- Session weight and % free-feed weight (no block variable) ---
f_mask = is_female  & bw_include;
m_mask = ~is_female & bw_include;

f_weight_d1 = session_weight_mat(f_mask, 1);
f_weight_d2 = session_weight_mat(f_mask, 2);
m_weight_d1 = session_weight_mat(m_mask, 1);
m_weight_d2 = session_weight_mat(m_mask, 2);

% Index into the session-level arrays computed at the top — no repeated division
f_pct_d1 = weight_pct_D1(f_mask);
f_pct_d2 = weight_pct_D2(f_mask);
m_pct_d1 = weight_pct_D1(m_mask);
m_pct_d2 = weight_pct_D2(m_mask);

plot_weight_bars(f_weight_d1, f_weight_d2, m_weight_d1, m_weight_d2, ...
    all_colors, scatter_colors, scatter_markers, 'Session weight (g)', ...
    [], [], 200, 350);
run_2way_anova_weight(f_weight_d1, f_weight_d2, m_weight_d1, m_weight_d2, 'Session weight');

plot_weight_bars(f_pct_d1, f_pct_d2, m_pct_d1, m_pct_d2, ...
    all_colors, scatter_colors, scatter_markers, '% free-feed weight', ...
    [], [], 200, 350);
run_2way_anova_weight(f_pct_d1, f_pct_d2, m_pct_d1, m_pct_d2, '% free-feed weight');

% --- Session completion donut charts (RDT D1 and D2) ---
% Colors in day-grouped order: F-D1, M-D1, F-D2, M-D2
donut_colors = {color_f_d1, color_m_d1, color_f_d2, color_m_d2};
plot_completion_donut(female_D1, female_D2, male_D1, male_D2, donut_colors, fig_width, fig_height);

%% Scatter plots with regression lines
% Each call produces two figures: one for RDT D1, one for RDT D2.
% Edit variable extractions and labels to add or swap pairs.
% Females = squares (dark shade of session color), Males = circles.

% Convenience: mean riskiness (blocks 2-3 large choice, %)
% Block 1 is always safe so is often excluded from a riskiness summary.
mean_risk_f_d1 = nanmean([female_D1.block_2_large, female_D1.block_3_large], 2) * pct;
mean_risk_m_d1 = nanmean([male_D1.block_2_large,   male_D1.block_3_large],   2) * pct;
mean_risk_f_d2 = nanmean([female_D2.block_2_large, female_D2.block_3_large], 2) * pct;
mean_risk_m_d2 = nanmean([male_D2.block_2_large,   male_D2.block_3_large],   2) * pct;

% ---- Mean riskiness vs session weight % ----
plot_scatter_regression( ...
    mean_risk_f_d1, mean_risk_m_d1, female_D1.weight_pct, male_D1.weight_pct, ...
    '% Large choice (blocks 2-3)', '% Free-feed weight', 'RDT D1', ...
    color_f_d1, color_m_d1, scatter_fig_size, scatter_fig_size, [20 90], [20:10:90], [60 90], [60:10:90]);
plot_scatter_regression( ...
    mean_risk_f_d2, mean_risk_m_d2, female_D2.weight_pct, male_D2.weight_pct, ...
    '% Large choice (blocks 2-3)', '% Free-feed weight', 'RDT D2', ...
    color_f_d2, color_m_d2, scatter_fig_size, scatter_fig_size, [20 90], [20:10:90], [60 90], [60:10:90]);

% ---- Mean riskiness vs session weight (raw grams) ----
plot_scatter_regression( ...
    mean_risk_f_d1, mean_risk_m_d1, female_D1.session_weight_g, male_D1.session_weight_g, ...
    '% Large choice (blocks 2-3)', 'Session weight (g)', 'RDT D1', ...
    color_f_d1, color_m_d1, scatter_fig_size, scatter_fig_size, [20 90], [20:10:90], [15 30], [15:5:30]);
plot_scatter_regression( ...
    mean_risk_f_d2, mean_risk_m_d2, female_D2.session_weight_g, male_D2.session_weight_g, ...
    '% Large choice (blocks 2-3)', 'Session weight (g)', 'RDT D2', ...
    color_f_d2, color_m_d2, scatter_fig_size, scatter_fig_size, [20 90], [20:10:90], [15 30], [15:5:30]);

% ---- Mean riskiness vs free-feed weight (raw grams) ----
% Free-feed weight is animal-level (same value in D1 and D2 tables).
plot_scatter_regression( ...
    mean_risk_f_d1, mean_risk_m_d1, female_D1.freefeed_weight_g, male_D1.freefeed_weight_g, ...
    '% Large choice (blocks 2-3)', 'Free-feed weight (g)', 'RDT D1', ...
    color_f_d1, color_m_d1, scatter_fig_size, scatter_fig_size, [20 90], [20:10:90], [20 40], [20:5:40]);
plot_scatter_regression( ...
    mean_risk_f_d2, mean_risk_m_d2, female_D2.freefeed_weight_g, male_D2.freefeed_weight_g, ...
    '% Large choice (blocks 2-3)', 'Free-feed weight (g)', 'RDT D2', ...
    color_f_d2, color_m_d2, scatter_fig_size, scatter_fig_size, [20 90], [20:10:90], [20 40], [20:5:40]);

% ---- Mean riskiness vs shock sensitivity ----
% shock_level is animal-level; identical in D1 and D2 tables
plot_scatter_regression( ...
    mean_risk_f_d1, mean_risk_m_d1, female_D1.shock_level, male_D1.shock_level, ...
    '% Large choice (blocks 2-3)', 'Shock sensitivity (mA)', 'RDT D1', ...
    color_f_d1, color_m_d1, scatter_fig_size, scatter_fig_size, [], [], [], []);
plot_scatter_regression( ...
    mean_risk_f_d2, mean_risk_m_d2, female_D2.shock_level, male_D2.shock_level, ...
    '% Large choice (blocks 2-3)', 'Shock sensitivity (mA)', 'RDT D2', ...
    color_f_d2, color_m_d2, scatter_fig_size, scatter_fig_size, [], [], [], []);

% ---- Mean riskiness vs sessions to criterion ----
plot_scatter_regression( ...
    mean_risk_f_d1, mean_risk_m_d1, female_D1.sessions_crit, male_D1.sessions_crit, ...
    '% Large choice (blocks 2-3)', 'Sessions to criterion', 'RDT D1', ...
    color_f_d1, color_m_d1, scatter_fig_size, scatter_fig_size, [], [], [], []);
plot_scatter_regression( ...
    mean_risk_f_d2, mean_risk_m_d2, female_D2.sessions_crit, male_D2.sessions_crit, ...
    '% Large choice (blocks 2-3)', 'Sessions to criterion', 'RDT D2', ...
    color_f_d2, color_m_d2, scatter_fig_size, scatter_fig_size, [], [], [], []);

% ---- D1 vs D2 riskiness (cross-session consistency) ----
% Animals missing D2 are excluded automatically via NaN filtering.
plot_scatter_regression( ...
    mean_risk_f_d1, mean_risk_m_d1, mean_risk_f_d2, mean_risk_m_d2, ...
    '% Large choice D1 (blocks 2-3)', '% Large choice D2 (blocks 2-3)', 'D1 vs D2', ...
    color_f_d2, color_m_d2, scatter_fig_size, scatter_fig_size,  [20 90], [20:10:90], [20 90], [20:10:90]);

%% ================================================================
%  Local functions
%  ================================================================

% ----- Grouped bar plot with scatter overlay -----
function plot_sex_session_bars(f_d1, f_d2, m_d1, m_d2, var_names, scale, ...
    all_colors, scatter_colors, scatter_markers, ylabel_str, ylim_vals, ytick_vals, fig_w, fig_h)
% Bar order per block: Female-D1 | Male-D1 | Female-D2 | Male-D2  (day-grouped)
% Females = squares, Males = circles
% Input color/marker arrays are in sex order {F-D1,F-D2,M-D1,M-D2}; reindex to day order here.
day_order = [1, 3, 2, 4];  % F-D1, M-D1, F-D2, M-D2
all_colors     = all_colors(day_order);
scatter_colors = scatter_colors(day_order);
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
if ~isempty(ylim_vals)
    ylim(ylim_vals);
end
if ~isempty(ytick_vals)
    yticks(ytick_vals);
end
legend(bar_h, {'Female D1','Male D1','Female D2','Male D2'}, ...
    'Location', 'best', 'FontSize', 9);
box off;

% Fix the inner axes area to a constant pixel size regardless of label widths.
% Adjust these margins (px) if labels get clipped.
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

data_fd1 = table2array(f_d1(:, var_names)) * scale;  % nF x 3
data_fd2 = table2array(f_d2(:, var_names)) * scale;
data_md1 = table2array(m_d1(:, var_names)) * scale;  % nM x 3
data_md2 = table2array(m_d2(:, var_names)) * scale;

nF = size(data_fd1, 1);
nM = size(data_md1, 1);

% Wide format: [D1B1 D1B2 D1B3 D2B1 D2B2 D2B3] — one row per subject
all_wide   = [data_fd1, data_fd2; data_md1, data_md2];
sex_labels = categorical([repmat({'Female'}, nF, 1); repmat({'Male'}, nM, 1)]);

% Apply listwise deletion before fitrm so both Sex levels are always explicit.
complete     = all(~isnan(all_wide), 2);
n_complete   = sum(complete);
n_excluded   = (nF + nM) - n_complete;
nF_c         = sum(complete(1:nF));
nM_c         = sum(complete(nF+1:end));

all_wide_c   = all_wide(complete, :);
sex_labels_c = sex_labels(complete);

tbl     = array2table(all_wide_c, 'VariableNames', {'D1B1','D1B2','D1B3','D2B1','D2B2','D2B3'});
tbl.Sex = sex_labels_c;

within = table( ...
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

    ws_rows = ws_tbl.Properties.RowNames;   % look up row indices by name
    bs_rows = bs_tbl.Properties.RowNames;

    % --- Print between-subjects (Sex) ---
    fprintf('\n  Between-subjects:\n');
    bs_vars = bs_tbl.Properties.VariableNames;
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

    % --- Print within-subjects and interactions ---
    % Each effect row is paired with its error row for df1 / df2
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
                         ws_tbl.DF(eff_idx), ws_tbl.DF(err_idx), ...
                         ws_tbl.pValue(eff_idx));
        end
    end

    % --- Post-hoc for significant effects ---
    p_3way       = get_p(ws_tbl, 'Sex:Day:Block');
    p_sex_block  = get_p(ws_tbl, 'Sex:Block');
    p_sex_day    = get_p(ws_tbl, 'Sex:Day');
    p_day_block  = get_p(ws_tbl, '(Intercept):Day:Block');
    p_block      = get_p(ws_tbl, '(Intercept):Block');
    p_day        = get_p(ws_tbl, '(Intercept):Day');

    ran_posthoc = false;

    if p_3way < 0.05
        fprintf('\n  Significant 3-way interaction — Sex x Day follow-ups within each Block (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        fprintf('\n  Day follow-ups within each Block (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Day', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end

    if p_sex_block < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex x Block — Sex comparisons within each Block (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end

    if p_sex_day < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex x Day — Sex comparisons within each Day (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end

    if p_day_block < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Day x Block — Day comparisons within each Block (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Day', 'By', 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end

    % Main effects with no significant interactions — simple pairwise
    if p_block < 0.05 && p_sex_block >= 0.05 && p_day_block >= 0.05 && p_3way >= 0.05
        fprintf('\n  Significant Block main effect — pairwise block comparisons (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Block', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
    end

    if p_day < 0.05 && p_sex_day >= 0.05 && p_day_block >= 0.05 && p_3way >= 0.05
        fprintf('\n  Significant Day main effect — D1 vs D2 (Tukey-Kramer):\n');
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
        if ~isempty(idx)
            p = tbl.pValue(idx(1));
        else
            p = NaN;
        end
    end
end

% ----- Helper: print one ANOVA effect line -----
function print_effect(label, F, df1, df2, p)
    fprintf('    %-26s F(%d,%d) = %7.3f,  p = %.4f%s\n', label, df1, df2, F, p, sig_stars(p));
end

% ----- Helper: print multcompare table -----
function print_multcomp(mc)
    if ~istable(mc)
        % older MATLAB returns a matrix [g1 g2 lower diff upper pval]
        for r = 1:size(mc, 1)
            p_val = mc(r, 6);
            fprintf('      Group %d vs %d: diff=%7.3f, 95%%CI[%6.3f,%6.3f], p=%.4f%s\n', ...
                mc(r,1), mc(r,2), mc(r,4), mc(r,3), mc(r,5), p_val, sig_stars(p_val));
        end
        return;
    end

    % Newer MATLAB returns a table. 'pValue' is a known column name.
    if ~ismember('pValue', mc.Properties.VariableNames)
        disp(mc);   % fallback if column name is unexpected
        return;
    end

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

        fprintf('      %-42s p = %.4f%s\n', ...
            [lv1 ' vs ' lv2 by_str], p_val, sig_stars(p_val));
    end
end

% ----- Weight bar plot (no block; x-axis = D1 / D2, bars = Female | Male) -----
function plot_weight_bars(f_d1, f_d2, m_d1, m_d2, all_colors, scatter_colors, scatter_markers, ylabel_str, ylim_vals, ytick_vals, fig_w, fig_h)
% X-axis: D1, D2  |  Within each day: Female | Male
% Colors match the block plots: D1 shades for D1 group, D2 shades for D2 group.
% all_colors input order: {F-D1, F-D2, M-D1, M-D2}

% Day-grouped order: F-D1, M-D1, F-D2, M-D2 → indices [1,3,2,4]
day_order      = [1, 3, 2, 4];
plot_colors    = all_colors(day_order);
plot_scatter   = scatter_colors(day_order);
plot_markers   = scatter_markers(day_order);

all_vals = {f_d1(:), m_d1(:), f_d2(:), m_d2(:)};

% Manual x-positions: D1 pair centred on 1, D2 pair centred on 2
bar_w = 0.35;
x_pos = [0.82, 1.18, 1.82, 2.18];   % F-D1, M-D1, F-D2, M-D2

means_vec = cellfun(@(v) nanmean(v),                              all_vals);
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

set(gca, 'XTick', [1, 2], 'XTickLabel', {'RDT D1', 'RDT D2'}, 'FontSize', 12);
xlim([0.4, 2.6]);
ylabel(ylabel_str, 'FontSize', 12);
if ~isempty(ylim_vals);  ylim(ylim_vals);       end
if ~isempty(ytick_vals); yticks(ytick_vals);    end
legend(bar_h, {'Female D1', 'Male D1', 'Female D2', 'Male D2'}, ...
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

% ----- Session completion donut charts -----
function plot_completion_donut(f_d1, f_d2, m_d1, m_d2, group_colors, fig_w, fig_h)
% One donut per group (F-D1, M-D1, F-D2, M-D2).
% Criterion: trials_completed >= 90 = completed session.
% group_colors: {color_f_d1, color_m_d1, color_f_d2, color_m_d2}

completion_thresh = 90;

groups       = {f_d1,        m_d1,      f_d2,        m_d2};
group_labels = {'Female D1', 'Male D1', 'Female D2', 'Male D2'};
not_complete_color = [0.82, 0.82, 0.82];   % light grey for incomplete slice

fprintf('\n--- Session Completion (threshold: >= %d trials) ---\n', completion_thresh);

for g = 1:4
    tc = groups{g}.trials_completed;
    tc = tc(~isnan(tc));
    n_total      = numel(tc);
    n_complete   = sum(tc >= completion_thresh);
    n_incomplete = n_total - n_complete;
    pct_complete   = n_complete   / n_total * 100;
    pct_incomplete = n_incomplete / n_total * 100;

    fprintf('  %s: %d / %d completed  (%.1f%% complete, %.1f%% not)\n', ...
        group_labels{g}, n_complete, n_total, pct_complete, pct_incomplete);

    figure;
    set(gcf, 'Position', [100, 100, fig_w, fig_h]);
    % Set a 2-entry colormap before creating the chart so each slice gets its own colour.
    % Slice order: [not completed, completed]
    colormap(gcf, [not_complete_color; group_colors{g}]);

    dc = donutchart([pct_incomplete, pct_complete], ...
        {'Not completed', 'Completed'}, ...
        'InnerRadius', 0.7, 'StartAngle', 270);

    dc.Title    = sprintf('%s  —  %d / %d completed', group_labels{g}, n_complete, n_total);
    dc.FontSize = 3;
end
end

% ----- 2-way mixed ANOVA: Sex (between) x Day (within) — weight variables -----
function run_2way_anova_weight(f_d1, f_d2, m_d1, m_d2, measure_name)

f_d1 = f_d1(:);  f_d2 = f_d2(:);
m_d1 = m_d1(:);  m_d2 = m_d2(:);
nF = numel(f_d1);
nM = numel(m_d1);

all_wide   = [f_d1, f_d2; m_d1, m_d2];
sex_labels = categorical([repmat({'Female'}, nF, 1); repmat({'Male'}, nM, 1)]);

complete     = all(~isnan(all_wide), 2);
n_complete   = sum(complete);
n_excluded   = (nF + nM) - n_complete;
nF_c         = sum(complete(1:nF));
nM_c         = sum(complete(nF+1:end));

all_wide_c   = all_wide(complete, :);
sex_labels_c = sex_labels(complete);

tbl     = array2table(all_wide_c, 'VariableNames', {'D1', 'D2'});
tbl.Sex = sex_labels_c;
within  = table(categorical({'D1'; 'D2'}), 'VariableNames', {'Day'});

fprintf('\n========== 2-Way Mixed ANOVA: %s ==========\n', measure_name);
fprintf('N = %d subjects (%d female, %d male)', n_complete, nF_c, nM_c);
if n_excluded > 0; fprintf('  [%d excluded — missing data]', n_excluded); end
fprintf('\n');

if nF_c == 0 || nM_c == 0
    fprintf('  [Only one sex has complete data — between-subjects test unavailable]\n');
    return;
end

try
    rm     = fitrm(tbl, 'D1,D2 ~ Sex', 'WithinDesign', within);
    bs_tbl = rm.anova();
    ws_tbl = ranova(rm, 'WithinModel', 'Day');

    bs_rows = bs_tbl.Properties.RowNames;
    ws_rows = ws_tbl.Properties.RowNames;

    % Between-subjects: Sex
    fprintf('\n  Between-subjects:\n');
    bs_vars = bs_tbl.Properties.VariableNames;
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
        end
    else
        fprintf('    [Sex row not found — bs_tbl rows: {%s}  cols: {%s}]\n', ...
            strjoin(bs_rows, ', '), strjoin(bs_vars, ', '));
    end

    % Within-subjects: Day, Sex x Day
    fprintf('  Within-subjects:\n');
    effects_ws = { ...
        '(Intercept):Day', 'Day',       'Error(Day)'; ...
        'Sex:Day',         'Sex x Day', 'Error(Day)'};
    for e = 1:size(effects_ws, 1)
        eff_idx = find(strcmp(ws_rows, effects_ws{e,1}));
        err_idx_w = find(strcmp(ws_rows, effects_ws{e,3}));
        if ~isempty(eff_idx) && ~isempty(err_idx_w)
            print_effect(effects_ws{e,2}, ws_tbl.F(eff_idx), ...
                         ws_tbl.DF(eff_idx), ws_tbl.DF(err_idx_w), ws_tbl.pValue(eff_idx));
        end
    end

    % Post-hoc (cascading: 3-way interaction → 2-way interactions → main effects)
    p_sex_day = get_p(ws_tbl, 'Sex:Day');
    p_day     = get_p(ws_tbl, '(Intercept):Day');
    p_sex     = NaN;
    if ~isempty(sex_idx)
        p_sex = bs_tbl.pValue(sex_idx);
    end

    ran_posthoc = false;
    if ~isnan(p_sex_day) && p_sex_day < 0.05
        fprintf('\n  Significant Sex x Day — Sex comparisons within each Day (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'By', 'Day', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if ~isnan(p_day) && p_day < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Day main effect — D1 vs D2 (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Day', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
        ran_posthoc = true;
    end
    if ~isnan(p_sex) && p_sex < 0.05 && ~ran_posthoc
        fprintf('\n  Significant Sex main effect — Female vs Male (Tukey-Kramer):\n');
        mc = multcompare(rm, 'Sex', 'ComparisonType', 'tukey-kramer');
        print_multcomp(mc);
    end

catch ME
    fprintf('  [ANOVA failed: %s]\n', ME.message);
end
end

% ----- Scatter plot with per-sex regression lines -----
function plot_scatter_regression(x_female, x_male, y_female, y_male, ...
    xlabel_str, ylabel_str, title_str, color_f, color_m, fig_w, fig_h, ...
    xlim_vals, xtick_vals, ylim_vals, ytick_vals)
% Plots scatter points and a linear regression line for each sex group.
% Females = filled squares, Males = filled circles.
% Prints r, R², and p to the command window.
% xlim_vals / ylim_vals : [lo hi] or [] for auto
% xtick_vals / ytick_vals: vector or [] for auto

% Strip NaN pairs
vf = ~isnan(x_female) & ~isnan(y_female);
vm = ~isnan(x_male)   & ~isnan(y_male);
xf = x_female(vf);  yf = y_female(vf);
xm = x_male(vm);    ym = y_male(vm);

% Regression statistics
[r_f, p_f, coef_f, r2_f] = reg_stats(xf, yf);
[r_m, p_m, coef_m, r2_m] = reg_stats(xm, ym);

% Figure
figure;
set(gcf, 'Position', [100, 100, fig_w, fig_h]);
hold on;

scatter(xf, yf, 50, color_f, 's', 'filled', ...
    'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');
scatter(xm, ym, 50, color_m, 'o', 'filled', ...
    'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');

% Regression lines (only if enough points)
if numel(xf) >= 2
    xfit = linspace(min(xf), max(xf), 100);
    plot(xfit, polyval(coef_f, xfit), 'Color', color_f, 'LineWidth', 2);
end
if numel(xm) >= 2
    xfit = linspace(min(xm), max(xm), 100);
    plot(xfit, polyval(coef_m, xfit), 'Color', color_m, 'LineWidth', 2);
end

xlabel(xlabel_str, 'FontSize', 12);
ylabel(ylabel_str, 'FontSize', 12);
title(title_str,   'FontSize', 12);
box off;

% Apply manual axis limits / ticks before reading them for annotation
ax = gca;
if ~isempty(xlim_vals);  xlim(xlim_vals);    end
if ~isempty(xtick_vals); xticks(xtick_vals); end
if ~isempty(ylim_vals);  ylim(ylim_vals);    end
if ~isempty(ytick_vals); yticks(ytick_vals); end

% Stat annotations (top-left, stacked)
xl = xlim(ax);  yl = ylim(ax);
text(xl(1) + 0.03*range(xl), yl(2) - 0.03*range(yl), ...
    sprintf('Female: r=%.2f, p=%.3f', r_f, p_f), ...
    'Color', color_f*0.8, 'FontSize', 9, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'top');
text(xl(1) + 0.03*range(xl), yl(2) - 0.18*range(yl), ...
    sprintf('Male:   r=%.2f, p=%.3f', r_m, p_m), ...
    'Color', max(color_m*0.7, 0), 'FontSize', 9, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'top');

% Fixed inner axes size
margin_left = 65;  margin_bottom = 55;
margin_right = 15; margin_top   = 30;
ax.Units = 'pixels';
ax.Position = [margin_left, margin_bottom, ...
               fig_w - margin_left - margin_right, ...
               fig_h - margin_bottom - margin_top];
hold off;

% Console output
fprintf('\n--- Scatter: %s  vs  %s  [%s] ---\n', xlabel_str, ylabel_str, title_str);
fprintf('  Female (n=%d): r=%.3f, R²=%.3f, p=%.4f\n', numel(xf), r_f, r2_f, p_f);
fprintf('  Male   (n=%d): r=%.3f, R²=%.3f, p=%.4f\n', numel(xm), r_m, r2_m, p_m);
end

% ----- Helper: linear regression + correlation statistics -----
function [r_val, p_val, coef, r2] = reg_stats(x, y)
if numel(x) < 2
    r_val = NaN;  p_val = NaN;  coef = [NaN NaN];  r2 = NaN;
    return;
end
coef   = polyfit(x(:), y(:), 1);
y_pred = polyval(coef, x(:));
ss_res = sum((y(:) - y_pred).^2);
ss_tot = sum((y(:) - mean(y(:))).^2);
r2     = 1 - ss_res / max(ss_tot, eps);
[rmat, pmat] = corrcoef(x(:), y(:));
r_val  = rmat(1,2);
p_val  = pmat(1,2);
end

% ----- Helper: significance stars -----
function s = sig_stars(p)
    if p < 0.001;    s = ' ***';
    elseif p < 0.01; s = ' **';
    elseif p < 0.05; s = ' *';
    else;            s = '';
    end
end

