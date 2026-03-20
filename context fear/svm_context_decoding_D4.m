%% SVM Context Decoding — D4 Session (Experimental vs One Context vs No Shock)
%
% D4 session structure (12 min total, 10 Hz) — inverse of D3:
%   Aversive : frames   1 – 1200   (0–2 min)
%   Neutral  : frames 1201 – 2400  (2–4 min)
%   Aversive : frames 2401 – 3600  (4–6 min)
%   Neutral  : frames 3601 – 4800  (6–8 min)
%   Aversive : frames 4801 – 6000  (8–10 min)
%   Neutral  : frames 6001 – 7200  (10–12 min)

%% ── Parameters ──────────────────────────────────────────────────────────

session_to_analyze = 'D4';
ca_data_type       = 'C_raw';
bin_size_s         = 1;
n_permutations     = 1000;
target_groups      = {'Experimental', 'One Context', 'No Shock'};
sampling_rate      = 10;

aversive_blocks = [1,    1200; 2401, 3600; 4801, 6000];
neutral_blocks  = [1201, 2400; 3601, 4800; 6001, 7200];

n_blocks       = size(neutral_blocks, 1);
bin_frames     = bin_size_s * sampling_rate;
bins_per_block = (neutral_blocks(1,2) - neutral_blocks(1,1) + 1) / bin_frames;

%% ── Run decoding for each group ──────────────────────────────────────────

all_results = struct();

for gg = 1:numel(target_groups)
    target_group = target_groups{gg};
    fprintf('\n##########  Group: %s  ##########\n', target_group);

    % Identify mice
    animalIDs = fieldnames(final);
    exp_mice  = {};
    for ii = 1:numel(animalIDs)
        mid = animalIDs{ii};
        if isfield(final.(mid), 'experimental_grp') && ...
           strcmp(final.(mid).experimental_grp, target_group) && ...
           isfield(final.(mid), session_to_analyze)
            exp_mice{end+1} = mid; %#ok<SAGROW>
        end
    end
    fprintf('Found %d mice with %s session.\n', numel(exp_mice), session_to_analyze);

    if isempty(exp_mice)
        fprintf('  No mice found — skipping group.\n');
        all_results(gg).group     = target_group;
        all_results(gg).mouse_id  = {};
        all_results(gg).accuracy  = [];
        continue
    end

    % Per-mouse storage
    r.mouse_id      = {};
    r.n_neurons     = [];
    r.accuracy      = [];
    r.chance_mean   = [];
    r.chance_std    = [];
    r.p_value       = [];
    r.fold_accuracy = {};

    for mm = 1:numel(exp_mice)
        mid = exp_mice{mm};
        fprintf('\n── Mouse %s ──\n', mid);

        ca_trace = final.(mid).(session_to_analyze).CNMFe_data.(ca_data_type);
        if issparse(ca_trace), ca_trace = full(ca_trace); end

        n_neurons    = size(ca_trace, 1);
        n_timepoints = size(ca_trace, 2);
        fprintf('  Neurons: %d,  Timepoints: %d\n', n_neurons, n_timepoints);

        if n_timepoints < neutral_blocks(end, 2)
            fprintf('  WARNING: Only %d frames; expected >= %d. Skipping.\n', ...
                    n_timepoints, neutral_blocks(end, 2));
            continue
        end

        X_neutral  = bin_ca_blocks(ca_trace, neutral_blocks,  bin_frames, bins_per_block);
        X_aversive = bin_ca_blocks(ca_trace, aversive_blocks, bin_frames, bins_per_block);

        block_id_neutral  = repelem((1:n_blocks)', bins_per_block);
        block_id_aversive = repelem((1:n_blocks)', bins_per_block);

        y_neutral  = zeros(size(X_neutral,  1), 1);
        y_aversive = ones( size(X_aversive, 1), 1);

        % Leave-one-block-out CV
        fold_acc = nan(n_blocks, 1);
        for fold = 1:n_blocks
            test_n  = block_id_neutral  == fold;
            test_av = block_id_aversive == fold;

            X_train = [X_neutral(~test_n, :);  X_aversive(~test_av, :)];
            y_train = [y_neutral(~test_n);      y_aversive(~test_av)];
            X_test  = [X_neutral(test_n, :);   X_aversive(test_av, :)];
            y_test  = [y_neutral(test_n);       y_aversive(test_av)];

            mu = mean(X_train, 1);
            sd = std(X_train,  0, 1);  sd(sd == 0) = 1;
            X_train_z = (X_train - mu) ./ sd;
            X_test_z  = (X_test  - mu) ./ sd;

            mdl = fitclinear(X_train_z, y_train, ...
                             'Learner', 'svm', 'Regularization', 'ridge');
            fold_acc(fold) = mean(predict(mdl, X_test_z) == y_test) * 100;
            fprintf('  Fold %d accuracy: %.1f%%\n', fold, fold_acc(fold));
        end

        mean_acc = mean(fold_acc);
        fprintf('  Mean accuracy: %.1f%%\n', mean_acc);

        % Permutation test
        X_all   = [X_neutral;  X_aversive];
        y_all   = [y_neutral;  y_aversive];
        blk_all = [block_id_neutral; block_id_aversive];
        cxt_all = [zeros(size(X_neutral,1),1); ones(size(X_aversive,1),1)];

        perm_acc = nan(n_permutations, 1);
        fprintf('  Running %d permutations', n_permutations);
        for pp = 1:n_permutations
            if mod(pp, 100) == 0, fprintf(' %d', pp); end
            y_shuf = y_all(randperm(numel(y_all)));
            pa = nan(n_blocks, 1);
            for fold = 1:n_blocks
                tm = (blk_all == fold & cxt_all == 0) | (blk_all == fold & cxt_all == 1);
                tr = ~(blk_all == fold);
                X_tr = X_all(tr, :);  y_tr = y_shuf(tr);
                X_te = X_all(tm, :);  y_te = y_shuf(tm);
                mu = mean(X_tr,1); sd = std(X_tr,0,1); sd(sd==0)=1;
                X_tr_z = (X_tr-mu)./sd;  X_te_z = (X_te-mu)./sd;
                m = fitclinear(X_tr_z, y_tr, 'Learner', 'svm', 'Regularization', 'ridge');
                pa(fold) = mean(predict(m, X_te_z) == y_te) * 100;
            end
            perm_acc(pp) = mean(pa);
        end
        fprintf('\n');

        p_val = mean(perm_acc >= mean_acc);
        fprintf('  Permutation p-value: %.4f  (chance = %.1f +/- %.1f%%)\n', ...
                p_val, mean(perm_acc), std(perm_acc));

        r.mouse_id{end+1}      = mid;
        r.n_neurons(end+1)     = n_neurons;
        r.accuracy(end+1)      = mean_acc;
        r.chance_mean(end+1)   = mean(perm_acc);
        r.chance_std(end+1)    = std(perm_acc);
        r.p_value(end+1)       = p_val;
        r.fold_accuracy{end+1} = fold_acc;
    end

    all_results(gg).group     = target_group;
    all_results(gg).mouse_id  = r.mouse_id;
    all_results(gg).n_neurons = r.n_neurons;
    all_results(gg).accuracy  = r.accuracy;
    all_results(gg).chance_mean = r.chance_mean;
    all_results(gg).chance_std  = r.chance_std;
    all_results(gg).p_value     = r.p_value;
    all_results(gg).fold_accuracy = r.fold_accuracy;

    % Per-group summary
    n_mice     = numel(r.mouse_id);
    group_mean = mean(r.accuracy);
    group_sem  = std(r.accuracy) / sqrt(n_mice);
    [~, p_tt, ~, st] = ttest(r.accuracy, 50);

    fprintf('\n--- %s summary ---\n', target_group);
    fprintf('Mean accuracy: %.1f +/- %.1f%% SEM  (n=%d)\n', group_mean, group_sem, n_mice);
    fprintf('t-test vs 50%%: t(%d)=%.2f, p=%.4f\n', st.df, st.tstat, p_tt);
    fprintf('\n%-20s  %6s  %10s  %15s  %8s\n', 'Mouse','Cells','Accuracy%','Chance+/-std%','p-val');
    fprintf('%s\n', repmat('-',1,65));
    for mm = 1:n_mice
        fprintf('%-20s  %6d  %10.1f  %5.1f +/- %-5.1f  %8.4f\n', ...
                r.mouse_id{mm}, r.n_neurons(mm), r.accuracy(mm), ...
                r.chance_mean(mm), r.chance_std(mm), r.p_value(mm));
    end

    all_results(gg).group_mean = group_mean;
    all_results(gg).group_sem  = group_sem;
    all_results(gg).p_ttest    = p_tt;
    all_results(gg).tstat      = st;
end

%% ── Between-group comparison (Kruskal-Wallis + post-hoc rank-sum) ────────

acc_all    = [];
grp_labels = [];
for gg = 1:numel(all_results)
    acc_all    = [acc_all,    all_results(gg).accuracy];
    grp_labels = [grp_labels, gg * ones(1, numel(all_results(gg).accuracy))];
end

[p_kw, ~, stats_kw] = kruskalwallis(acc_all, grp_labels, 'off');
fprintf('\nKruskal-Wallis across groups: p = %.4f\n', p_kw);

% Pairwise rank-sum post-hoc
pairs = nchoosek(1:numel(target_groups), 2);
p_posthoc = nan(size(pairs, 1), 1);
for pp = 1:size(pairs, 1)
    a = all_results(pairs(pp,1)).accuracy;
    b = all_results(pairs(pp,2)).accuracy;
    if ~isempty(a) && ~isempty(b)
        p_posthoc(pp) = ranksum(a, b);
        fprintf('  %s vs %s: p = %.4f\n', ...
                target_groups{pairs(pp,1)}, target_groups{pairs(pp,2)}, p_posthoc(pp));
    end
end

%% ── Figure ───────────────────────────────────────────────────────────────

group_colors = {[0.2 0.5 0.9], [0.9 0.4 0.2], [0.3 0.75 0.4]};   % blue, orange, green

fig = figure('Position', [100, 100, 1000, 480], 'Color', 'w', 'Visible', 'off');
tl  = tiledlayout(1, numel(target_groups) + 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Per-group per-mouse tiles
for gg = 1:numel(target_groups)
    ax = nexttile;
    hold(ax, 'on');

    r   = all_results(gg);
    n_m = numel(r.mouse_id);
    if n_m == 0, title(ax, r.group); continue; end

    bar(ax, 1:n_m, r.accuracy, 0.6, 'FaceColor', group_colors{gg}, 'EdgeColor', 'none');
    errorbar(ax, 1:n_m, r.chance_mean, r.chance_std, 'o', ...
             'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], ...
             'MarkerSize', 4, 'LineWidth', 1, 'CapSize', 5);
    yline(ax, 50, '--k', 'LineWidth', 1);

    for mm = 1:n_m
        if r.p_value(mm) < 0.05
            text(ax, mm, r.accuracy(mm)+2, '*', 'HorizontalAlignment','center','FontSize',13);
        end
    end

    set(ax, 'XTick', 1:n_m, 'XTickLabel', r.mouse_id, ...
            'XTickLabelRotation', 45, 'TickDir', 'out', 'Box', 'off');
    ylabel(ax, 'Decoding accuracy (%)');
    title(ax, r.group, 'Interpreter', 'none');
    ylim(ax, [0, 112]);
end

% Group comparison tile
ax3 = nexttile;
hold(ax3, 'on');
yline(ax3, 50, '--k', 'LineWidth', 1);

for gg = 1:numel(target_groups)
    r   = all_results(gg);
    n_m = numel(r.mouse_id);
    if n_m == 0, continue; end

    scatter(ax3, gg*ones(1,n_m), r.accuracy, 55, group_colors{gg}, 'filled', ...
            'jitter', 'on', 'JitterAmount', 0.06);
    errorbar(ax3, gg, r.group_mean, r.group_sem, 'o', ...
             'Color', 'k', 'MarkerFaceColor', 'k', ...
             'MarkerSize', 7, 'LineWidth', 2, 'CapSize', 8);

    % t-test annotation
    if r.p_ttest < 0.001, ps = 'p<0.001';
    else, ps = sprintf('p=%.3f', r.p_ttest); end
    text(ax3, gg, max(r.accuracy)+4, ps, 'HorizontalAlignment','center','FontSize',9);
end

% Between-group significance brackets (stacked)
y_top = max([all_results.accuracy]) + 8;
for pp = 1:size(pairs, 1)
    if isnan(p_posthoc(pp)), continue; end
    x1 = pairs(pp,1);  x2 = pairs(pp,2);
    yb = y_top + (pp-1) * 7;
    plot(ax3, [x1 x1 x2 x2], [yb-1.5 yb yb yb-1.5], '-k', 'LineWidth', 0.8);
    if p_posthoc(pp) < 0.001, bstr = 'p<0.001';
    elseif p_posthoc(pp) < 0.05, bstr = sprintf('p=%.3f', p_posthoc(pp));
    else, bstr = 'ns'; end
    text(ax3, (x1+x2)/2, yb+0.5, bstr, 'HorizontalAlignment','center','FontSize',8);
end

set(ax3, 'XTick', 1:numel(target_groups), 'XTickLabel', target_groups, ...
         'TickDir', 'out', 'Box', 'off');
ylabel(ax3, 'Decoding accuracy (%)');
title(ax3, 'Group comparison');
xlim(ax3, [0.4, numel(target_groups)+0.6]);
ylim(ax3, [0, 135]);

sgtitle('SVM Context Decoding — D4', 'FontWeight', 'bold');

out_fig = 'I:\MATLAB\my_repo\context fear\svm_context_decoding_D4.png';
exportgraphics(fig, out_fig, 'Resolution', 150);
fprintf('\nFigure saved to: %s\n', out_fig);

%% ── Local function ───────────────────────────────────────────────────────

function X = bin_ca_blocks(ca_trace, block_edges, bin_frames, bins_per_block)
    n_b = size(block_edges, 1);
    X   = zeros(n_b * bins_per_block, size(ca_trace, 1));
    row = 1;
    for b = 1:n_b
        block_data = ca_trace(:, block_edges(b,1):block_edges(b,2));
        for k = 1:bins_per_block
            idx = (k-1)*bin_frames + 1 : k*bin_frames;
            X(row, :) = mean(block_data(:, idx), 2)';
            row = row + 1;
        end
    end
end
