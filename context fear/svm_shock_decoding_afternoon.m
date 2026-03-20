%% SVM Shock Decoding — Morning & Afternoon Sessions (D1 & D2)
%
% The same 5s-window decoding framework is applied to both session types:
%
%   Morning   (D1_Morning, D2_Morning):   NO shocks — negative control
%   Afternoon (D1_Afternoon, D2_Afternoon): shocks at minutes 4,5,6,7,8,9
%
%   "Target"   class: 5s window starting at each of 6 "shock" timepoints
%                     (minutes 4-9 → frames 2401,3001,3601,4201,4801,5401)
%   "Baseline" class: 6 matched 5s windows evenly spaced in pre-period (min 0-4)
%   Bins: 1s → 5 bins/window × 6 windows = 30 bins per class
%   CV:   leave-one-window-out (6-fold, paired target/baseline hold-out)
%
% Per-mouse accuracy = mean across available sessions of that type (D1 and/or D2).
% Morning sessions serve as a within-design negative control: if decoding is
% shock-driven, all groups should perform at chance in the morning.

%% ── Parameters ──────────────────────────────────────────────────────────

ca_data_type   = 'C_raw';
bin_size_s     = 1;
sampling_rate  = 10;
n_permutations = 1000;
target_groups  = {'Experimental', 'One Context', 'No Shock'};

session_sets   = { {'D1_Morning','D2_Morning'}, {'D1_Afternoon','D2_Afternoon'} };
session_labels = {'Morning', 'Afternoon'};
n_sess_types   = numel(session_sets);

% "Shock" window definitions (identical timepoints applied to both session types)
shock_onset_min    = [4, 5, 6, 7, 8, 9];
shock_onset_frames = shock_onset_min * 60 * sampling_rate + 1;   % [2401,3001,...,5401]
window_s           = 5;
window_frames      = window_s * sampling_rate;                    % 50 frames
n_windows          = numel(shock_onset_frames);                   % 6
bin_frames         = bin_size_s * sampling_rate;                  % 10 frames/bin
bins_per_window    = window_s / bin_size_s;                       % 5 bins/window

% Evenly-spaced baseline windows in pre-period (frames 1–2400)
pre_shock_end  = 4 * 60 * sampling_rate;                          % 2400
spacing        = floor((pre_shock_end - window_frames) / (n_windows - 1));  % 470
noshock_onsets = 1 : spacing : (1 + spacing*(n_windows-1));       % [1,471,941,1411,1881,2351]

min_frames_needed = shock_onset_frames(end) + window_frames - 1;  % last frame required

%% ── Run decoding: groups × session types ────────────────────────────────

% Pre-allocate with empty struct fields so indexing into (gg,st) works
empty_r = struct('group','','session_type','','mouse_id',{{}},'n_neurons',[], ...
                 'accuracy',[],'chance_mean',[],'chance_std',[],'p_value',[], ...
                 'group_mean',NaN,'group_sem',NaN,'p_ttest',NaN,'tstat',struct('df',0,'tstat',0));
all_results = repmat(empty_r, numel(target_groups), n_sess_types);

for gg = 1:numel(target_groups)
    target_group = target_groups{gg};

    % Find mice in this group that have data in at least one session set
    all_sessions = [session_sets{1}, session_sets{2}];
    animalIDs    = fieldnames(final);
    exp_mice     = {};
    for ii = 1:numel(animalIDs)
        mid = animalIDs{ii};
        if isfield(final.(mid), 'experimental_grp') && ...
           strcmp(final.(mid).experimental_grp, target_group) && ...
           any(cellfun(@(s) isfield(final.(mid), s), all_sessions))
            exp_mice{end+1} = mid; %#ok<SAGROW>
        end
    end

    fprintf('\n##########  Group: %s  (%d mice)  ##########\n', target_group, numel(exp_mice));

    for st = 1:n_sess_types
        sess_type          = session_labels{st};
        sessions_to_analyze = session_sets{st};
        fprintf('\n  ── Session type: %s ──\n', sess_type);

        r = empty_r;
        r.group        = target_group;
        r.session_type = sess_type;

        for mm = 1:numel(exp_mice)
            mid = exp_mice{mm};

            % Check mouse has any sessions of this type
            has_sess = any(cellfun(@(s) isfield(final.(mid), s), sessions_to_analyze));
            if ~has_sess, continue; end

            fprintf('\n    Mouse %s\n', mid);
            sess_accs    = [];
            sess_perm    = [];
            sess_neurons = [];

            for ss = 1:numel(sessions_to_analyze)
                sess = sessions_to_analyze{ss};
                if ~isfield(final.(mid), sess), continue; end

                ca_trace = final.(mid).(sess).CNMFe_data.(ca_data_type);
                if issparse(ca_trace), ca_trace = full(ca_trace); end

                n_neurons    = size(ca_trace, 1);
                n_timepoints = size(ca_trace, 2);
                fprintf('      [%s] Neurons: %d,  Timepoints: %d\n', sess, n_neurons, n_timepoints);

                if n_timepoints < min_frames_needed
                    fprintf('      WARNING: only %d frames (need %d). Skipping.\n', ...
                            n_timepoints, min_frames_needed);
                    continue
                end

                X_target   = bin_windows(ca_trace, shock_onset_frames, window_frames, bin_frames, bins_per_window);
                X_baseline = bin_windows(ca_trace, noshock_onsets,     window_frames, bin_frames, bins_per_window);

                win_id   = repelem((1:n_windows)', bins_per_window);
                y_target = ones( size(X_target,   1), 1);
                y_base   = zeros(size(X_baseline, 1), 1);

                % ── Leave-one-window-out CV ──────────────────────────────
                fold_acc = nan(n_windows, 1);
                for fold = 1:n_windows
                    tm = win_id == fold;
                    X_train = [X_target(~tm,:); X_baseline(~tm,:)];
                    y_train = [y_target(~tm);   y_base(~tm)];
                    X_test  = [X_target( tm,:); X_baseline( tm,:)];
                    y_test  = [y_target( tm);   y_base( tm)];
                    mu = mean(X_train,1); sd = std(X_train,0,1); sd(sd==0)=1;
                    mdl = fitclinear((X_train-mu)./sd, y_train, ...
                                     'Learner','svm','Regularization','ridge');
                    fold_acc(fold) = mean(predict(mdl,(X_test-mu)./sd)==y_test)*100;
                end

                mean_acc = mean(fold_acc);
                fprintf('      [%s] Accuracy: %.1f%%\n', sess, mean_acc);
                sess_accs(end+1)    = mean_acc; %#ok<AGROW>
                sess_neurons(end+1) = n_neurons;

                % ── Permutation test ─────────────────────────────────────
                X_all   = [X_target;   X_baseline];
                y_all   = [y_target;   y_base];
                blk_all = [win_id;     win_id];

                perm_dist = nan(n_permutations, 1);
                fprintf('      Permutations:');
                for pp = 1:n_permutations
                    if mod(pp,200)==0, fprintf(' %d',pp); end
                    y_shuf = y_all(randperm(numel(y_all)));
                    pa = nan(n_windows,1);
                    for fold = 1:n_windows
                        tm = blk_all==fold;
                        X_tr=X_all(~tm,:); y_tr=y_shuf(~tm);
                        X_te=X_all( tm,:); y_te=y_shuf(tm);
                        mu=mean(X_tr,1); sd=std(X_tr,0,1); sd(sd==0)=1;
                        m=fitclinear((X_tr-mu)./sd,y_tr,'Learner','svm','Regularization','ridge');
                        pa(fold)=mean(predict(m,(X_te-mu)./sd)==y_te)*100;
                    end
                    perm_dist(pp) = mean(pa);
                end
                fprintf('\n');
                sess_perm = [sess_perm, perm_dist]; %#ok<AGROW>
            end % sessions loop

            if isempty(sess_accs), continue; end

            mouse_acc   = mean(sess_accs);
            pooled_perm = mean(sess_perm, 2);
            p_val       = mean(pooled_perm >= mouse_acc);

            fprintf('      Mean accuracy: %.1f%%  chance=%.1f±%.1f%%  p=%.4f\n', ...
                    mouse_acc, mean(pooled_perm), std(pooled_perm), p_val);

            r.mouse_id{end+1}    = mid;
            r.n_neurons(end+1)   = round(mean(sess_neurons));
            r.accuracy(end+1)    = mouse_acc;
            r.chance_mean(end+1) = mean(pooled_perm);
            r.chance_std(end+1)  = std(pooled_perm);
            r.p_value(end+1)     = p_val;
        end % mice loop

        % Group-level stats
        n_mice     = numel(r.mouse_id);
        if n_mice >= 1
            r.group_mean = mean(r.accuracy);
            r.group_sem  = std(r.accuracy) / sqrt(n_mice);
        end
        if n_mice >= 2
            [~, p_tt, ~, tst] = ttest(r.accuracy, 50);
            r.p_ttest = p_tt;
            r.tstat   = tst;
        end

        all_results(gg, st) = r;

        fprintf('\n    %s | %s: %.1f ± %.1f%% SEM  (n=%d)', ...
                target_group, sess_type, r.group_mean, r.group_sem, n_mice);
        if ~isnan(r.p_ttest)
            fprintf('  t(%d)=%.2f, p=%.4f', r.tstat.df, r.tstat.tstat, r.p_ttest);
        end
        fprintf('\n');
    end % session type loop
end % group loop

%% ── Between-group comparisons (afternoon only) ───────────────────────────

fprintf('\n=== Between-group comparison (Afternoon) ===\n');
aft_idx = 2;
acc_aft    = [];
grp_labels = [];
for gg = 1:numel(target_groups)
    a = all_results(gg, aft_idx).accuracy;
    acc_aft    = [acc_aft,    a];
    grp_labels = [grp_labels, gg*ones(1,numel(a))];
end

[p_kw, ~, ~] = kruskalwallis(acc_aft, grp_labels, 'off');
fprintf('Kruskal-Wallis (afternoon): p = %.4f\n', p_kw);

pairs = nchoosek(1:numel(target_groups), 2);
p_posthoc = nan(size(pairs,1), 1);
for pp = 1:size(pairs,1)
    a = all_results(pairs(pp,1), aft_idx).accuracy;
    b = all_results(pairs(pp,2), aft_idx).accuracy;
    if ~isempty(a) && ~isempty(b)
        p_posthoc(pp) = ranksum(a, b);
        fprintf('  %s vs %s: p=%.4f\n', ...
                target_groups{pairs(pp,1)}, target_groups{pairs(pp,2)}, p_posthoc(pp));
    end
end

%% ── Figure ───────────────────────────────────────────────────────────────
%
% Layout: (n_groups + 1) columns
%   Columns 1–n_groups : per-group paired morning/afternoon plot
%   Last column        : group comparison panel (afternoon, with morning overlay)

group_colors  = {[0.2 0.5 0.9], [0.9 0.4 0.2], [0.3 0.75 0.4]};
morn_color    = [0.75 0.75 0.75];   % light grey for morning
n_groups      = numel(target_groups);

fig = figure('Position',[50,50,1150,500],'Color','w','Visible','off');
tlo = tiledlayout(1, n_groups+1, 'TileSpacing','compact','Padding','compact');

% ── Per-group paired panels ───────────────────────────────────────────────
for gg = 1:n_groups
    ax = nexttile; hold(ax,'on');
    r_morn = all_results(gg, 1);
    r_aft  = all_results(gg, 2);

    yline(ax, 50, '--', 'Color',[0.6 0.6 0.6], 'LineWidth',1);

    % Build a shared mouse list to enable pairing
    all_mids  = union(r_morn.mouse_id, r_aft.mouse_id);
    n_mice    = numel(all_mids);
    morn_acc  = nan(1, n_mice);
    aft_acc   = nan(1, n_mice);
    for mm = 1:n_mice
        im = find(strcmp(r_morn.mouse_id, all_mids{mm}), 1);
        ia = find(strcmp(r_aft.mouse_id,  all_mids{mm}), 1);
        if ~isempty(im), morn_acc(mm) = r_morn.accuracy(im); end
        if ~isempty(ia), aft_acc(mm)  = r_aft.accuracy(ia);  end
    end

    % Connecting lines for paired mice (both sessions available)
    paired = ~isnan(morn_acc) & ~isnan(aft_acc);
    for mm = 1:n_mice
        if paired(mm)
            plot(ax, [1 2], [morn_acc(mm) aft_acc(mm)], '-', ...
                 'Color',[0.7 0.7 0.7], 'LineWidth',0.8);
        end
    end

    % Individual points
    scatter(ax, ones(1,n_mice),  morn_acc, 45, morn_color, 'filled', ...
            'jitter','on','JitterAmount',0.06, 'MarkerEdgeColor','none');
    scatter(ax, 2*ones(1,n_mice), aft_acc, 45, group_colors{gg}, 'filled', ...
            'jitter','on','JitterAmount',0.06, 'MarkerEdgeColor','none');

    % Group means ± SEM
    for st_idx = 1:2
        r_st = all_results(gg, st_idx);
        if isnan(r_st.group_mean), continue; end
        col  = [morn_color; group_colors{gg}];
        errorbar(ax, st_idx, r_st.group_mean, r_st.group_sem, 'o', ...
                 'Color','k','MarkerFaceColor',col(st_idx,:), ...
                 'MarkerSize',9,'LineWidth',2,'CapSize',8);

        % t-test vs 50% annotation
        if isnan(r_st.p_ttest),        ps = '';
        elseif r_st.p_ttest < 0.001,   ps = 'p<0.001';
        elseif r_st.p_ttest < 0.05,    ps = sprintf('p=%.3f',r_st.p_ttest);
        else,                           ps = sprintf('ns (p=%.2f)',r_st.p_ttest);
        end
        text(ax, st_idx, 105, ps, 'HorizontalAlignment','center','FontSize',8);
    end

    % Morning vs Afternoon within-group significance (Wilcoxon signed-rank on paired mice)
    if sum(paired) >= 3
        [p_sr, ~] = signrank(morn_acc(paired), aft_acc(paired));
        y_brk = 112;
        plot(ax,[1 1 2 2],[y_brk-2 y_brk y_brk y_brk-2],'-k','LineWidth',0.8);
        if p_sr < 0.001, bstr = 'p<0.001';
        elseif p_sr < 0.05, bstr = sprintf('p=%.3f',p_sr);
        else, bstr = 'ns'; end
        text(ax, 1.5, y_brk+1, bstr,'HorizontalAlignment','center','FontSize',8);
    end

    set(ax,'XTick',[1 2],'XTickLabel',{'Morning','Afternoon'},'TickDir','out','Box','off');
    ylabel(ax,'Decoding accuracy (%)');
    title(ax, target_groups{gg},'Interpreter','none');
    ylim(ax,[0,122]);
    xlim(ax,[0.5 2.5]);
end

% ── Group comparison panel (afternoon foreground, morning background) ─────
ax_cmp = nexttile; hold(ax_cmp,'on');
yline(ax_cmp, 50, '--', 'Color',[0.6 0.6 0.6], 'LineWidth',1);

x_morn = (1:n_groups) - 0.2;
x_aft  = (1:n_groups) + 0.2;

for gg = 1:n_groups
    r_m = all_results(gg, 1);   % morning
    r_a = all_results(gg, 2);   % afternoon

    % Morning (open circles, light)
    if ~isempty(r_m.accuracy)
        scatter(ax_cmp, x_morn(gg)*ones(1,numel(r_m.accuracy)), r_m.accuracy, 40, ...
                morn_color, 'filled','jitter','on','JitterAmount',0.04);
        errorbar(ax_cmp, x_morn(gg), r_m.group_mean, r_m.group_sem, 'o', ...
                 'Color',[0.4 0.4 0.4],'MarkerFaceColor',morn_color, ...
                 'MarkerSize',8,'LineWidth',1.5,'CapSize',6);
    end

    % Afternoon (filled circles, group color)
    if ~isempty(r_a.accuracy)
        scatter(ax_cmp, x_aft(gg)*ones(1,numel(r_a.accuracy)), r_a.accuracy, 40, ...
                group_colors{gg}, 'filled','jitter','on','JitterAmount',0.04);
        errorbar(ax_cmp, x_aft(gg), r_a.group_mean, r_a.group_sem, 'o', ...
                 'Color','k','MarkerFaceColor',group_colors{gg}, ...
                 'MarkerSize',8,'LineWidth',1.5,'CapSize',6);
    end
end

% Between-group post-hoc brackets (afternoon accuracy)
valid_accs = cellfun(@(g) all_results(find(strcmp(target_groups,g),1), aft_idx).accuracy, ...
                     target_groups, 'UniformOutput', false);
y_top = max(cellfun(@(a) max([a, 50]), valid_accs)) + 8;

for pp = 1:size(pairs,1)
    if isnan(p_posthoc(pp)), continue; end
    x1 = x_aft(pairs(pp,1)); x2 = x_aft(pairs(pp,2));
    yb = y_top + (pp-1)*8;
    plot(ax_cmp,[x1 x1 x2 x2],[yb-1.5 yb yb yb-1.5],'-k','LineWidth',0.8);
    if p_posthoc(pp) < 0.001, bstr='p<0.001';
    elseif p_posthoc(pp) < 0.05, bstr=sprintf('p=%.3f',p_posthoc(pp));
    else, bstr='ns'; end
    text(ax_cmp,(x1+x2)/2,yb+1,bstr,'HorizontalAlignment','center','FontSize',8);
end

% Legend proxies
h_morn = scatter(ax_cmp, nan, nan, 40, morn_color,       'filled');
h_aft  = scatter(ax_cmp, nan, nan, 40, [0.5 0.5 0.5],   'filled');
legend(ax_cmp,[h_morn, h_aft],{'Morning','Afternoon'},'Location','southeast','Box','off','FontSize',8);

set(ax_cmp,'XTick',1:n_groups,'XTickLabel',target_groups,'TickDir','out','Box','off');
ylabel(ax_cmp,'Decoding accuracy (%)');
title(ax_cmp,'Group comparison (Aft.)');
xlim(ax_cmp,[0.4, n_groups+0.6]);
ylim(ax_cmp,[0, y_top + size(pairs,1)*8 + 8]);

sgtitle(tlo,'SVM Shock Decoding — Morning (control) vs Afternoon (shock)','FontWeight','bold');

out_fig = 'I:\MATLAB\my_repo\context fear\svm_shock_decoding_morn_vs_aft.png';
exportgraphics(fig, out_fig, 'Resolution',150);
fprintf('\nFigure saved to: %s\n', out_fig);

%% ── Local function ───────────────────────────────────────────────────────

function X = bin_windows(ca_trace, onsets, window_frames, bin_frames, bins_per_window)
% Extract and 1s-bin fixed-length windows from ca_trace.
% onsets : 1-indexed start frames (row vec)
% Output : [n_windows*bins_per_window  ×  n_neurons]
    n_win = numel(onsets);
    X     = zeros(n_win * bins_per_window, size(ca_trace,1));
    row   = 1;
    for w = 1:n_win
        win_data = ca_trace(:, onsets(w) : onsets(w)+window_frames-1);
        for k = 1:bins_per_window
            idx = (k-1)*bin_frames+1 : k*bin_frames;
            X(row,:) = mean(win_data(:,idx), 2)';
            row = row + 1;
        end
    end
end
