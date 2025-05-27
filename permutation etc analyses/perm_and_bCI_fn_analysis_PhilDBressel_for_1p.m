function [comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, x_limits, y_limits, input_threshold) 

% data_array is expected to be a 1x# cell array. each cell should the
% following orientation:
% rows = neurons
% columns = time samples

%% Run access_behav_struct_v# and then run the below scripts (choose which CI or permutation-based approach you want)


% for now just add accuracy_per_iteration from decoding manually
% ZallMean_for_perm_test = {accuracy_per_iteration{1, 1}{1, 1}', accuracy_per_iteration{1, 2}{1, 1}'};  

% ZallMean_for_perm_test = {aa_large, aa_small}; 




ZallMean_for_perm_test = mean_data_array; 
ZallSEM_for_perm_test = sem_data_array;

bCI_tCI_CI_threshold = 0; % use for things that hover around 0 (dF/F for example)
% bCI_tCI_CI_threshold = 0.5; % use for things like decoding


%% ERT example
sig = 0.001;
consec_thresh = input_threshold; %formerly hardcoded to 10 (for 30 Hz data) or 3 (for 10 Hz data) % 1017.3Hz sample rate / 3Hz filter %340 PRD used 340 because his data WERE NOT DOWNSAMPLED (e.g., they were 1018 samples per sec) our data are 30 samples per sec or so.
% this means if his threshold is 1017/340 = ~3, ours should be 30/x = 3,
% which is 10

% Graphing parameters
ylims = y_limits;
% xlims = [ts1(1) round(ts1(end))];
xlims = x_limits;
% sig_plot_level = linspace(4,3.2,7);

ind_2 = ts1;

comparison = struct;

for ii = 1:size(ZallMean_for_perm_test, 2)
    comparison(ii).data = ZallMean_for_perm_test{1,ii};
    comparison(ii).sem_data = ZallSEM_for_perm_test{1,ii};
%     comparison(ii).mean_Cp = mean(comparison(ii).data,1, 'omitnan'); %
    comparison(ii).mean_Cp = mean(comparison(ii).data,1); %
%     comparison(ii).sem_Cp = nansem(comparison(ii).data); %
    comparison(ii).sem_Cp = mean(comparison(ii).sem_data, 1); %
    adjust_labels(ii) = max(comparison(ii).mean_Cp)+2*max(comparison(ii).sem_Cp);
    max_mean(ii) = max(comparison(ii).mean_Cp);
    max_SEM(ii) = max(comparison(ii).sem_Cp);
end

max_adjustment = max(max_mean)+2*max(max_SEM);

sig_plot_level = linspace(max_adjustment+2*max(max_SEM), max_adjustment-max(max_SEM), 7);

% sig_plot_level_v2 = linspace(max_adjustment+0.5, max_adjustment, 6);

% sig_plot_level_v2 = linspace(max_adjustment, max_adjustment, 6);

sig_plot_level_v2 = 0;

% arg_string = string(varargin_array);

% arg_string_combine = join(arg_string, '=');
% 
%% Get labels for each comparison based on the events that the data are filtered on

% 
% for ii = 1:size(arg_string, 2)
%     for qq = 1:size(arg_string, 1)
%         if arg_string(qq,ii) == 'BLOCK';
%             arg_string_block(qq,:) = join([arg_string(qq,ii) arg_string(qq, ii+1)], '=');
%         end
%     end
% end
% 
% cc = 1;
% 
% for ii = 1:size(varargin_array, 2)
%     for qq = 1:size(varargin_array, 1)-1
%         if ~isequal(varargin_array(qq,ii), varargin_array(qq+1,ii))
%             if isstring(varargin_array(qq,ii))
%                 arg_string_other = join([arg_string(:,ii) arg_string(:, ii+1)], '=')
%                 cc = cc+1;
%             elseif ~isstring(varargin_array(qq,ii))
%                 arg_string_other = join([arg_string(:,ii-1) arg_string(:, ii)], '=')
%                 cc = cc+1;
%             end
%         end
%     end
% end

%%
clear ii

for ii = 1:size(comparison, 2)
    [comparison(ii).uv.n_Cp, comparison(ii).uv.ev_win] = size(comparison(ii).data);
%     [n_Cm,~] = size(ZallMean_small_trunc);
    timeline = linspace(ts1(1),ts1(end),comparison(1).uv.ev_win);

    comparison(ii).uv.Cp_t_crit = tinv(1-sig/2,comparison(ii).uv.n_Cp-1);
%     Cm_t_crit = tinv(1-sig/2,n_Cm-1);


    comparison(ii).Cp_bCI = boot_CI(comparison(ii).data,1000,sig);
    [comparison(ii).adjLCI,comparison(ii).adjUCI] = CIadjust(comparison(ii).Cp_bCI(1,:),comparison(ii).Cp_bCI(2,:),[],comparison(ii).uv.n_Cp,2);
    comparison(ii).Cp_bCIexp = [comparison(ii).adjLCI;comparison(ii).adjUCI];
    comparison(ii).Cp_tCI = [comparison(ii).mean_Cp - comparison(ii).sem_Cp*comparison(ii).uv.Cp_t_crit ; comparison(ii).mean_Cp + comparison(ii).sem_Cp*comparison(ii).uv.Cp_t_crit];
    
    %tCI
    comparison(ii).Cp_tCI_sig = NaN(1,comparison(ii).uv.ev_win);
    comparison(ii).sig_idx_tCI = find((comparison(ii).Cp_tCI(1,:) > bCI_tCI_CI_threshold) | (comparison(ii).Cp_tCI(2,:) < bCI_tCI_CI_threshold));
    comparison(ii).consec_tCI = consec_idx(comparison(ii).sig_idx_tCI, consec_thresh);
    comparison(ii).Cp_tCI_sig(comparison(ii).sig_idx_tCI(comparison(ii).consec_tCI)) = sig_plot_level_v2;

    %bCI
    comparison(ii).Cp_bCIexp_sig = NaN(1,comparison(ii).uv.ev_win);
    comparison(ii).sig_idx_bCI = find((comparison(ii).Cp_bCIexp(1,:) > bCI_tCI_CI_threshold) | (comparison(ii).Cp_bCIexp(2,:) < bCI_tCI_CI_threshold));
    comparison(ii).consec_bCI = consec_idx(comparison(ii).sig_idx_bCI,consec_thresh);
    comparison(ii).Cp_bCIexp_sig(comparison(ii).sig_idx_bCI(comparison(ii).consec_bCI)) = sig_plot_level_v2;



%     mean_Cm = mean(ZallMean_small_trunc,1);
%     sem_Cm = sem(ZallMean_small_trunc);
%     Cm_bCI = boot_CI(ZallMean_small_trunc,1000,sig);
%     [adjLCI,adjUCI] = CIadjust(Cm_bCI(1,:),Cm_bCI(2,:),[],n_Cm,2);
%     Cm_bCIexp = [adjLCI;adjUCI];
%     Cm_tCI = [mean_Cm - sem_Cm*Cm_t_crit ; mean_Cm + sem_Cm*Cm_t_crit];

%     perm_p = permTest_array(ZallMean_large_trunc,ZallMean_small_trunc,1000);% permTest_array(ERT_test.Cp_off1,ERT_test.Cm_off3,1000);
%     diff_bCI = boot_diffCI(ZallMean_large_trunc,ZallMean_small_trunc,1000,sig);
%     [adjLCI,adjUCI] = CIadjust(diff_bCI(1,:),diff_bCI(2,:),[],n_Cm,2);
%     diff_bCIexp = [adjLCI;adjUCI];
end



%% 

num_comparisons = 1:size(comparison, 2);

% Get all possible comparisons given the input data
if size(comparison,2) == 2
    pairwise_comps = num_comparisons;
elseif size(comparison, 2) > 2
    pairwise_comps = nchoosek(num_comparisons, (length(num_comparisons)-1));

end

for qq = 1:size(pairwise_comps, 1)
    perm_p(qq,:) = permTest_array(comparison(pairwise_comps(qq, 1)).data,comparison(pairwise_comps(qq, 2)).data,1000);% permTest_array(ERT_test.Cp_off1,ERT_test.Cm_off3,1000);
    diff_bCI{qq} = boot_diffCI(comparison(pairwise_comps(qq, 1)).data,comparison(pairwise_comps(qq, 2)).data,1000,sig);
    [adjLCI,adjUCI] = CIadjust(diff_bCI{1, qq}(1,:),diff_bCI{1, qq}(2,:),[],comparison(pairwise_comps(qq, 1)).uv.n_Cp,2);
    diff_bCIexp{qq} = [adjLCI;adjUCI];
    clear adjLCI adjUCI

    diff_tCI_sig = NaN(1,comparison(pairwise_comps(qq, 1)).uv.ev_win);
    diff_tCI_sig_idx = ttest2(comparison(pairwise_comps(qq, 1)).mean_Cp,comparison(pairwise_comps(qq, 2)).mean_Cp);
    diff_tCI_sig_idx = find(diff_tCI_sig_idx == 1);
    diff_tCI_consec = consec_idx(diff_tCI_sig_idx,consec_thresh);
    diff_tCI_sig(qq, diff_tCI_sig_idx(diff_tCI_consec)) = sig_plot_level_v2;
    clear diff_tCI_sig_idx diff_tCI_consec


    diff_bCIexp_sig(qq, :) = NaN(1,comparison(pairwise_comps(qq, 1)).uv.ev_win);
    diff_bCIexp_sig_idx = find((diff_bCIexp{1, qq}(1,:) > 0) | (diff_bCIexp{1, qq}(2,:) < 0));
    diff_bCIexp_consec = consec_idx(diff_bCIexp_sig_idx,consec_thresh);
    diff_bCIexp_sig(qq, diff_bCIexp_sig_idx(diff_bCIexp_consec)) = sig_plot_level_v2;
    clear diff_bCIexp_sig_idx diff_bCIexp_consec

    %Permutation test
    perm_p_sig(qq,:) = NaN(1,comparison(pairwise_comps(qq, 1)).uv.ev_win);
    perm_p_sig_idx = find(perm_p(qq, :) < sig);
    perm_p_consec = consec_idx(perm_p_sig_idx,consec_thresh);
    perm_p_sig(qq, perm_p_sig_idx(perm_p_consec)) = sig_plot_level_v2;
    clear perm_p_sig_idx perm_p_consec
end




%% Plot permutation test data

% create a label for each comparison to be made
% if size(arg_string_other, 1) > 2
%     if contains(arg_string_other, 'block', 'IgnoreCase', true)
%         for zz = 1:size(arg_string_block, 1)
%             for hh = 1:size(pairwise_comps, 2)
%                 comparison_labels(zz, hh) = arg_string_block(contains(arg_string_block, string(pairwise_comps(zz, hh))));
%             end
%         end
%     elseif ~contains(arg_string_other, 'block', 'IgnoreCase', true)
%         for zz = 1:size(arg_string_other, 1)
%             for hh = 1:size(pairwise_comps, 2)-1
%                 comparison_labels(hh, zz) = arg_string_other(zz);
%             end
%         end
%     end
% elseif size(arg_string_other,1) <= 2
%     for zz = 1:size(arg_string_other, 1)
%         for hh = 1:size(pairwise_comps, 2)-1
%             comparison_labels(hh, zz) = arg_string_other(zz);
%         end
%     end
% 
% end
% % combine the comparisons made above to directly display what is being
% % % compared
% comparison_labels_join = join(comparison_labels, ' vs ');
% 
% num_trials_sum_string = string(num_trials_sum)';
% 
% num_sessions_sum_strong = string(num_sessions_sum)';
% cmb_strings = append(arg_string_other, ' ','num events=', num_trials_sum_string, ' ', 'num sessions=', num_sessions_sum_strong);
% 

% figure; hold on
% for vv = 1:size(pairwise_comps, 1)
% 
%     if size(comparison,2) == 2
%         for qq = 1:size(comparison,2)
% 
%             comp_mean(qq,:) = comparison(qq).mean_Cp;
%             comp_sem(qq,:) = comparison(qq).sem_Cp;
% 
%             plot(timeline,comp_mean(qq,:),'Color',col_rep(qq))
%             
%             errorplot3(comp_mean(qq,:)-comp_sem(qq,:),comp_mean(qq,:)+comp_sem(qq,:),[-8 8],col_rep(qq),.15)
% 
%             plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.')
%             text(xlims(1),sig_plot_level_v2(vv), comparison_labels_join(vv),'Color',col_rep(vv+1),'FontSize', 6);
%             lgd = legend(comparison_labels)
%             
%             
%         end
%     elseif size(comparison,2) > 2
%             comp_mean(vv,:) = comparison(vv).mean_Cp;
%             comp_sem(vv,:) = comparison(vv).sem_Cp;
% 
%             plot(timeline,comp_mean(vv,:),'Color',col_rep(vv))
%             errorplot3(comp_mean(vv,:)-comp_sem(vv,:),comp_mean(vv,:)+comp_sem(vv,:),[-8 8],col_rep(vv),.15)
% 
%             plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.')
%             text(xlims(1),sig_plot_level_v2(vv),comparison_labels_join(vv),'Color',col_rep(vv+1), 'FontSize', 6);
%             
%     end
%     plot([0 0],ylim,'k:')
%     plot(xlim,[0 0],'k--')
%     
%     xlim(xlims);
% 
% end
% 
% figure; hold on
% width = 300; % Width of the figure
% height = 600; % Height of the figure (width is half of height)
% set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
% xlim(xlims);
% % Set X-axis ticks
% set(gca);
% ytickformat('%.2f');
% 
% for vv = 1:size(pairwise_comps, 1)
% 
%     if size(comparison,2) == 2
%         for qq = 1:size(comparison,2)
% 
%             comp_mean(qq,:) = comparison(qq).mean_Cp;
%             comp_sem(qq,:) = comparison(qq).sem_Cp;
%             % z = shadedErrorBar(timeline, comp_mean(qq,:), comp_sem(qq,:), 'lineProps', {'color', col_rep(qq)});
%             z = shadedErrorBar(timeline, comp_mean(qq,:), comp_sem(qq,:), 'lineProps', {'color', col_rep(qq)});
%             % set(z.edge,'LineStyle','none')
%             % z.mainLine.LineWidth = 3;            
%             % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'westoutside')
% %             plot(timeline,comp_mean(qq,:),'Color',col_rep(qq))
% 
% %             errorplot3(comp_mean(qq,:)-comp_sem(qq,:),comp_mean(qq,:)+comp_sem(qq,:),[-8 8],col_rep(qq),.15)
% 
%             f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'linestyle','-');
%             % f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.');
%             f.Annotation.LegendInformation.IconDisplayStyle = 'off';
%             % text(xlims(1),sig_plot_level_v2(vv), comparison_labels_join(vv),'Color',col_rep(vv+1),'FontSize', 6);
% 
% 
% 
%         end
% %     elseif size(comparison,2) > 2
% %             comp_mean(vv,:) = comparison(vv).mean_Cp;
% %             comp_sem(vv,:) = comparison(vv).sem_Cp;
% %             z = shadedErrorBar(timeline, comp_mean(vv,:), comp_sem(vv,:), 'lineProps', {'color', col_rep(qq+1)});
% %             % set(z.edge,'LineStyle','none')
% %             % z.mainLine.LineWidth = 3;
% % 
% %             % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'westoutside')
% % %             legend(arg_string_other, 'AutoUpdate','off')
% % %             plot(timeline,comp_mean(vv,:),'Color',col_rep(vv))
% % %             errorplot3(comp_mean(vv,:)-comp_sem(vv,:),comp_mean(vv,:)+comp_sem(vv,:),[-8 8],col_rep(vv),.15)
% % 
% %             f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'linestyle','-');
% %             % f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.');
% %             f.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % 
% %             % text(xlims(1),sig_plot_level_v2(vv),comparison_labels_join(vv),'Color',col_rep(vv+1), 'FontSize', 6);
% % 
% %     end
% %
%     elseif size(comparison,2) > 2
%         for qq = 1:size(comparison,2)
% 
%             comp_mean(qq,:) = comparison(qq).mean_Cp;
%             comp_sem(qq,:) = comparison(qq).sem_Cp;
%             % z = shadedErrorBar(timeline, comp_mean(qq,:), comp_sem(qq,:), 'lineProps', {'color', col_rep(qq)});
%             z = shadedErrorBar(timeline, comp_mean(qq,:), comp_sem(qq,:), 'lineProps', {'color', col_rep(qq)});
%             % set(z.edge,'LineStyle','none')
%             % z.mainLine.LineWidth = 3;
%             % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'westoutside')
%             %             plot(timeline,comp_mean(qq,:),'Color',col_rep(qq))
% 
%             %             errorplot3(comp_mean(qq,:)-comp_sem(qq,:),comp_mean(qq,:)+comp_sem(qq,:),[-8 8],col_rep(qq),.15)
% 
%             % f = plot(timeline,perm_p_sig(qq,:),'Color',col_rep(qq),'linestyle','-');
%             % f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.');
%             % f.Annotation.LegendInformation.IconDisplayStyle = 'off';
%             % text(xlims(1),sig_plot_level_v2(vv), comparison_labels_join(vv),'Color',col_rep(vv+1),'FontSize', 6);
%             f = plot(timeline,perm_p_sig(qq,:),'Color',col_rep(qq),'linestyle','-');
%             % f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.');
%             f.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% 
%         end
% 
%     end
% 
% 
% end
% 
%     xlim(xlims);
%     ylim(ylims)
%     yticks(unique([yticks, ylims(2)]));
%     % set(gcf, 'position', [10, 10, 900, 600]);
%     % fontsize(16, "pixels")
%     title('perm test')
%     % set(gcf, 'resize', 'off' );
% 
% 
% 
    %% Plot permutation test data

% create a label for each comparison to be made
% if size(arg_string_other, 1) > 2
%     if contains(arg_string_other, 'block', 'IgnoreCase', true)
%         for zz = 1:size(arg_string_block, 1)
%             for hh = 1:size(pairwise_comps, 2)
%                 comparison_labels(zz, hh) = arg_string_block(contains(arg_string_block, string(pairwise_comps(zz, hh))));
%             end
%         end
%     elseif ~contains(arg_string_other, 'block', 'IgnoreCase', true)
%         for zz = 1:size(arg_string_other, 1)
%             for hh = 1:size(pairwise_comps, 2)-1
%                 comparison_labels(hh, zz) = arg_string_other(zz);
%             end
%         end
%     end
% elseif size(arg_string_other,1) <= 2
%     for zz = 1:size(arg_string_other, 1)
%         for hh = 1:size(pairwise_comps, 2)-1
%             comparison_labels(hh, zz) = arg_string_other(zz);
%         end
%     end
% 
% end
% % combine the comparisons made above to directly display what is being
% % % compared
% comparison_labels_join = join(comparison_labels, ' vs ');
% 
% num_trials_sum_string = string(num_trials_sum)';
% 
% num_sessions_sum_strong = string(num_sessions_sum)';
% cmb_strings = append(arg_string_other, ' ','num events=', num_trials_sum_string, ' ', 'num sessions=', num_sessions_sum_strong);
% 

% figure; hold on
% for vv = 1:size(pairwise_comps, 1)
% 
%     if size(comparison,2) == 2
%         for qq = 1:size(comparison,2)
% 
%             comp_mean(qq,:) = comparison(qq).mean_Cp;
%             comp_sem(qq,:) = comparison(qq).sem_Cp;
% 
%             plot(timeline,comp_mean(qq,:),'Color',col_rep(qq))
%             
%             errorplot3(comp_mean(qq,:)-comp_sem(qq,:),comp_mean(qq,:)+comp_sem(qq,:),[-8 8],col_rep(qq),.15)
% 
%             plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.')
%             text(xlims(1),sig_plot_level_v2(vv), comparison_labels_join(vv),'Color',col_rep(vv+1),'FontSize', 6);
%             lgd = legend(comparison_labels)
%             
%             
%         end
%     elseif size(comparison,2) > 2
%             comp_mean(vv,:) = comparison(vv).mean_Cp;
%             comp_sem(vv,:) = comparison(vv).sem_Cp;
% 
%             plot(timeline,comp_mean(vv,:),'Color',col_rep(vv))
%             errorplot3(comp_mean(vv,:)-comp_sem(vv,:),comp_mean(vv,:)+comp_sem(vv,:),[-8 8],col_rep(vv),.15)
% 
%             plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.')
%             text(xlims(1),sig_plot_level_v2(vv),comparison_labels_join(vv),'Color',col_rep(vv+1), 'FontSize', 6);
%             
%     end
%     plot([0 0],ylim,'k:')
%     plot(xlim,[0 0],'k--')
%     
%     xlim(xlims);
% 
% end

figure; hold on
width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim(xlims);
% Set X-axis ticks
set(gca);
ytickformat('%.2f');

for vv = 1:size(pairwise_comps, 2)

    if size(comparison,2) == 2


            comp_mean(vv,:) = comparison(vv).mean_Cp;
            comp_sem(vv,:) = comparison(vv).sem_Cp;
            % z = shadedErrorBar(timeline, comp_mean(qq,:), comp_sem(qq,:), 'lineProps', {'color', col_rep(qq)});
            z = shadedErrorBar(timeline, comp_mean(vv,:), comp_sem(vv,:), 'lineProps', {'color', col_rep(vv)});
            % set(z.edge,'LineStyle','none')
            % z.mainLine.LineWidth = 3;            
            % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'westoutside')
%             plot(timeline,comp_mean(qq,:),'Color',col_rep(qq))
            
%             errorplot3(comp_mean(qq,:)-comp_sem(qq,:),comp_mean(qq,:)+comp_sem(qq,:),[-8 8],col_rep(qq),.15)

            f = plot(timeline,perm_p_sig,'Color',col_rep(vv),'linestyle','-');
            % f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.');
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % text(xlims(1),sig_plot_level_v2(vv), comparison_labels_join(vv),'Color',col_rep(vv+1),'FontSize', 6);
            

%     elseif size(comparison,2) > 2
%             comp_mean(vv,:) = comparison(vv).mean_Cp;
%             comp_sem(vv,:) = comparison(vv).sem_Cp;
%             z = shadedErrorBar(timeline, comp_mean(vv,:), comp_sem(vv,:), 'lineProps', {'color', col_rep(qq+1)});
%             % set(z.edge,'LineStyle','none')
%             % z.mainLine.LineWidth = 3;
% 
%             % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'westoutside')
% %             legend(arg_string_other, 'AutoUpdate','off')
% %             plot(timeline,comp_mean(vv,:),'Color',col_rep(vv))
% %             errorplot3(comp_mean(vv,:)-comp_sem(vv,:),comp_mean(vv,:)+comp_sem(vv,:),[-8 8],col_rep(vv),.15)
% 
%             f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'linestyle','-');
%             % f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.');
%             f.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
%             % text(xlims(1),sig_plot_level_v2(vv),comparison_labels_join(vv),'Color',col_rep(vv+1), 'FontSize', 6);
% 
%     end
%
    elseif size(comparison,2) > 2


            comp_mean(vv,:) = comparison(vv).mean_Cp;
            comp_sem(vv,:) = comparison(vv).sem_Cp;
            % z = shadedErrorBar(timeline, comp_mean(qq,:), comp_sem(qq,:), 'lineProps', {'color', col_rep(qq)});
            z = shadedErrorBar(timeline, comp_mean(vv,:), comp_sem(vv,:), 'lineProps', {'color', col_rep(vv)});
            % set(z.edge,'LineStyle','none')
            % z.mainLine.LineWidth = 3;
            % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'westoutside')
            %             plot(timeline,comp_mean(qq,:),'Color',col_rep(qq))

            %             errorplot3(comp_mean(qq,:)-comp_sem(qq,:),comp_mean(qq,:)+comp_sem(qq,:),[-8 8],col_rep(qq),.15)

            % f = plot(timeline,perm_p_sig(qq,:),'Color',col_rep(qq),'linestyle','-');
            % f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.');
            % f.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % text(xlims(1),sig_plot_level_v2(vv), comparison_labels_join(vv),'Color',col_rep(vv+1),'FontSize', 6);
            f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv),'linestyle','-');
            % f = plot(timeline,perm_p_sig(vv,:),'Color',col_rep(vv+1),'Marker','.');
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';




    end


end

    xlim(xlims);
    ylim(ylims)
    yticks(unique([yticks, ylims(2)]));
    %yticks(ylims/2);
    % set(gcf, 'position', [10, 10, 900, 600]);
    % fontsize(16, "pixels")
    title('perm test')
    % set(gcf, 'resize', 'off' );
%% Plot diff bCI


figure; hold on


for vv = 1:size(pairwise_comps, 1)

    if size(comparison,2) == 2
        for qq = 1:size(comparison,2)

            comp_mean(qq,:) = comparison(qq).mean_Cp;
            comp_sem(qq,:) = comparison(qq).sem_Cp;
            z = shadedErrorBar(timeline, comp_mean(qq,:), comp_sem(qq,:), 'lineProps', {'color', col_rep(qq)});
            % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'northoutside')            

%             plot(timeline,comp_mean(qq,:),'Color',col_rep(qq))
%             errorplot3(comp_mean(qq,:)-comp_sem(qq,:),comp_mean(qq,:)+comp_sem(qq,:),[-8 8],col_rep(qq),.15)

            f = plot(timeline, diff_bCIexp_sig(vv,:),'Color',col_rep(vv),'linestyle','-');
            % f = plot(timeline, diff_bCIexp_sig(vv,:),'Color',col_rep(vv),'Marker','.');
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % text(xlims(1),sig_plot_level_v2(vv), comparison_labels_join(vv),'Color',col_rep(vv),'FontSize', 6);
            
        end
    elseif size(comparison,2) > 2
            comp_mean(vv,:) = comparison(vv).mean_Cp;
            comp_sem(vv,:) = comparison(vv).sem_Cp;
            z = shadedErrorBar(timeline, comp_mean(vv,:), comp_sem(vv,:), 'lineProps', {'color', col_rep(vv)});
            % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'northoutside')  
%             plot(timeline,comp_mean(vv,:),'Color',col_rep(vv))
%             errorplot3(comp_mean(vv,:)-comp_sem(vv,:),comp_mean(vv,:)+comp_sem(vv,:),[-8 8],col_rep(vv),.15)

            f = plot(timeline,diff_bCIexp_sig(vv,:),'Color',col_rep(vv),'linestyle','-');
            % f = plot(timeline,diff_bCIexp_sig(vv,:),'Color',col_rep(vv),'Marker','.');
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % text(xlims(1),sig_plot_level_v2(vv),comparison_labels_join(vv),'Color',col_rep(vv), 'FontSize', 6);

    end

    ylabel('z-scored dF/F', 'FontSize', 12);
    xlabel('Time from choice (s)');
    xlim(xlims);
    ylim(ylims)
    yticks(unique([yticks, ylims(2)]));
    set(gcf, 'position', [10, 10, 400, 800]);
    title('diff bCI')
    
end
    
%% Plot bCI



figure; hold on

width = 300; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
% xtickformat('%.2f');
ytickformat('%.2f');


for vv = 1:size(comparison, 2)

    if size(comparison,2) == 2
        for qq = 1:size(comparison,2)

            comp_mean(qq,:) = comparison(qq).mean_Cp;
            comp_sem(qq,:) = comparison(qq).sem_Cp;
            z = shadedErrorBar(timeline, comp_mean(qq,:), comp_sem(qq,:), 'lineProps', {'color', col_rep(qq)});
            % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'northoutside')  

%             plot(timeline,comp_mean(qq,:),'Color',col_rep(qq))
%             errorplot3(comp_mean(qq,:)-comp_sem(qq,:),comp_mean(qq,:)+comp_sem(qq,:),[-8 8],col_rep(qq),.15)

            f = plot(timeline,comparison(vv).Cp_bCIexp_sig,'Color',col_rep(vv),'linestyle','-');
            % f = plot(timeline,comparison(vv).Cp_bCIexp_sig,'Color',col_rep(vv),'Marker','.');
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';
%             text(xlims(1),sig_plot_level_v2(vv), comparison_labels(vv),'Color',col_rep(vv),'FontSize', 6);
            % text(xlims(1),sig_plot_level_v2(vv), arg_string_other(vv),'Color',col_rep(vv),'FontSize', 6);

        end
    elseif size(comparison,2) > 2

            comp_mean(vv,:) = comparison(vv).mean_Cp;
            comp_sem(vv,:) = comparison(vv).sem_Cp;
            z = shadedErrorBar(timeline, comp_mean(vv,:), comp_sem(vv,:), 'lineProps', {'color', col_rep(vv)});
            % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'northoutside')

%             plot(timeline,comp_mean(vv,:),'Color',col_rep(vv))
%             errorplot3(comp_mean(vv,:)-comp_sem(vv,:),comp_mean(vv,:)+comp_sem(vv,:),[-8 8],col_rep(vv),.15)

            f = plot(timeline,comparison(vv).Cp_bCIexp_sig,'Color',col_rep(vv),'linestyle','-');
            % f = plot(timeline,comparison(vv).Cp_bCIexp_sig,'Color',col_rep(vv),'Marker','.');
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % text(xlims(1),sig_plot_level_v2(vv), arg_string_other(vv),'Color',col_rep(vv), 'FontSize', 6);

    end
    % plot([0 0],ylim,'k:')
    % plot(xlim,[0 0],'k')
    % title('bootstrapped CI: 95% CI does not include 0 dF/F')
    % ylabel('z-scored dF/F', 'FontSize', 12);
    % xlabel('Time from choice (s)');
    xlim(xlims);
    ylim(ylims)
    yticks(unique([yticks, ylims(2)]));

    % Set X-axis ticks

    title('bCI')
    
end


%% Plot tCI




figure; hold on


for vv = 1:size(comparison, 2)

    if size(comparison,2) == 2
        for qq = 1:size(comparison,2)

            comp_mean(qq,:) = comparison(qq).mean_Cp;
            comp_sem(qq,:) = comparison(qq).sem_Cp;
            z = shadedErrorBar(timeline, comp_mean(qq,:), comp_sem(qq,:), 'lineProps', {'color', col_rep(qq)});
            % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'northoutside')  

%             plot(timeline,comp_mean(qq,:),'Color',col_rep(qq))
%             errorplot3(comp_mean(qq,:)-comp_sem(qq,:),comp_mean(qq,:)+comp_sem(qq,:),[-8 8],col_rep(qq),.15)

            % f = plot(timeline,comparison(vv).Cp_tCI_sig,'Color',col_rep(vv),'Marker','.');
            f = plot(timeline,comparison(vv).Cp_tCI_sig,'Color',col_rep(vv),'linestyle','-');
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % text(xlims(1),sig_plot_level_v2(vv), arg_string_other(vv),'Color',col_rep(vv),'FontSize', 6);

        end
    elseif size(comparison,2) > 2
            comp_mean(vv,:) = comparison(vv).mean_Cp;
            comp_sem(vv,:) = comparison(vv).sem_Cp;

            z = shadedErrorBar(timeline, comp_mean(vv,:), comp_sem(vv,:), 'lineProps', {'color', col_rep(vv)});
            % legend(cmb_strings, 'AutoUpdate','off', 'Location', 'northoutside')

%             plot(timeline,comp_mean(vv,:),'Color',col_rep(vv))
%             errorplot3(comp_mean(vv,:)-comp_sem(vv,:),comp_mean(vv,:)+comp_sem(vv,:),[-8 8],col_rep(vv),.15)

            % f = plot(timeline,comparison(vv).Cp_tCI_sig,'Color',col_rep(vv),'Marker','.');
            f = plot(timeline,comparison(vv).Cp_tCI_sig,'Color',col_rep(vv),'linestyle','-');
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % text(xlims(1),sig_plot_level_v2(vv), arg_string_other(vv),'Color',col_rep(vv), 'FontSize', 6);

    end

    title('Parametric t interval CI: 95% CI does not include 0 dF/F')
    ylabel('z-scored dF/F', 'FontSize', 12);
    xlabel('Time from choice (s)');
    xlim(xlims);
    ylim(ylims)
    yticks(unique([yticks, ylims(2)]));
    title('tCI')

end

%%
% mean(vv) = comparison(1).mean_Cp;
% mean_Cp = comparison(2).mean_Cp;
% sem_Cm = comparison(1).sem_Cp;
% sem_Cp = comparison(2).sem_Cp;
% 
% figure; hold on
% plot(timeline,mean_Cm,'Color',col_rep(3))
% errorplot3(mean_Cm-sem_Cm,mean_Cm+sem_Cm,[-8 8],col_rep(3),.15)
% 
% plot(timeline,mean_Cp,'Color',col_rep(2))
% errorplot3(mean_Cp-sem_Cp,mean_Cp+sem_Cp,[-8 8],col_rep(2),.15)

% %Plor tCI sig
% plot(timeline,comparison(1).Cp_tCI_sig,'Color',col_rep(2),'Marker','.')
% text(xlims(1),sig_plot_level(1),'\bf CS+ tCI','Color',col_rep(2));
% plot(timeline,comparison(2).Cp_tCI_sig,'Color',col_rep(3),'Marker','.')
% text(xlims(1),sig_plot_level(2),'\bf CS- tCI','Color',col_rep(3));
% plot(timeline,diff_tCI_sig,'Color',col_rep(4),'Marker','.')
% text(xlims(1),sig_plot_level(3),'\bf Diff tCI','Color',col_rep(4));
% 
% %Plot bCI sig
% plot(timeline,comparison(1).Cp_bCIexp_sig,'Color',col_rep(2),'Marker','.')
% text(xlims(1),sig_plot_level(4),'\bf CS+ bCI','Color',col_rep(2));
% plot(timeline,comparison(2).Cp_bCIexp_sig,'Color',col_rep(3),'Marker','.')
% text(xlims(1),sig_plot_level(5),'\bf CS- bCI','Color',col_rep(3));
% plot(timeline,diff_bCIexp_sig,'Color',col_rep(4),'Marker','.')
% text(xlims(1),sig_plot_level(6),'\bf Diff bCI','Color',col_rep(4));
% 
% %Plot permutation test sig
% % plot(timeline,perm_p_sig,'Color',col_rep(1),'Marker','.')
% % text(xlims(1),sig_plot_level(7),'\bf Perm','Color',col_rep(1));
% 
% plot([-0.5 -0.5],ylim,'k:')
% plot(xlim,[0 0],'k--')
% 
% xlim(xlims);

