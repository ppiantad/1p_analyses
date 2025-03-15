function [comparison] = bCI_and_tCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1, x_limits, y_limits)

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
consec_thresh = 10; % 1017.3Hz sample rate / 3Hz filter %340 PRD used 340 because his data WERE NOT DOWNSAMPLED (e.g., they were 1018 samples per sec) our data are 30 samples per sec or so.
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
    comparison(ii).Cp_tCI_sig(comparison(ii).sig_idx_tCI(comparison(ii).consec_tCI)) = sig_plot_level_v2(ii);

    %bCI
    comparison(ii).Cp_bCIexp_sig = NaN(1,comparison(ii).uv.ev_win);
    comparison(ii).sig_idx_bCI = find((comparison(ii).Cp_bCIexp(1,:) > bCI_tCI_CI_threshold) | (comparison(ii).Cp_bCIexp(2,:) < bCI_tCI_CI_threshold));
    comparison(ii).consec_bCI = consec_idx(comparison(ii).sig_idx_bCI,consec_thresh);
    comparison(ii).Cp_bCIexp_sig(comparison(ii).sig_idx_bCI(comparison(ii).consec_bCI)) = sig_plot_level_v2(ii);



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





    
%% Plot bCI



figure; hold on

width = 300; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
% xtickformat('%.2f');
ytickformat('%.2f');


for vv = 1:size(comparison, 2)

    if size(comparison,2) == 1
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

    if size(comparison,2) == 1
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
    plot([0 0],ylim,'k:')
    plot(xlim,[0 0],'k--')
    title('Parametric t interval CI: 95% CI does not include 0 dF/F')
    ylabel('z-scored dF/F', 'FontSize', 12);
    xlabel('Time from choice (s)');
    xlim(xlims);
    ylim(ylims);
    title('tCI')
    yticks(unique([yticks, ylims(2)]));
end

