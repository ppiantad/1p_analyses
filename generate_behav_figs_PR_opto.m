% iter = iter+1;
% neuron_num = 0;
% these are mice that did not complete the entire session - kinda have to
% toss them to do some comparisons during RDT

% final_behavior = final_SLEAP; % for hM4Di data;

session_to_analyze = {'PR_D1', 'PR_D2'};

if strcmp('RM_D1', session_to_analyze)| strcmp('RDT_D1', session_to_analyze) | strcmp('Pre_RDT_RM', session_to_analyze)
    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_29', 'BLA_Insc_38', 'BLA_Insc_39', 'BLA_Insc_13'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_behavior, fieldsToRemove{i})
            final_behavior = rmfield(final_behavior, fieldsToRemove{i});
        end
    end
elseif strcmp('RDT_D2', session_to_analyze)

    fieldsToRemove = {'BLA_Insc_28', 'BLA_Insc_39'};

    for i = 1:length(fieldsToRemove)
        if isfield(final_behavior, fieldsToRemove{i})
            final_behavior = rmfield(final_behavior, fieldsToRemove{i});
        end
    end
end


%%


animalIDs = (fieldnames(final_behavior));



possible_effectors = {'stGtACR', 'PdCO', 'ChrimsonR', 'hM4Di'}
possible_controls = {'mCherry', 'EGFP'}


risk_table = table;


if exist('hM4Di_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) all(ismember(session_to_analyze, sessions)), valid_sessions);
    valid_animalIDs = hM4Di_IDs(valid_mice);
    treatment_group_all = hM4Di_treatment_groups;
    
elseif exist('stGtACR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = stGtACR_IDs(valid_mice);
    valid_animalIDs = valid_animalIDs(~strcmp(valid_animalIDs, 'RRD391'));


elseif exist('PdCO_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = PdCO_IDs(valid_mice);

elseif exist('ChrimsonR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = ChrimsonR_IDs(valid_mice);

elseif exist('treatment_groups', 'var') == 1
    valid_mice = cellfun(@(sessions) any(strcmp(sessions, session_to_analyze)), valid_sessions);
    valid_animalIDs = mouse_IDs(valid_mice);

elseif exist('NAcSh_eNpHR_treatment_groups', 'var') == 1
    valid_mice = cellfun(@(mice) any(strcmp(mice, animalIDs)), mouse_IDs);
    valid_animalIDs = mouse_IDs(valid_mice);
    treatment_group_all = NAcSh_eNpHR_treatment_groups;
else
    valid_animalIDs = animalIDs;

end

valid_session_to_add = 0;

for ii = 1:size(valid_animalIDs,1) % 1:size(fieldnames(final),1)
    currentanimal = char(valid_animalIDs(ii));
    treatment_group = treatment_group_all{ii};
    for hh = 1:size(session_to_analyze, 2)


        if isfield(final_behavior.(currentanimal), session_to_analyze{hh})
            % if contains(session_to_analyze, 'CNO')
            %     BehavData = final_behavior.(currentanimal).(session_to_analyze).BehavData;
            % else
            %     BehavData = final_behavior.(currentanimal).(session_to_analyze).uv.BehavData;
            % end
            valid_session_to_add = valid_session_to_add + 1;

            animalID_list(valid_session_to_add) = {currentanimal};
            laser_treatment_current = laser_treatment_order_for_PR{ii}{1, hh};
            BehavData = final_behavior.(currentanimal).(session_to_analyze{hh}).uv.BehavData;
            % trim BehavData down to the max session length (3600s), mice
            % occasionally get 1.5 hr session due to user error (rare)
            max_session_time_seconds = 3600;
            valid_press_times = BehavData.PressTime <= 3600;
            BehavData = BehavData(valid_press_times, :);
            total_presses = BehavData.Trial(end);
            rewards_obtained = sum(BehavData.collectTrial);

            final_ratio_attained = max(BehavData.Ratio - 1);
            final_ratio_reached = max(BehavData.Ratio);
            BehavData.collect_latency = BehavData.collectionTime - BehavData.rewDeliveryTime;
            mean_collect_latency = mean(BehavData.collect_latency(BehavData.collectTrial == 1));


            % block_1_choice_latency_all = mean(BehavData.choice_latency(BehavData.Block == 1 & (BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3)));
            % block_2_choice_latency_all = mean(BehavData.choice_latency(BehavData.Block == 2 & (BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3)));
            % block_3_choice_latency_all = mean(BehavData.choice_latency(BehavData.Block == 3 & (BehavData.bigSmall == 1.2 | BehavData.bigSmall == 0.3)));
            %
            % block_1_large_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 1));
            % block_2_large_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 2));
            % block_3_large_collect_latency = mean(BehavData.collect_latency(BehavData.bigSmall == 1.2 & BehavData.Block == 3));

            session_length = BehavData.PressTime(end);


            risk_table(valid_session_to_add,:) = table(...
                [final_ratio_attained],...
                [final_ratio_reached],...
                [total_presses],...
                [rewards_obtained],...
                [session_length],...
                string(currentanimal),...
                string(treatment_group),...
                string(laser_treatment_current),...
                [mean_collect_latency]...
                );



        elseif ~isfield(final_behavior.(currentanimal), session_to_analyze)
            valid_session_to_add = valid_session_to_add + 1;
            animalID_list(valid_session_to_add) = {currentanimal};
            risk_table{valid_session_to_add, :} = NaN;
           
        end

    end
end

risk_table = renamevars(risk_table,["Var6", "Var7", "Var8"], ["animalID", "TreatmentCondition", "LaserTreatment"]);
% animal_id_table =  array2table(animalID_list');
% risk_table = [animal_id_table risk_table];








%% for D1 vs D2

%% 06/27/2025 data pulled from RRD Photometry and Opto Progressive Ratio v2 - need to get raw data from PCs in FLAC

PR_presses_D1_cre_laser_ON = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'D1') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_presses_D1_cre_laser_OFF = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'D1') & strcmp(risk_table.LaserTreatment, 'Laser Off'));

PR_presses_A2A_cre_laser_ON = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'A2A') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_presses_A2A_cre_laser_OFF = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'A2A') & strcmp(risk_table.LaserTreatment, 'Laser Off'));

PR_presses_Control_laser_ON = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'Control') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_presses_Control_laser_OFF = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'Control') & strcmp(risk_table.LaserTreatment, 'Laser Off'));

PR_ratios_D1_cre_laser_ON = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'D1') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_ratios_D1_cre_laser_OFF = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'D1') & strcmp(risk_table.LaserTreatment, 'Laser Off'));

PR_ratios_A2A_cre_laser_ON = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'A2A') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_ratios_A2A_cre_laser_OFF = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'A2A') & strcmp(risk_table.LaserTreatment, 'Laser Off'));

PR_ratios_Control_laser_ON = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'Control') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_ratios_Control_laser_OFF = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'Control') & strcmp(risk_table.LaserTreatment, 'Laser Off'));




% Compute means and SEMs for each group
% PR Presses
mean_presses = [mean(PR_presses_D1_cre_laser_ON),...
                mean(PR_presses_D1_cre_laser_OFF),...
                mean(PR_presses_A2A_cre_laser_ON),...
                mean(PR_presses_A2A_cre_laser_OFF),...
                mean(PR_presses_Control_laser_ON),...
                mean(PR_presses_Control_laser_OFF),...
                ];
sem_presses = [std(PR_presses_D1_cre_laser_ON)/sqrt(length(PR_presses_D1_cre_laser_ON)), ...
               std(PR_presses_D1_cre_laser_OFF)/sqrt(length(PR_presses_D1_cre_laser_OFF)), ...
               std(PR_presses_A2A_cre_laser_ON)/sqrt(length(PR_presses_A2A_cre_laser_ON)), ...
               std(PR_presses_A2A_cre_laser_OFF)/sqrt(length(PR_presses_A2A_cre_laser_OFF)), ...
               std(PR_presses_Control_laser_ON)/sqrt(length(PR_presses_Control_laser_ON)), ...
               std(PR_presses_Control_laser_OFF)/sqrt(length(PR_presses_Control_laser_OFF)), ...
               ];

% PR Ratios
mean_ratios = [mean(PR_ratios_D1_cre_laser_ON),...
                mean(PR_ratios_D1_cre_laser_OFF),...
                mean(PR_ratios_A2A_cre_laser_ON),...
                mean(PR_ratios_A2A_cre_laser_OFF),...
                mean(PR_ratios_Control_laser_ON),...
                mean(PR_ratios_Control_laser_OFF),...
                ];
sem_ratios = [std(PR_ratios_D1_cre_laser_ON)/sqrt(length(PR_ratios_D1_cre_laser_ON)), ...
               std(PR_ratios_D1_cre_laser_OFF)/sqrt(length(PR_ratios_D1_cre_laser_OFF)), ...
               std(PR_ratios_A2A_cre_laser_ON)/sqrt(length(PR_ratios_A2A_cre_laser_ON)), ...
               std(PR_ratios_A2A_cre_laser_OFF)/sqrt(length(PR_ratios_A2A_cre_laser_OFF)), ...
               std(PR_ratios_Control_laser_ON)/sqrt(length(PR_ratios_Control_laser_ON)), ...
               std(PR_ratios_Control_laser_OFF)/sqrt(length(PR_ratios_Control_laser_OFF)), ...
               ];

% Define colors for each group
colors = [0.6, 0.8, 1.0;    % Light blue for Control Laser On
          0.2, 0.4, 0.8;    % Dark blue for Control Laser Off
          1.0, 0.7, 0.7;    % Light red for ChrimsonR Laser On
          0.8, 0.2, 0.2;    % Dark red for ChrimsonR Laser Off
          1.0, 0.7, 0.7;    % Light red for ChrimsonR Laser On  
          0.8, 0.2, 0.2;    % Dark red for ChrimsonR Laser Off
          ];   

% Create figure
figure;

% First subplot - Final PR Ratio Attained
subplot(1,2,1);
b1 = bar(1:6, mean_ratios, 0.6, 'FaceColor', 'flat');
% Set colors for each bar
for i = 1:6
    b1.CData(i,:) = colors(i,:);
end
hold on;

% Add error bars
errorbar(1:6, mean_ratios, sem_ratios, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Add individual data points with jitter
all_ratio_data = {PR_ratios_D1_cre_laser_ON,...
                    PR_ratios_D1_cre_laser_OFF,...
                    PR_ratios_A2A_cre_laser_ON,...
                    PR_ratios_A2A_cre_laser_OFF,...
                    PR_ratios_Control_laser_ON,...
                    PR_ratios_Control_laser_OFF,...                    
                   };
for i = 1:6
    x_pos = i + 0.1 * (rand(size(all_ratio_data{i})) - 0.5); % Add jitter
    scatter(x_pos, all_ratio_data{i}, 40, 'k', 'filled');
end

hold off;
xlim([0.5 6.5]);
ylim([0 20]);
xticks(1:3);
xticklabels({'Control\nLaser On', 'Control\nLaser Off', 'ChrimsonR\nLaser On', 'ChrimsonR\nLaser Off'});
yticks([0:5:20]);
ylabel('Final PR ratio attained');

% Second subplot - Total Presses
subplot(1,2,2);
b2 = bar(1:6, mean_presses, 0.6, 'FaceColor', 'flat');
% Set colors for each bar
for i = 1:6
    b2.CData(i,:) = colors(i,:);
end
hold on;

% Add error bars
errorbar(1:6, mean_presses, sem_presses, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);


% Add individual data points with jitter
all_presses_data = {PR_presses_D1_cre_laser_ON,...
                    PR_presses_D1_cre_laser_OFF,...
                    PR_presses_A2A_cre_laser_ON,...
                    PR_presses_A2A_cre_laser_OFF,...
                    PR_presses_Control_laser_ON,...
                    PR_presses_Control_laser_OFF,...                    
                   };
for i = 1:6
    x_pos = i + 0.1 * (rand(size(all_presses_data{i})) - 0.5); % Add jitter
    scatter(x_pos, all_presses_data{i}, 40, 'k', 'filled');
end

hold off;
xlim([0.5 6.5]);
ylim([0 600]);
xticks(1:3);
xticklabels({'Control\nLaser On', 'Control\nLaser Off', 'ChrimsonR\nLaser On', 'ChrimsonR\nLaser Off'});
yticks([0 200 400 600]);
ylabel('Total presses');

% Adjust layout
set(gcf, 'Position', [100, 100, 500, 350]); % Increased width for 4 bars



%% ANOVA Analysis

% Combine all PR ratio data
ratios = [PR_ratio_D1_cre, ...
          PR_ratio_A2A_cre, ...
          mean_PR_ratio_ChrimsonR_laser_on, ...
          mean_PR_ratio_ChrimsonR_laser_off]';

% Define group labels
group_treatment = [repmat({'Control'}, length(PR_ratio_D1_cre), 1); 
                   repmat({'Control'}, length(PR_ratio_A2A_cre), 1); 
                   repmat({'ChrimsonR'}, length(mean_PR_ratio_ChrimsonR_laser_on), 1); 
                   repmat({'ChrimsonR'}, length(mean_PR_ratio_ChrimsonR_laser_off), 1)];

group_laser = [repmat({'On'}, length(PR_ratio_D1_cre), 1); 
               repmat({'Off'}, length(PR_ratio_A2A_cre), 1); 
               repmat({'On'}, length(mean_PR_ratio_ChrimsonR_laser_on), 1); 
               repmat({'Off'}, length(mean_PR_ratio_ChrimsonR_laser_off), 1)];

% Create a table
T_ratio = table(ratios, group_treatment, group_laser, ...
    'VariableNames', {'PR_Ratio', 'Treatment', 'Laser'});

% Run 2-way ANOVA
[p_ratio, tbl_ratio, stats_ratio] = anovan(T_ratio.PR_Ratio, ...
    {T_ratio.Treatment, T_ratio.Laser}, ...
    'model', 'interaction', ...
    'varnames', {'Treatment', 'Laser'});

% Posthoc if interaction is significant
if p_ratio(3) < 0.05
    fprintf('\nSignificant interaction detected (p = %.4f), running Tukey posthoc:\n', p_ratio(3));
    figure;
    c_ratio = multcompare(stats_ratio, 'Dimension', [1 2]);
end

% Run multcompare with group indices [1 2] (Treatment × Laser)
[c_ratio, m_ratio, h_ratio, gnames_ratio] = multcompare(stats_ratio, 'Dimension', [1 2]);

% Display comparisons in readable format
fprintf('\nPosthoc comparisons for PR Ratio:\n');
for i = 1:size(c_ratio,1)
    group1 = gnames_ratio{c_ratio(i,1)};
    group2 = gnames_ratio{c_ratio(i,2)};
    pval   = c_ratio(i,6);  % p-value is in column 6
    fprintf('  %s vs %s: p = %.4f\n', group1, group2, pval);
end

%%
% Combine total presses data
presses = [PR_presses_D1_cre, ...
           PR_presses_A2A_cre, ...
           mean_PR_presses_ChrimsonR_laser_on, ...
           mean_PR_presses_ChrimsonR_laser_off]';

% Define same group labels
group_treatment = [repmat({'Control'}, length(PR_presses_D1_cre), 1); 
                   repmat({'Control'}, length(PR_presses_A2A_cre), 1); 
                   repmat({'ChrimsonR'}, length(mean_PR_presses_ChrimsonR_laser_on), 1); 
                   repmat({'ChrimsonR'}, length(mean_PR_presses_ChrimsonR_laser_off), 1)];

group_laser = [repmat({'On'}, length(PR_presses_D1_cre), 1); 
               repmat({'Off'}, length(PR_presses_A2A_cre), 1); 
               repmat({'On'}, length(mean_PR_presses_ChrimsonR_laser_on), 1); 
               repmat({'Off'}, length(mean_PR_presses_ChrimsonR_laser_off), 1)];

% Create table
T_presses = table(presses, group_treatment, group_laser, ...
    'VariableNames', {'Presses', 'Treatment', 'Laser'});

% Run 2-way ANOVA
[p_press, tbl_press, stats_press] = anovan(T_presses.Presses, ...
    {T_presses.Treatment, T_presses.Laser}, ...
    'model', 'interaction', ...
    'varnames', {'Treatment', 'Laser'});

% Posthoc if interaction is significant
if p_press(3) < 0.05
    fprintf('\nSignificant interaction detected (p = %.4f), running Tukey posthoc:\n', p_press(3));
    figure;
    c_press = multcompare(stats_press, 'Dimension', [1 2]);
end


[c_press, m_press, h_press, gnames_press] = multcompare(stats_press, 'Dimension', [1 2]);

fprintf('\nPosthoc comparisons for Total Presses:\n');
for i = 1:size(c_press,1)
    group1 = gnames_press{c_press(i,1)};
    group2 = gnames_press{c_press(i,2)};
    pval   = c_press(i,6);
    fprintf('  %s vs %s: p = %.4f\n', group1, group2, pval);
end

%%
% Example 1: Basic usage for final_ratio_attained

colors_struct.Control = Control_color;
colors_struct.A2A = A2A_color;
colors_struct.D1 = D1_color;

symbols_struct.Control = Control_symbol;  % e.g., 'o'
symbols_struct.A2A = A2A_symbol;          % e.g., 's'
symbols_struct.D1 = D1_symbol;            % e.g., '^'

[diff_data, fig] = plot_within_subject_laser_difference(risk_table, 'final_ratio_attained', ...
                                                        'Colors', colors_struct, ...
                                                        'Symbols', symbols_struct, ...
                                                        'YLabel', 'Δ Final PR Ratio', ...
                                                        'YLim', [-10 10], ...
                                                        'FigureWidth', 200, ...
                                                        'FigureHeight', 500);

% Prepare for ANOVA
control_diffs = diff_data.Control.differences;
a2a_diffs = diff_data.A2A.differences;
d1_diffs = diff_data.D1.differences;

% Run one-way ANOVA on the difference scores
all_diffs = [control_diffs; a2a_diffs; d1_diffs];
groups = [ones(size(control_diffs)); 2*ones(size(a2a_diffs)); 3*ones(size(d1_diffs))];
[p, tbl, stats] = anova1(all_diffs, groups, 'off');

% Post-hoc tests if significant
if p < 0.05
    [c, m, h, gnames] = multcompare(stats, 'Display', 'off');
    disp(c);
end
%%


[diff_data, fig] = plot_within_subject_laser_difference(risk_table, 'total_presses', ...
                                                        'Colors', colors_struct, ...
                                                        'Symbols', symbols_struct, ...
                                                        'YLabel', 'Δ total_presses', ...
                                                        'YLim', [-250 250], ...
                                                        'FigureWidth', 300, ...
                                                        'FigureHeight', 300);

% Prepare for ANOVA
control_diffs = diff_data.Control.differences;
a2a_diffs = diff_data.A2A.differences;
d1_diffs = diff_data.D1.differences;

% Run one-way ANOVA on the difference scores
all_diffs = [control_diffs; a2a_diffs; d1_diffs];
groups = [ones(size(control_diffs)); 2*ones(size(a2a_diffs)); 3*ones(size(d1_diffs))];
[p, tbl, stats] = anova1(all_diffs, groups, 'off');

% Post-hoc tests if significant
if p < 0.05
    [c, m, h, gnames] = multcompare(stats, 'Display', 'off');
    disp(c);
end



%% 06/27/2025 data pulled from RRD Photometry and Opto Progressive Ratio v2 - need to get raw data from PCs in FLAC

PR_presses_hM4Di_laser_ON = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'hM4Di') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_presses_hM4Di_laser_OFF = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'hM4Di') & strcmp(risk_table.LaserTreatment, 'Laser Off'));

PR_presses_mCherry_laser_ON = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'mCherry') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_presses_mCherry_laser_OFF = risk_table.total_presses(strcmp(risk_table.TreatmentCondition, 'mCherry') & strcmp(risk_table.LaserTreatment, 'Laser Off'));


PR_ratios_hM4Di_laser_ON = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'hM4Di') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_ratios_hM4Di_laser_OFF = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'hM4Di') & strcmp(risk_table.LaserTreatment, 'Laser Off'));

PR_ratios_mCherry_laser_ON = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'mCherry') & strcmp(risk_table.LaserTreatment, 'Laser On'));
PR_ratios_mCherry_laser_OFF = risk_table.final_ratio_attained(strcmp(risk_table.TreatmentCondition, 'mCherry') & strcmp(risk_table.LaserTreatment, 'Laser Off'));





% Compute means and SEMs for each group
% PR Presses
mean_presses = [mean(PR_presses_hM4Di_laser_ON),...
                mean(PR_presses_hM4Di_laser_OFF),...
                mean(PR_presses_mCherry_laser_ON),...
                mean(PR_presses_mCherry_laser_OFF),...

                ];
sem_presses = [std(PR_presses_hM4Di_laser_ON)/sqrt(length(PR_presses_hM4Di_laser_ON)), ...
               std(PR_presses_hM4Di_laser_OFF)/sqrt(length(PR_presses_hM4Di_laser_OFF)), ...
               std(PR_presses_mCherry_laser_ON)/sqrt(length(PR_presses_mCherry_laser_ON)), ...
               std(PR_presses_mCherry_laser_OFF)/sqrt(length(PR_presses_mCherry_laser_OFF)), ...

               ];

% PR Ratios
mean_ratios = [mean(PR_ratios_hM4Di_laser_ON),...
                mean(PR_ratios_hM4Di_laser_OFF),...
                mean(PR_ratios_mCherry_laser_ON),...
                mean(PR_ratios_mCherry_laser_OFF),...
                ];
sem_ratios = [std(PR_ratios_hM4Di_laser_ON)/sqrt(length(PR_ratios_hM4Di_laser_ON)), ...
               std(PR_ratios_hM4Di_laser_OFF)/sqrt(length(PR_ratios_hM4Di_laser_OFF)), ...
               std(PR_ratios_mCherry_laser_ON)/sqrt(length(PR_ratios_mCherry_laser_ON)), ...
               std(PR_ratios_mCherry_laser_OFF)/sqrt(length(PR_ratios_mCherry_laser_OFF)), ...
               ];

% Define colors for each group
colors = [0.6, 0.8, 1.0;    % Light blue for Control Laser On
          0.2, 0.4, 0.8;    % Dark blue for Control Laser Off
          1.0, 0.7, 0.7;    % Light red for ChrimsonR Laser On
          0.8, 0.2, 0.2;    % Dark red for ChrimsonR Laser Off
          ];   

% Create figure
figure;

% First subplot - Final PR Ratio Attained
subplot(1,2,1);
b1 = bar(1:4, mean_ratios, 0.6, 'FaceColor', 'flat');
% Set colors for each bar
for i = 1:4
    b1.CData(i,:) = colors(i,:);
end
hold on;

% Add error bars
errorbar(1:4, mean_ratios, sem_ratios, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Add individual data points with jitter
all_ratio_data = {PR_ratios_hM4Di_laser_ON,...
                    PR_ratios_hM4Di_laser_OFF,...
                    PR_ratios_mCherry_laser_ON,...
                    PR_ratios_mCherry_laser_OFF,...              
                   };
for i = 1:4
    x_pos = i + 0.1 * (rand(size(all_ratio_data{i})) - 0.5); % Add jitter
    scatter(x_pos, all_ratio_data{i}, 40, 'k', 'filled');
end

hold off;
xlim([0.5 6.5]);
ylim([0 20]);
xticks(1:2);
xticklabels({'Control\nLaser On', 'Control\nLaser Off', 'ChrimsonR\nLaser On', 'ChrimsonR\nLaser Off'});
yticks([0:5:20]);
ylabel('Final PR ratio attained');

% Second subplot - Total Presses
subplot(1,2,2);
b2 = bar(1:4, mean_presses, 0.6, 'FaceColor', 'flat');
% Set colors for each bar
for i = 1:4
    b2.CData(i,:) = colors(i,:);
end
hold on;

% Add error bars
errorbar(1:4, mean_presses, sem_presses, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);


% Add individual data points with jitter
all_presses_data = {PR_presses_hM4Di_laser_ON,...
                    PR_presses_hM4Di_laser_OFF,...
                    PR_presses_mCherry_laser_ON,...
                    PR_presses_mCherry_laser_OFF,...            
                   };
for i = 1:4
    x_pos = i + 0.1 * (rand(size(all_presses_data{i})) - 0.5); % Add jitter
    scatter(x_pos, all_presses_data{i}, 40, 'k', 'filled');
end

hold off;
xlim([0.5 6.5]);
ylim([0 600]);
xticks(1:2);
xticklabels({'Control\nLaser On', 'Control\nLaser Off', 'ChrimsonR\nLaser On', 'ChrimsonR\nLaser Off'});
yticks([0 200 400 600]);
ylabel('Total presses');

% Adjust layout
set(gcf, 'Position', [100, 100, 500, 350]); % Increased width for 4 bars