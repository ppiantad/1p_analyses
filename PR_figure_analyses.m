% Compute mean and SEM
mean_PR_presses = mean(total_PR_presses);
sem_PR_presses = std(total_PR_presses) / sqrt(length(total_PR_presses));

mean_PR_ratio = mean(PR_final_ratio_reached);
sem_PR_ratio = std(PR_final_ratio_reached) / sqrt(length(PR_final_ratio_reached));

% Create figure
figure;

% First subplot - Final PR Ratio Attained
subplot(1,2,1);
bar(1, mean_PR_ratio, 0.5, 'FaceColor', [0.8, 0.4, 0.4]); % Narrow bar
hold on;
errorbar(1, mean_PR_ratio, sem_PR_ratio, 'k', 'LineWidth', 1.5);

% Scatter plot of individual data points
scatter(ones(size(PR_final_ratio_reached)), PR_final_ratio_reached, 40, 'k', 'filled', 'jitter', 'on', 'jitterAmount', 0.1);

hold off;
xlim([0.5 1.5]);
ylim([0 16])
xticks([]);
yticks([0 5 10 15])
ylabel('Final PR ratio attained');

% Second subplot - Total Presses
subplot(1,2,2);
bar(1, mean_PR_presses, 0.5, 'FaceColor', [0.3, 0.6, 0.8]); % Narrow bar
hold on;
errorbar(1, mean_PR_presses, sem_PR_presses, 'k', 'LineWidth', 1.5);

% Scatter plot of individual data points
scatter(ones(size(total_PR_presses)), total_PR_presses, 40, 'k', 'filled', 'jitter', 'on', 'jitterAmount', 0.1);

hold off;
xlim([0.5 1.5]);
ylim([0 320])
xticks([]);
yticks([0 100 200 300])
ylabel('Total presses');

% Adjust layout
set(gcf, 'Position', [100, 100, 200, 350]); % Adjust figure size


%%

% all_aborts = risk_table.large_abort + risk_table.small_aborts;

latency_change = risk_table.block_1_large_choice_latency - (risk_table.block_2_large_choice_latency + risk_table.block_2_small_choice_latency);

% x = remapped_prechoice_ratio(risk_table.risky == 1)';
% y = latency_change(risk_table.risky == 1);
x = total_PR_presses';
y = risk_table.Mean_1_to_3;
% x = remapped_prechoice_ratio(num_cells_mouse > 30)';
% y = risk_table.Mean_1_to_3(num_cells_mouse > 30);
% Create a new figure with specific dimensions
figure;
width = 250; % Width of the figure
height = 350; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]

hold on;
scatter(x, y, 100, 'filled')
% Set the axis labels to have 2 decimal places
xtickformat('%.2f');
ytickformat('%.2f');
% Add a regression line (You can keep this part unchanged)
coefficients = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(coefficients, x_fit);
plot(x_fit, y_fit, 'r');

% Calculate R-squared value (You can keep this part unchanged)
y_pred = polyval(coefficients, x);
ssr = sum((y_pred - mean(y)).^2);
sst = sum((y - mean(y)).^2);
r_squared = ssr / sst;

% can also calculate the r-squared this way
% Calculate the R^2 value
[r, pval] = corrcoef(x, y); % Compute correlation coefficient matrix
rsq = r(1, 2)^2; % Extract and square the correlation coefficient
% Get axes limits
ax = gca;
xLimits = xlim(ax);
yLimits = ylim(ax);

% Set text position relative to axes limits
xPos = xLimits(2) - 0.05 * range(xLimits); % Slightly inside the top-right
yPos = yLimits(2) - 0.05 * range(yLimits); % Slightly inside the top-right

% Add R-squared value to the plot (You can keep this part unchanged)
text(xPos, yPos, ...
    {['R^2 = ' num2str(r_squared, '%.2f')], ...
     ['p = ' num2str(pval(2), '%.2f')]}, ...
    'FontSize', 12, ...
    'Color', 'blue', ...
    'HorizontalAlignment', 'right', ... % Align text to the right
    'VerticalAlignment', 'top');       % Align text to the top

hold off;

%%
consum_Late_RM = respClass_all_array{1, 1} == 1;
consum_PR_D1 = respClass_all_array{1, 2} == 1;
consum_both = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 1;

consum_both_sum = sum(consum_Late_RM & consum_PR_D1);

consum_both_sum/sum(consum_Late_RM)
consum_both_sum/sum(consum_PR_D1)

consum_both_sum/(sum(consum_Late_RM) + sum(consum_PR_D1))

sum(consum_both)/sum(consum_Late_RM(consum_both ~=1)) + sum(consum_PR_D1(consum_both ~=1))

total_modulated = [(sum(consum_Late_RM)/neuron_num)*100 (sum(consum_PR_D1)/neuron_num)*100];
A = total_modulated;
I = (sum(consum_both)/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end


%%
figure;
width = 200; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
h(1) = shadedErrorBar(ts1, mean(neuron_mean_array{1,1}(consum_both == 1, :)), mean(neuron_sem_array{1,1}(consum_both == 1, :)), 'lineProps', {'color', 'r'});
h(2) = shadedErrorBar(ts1, mean(neuron_mean_array{1,2}(consum_both == 1, :)), mean(neuron_sem_array{1,2}(consum_both == 1, :)), 'lineProps', {'color', 'b'});
% h(3) = shadedErrorBar(ts1, mean(neuron_mean_array{1,4}(conserved_all == 1, :)), mean(neuron_sem_array{1,4}(conserved_all == 1, :)), 'lineProps', {'color', 'g'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(lost_all == 1, :)), nanmean(neuron_sem_array{1,4}(lost_all == 1, :)), 'lineProps', {'color', 'r'});
% h(3) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(remapped_all == 1, :)), nanmean(neuron_sem_array{1,4}(remapped_all == 1, :)), 'lineProps', {'color', 'r'});
% h(2) = shadedErrorBar(ts1, nanmean(neuron_mean_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), nanmean(neuron_sem_array{1,arrays_to_examine(2)}(eval(vars_to_use{1, 2}) == event_for_figures, :)), 'lineProps', {'color', 'b'});
% legend([h(1).mainLine h(2).mainLine], '1st block', '2nd and 3rd block')
xlim([-2 4]);
ylim([-0.5 0.6])
set(gca, 'YTick', [-.50 -.25 0 .25 0.5]);
ytickformat('%.2f');
% hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 3}(respClass_all_array{1,1} == 1,:)), nanmean(neuron_sem_array{1, 3}(respClass_all_array{1,1} == 1,:)), 'lineProps', {'color', batlowW(iter,:)});
xline(0);
% xline(median_start_time_all, 'g', {'Median', 'start', 'time'})
% xline(median_collect_times_all, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice (s)');

%% 
% assuming you have done: 
% for RDT_D1: choiceTime.postchoice_0to2.SHK_1
% for PR_D1: collectionTime.reward_collection_1to3.COLLECT_1

shk_event = respClass_all_array{1,1} == 1;

consumption_event = respClass_all_array{1, 2} == 1;

total_modulated = [(sum(shk_event)/neuron_num)*100 (sum(consumption_event)/neuron_num)*100];

shk_and_consum_both_excited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1;
% this is the start of checking if neurons are MORE active than during
% other events, i.e. if you wanted to check if REW activated neurons are
% significantly differentially activated by SHK. preliminary poking around
% seems to suggest that few large reward active neurons have their activity
% increase in response to SHK
co_activated_indices = find(shk_and_consum_both_excited(1,:) == 1);
co_activated_indices_sum = numel(co_activated_indices);


A = total_modulated;
I = (co_activated_indices_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end


%%
shk_consistent = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1;
sum(shk_consistent)
shk_first_only = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} ~= 1;
sum(shk_first_only)
shk_second_only = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} == 1;
sum(shk_second_only)

% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav
% this outputs a ever so slightly wonky diagram. a few nodes that do not
% actually overlap minimally overlap (but intersections are 0), and 1 node
% that has 1 overlap does not overlap at all. 
% shk_consistent_ind = find(respClass_all_array{1,1} == 1);
% pre_choice_active_ind = find(respClass_all_array{1,1} == 1);
shk_first_only_ind = find(respClass_all_array{1,1} == 1);
shk_second_only_ind = find(respClass_all_array{1,2} == 1);
% consum_inhibited_ind = find(all_consum_inhibited == 1);
% setListData = {shk_first_only_ind, shk_second_only_ind};
setListData = {shk_first_only_ind, shk_second_only_ind};
setLabels = ["Shock Test", "RDT D1"];
figure;
ve_diagram = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

ve_diagram.ShowIntersectionCounts = true;
ve_diagram.ShowIntersectionAreas = true;
% h.SetLabels = [];

%%

shk_and_consum_both_excited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1;
% this is the start of checking if neurons are MORE active than during
% other events, i.e. if you wanted to check if REW activated neurons are
% significantly differentially activated by SHK. preliminary poking around
% seems to suggest that few large reward active neurons have their activity
% increase in response to SHK
co_activated_indices = find(shk_and_consum_both_excited(1,:) == 1);
co_activated_indices_sum = numel(co_activated_indices);

% Calculate different groups as percentages
only_shk = (sum(respClass_all_array{1,1} == 1)/neuron_num)*100 - (co_activated_indices_sum/neuron_num)*100; % SHK only
only_consumption = (sum(respClass_all_array{1,2} == 1)/neuron_num)*100 - (co_activated_indices_sum/neuron_num)*100; % Consumption only
both = (co_activated_indices_sum/neuron_num)*100; % Overlap between SHK and Consumption
not_modulated = 100 - (only_shk + only_consumption + both); % Unmodulated neurons

% Data for the stacked bar (bottom to top: Not modulated, Consumption only, Both, SHK only)
data_for_bar_plot = [not_modulated, only_consumption, both, only_shk];

% Create the stacked bar plot
figure;
bar(1, data_for_bar_plot, 'stacked'); % Single bar at x = 1
colormap([0.7 0.7 0.7; 0 0 1; 0.8 0.8 0; 1 0 0]); % Grey (Unmodulated), Blue (Consumption), Yellow (Both), Red (SHK)

% Formatting
xticks(1); % Single bar on x-axis
xticklabels({'Neuron Modulation'});
ylabel('Percentage of Neurons (%)');
title('Neuron Modulation by Events');
ylim([0 100]); % Ensure the bar always reaches 100%
legend({'Not Modulated', 'Consumption Only', 'Both', 'SHK Only'}, 'Location', 'eastoutside');

% Adding percentage labels
y_offset = 0; % Start at the base of the bar
for i = 1:length(data_for_bar_plot)
    if data_for_bar_plot(i) > 0  % Only add labels for non-zero sections
        text(1, y_offset + data_for_bar_plot(i)/2, sprintf('%.1f%%', data_for_bar_plot(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', 'white', 'FontWeight', 'bold');
        y_offset = y_offset + data_for_bar_plot(i); % Move up for next section
    end
end



%% 06/27/2025 data pulled from RRD Photometry and Opto Progressive Ratio v2 - need to get raw data from PCs in FLAC

mean_PR_presses_control_laser_on = [419 417 275 461 186 346 548]
mean_PR_presses_control_laser_off = [302 503 213 229 519]

mean_PR_presses_ChrimsonR_laser_on = [120 231 269 184 217 378]
mean_PR_presses_ChrimsonR_laser_off = [511 462 363 379 503 294]

mean_PR_ratio_control_laser_on = [53 53 41 57 33 49 61]
mean_PR_ratio_control_laser_off = [45 61 37 37 61]

mean_PR_ratio_ChrimsonR_laser_on = [25 37 41 33 37 49]
mean_PR_ratio_ChrimsonR_laser_off = [61 57 49 53 61 45]




% Compute means and SEMs for each group
% PR Presses
mean_presses = [mean(mean_PR_presses_control_laser_on), mean(mean_PR_presses_control_laser_off), ...
                mean(mean_PR_presses_ChrimsonR_laser_on), mean(mean_PR_presses_ChrimsonR_laser_off)];
sem_presses = [std(mean_PR_presses_control_laser_on)/sqrt(length(mean_PR_presses_control_laser_on)), ...
               std(mean_PR_presses_control_laser_off)/sqrt(length(mean_PR_presses_control_laser_off)), ...
               std(mean_PR_presses_ChrimsonR_laser_on)/sqrt(length(mean_PR_presses_ChrimsonR_laser_on)), ...
               std(mean_PR_presses_ChrimsonR_laser_off)/sqrt(length(mean_PR_presses_ChrimsonR_laser_off))];

% PR Ratios
mean_ratios = [mean(mean_PR_ratio_control_laser_on), mean(mean_PR_ratio_control_laser_off), ...
               mean(mean_PR_ratio_ChrimsonR_laser_on), mean(mean_PR_ratio_ChrimsonR_laser_off)];
sem_ratios = [std(mean_PR_ratio_control_laser_on)/sqrt(length(mean_PR_ratio_control_laser_on)), ...
              std(mean_PR_ratio_control_laser_off)/sqrt(length(mean_PR_ratio_control_laser_off)), ...
              std(mean_PR_ratio_ChrimsonR_laser_on)/sqrt(length(mean_PR_ratio_ChrimsonR_laser_on)), ...
              std(mean_PR_ratio_ChrimsonR_laser_off)/sqrt(length(mean_PR_ratio_ChrimsonR_laser_off))];

% Define colors for each group
colors = [0.6, 0.8, 1.0;    % Light blue for Control Laser On
          0.2, 0.4, 0.8;    % Dark blue for Control Laser Off
          1.0, 0.7, 0.7;    % Light red for ChrimsonR Laser On
          0.8, 0.2, 0.2];   % Dark red for ChrimsonR Laser Off

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
all_ratio_data = {mean_PR_ratio_control_laser_on, mean_PR_ratio_control_laser_off, ...
                  mean_PR_ratio_ChrimsonR_laser_on, mean_PR_ratio_ChrimsonR_laser_off};
for i = 1:4
    x_pos = i + 0.1 * (rand(size(all_ratio_data{i})) - 0.5); % Add jitter
    scatter(x_pos, all_ratio_data{i}, 40, 'k', 'filled');
end

hold off;
xlim([0.5 4.5]);
ylim([0 70]);
xticks(1:4);
xticklabels({'Control\nLaser On', 'Control\nLaser Off', 'ChrimsonR\nLaser On', 'ChrimsonR\nLaser Off'});
yticks([0 20 40 60]);
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
all_presses_data = {mean_PR_presses_control_laser_on, mean_PR_presses_control_laser_off, ...
                    mean_PR_presses_ChrimsonR_laser_on, mean_PR_presses_ChrimsonR_laser_off};
for i = 1:4
    x_pos = i + 0.1 * (rand(size(all_presses_data{i})) - 0.5); % Add jitter
    scatter(x_pos, all_presses_data{i}, 40, 'k', 'filled');
end

hold off;
xlim([0.5 4.5]);
ylim([0 600]);
xticks(1:4);
xticklabels({'Control\nLaser On', 'Control\nLaser Off', 'ChrimsonR\nLaser On', 'ChrimsonR\nLaser Off'});
yticks([0 200 400 600]);
ylabel('Total presses');

% Adjust layout
set(gcf, 'Position', [100, 100, 500, 350]); % Increased width for 4 bars



%% ANOVA Analysis

% Combine all PR ratio data
ratios = [mean_PR_ratio_control_laser_on, ...
          mean_PR_ratio_control_laser_off, ...
          mean_PR_ratio_ChrimsonR_laser_on, ...
          mean_PR_ratio_ChrimsonR_laser_off]';

% Define group labels
group_treatment = [repmat({'Control'}, length(mean_PR_ratio_control_laser_on), 1); 
                   repmat({'Control'}, length(mean_PR_ratio_control_laser_off), 1); 
                   repmat({'ChrimsonR'}, length(mean_PR_ratio_ChrimsonR_laser_on), 1); 
                   repmat({'ChrimsonR'}, length(mean_PR_ratio_ChrimsonR_laser_off), 1)];

group_laser = [repmat({'On'}, length(mean_PR_ratio_control_laser_on), 1); 
               repmat({'Off'}, length(mean_PR_ratio_control_laser_off), 1); 
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
presses = [mean_PR_presses_control_laser_on, ...
           mean_PR_presses_control_laser_off, ...
           mean_PR_presses_ChrimsonR_laser_on, ...
           mean_PR_presses_ChrimsonR_laser_off]';

% Define same group labels
group_treatment = [repmat({'Control'}, length(mean_PR_presses_control_laser_on), 1); 
                   repmat({'Control'}, length(mean_PR_presses_control_laser_off), 1); 
                   repmat({'ChrimsonR'}, length(mean_PR_presses_ChrimsonR_laser_on), 1); 
                   repmat({'ChrimsonR'}, length(mean_PR_presses_ChrimsonR_laser_off), 1)];

group_laser = [repmat({'On'}, length(mean_PR_presses_control_laser_on), 1); 
               repmat({'Off'}, length(mean_PR_presses_control_laser_off), 1); 
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
