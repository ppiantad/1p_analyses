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
xticks([]);
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
xticks([]);
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