%% run eventRelatedActivity and then run data_loop to get AA for Block 1 and Block 2. should continue to update this & do a real probe of AAs
% actually just load the x10 dataset, then run data_loop with AA, 1, Block,
% 2 and AA, 1, Block, 3





%% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav
% this outputs a ever so slightly wonky diagram. a few nodes that do not
% actually overlap minimally overlap (but intersections are 0), and 1 node
% that has 1 overlap does not overlap at all. 
shk_ind = find(respClass_all_array{1,4} == 1);
% pre_choice_active_ind = find(respClass_all_array{1,1} == 1);
consum_active_ind = find(respClass_all_array{1,3} == 1);
consum_active_block_2_3 = find(respClass_all_array{1,10} == 1);
post_choice_active_ind = find(respClass_all_array{1,2} == 1);
aa_active_ind = find(respClass_all_array{1,11} == 1);
% consum_inhibited_ind = find(all_consum_inhibited == 1);
setListData = {shk_ind, consum_active_ind, aa_active_ind};
setLabels = ["Shk excited", "Consumption excited", "Approach-Abort excited"];
figure;
ve_diagram = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

ve_diagram.ShowIntersectionCounts = true;
ve_diagram.ShowIntersectionAreas = true;
% h.SetLabels = [];

shk_alone = respClass_all_array{1,4} == 1  & prechoice_block_1 ~=1 & postchoice_reward_block_1 ~= 1 & collect_block_1 ~= 1 & respClass_all_array{1,11} ~= 1;
%% chi-square test of proportions to see if more the proprtion of AA that are also shock cells is greater than consumption that are also shock
% based on https://www.mathworks.com/matlabcentral/answers/96572-how-can-i-perform-a-chi-square-test-to-determine-how-statistically-different-two-proportions-are-in
aa_and_shk = respClass_all_array{1,11} == 1 & respClass_all_array{1,4} == 1 & respClass_all_array{1,3} ~= 1;
aa_and_shk_sum = sum(aa_and_shk)
aa_not_shk = respClass_all_array{1,11} == 1 & respClass_all_array{1,4} ~= 1 & respClass_all_array{1,3} ~= 1;
aa_not_shk_sum = sum(aa_not_shk)

% using block 1 consumption neurons
consumption_and_shk = respClass_all_array{1,3} == 1 & respClass_all_array{1,4} == 1 & respClass_all_array{1,11} ~= 1;
consumption_and_shk_sum = sum(consumption_and_shk)
consumption_not_shk = respClass_all_array{1,3} == 1 & respClass_all_array{1,4} ~= 1 & respClass_all_array{1,11} ~= 1;
consumption_not_shk_sum = sum(consumption_not_shk)


postchoice_and_shk = respClass_all_array{1,2} == 1 & respClass_all_array{1,4} == 1 & respClass_all_array{1,11} ~= 1;
postchoice_and_shk_sum = sum(postchoice_and_shk)
postchoice_not_shk = respClass_all_array{1,2} == 1 & respClass_all_array{1,4} ~= 1 & respClass_all_array{1,11} ~= 1;
postchoice_not_shk_sum = sum(postchoice_not_shk)



n1 = aa_and_shk_sum; N1 = aa_and_shk_sum+aa_not_shk_sum;

n2 = consumption_and_shk_sum; N2 = consumption_and_shk_sum+consumption_not_shk_sum;

% n2 = postchoice_and_shk_sum; N2 = postchoice_and_shk_sum+postchoice_not_shk_sum;

% Pooled estimate of proportion

p0 = (n1+n2) / (N1+N2)

% Expected counts under H0 (null hypothesis)

n10 = N1 * p0;

n20 = N2 * p0;

% Chi-square test, by hand

observed = [n1 N1-n1 n2 N2-n2];

expected = [n10 N1-n10 n20 N2-n20];

chi2stat = sum((observed-expected).^2 ./ expected)

p = 1 - chi2cdf(chi2stat,1)


%%



abort_mean = mean(neuron_mean_array{1, 11}(:, ts1 > -1 & ts1 < 1),  2);

% [peak_values, time_of_peak_activity] = max(neuron_mean_array{1, 1}, [], 2);
[~, sort_indices] = sort(abort_mean);
neuron_mean_sorted = neuron_mean_array{1, 11}(sort_indices, :);


% Sort the rows of activated_neuron_mean based on peak_times.
% [~, sort_indices] = sort(time_of_peak_activity);
% activated_neuron_mean_sorted = activated_rows(sort_indices, :);

% Now, activated_neuron_mean_sorted contains the rows of neuron_mean filtered by respClass_all == 1
% and sorted by the time of peak activity.

figure;
% Generate the heatmap
imagesc(ts1, 1, neuron_mean_sorted);

% Add a colorbar and axis labels
colorbar;
xlabel('Time (s)');
ylabel('Neuron');

%%
custom_colormap = [
    1, 1, 1; % white
    1, 0.9, 0.9;
    1, 0.8, 0.8;
    1, 0.7, 0.7;
    1, 0.6, 0.6;
    1, 0.5, 0.5;
    1, 0.4, 0.4;
    1, 0.3, 0.3;
    1, 0.2, 0.2;
    1, 0.1, 0.1;
    1, 0, 0;   % red
];


% Generate more intermediate colors for a smoother transition
n = 256; % Number of colors
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, n));

% Create a figure with a narrow width and taller height
figure('Position', [100, 100, 350, 600]); % [left, bottom, width, height]
hold on
% Create a tiled layout with 2 rows and 1 column
% tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile (heatmap)
% ax1 = nexttile;
% hold on;

% Plot the heatmap
imagesc(ts1, 1, neuron_mean_sorted);

% Apply the custom colormap
colormap(custom_colormap);

% Restrict the color axis range to [-1, 1]
clim([-.5 .5]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 
ylim([1, neuron_num]);
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8]);
set(gca, 'YTick', [1, neuron_num]);
xline(0)
% scatter(time2Collect, Tris               , 'Marker', 'p')
% scatter(trialStartTime, Tris, 'Marker', 's')
fontsize(18, 'points')
hold off;



%%

figure;
hold on
% Create a histogram for allCorrelations

width = 300; % Width of the figure
height = 600; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 25, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
ylim([-0.5 0.5]);
% Set X-axis ticks
set(gca, 'XTick', [-8, 0, 8], 'YTick', [-0.5 0 0.5]);
shadedErrorBar(ts1, nanmean(zall_mean_all_array{1, 12}(respClass_all_array{1, 11}==1, :)), nanmean(sem_all_array{1, 12}(respClass_all_array{1, 11}==1, :)), 'lineProps', {'color', 'r'});
hold on;shadedErrorBar(ts1, nanmean(zall_mean_all_array{1, 13}(respClass_all_array{1, 11}==1, :)), nanmean(sem_all_array{1, 13}(respClass_all_array{1, 11}==1, :)), 'lineProps', {'color', 'k'});

xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption'}, 'Location','northwest')

hold off


%%
% Get AUCs for the relevant periods for the 3 defined events
% Define time windows
pre_choice_window = [-4 0];     % Pre-choice period: -4 to 0 s
post_choice_window = [0 2];     % Post-choice period: 0 to 2 s
consumption_window = [1 3];     % Consumption period: 1 to 3 s if using data aligned to collect, do 0 to 2 to keep things consistent

% Initialize arrays to store AUCs
% auc_pre_choice = zeros(size(neuron_mean_array));
% auc_post_choice = zeros(size(neuron_mean_array));
% auc_consumption = zeros(size(neuron_mean_array));
pre_choice_neuron_count = 0;
% Iterate over each element of neuron_mean_array
for i = 1:size(neuron_mean_array{1,1}, 1)
    % Select data where exclusive_activated_session_1 is 1
    if prechoice_block_1(i) == 1
        pre_choice_neuron_count = pre_choice_neuron_count+1;
        selected_data = neuron_mean_array{1, 1}(i, :);
        % % Extract time variable (assuming it's named 'ts1')
        % ts1_data = ts1{i}(exclusive_activated_session_1{i} == 1);

        % Find indices corresponding to each time window
        pre_choice_indices = ts1 >= pre_choice_window(1) & ts1 <= pre_choice_window(2);
        post_choice_indices = ts1 >= post_choice_window(1) & ts1 <= post_choice_window(2);
        consumption_indices = ts1 >= consumption_window(1) & ts1 <= consumption_window(2);

        % Compute AUC for each time window
        % AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
        action_auc_pre_choice(pre_choice_neuron_count) = trapz(selected_data(pre_choice_indices));
        action_auc_post_choice(pre_choice_neuron_count) = trapz(selected_data(post_choice_indices));
        action_auc_consumption(pre_choice_neuron_count) = trapz(selected_data(consumption_indices));
    else

    end
end

%%

figure; imagesc(ts1, [], zall_mouse{4, 11}{1, 75})
figure; imagesc(ts1, [], zall_mouse{4, 4}{1, 75})

figure; imagesc(ts1, [], zall_mouse{8, 11}{1, 7})
figure; imagesc(ts1, [], zall_mouse{8, 4}{1, 7})

figure; imagesc(ts1, [], zall_mouse{8, 11}{1, 41})
figure; imagesc(ts1, [], zall_mouse{8, 4}{1, 41})


aa_large = zall_mean_all_array{1, 12}(respClass_all_array{1, 11}==1, :);
aa_small = zall_mean_all_array{1, 13}(respClass_all_array{1, 11}==1, :);
aa_large = aa_large(~any(isnan(aa_large), 2), :);
aa_small = aa_small(~any(isnan(aa_small), 2), :);


%%
% run data_loop_SLEAP
% [BehavData,trials, varargin_identity_class]=TrialFilter_test(BehavData, 'REW', 1.2, 'BLOCK', 2, 'BLOCK', 3, 'SHK', 0);
% [BehavData,trials, varargin_identity_class]=TrialFilter_test(BehavData, 'AA', 1);


mean_data_array = {zall_mean_all_array{1, 1}, zall_mean_all_array{1, 2}};
sem_data_array = {sem_all_array{1, 1}  , sem_all_array{1, 2}};

[comparison, perm_p_sig] = perm_and_bCI_fn_analysis_PhilDBressel_for_1p(mean_data_array, sem_data_array, ts1)