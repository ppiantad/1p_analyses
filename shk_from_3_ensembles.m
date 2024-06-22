% assuming input data are:
%pre-choice 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1 choiceTime -4 to 0
%post-choice 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1 choiceTime 0 to 2
%consumption 'OMITALL', 0, 'BLANK_TOUCH', 0, 'BLOCK', 1 collectionTime 1 to 3
%shk 'SHK', 1 choiceTime 0 to 2


%%
%These data can be used to plot the median or mean choice
% time on a PCA graph, for example

behav_tbl_iter_single = behav_tbl_iter(1);

% Initialize the concatenated table
concatenatedTable = table();

% Iterate through the 3x1 cell array
for i = 1:numel(behav_tbl_iter_single)
    % Assuming each cell contains a 12x1 cell array of tables
    twelveByOneCellArray = behav_tbl_iter_single{i};
    
    % Initialize a temporary table to store the concatenated tables for this cell
    tempTable = table();
    
    % Iterate through the 12x1 cell array
    for j = 1:numel(twelveByOneCellArray)
        % Assuming each cell in the 12x1 cell array contains a table
        currentTable = twelveByOneCellArray{j};
        
        % Concatenate the current table to the temporary table vertically
        tempTable = vertcat(tempTable, currentTable);
    end
    
    % Concatenate the temporary table to the overall concatenated table vertically
    concatenatedTable = vertcat(concatenatedTable, tempTable);
end

median_choice_time_block_1 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_choice_time_block_2 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_choice_time_block_3 = median(concatenatedTable.choiceTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));

median_collect_time_block_1 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 1) - concatenatedTable.stTime(concatenatedTable.Block == 1));
median_collect_time_block_2 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 2) - concatenatedTable.stTime(concatenatedTable.Block == 2));
median_collect_time_block_3 = median(concatenatedTable.collectionTime(concatenatedTable.Block == 3) - concatenatedTable.stTime(concatenatedTable.Block == 3));





median_start_time_from_choice = median(concatenatedTable.stTime - concatenatedTable.choiceTime);
median_collect_time_from_choice = median(concatenatedTable.collectionTime - concatenatedTable.choiceTime);



%%

excited_to_excited_all = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} == 1;
excited_to_excited_sum = sum(excited_to_excited_all);

excited_to_excited_1_to_2 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} ~= 1;
excited_to_excited_1_to_2_sum = sum(excited_to_excited_1_to_2);


excited_to_excited_1_to_3 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} == 1;
excited_to_excited_1_to_3_sum = sum(excited_to_excited_1_to_3);

excited_to_excited_2_to_3 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} == 1;
excited_to_excited_2_to_3_sum = sum(excited_to_excited_2_to_3);



exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} ~= 1 & respClass_all_array{1,4} ~= 1;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1)
exclusive_activated_session_2 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} ~= 1 & respClass_all_array{1,4} ~= 1;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2)
exclusive_activated_session_3 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} == 1 & respClass_all_array{1,4} ~= 1;
exclusive_activated_session_3_sum = sum(exclusive_activated_session_3)
exclusive_activated_session_4 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} ~= 1 & respClass_all_array{1,4} == 1;
exclusive_activated_session_4_sum = sum(exclusive_activated_session_4)


not_active = neuron_num - (exclusive_activated_session_1_sum + exclusive_activated_session_2_sum + exclusive_activated_session_3_sum + exclusive_activated_session_4_sum);
sum(not_active)

not_active_ind = ~(exclusive_activated_session_1 + exclusive_activated_session_2 + exclusive_activated_session_3+ exclusive_activated_session_4);

%%
% Example 2: Nested pie chart with custom colors for each wedge

% Initialize data points
inner_pie = [exclusive_activated_session_1_sum/neuron_num,...
            
            exclusive_activated_session_2_sum/neuron_num,...
            
            exclusive_activated_session_3_sum/neuron_num,...
           
            exclusive_activated_session_4_sum/neuron_num,...


            
            
            not_active/neuron_num];

figure; donutchart(inner_pie, 'InnerRadius', 0.5)
figure; pie(inner_pie)





figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]


shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(exclusive_activated_session_1==1, :)), nanmean(neuron_sem_array{1, 1}(exclusive_activated_session_1==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(exclusive_activated_session_2==1, :)), nanmean(neuron_sem_array{1, 1}(exclusive_activated_session_2==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(exclusive_activated_session_3==1, :)), nanmean(neuron_sem_array{1, 1}(exclusive_activated_session_3==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(exclusive_activated_session_4==1, :)), nanmean(neuron_sem_array{1, 4}(exclusive_activated_session_4==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(not_active_ind==1, :)), nanmean(neuron_sem_array{1, 1}(not_active_ind==1, :)), 'lineProps', {'color', batlowW(iter,:)});
% xtickformat('%.2f');
ytickformat('%.2f');
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'pre-choice active', 'post-choice reward active', 'consumption', 'not active'}, 'Location','northwest')


true_neutral = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 3; %omit shock == 3 from classification
sum(true_neutral)


figure;
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 1}  ==1, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array{1, 1}  ==1, :)), 'lineProps', {'color', batlowW(5,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 1}  ==2, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array{1, 1}  ==2, :)), 'lineProps', {'color', batlowW(50,:)});


figure;
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 2}  ==1, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array{1, 2}  ==1, :)), 'lineProps', {'color', batlowW(5,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 2}  ==2, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array{1, 2}  ==2, :)), 'lineProps', {'color', batlowW(50,:)});



figure;
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 3}  ==1, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array{1, 3}  ==1, :)), 'lineProps', {'color', batlowW(5,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(respClass_all_array{1, 3}  ==2, :)), nanmean(neuron_sem_array{1, 1}(respClass_all_array{1, 3}  ==2, :)), 'lineProps', {'color', batlowW(50,:)});



figure;
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(respClass_all_array{1, 4}  ==1, :)), nanmean(neuron_sem_array{1, 4}(respClass_all_array{1, 4}  ==1, :)), 'lineProps', {'color', batlowW(5,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 4}(respClass_all_array{1, 4}  ==2, :)), nanmean(neuron_sem_array{1, 4}(respClass_all_array{1, 4}  ==2, :)), 'lineProps', {'color', batlowW(50,:)});





% true_neutral = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} ~= 1 & respClass_all_array{1,3} ~= 1 & respClass_all_array{1,4} ~= 1;
% sum(true_neutral)

%%

exclusive_shk_activated = respClass_all_array{1,4} == 1 & respClass_all_array{1,1} == 3 &respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 3;
exclusive_collection_activated = respClass_all_array{1,4} == 3 & respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 1;
shk_event = respClass_all_array{1,4} == 1;
post_choice_rew_event = respClass_all_array{1,2} == 1;
post_choice_both_excited = respClass_all_array{1,2} == 1 & respClass_all_array{1,4} == 1;
% this is the start of checking if neurons are MORE active than during
% other events, i.e. if you wanted to check if REW activated neurons are
% significantly differentially activated by SHK. preliminary poking around
% seems to suggest that few large reward active neurons have their activity
% increase in response to SHK
co_activated_indices = find(post_choice_both_excited(1,:) == 1);
co_activated_indices_sum = numel(co_activated_indices);
clear h p ci stats
for qq = 1:size(co_activated_indices, 2)
    [h(qq),p(qq),ci{qq},stats{qq}] = ttest(neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx),neuron_mean_array{1, 4}(co_activated_indices(qq),sub_window_idx));
    mean_diff(qq) = mean(neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx) - mean(neuron_mean_array{1, 4}(co_activated_indices(qq),sub_window_idx)));
end

sig_increase_shk_from_large = co_activated_indices(p < 0.05 & mean_diff < 0);
sig_increase_shk_from_large_sum = numel(sig_increase_shk_from_large);
sig_increase_shk_from_large_ind = zeros(1, size(respClass_all_array{1,4}, 2));

sig_increase_shk_from_large_ind(:, sig_increase_shk_from_large) = 1;

no_sig_increase_shk_from_large = co_activated_indices(p > 0.05);
no_sig_increase_shk_from_large_sum = numel(no_sig_increase_shk_from_large);
no_sig_increase_shk_from_large_ind = zeros(1, size(respClass_all_array{1,4}, 2));
no_sig_increase_shk_from_large_ind(:, no_sig_increase_shk_from_large ) = 1;

shk_activated = respClass_all_array{1,4} == 1 & no_sig_increase_shk_from_large_ind ~= 1
shk_activated_sum = sum(shk_activated);

figure; plot(ts1, nanmean(neuron_mean_array{1, 4}(shk_activated,:))); hold on; plot(ts1,  nanmean(neuron_mean_array{1, 2}(respClass_all_array{1,2} == 1,:)));


total_modulated = [(sum(post_choice_rew_event)/neuron_num)*100 (sum(shk_event)/neuron_num)*100];

A = total_modulated;
I = (co_activated_indices_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end

total_modulated = [(sum(shk_event)/neuron_num)*100 (sum(post_choice_rew_event)/neuron_num)*100];
A = total_modulated;
I = (co_activated_indices_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end


for zz = 1:size(respClass_all_array_mouse, 1)
    exclusive_shk_activated_mouse{zz} = respClass_all_array_mouse{zz,4} == 1 & respClass_all_array_mouse{zz,1} == 3 & respClass_all_array_mouse{zz,2} == 3 & respClass_all_array_mouse{zz,3} == 3;
    sum_exclusive_shk_activated_mouse(zz) = sum(exclusive_shk_activated_mouse{zz});
    percent_exclusive_shk_activated_mouse(zz) = sum_exclusive_shk_activated_mouse(zz)/size(exclusive_shk_activated_mouse{zz}, 2);
    max_activity_exclusive_shk(zz) = max(mean(neuron_mean_mouse{zz, 4}(exclusive_shk_activated_mouse{1, zz}   == 1, :)));
end

%%

pre_active_choice_to_shk = respClass_all_array{1,1} == 1 & shk_activated == 1;
sum(pre_active_choice_to_shk)

post_choice_reward_active_to_shk = respClass_all_array{1,2} == 1 & shk_activated == 1;
sum(post_choice_reward_active_to_shk)

consumption_active_to_shk = respClass_all_array{1,3} == 1 & shk_activated == 1;
sum(consumption_active_to_shk)

neutral_to_shk = true_neutral == 1 & shk_activated == 1;
sum(neutral_to_shk)

pre_inhibited_choice_to_shk = respClass_all_array{1,1} == 2 & shk_activated == 1;
sum(pre_inhibited_choice_to_shk)

post_choice_reward_inhibited_to_shk = respClass_all_array{1,2} == 2 & shk_activated == 1;
sum(post_choice_reward_inhibited_to_shk)

consumption_inhibited_to_shk = respClass_all_array{1,3} == 2 & shk_activated == 1;
sum(consumption_active_to_shk)


pie_slices = [sum(pre_active_choice_to_shk) sum(post_choice_reward_active_to_shk) sum(consumption_active_to_shk) sum(neutral_to_shk)]
figure;
pie(pie_slices)


ylimits_for_figs = [-0.8 1.5]


figure; 

shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(pre_active_choice_to_shk , :)), nanmean(neuron_sem_array{1, 1}(pre_active_choice_to_shk , :)))
hold on; shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(pre_active_choice_to_shk , :)), nanmean(neuron_sem_array{1, 4}(pre_active_choice_to_shk , :)))
ylim(ylimits_for_figs)

figure; 
ylim(ylimits_for_figs)
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(post_choice_reward_active_to_shk , :)), nanmean(neuron_sem_array{1, 1}(post_choice_reward_active_to_shk , :)))
hold on; shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(post_choice_reward_active_to_shk , :)), nanmean(neuron_sem_array{1, 4}(post_choice_reward_active_to_shk , :)))
ylim(ylimits_for_figs)

figure; 
ylim(ylimits_for_figs)
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(consumption_active_to_shk , :)), nanmean(neuron_sem_array{1, 1}(consumption_active_to_shk , :)))
hold on; shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(consumption_active_to_shk , :)), nanmean(neuron_sem_array{1, 4}(consumption_active_to_shk , :)))
ylim(ylimits_for_figs)

figure; 
ylim(ylimits_for_figs)
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(neutral_to_shk , :)), nanmean(neuron_sem_array{1, 1}(neutral_to_shk , :)))
hold on; shadedErrorBar(ts1, nanmean(neuron_mean_array{1,4}(neutral_to_shk , :)), nanmean(neuron_sem_array{1, 4}(neutral_to_shk , :)))
ylim(ylimits_for_figs)

figure; plot(ts1, mean(neuron_mean_array{1, 1}(consumption_active_to_shk , :)))
hold on; plot(ts1, mean(neuron_mean_array{1,4}(consumption_active_to_shk , :)))

figure; 

plot(ts1, mean(neuron_mean_array{1, 1}(neutral_to_shk , :)))
hold on; plot(ts1, mean(neuron_mean_array{1,4}(neutral_to_shk , :)))


hold off;



%%
% Example 2: Nested pie chart with custom colors for each wedge

% Initialize data points
inner_pie = [exclusive_activated_session_1_sum/neuron_num,...
            
            exclusive_activated_session_2_sum/neuron_num,...
            
            exclusive_activated_session_3_sum/neuron_num,...
           
            excited_to_excited_sum/neuron_num,...

            excited_to_excited_1_to_2_sum/neuron_num,...

            excited_to_excited_1_to_3_sum/neuron_num,...

            excited_to_excited_2_to_3_sum/neuron_num,...
            
            
            not_active/neuron_num];

figure; donutchart(inner_pie, 'InnerRadius', 0.5)
figure; pie(inner_pie)

%%
test = [neuron_mean_array{1, 1}(exclusive_activated_session_1==1, :)];
test = [test; neuron_mean_array{1, 1}(exclusive_activated_session_2==1, :)];
test = [test; neuron_mean_array{1, 1}(exclusive_activated_session_3==1, :)];
test = [test; neuron_mean_array{1, 4}(exclusive_activated_session_4==1, :)];


action_index = [1:sum(respClass_all_array{1, 1}==1)];
consumption_index = [action_index(end)+1:action_index(end)+sum(respClass_all_array{1, 2}==1)];
neutral_index = [consumption_index(end)+1:consumption_index(end)+sum(respClass_all_array{1, 2}==3 & respClass_all_array{1, 1}==3)];


data = test;

alpha = 0.0001;

% Initialize matrices to store correlation coefficients and p-values
correlation_matrix = zeros(size(data, 1));
p_value_matrix = zeros(size(data, 1));

% Calculate correlation coefficients and p-values between rows
for i = 1:size(data, 1)
    for j = 1:size(data, 1)
        [corr_coeff, p_value] = corrcoef(data(i, :)', data(j, :)');
        correlation_matrix(i, j) = corr_coeff(1, 2); % Store correlation coefficient
        p_value_matrix(i, j) = p_value(1, 2); % Store p-value
    end
end

% Plot the correlation matrix
figure;
imagesc(correlation_matrix);
colorbar; % Add a colorbar to the plot
axis square; % Make the plot square for better visualization
title('Correlation Matrix');
xlabel('Neuron Number');
ylabel('Neuron Number');

% Show row and column indices on the plot
xticks(0:50:size(data, 1));
yticks(0:50:size(data, 1));

% If you want to customize the color map, you can use colormap function
% For example, using a blue-white-red colormap:
colormap(bluewhitered);

% If you want to limit the color scale to the range [0, 1]
caxis([-1 1]); % Assuming correlations range from -1 to 1

% % Display p-values as text on the plot
% for i = 1:size(data, 1)
%     for j = 1:size(data, 1)
%         text(j, i, sprintf('p = %.4f', p_value_matrix(i, j)), ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     end
% end
%% PRECHOICE
%run data_loop with choiceTime large rew [-10 to 5]
large_pre_choice_ensemble_block_1_zall = zall_mean_all_array{1, 1}(exclusive_activated_session_1==1, :);
large_pre_choice_ensemble_block_1_sem = sem_all_array{1, 1}(exclusive_activated_session_1==1, :);

large_pre_choice_ensemble_block_2_zall = zall_mean_all_array{1, 2}(exclusive_activated_session_1==1, :);
large_pre_choice_ensemble_block_2_sem = sem_all_array{1, 2}(exclusive_activated_session_1==1, :);

large_pre_choice_ensemble_block_3_zall = zall_mean_all_array{1, 3}(exclusive_activated_session_1==1, :);
large_pre_choice_ensemble_block_3_sem = sem_all_array{1, 3}(exclusive_activated_session_1==1, :);




figure;
shadedErrorBar(ts1, nanmean(large_pre_choice_ensemble_block_1_zall), nanmean(large_pre_choice_ensemble_block_1_sem), 'lineProps', {'color', 'b'});
hold on;shadedErrorBar(ts1, nanmean(large_pre_choice_ensemble_block_2_zall), nanmean(large_pre_choice_ensemble_block_2_sem), 'lineProps', {'color','g'});
hold on;shadedErrorBar(ts1, nanmean(large_pre_choice_ensemble_block_3_zall), nanmean(large_pre_choice_ensemble_block_3_sem), 'lineProps', {'color', 'r'});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'large_pre_choice_ensemble_block_1', 'consumption active', 'neutral'}, 'Location','northwest')

%% POST-CHOICE REWARD
%run data_loop with choiceTime large rew [-10 to 5]
large_pre_choice_ensemble_block_1_zall = zall_mean_all_array{1, 1}(exclusive_activated_session_2==1, :);
large_pre_choice_ensemble_block_1_sem = sem_all_array{1, 1}(exclusive_activated_session_2==1, :);

large_pre_choice_ensemble_block_2_zall = zall_mean_all_array{1, 2}(exclusive_activated_session_2==1, :);
large_pre_choice_ensemble_block_2_sem = sem_all_array{1, 2}(exclusive_activated_session_2==1, :);

large_pre_choice_ensemble_block_3_zall = zall_mean_all_array{1, 3}(exclusive_activated_session_2==1, :);
large_pre_choice_ensemble_block_3_sem = sem_all_array{1, 3}(exclusive_activated_session_2==1, :);




figure;
shadedErrorBar(ts1, nanmean(large_pre_choice_ensemble_block_1_zall), nanmean(large_pre_choice_ensemble_block_1_sem), 'lineProps', {'color', 'b'});
hold on;shadedErrorBar(ts1, nanmean(large_pre_choice_ensemble_block_2_zall), nanmean(large_pre_choice_ensemble_block_2_sem), 'lineProps', {'color','g'});
hold on;shadedErrorBar(ts1, nanmean(large_pre_choice_ensemble_block_3_zall), nanmean(large_pre_choice_ensemble_block_3_sem), 'lineProps', {'color', 'r'});
xline(0);
xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
xline(median_collect_time_from_choice, '--r', {'Median', 'collect', 'latency'})
xlabel('Time from Large Rew Choice (s)');
legend({'large_pre_choice_ensemble_block_1', 'consumption active', 'neutral'}, 'Location','northwest')


%% attempts to correlate shk-identified neurons w/ behavior 
% need to run get_risk_percent_from_BehavData
% load a RDT dataset with the following variables filtered (in order):
    % "choiceTime.Outcome_Minus_4to0.OMITALL_0_BLANK_TOUCH_0_BLOCK_1"
    % "choiceTime.Outcome_0to2.OMITALL_0_BLANK_TOUCH_0_BLOCK_1"
    % "collectionTime.Outcome_1to3.OMITALL_0_BLANK_TOUCH_0_BLOCK_1"
    % "choiceTime.Outcome_0to2.SHK_1"
    % "choiceTime.Outcome_0to2.LOSS_PLUS_ONE_1"


epoch_to_analyze = 'choiceTime'; % choiceTime collectionTime
window_to_analyze = 'Outcome_0to2'; % Minus_4to0  Outcome_0to2 Outcome_1to3
outcome_variable = 'SHK_1'; % OMITALL_0_BLANK_TOUCH_0_BLOCK_1 SHK_1
array_to_select = 4; 

for ff = 1:size(neuron_mean_mouse, 1)
    max_response(ff) = max(mean(neuron_mean_mouse{ff, (array_to_select)}(respClass_mouse.(animalIDs{ff}).(session_to_analyze).(epoch_to_analyze).(window_to_analyze).(outcome_variable)==1, :)));
    min_response(ff) = min(mean(neuron_mean_mouse{ff, (array_to_select)}(respClass_mouse.(animalIDs{ff}).(session_to_analyze).(epoch_to_analyze).(window_to_analyze).(outcome_variable)==2, :)));
    num_shk_cells(ff) = sum(respClass_mouse.(animalIDs{ff}).(session_to_analyze).(epoch_to_analyze).(window_to_analyze).(outcome_variable)==1);
    num_consumption_cells(ff) = sum(respClass_mouse.(animalIDs{ff}).(session_to_analyze).collectionTime.Outcome_1to3.OMITALL_0_BLANK_TOUCH_0_BLOCK_1==1);
    percent_shk_ensemble(ff) = num_shk_cells(ff)/size(respClass_mouse.(animalIDs{ff}).(session_to_analyze).(epoch_to_analyze).(window_to_analyze).(outcome_variable), 1);
end

consum_to_shk_ratio = num_consumption_cells./num_shk_cells;

x = max_response;
y = riskiness;
% Create a new figure with specific dimensions
figure;
width = 650; % Width of the figure
height = 350; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]

hold on;
scatter(x, y, 100)
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

% Add R-squared value to the plot (You can keep this part unchanged)
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared, '%.2f')], 'FontSize', 12);
text(min(x), max(y) - (max(y)-min(y))/10, ['p = ' num2str(pval(2), '%.2f')], 'FontSize', 12, 'Color', 'blue');
hold off; 

x = num_shk_cells;
y = riskiness;
figure;
width = 500; % Width of the figure
height = 1000; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]

hold on;
scatter(x, y, 100)
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

% Add R-squared value to the plot (You can keep this part unchanged)
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared, '%.2f')], 'FontSize', 12);
text(min(x), max(y) - (max(y)-min(y))/10, ['p = ' num2str(pval(2), '%.2f')], 'FontSize', 12, 'Color', 'blue');
hold off;  

x = percent_shk_ensemble;
y = riskiness;
figure;
width = 650; % Width of the figure
height = 350; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
hold on;
scatter(x, y, 100)
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

% Add R-squared value to the plot (You can keep this part unchanged)
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared, '%.2f')], 'FontSize', 12);
text(min(x), max(y) - (max(y)-min(y))/10, ['p = ' num2str(pval(2), '%.2f')], 'FontSize', 12, 'Color', 'blue');
hold off; 


x = percent_exclusive_shk_activated_mouse;
y = riskiness;
figure;
width = 500; % Width of the figure
height = 1000; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]

hold on;
scatter(x, y, 100)
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

% Add R-squared value to the plot (You can keep this part unchanged)
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared, '%.2f')], 'FontSize', 12);
text(min(x), max(y) - (max(y)-min(y))/10, ['p = ' num2str(pval(2), '%.2f')], 'FontSize', 12, 'Color', 'blue');
hold off;  


x = max_activity_exclusive_shk;
y = riskiness;
figure;
width = 500; % Width of the figure
height = 1000; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]

hold on;
scatter(x, y, 100)
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

% Add R-squared value to the plot (You can keep this part unchanged)
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared, '%.2f')], 'FontSize', 12);
text(min(x), max(y) - (max(y)-min(y))/10, ['p = ' num2str(pval(2), '%.2f')], 'FontSize', 12, 'Color', 'blue');
hold off; 

x = max_response;
y = aborts;
figure;
width = 250; % Width of the figure
height = 500; % Height of the figure (width is half of height)
set(gcf, 'Position', [50, 50, width, height]); % Set position and size [left, bottom, width, height]

hold on;
scatter(x, y, 100)
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

% Add R-squared value to the plot (You can keep this part unchanged)
text(min(x) + 0.1, max(y) - 0.1, ['R^2 = ' num2str(r_squared, '%.2f')], 'FontSize', 12);
text(min(x), max(y) - (max(y)-min(y))/10, ['p = ' num2str(pval(2), '%.2f')], 'FontSize', 12, 'Color', 'blue');
hold off; 