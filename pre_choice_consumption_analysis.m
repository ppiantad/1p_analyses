all_pre_choice_activated = respClass_all_array{1,1} == 1;
all_pre_choice_inhibited = respClass_all_array{1,1} == 2;
pre_choice_activated_exclusive = respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 3; 
pre_choice_inhibited_exclusive = respClass_all_array{1,1} == 2 & respClass_all_array{1,3} == 3; 
pre_choice_activated_consum_inhibited = respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 2; 
pre_choice_activated_consum_activated = respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 1; 
pre_choice_inhibited_consum_inhibited = respClass_all_array{1,1} == 2 & respClass_all_array{1,3} == 2; 
pre_choice_inhibited_consum_activated = respClass_all_array{1,1} == 2 & respClass_all_array{1,3} == 1; 



all_consum_activated = respClass_all_array{1,3} == 1;
all_consum_inhibited = respClass_all_array{1,3} == 2;
consum_activated_exclusive = respClass_all_array{1,3} == 1 & respClass_all_array{1,1} == 3; 
consum_inhibited_exclusive = respClass_all_array{1,3} == 2 & respClass_all_array{1,1} == 3; 
consum_activated_pre_choice_inhibited = respClass_all_array{1,3} == 1 & respClass_all_array{1,1} == 2; 
consum_activated_pre_choice_activated = respClass_all_array{1,3} == 1 & respClass_all_array{1,1} == 1; 
consum_inhibited_pre_choice_inhibited = respClass_all_array{1,3} == 2 & respClass_all_array{1,1} == 2; 
consum_inhibited_pre_choice_activated = respClass_all_array{1,3} == 2 & respClass_all_array{1,1} == 1; 


not_active = neuron_num - (sum(all_pre_choice_activated) + sum(all_pre_choice_inhibited) + sum(all_consum_activated) + sum(all_consum_inhibited));
sum(not_active)

% Initialize data points
inner_pie = [sum(all_pre_choice_activated)/neuron_num,...

            sum(all_pre_choice_inhibited)/neuron_num,...
            
            sum(all_consum_activated)/neuron_num,...

            sum(all_consum_inhibited)/neuron_num,...
            
            not_active/neuron_num];


figure; pie(inner_pie)

outer_donut = [ sum(pre_choice_activated_exclusive)/neuron_num...
                
                sum(pre_choice_activated_consum_inhibited)/neuron_num...

                sum(pre_choice_activated_consum_activated)/neuron_num...

                sum(pre_choice_inhibited_exclusive)/neuron_num...
                
                sum(pre_choice_inhibited_consum_inhibited)/neuron_num...

                sum(pre_choice_inhibited_consum_activated)/neuron_num...
                %
                sum(consum_activated_exclusive)/neuron_num...
                
                sum(consum_activated_pre_choice_inhibited)/neuron_num...

                sum(consum_activated_pre_choice_activated)/neuron_num...

                sum(consum_inhibited_exclusive)/neuron_num...
                
                sum(consum_inhibited_pre_choice_inhibited)/neuron_num...

                sum(consum_inhibited_pre_choice_activated)/neuron_num...
                
                not_active/neuron_num...
                



];



figure; donutchart(inner_pie, 'InnerRadius', 0)
figure; donutchart(inner_pie, 'InnerRadius', .3)
figure; donutchart(outer_donut, 'InnerRadius', 0.7)


mean_data_pre_choice_activated_exclusive = neuron_mean_array{1,1}(pre_choice_activated_exclusive, :);
mean_data_consum_activated_exclusive = neuron_mean_array{1,1}(consum_activated_exclusive, :);
mean_data_pre_choice_inhibited_exclusive = neuron_mean_array{1,1}(pre_choice_inhibited_exclusive, :);
mean_data_consum_inhibited_exclusive = neuron_mean_array{1,1}(consum_inhibited_exclusive, :);


data_for_big_heatmap = [mean_data_pre_choice_activated_exclusive; mean_data_consum_activated_exclusive];

figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
imagesc(ts1, [], data_for_big_heatmap)
% Apply the custom colormap
colormap(flipud(gray));

% Restrict the color axis range to [-1, 1]
clim([-.5 .5]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 

ylim([1  size(zall_array{1, plot_num}, 1)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(zall_array{1, plot_num}, 1)]);
xline(0)
hold off;

%%

% Calculate the number of rows for each group
num_rows_pre_choice_activated_exclusive = size(mean_data_pre_choice_activated_exclusive, 1);
num_rows_pre_choice_activated_consum_inhibited = size(mean_data_pre_choice_activated_consum_inhibited, 1);
num_rows_pre_choice_activated_consum_activated = size(mean_data_pre_choice_activated_consum_activated, 1);
num_rows_consum_inhibited_exclusive = size(mean_data_consum_inhibited_exclusive, 1);
num_rows_consum_activated_exclusive = size(mean_data_consum_activated_exclusive, 1);
num_rows_pre_choice_inhibited_consum_inhibited = size(mean_data_pre_choice_inhibited_consum_inhibited, 1);
num_rows_pre_choice_inhibited_consum_activated = size(mean_data_pre_choice_inhibited_consum_activated, 1);
num_rows_pre_choice_inhibited_exclusive = size(mean_data_pre_choice_inhibited_exclusive, 1);

pre_choice_activated_consum_inhibited = respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 2; 
pre_choice_activated_consum_activated = respClass_all_array{1,1} == 1 & respClass_all_array{1,3} == 1; 
pre_choice_inhibited_consum_inhibited = respClass_all_array{1,1} == 2 & respClass_all_array{1,3} == 2; 
pre_choice_inhibited_consum_activated = respClass_all_array{1,1} == 2 & respClass_all_array{1,3} == 1; 

mean_data_pre_choice_activated_consum_inhibited = neuron_mean_array{1,1}(pre_choice_activated_consum_inhibited, :);
mean_data_pre_choice_activated_consum_activated = neuron_mean_array{1,1}(pre_choice_activated_consum_activated, :);
mean_data_pre_choice_inhibited_consum_inhibited = neuron_mean_array{1,1}(pre_choice_inhibited_consum_inhibited, :);
mean_data_pre_choice_inhibited_consum_activated = neuron_mean_array{1,1}(pre_choice_inhibited_consum_activated, :);


data_for_big_heatmap = [    mean_data_pre_choice_activated_exclusive;...
                            mean_data_pre_choice_activated_consum_inhibited;...
                            mean_data_pre_choice_activated_consum_activated;...
                            mean_data_consum_inhibited_exclusive;...
                            mean_data_consum_activated_exclusive;...
                            mean_data_pre_choice_inhibited_consum_inhibited;...
                            mean_data_pre_choice_inhibited_consum_activated;...
                            mean_data_pre_choice_inhibited_exclusive;...
                            ];

% Calculate the positions where to add the horizontal lines
yline_positions = cumsum([num_rows_pre_choice_activated_exclusive, ...
                          num_rows_pre_choice_activated_consum_inhibited, ...
                          num_rows_pre_choice_activated_consum_activated, ...
                          num_rows_consum_inhibited_exclusive, ...
                          num_rows_consum_activated_exclusive, ...
                          num_rows_pre_choice_inhibited_consum_inhibited, ...
                          num_rows_pre_choice_inhibited_consum_activated, ...
                          num_rows_pre_choice_inhibited_exclusive]);



figure('Position', [100, 100, 300, 600]); % [left, bottom, width, height]
imagesc(ts1, [], data_for_big_heatmap)
% Apply the custom colormap
colormap(flipud(gray));

% Restrict the color axis range to [-1, 1]
clim([-.5 .5]);

% Add a separate axes for the colorbar to associate it only with the upper tile
c = colorbar('eastoutside');
set(c, 'YTick', clim); % 


% Add ylines at the end of each set of data
hold on;
for i = 1:length(yline_positions)
    yline(yline_positions(i) + 0.5, 'r--', 'LineWidth', 1); % Adjust +0.5 to place the line between rows
end
hold off;

ylim([1  size(data_for_big_heatmap, 1)])
xlim([-8 8]);
% Set X-axis ticks
set(gca, 'XTick', [-8, -4, 0, 4, 8]);
set(gca, 'YTick', [1, size(data_for_big_heatmap, 1)]);
xline(0)
hold off;


%%
figure;
width = 450; % Width of the figure
height = 650; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]

shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(all_pre_choice_activated==1, :)), nanmean(neuron_sem_array{1, 1}(all_pre_choice_activated==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(all_pre_choice_inhibited==1, :)), nanmean(neuron_sem_array{1, 1}(all_pre_choice_inhibited==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(all_consum_activated==1, :)), nanmean(neuron_sem_array{1, 1}(all_consum_activated==1, :)), 'lineProps', {'color', batlowW(iter,:)});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(all_consum_inhibited==1, :)), nanmean(neuron_sem_array{1, 1}(all_consum_inhibited==1, :)), 'lineProps', {'color', batlowW(iter,:)});
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

%%
% Initialize data points

activated_inner_pie = [sum(all_pre_choice_activated)/neuron_num,...


            
            sum(all_consum_activated)/neuron_num,...


            
            (neuron_num - (sum(all_pre_choice_activated) + sum(all_consum_activated))) / neuron_num];







activated_outer_donut = [ sum(pre_choice_activated_exclusive)/neuron_num...
                
                sum(pre_choice_activated_consum_inhibited)/neuron_num...

                sum(pre_choice_activated_consum_activated)/neuron_num...
                %
                sum(consum_activated_exclusive)/neuron_num...
                
                sum(consum_activated_pre_choice_inhibited)/neuron_num...

                sum(consum_activated_pre_choice_activated)/neuron_num...
                
                (neuron_num - (sum(pre_choice_activated_exclusive) + sum(pre_choice_activated_consum_inhibited) + sum(pre_choice_activated_consum_activated) + sum(consum_activated_exclusive) + sum(consum_activated_pre_choice_inhibited) + sum(consum_activated_pre_choice_activated)))/neuron_num...
                
];

figure; donutchart(activated_inner_pie, 'InnerRadius', .5)
figure; donutchart(activated_outer_donut, 'InnerRadius', 0.7)




inhibited_inner_pie = [sum(all_pre_choice_inhibited)/neuron_num,...

            sum(all_consum_inhibited)/neuron_num,...
            
            (neuron_num - (sum(all_pre_choice_inhibited) + sum(all_consum_inhibited))) / neuron_num];


inhibited_outer_donut = [ sum(pre_choice_inhibited_exclusive)/neuron_num...
                
                sum(pre_choice_inhibited_consum_inhibited)/neuron_num...

                sum(pre_choice_inhibited_consum_activated)/neuron_num...
                %
                sum(consum_inhibited_exclusive)/neuron_num...
                
                sum(consum_inhibited_pre_choice_inhibited)/neuron_num...

                sum(consum_inhibited_pre_choice_activated)/neuron_num...
                
                (neuron_num - (sum(pre_choice_inhibited_exclusive) + sum(pre_choice_inhibited_consum_inhibited) + sum(pre_choice_inhibited_consum_activated) + sum(consum_inhibited_exclusive) + sum(consum_inhibited_pre_choice_inhibited) + sum(consum_inhibited_pre_choice_activated)))/neuron_num...
                
];




figure; donutchart(inhibited_inner_pie, 'InnerRadius', .5)
figure; donutchart(inhibited_outer_donut, 'InnerRadius', 0.7)

%%
percent_pre_choice_activated_exclusive = sum(pre_choice_activated_exclusive)/sum(all_pre_choice_activated)*100;
percent_pre_choice_activated_consum_inhibited = sum(pre_choice_activated_consum_inhibited)/sum(all_pre_choice_activated)*100;
percent_pre_choice_activated_consum_activated = sum(pre_choice_activated_consum_activated)/sum(all_pre_choice_activated)*100;

figure;
bar(1, [percent_pre_choice_activated_exclusive; percent_pre_choice_activated_consum_inhibited; percent_pre_choice_activated_consum_activated], 'stacked')

%%
percent_consum_activated_exclusive = sum(consum_activated_exclusive)/sum(all_consum_activated)*100;
percent_consum_activated_pre_choice_inhibited = sum(consum_activated_pre_choice_inhibited)/sum(all_consum_activated)*100;
percent_consum_activated_pre_choice_activated = sum(consum_activated_pre_choice_activated)/sum(all_consum_activated)*100;

figure;
bar(1, [percent_consum_activated_exclusive; percent_consum_activated_pre_choice_inhibited; percent_consum_activated_pre_choice_activated], 'stacked')

%%
percent_pre_choice_active = (sum(all_pre_choice_activated)/neuron_num)*100
percent_consum_inhib = (sum(all_consum_inhibited)/neuron_num)*100
percent_overlap_from_total_neurons =  (sum(pre_choice_activated_consum_inhibited)/neuron_num)*100
total_modulated = [percent_pre_choice_active percent_consum_inhib];


A = total_modulated;
I = percent_overlap_from_total_neurons;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end

%%
percent_pre_choice_inhib = (sum(all_pre_choice_inhibited)/neuron_num)*100
percent_consum_active = (sum(all_consum_activated)/neuron_num)*100
percent_overlap_from_total_neurons =  (sum(pre_choice_inhibited_consum_activated)/neuron_num)*100
total_modulated = [percent_pre_choice_inhib percent_consum_active];


A = total_modulated;
I = percent_overlap_from_total_neurons;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end
%% requires https://www.mathworks.com/matlabcentral/fileexchange/98974-venn-euler-diagram?s_tid=FX_rc3_behav
% this outputs a ever so slightly wonky diagram. a few nodes that do not
% actually overlap minimally overlap (but intersections are 0), and 1 node
% that has 1 overlap does not overlap at all. 
pre_choice_active_ind = find(all_pre_choice_activated == 1);
consum_active_ind = find(all_consum_activated == 1);
pre_choice_inhibited_ind = find(all_pre_choice_inhibited == 1);
consum_inhibited_ind = find(all_consum_inhibited == 1);
setListData = {pre_choice_active_ind, consum_active_ind, pre_choice_inhibited_ind, consum_inhibited_ind};
setLabels = ["Pre-choice excited", "Consumption excited", "Pre-choice inhibited", "Consumption inhibited"];

h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

h.ShowIntersectionCounts = true;
h.ShowIntersectionAreas = true;
% h.SetLabels = [];