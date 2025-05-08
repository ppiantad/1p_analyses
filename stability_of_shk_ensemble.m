%% load matched RDT_D1 & RDT_D2 data / RDT_D1 & SHOCK_TEST data
% filter on SHK, 1

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

% String to compare
targetAnimal = 'BLA_Insc_27';

% Perform element-wise comparison
animal_index_to_plot = find(strcmp(animalIDs, targetAnimal));

session_to_analyze_1 = 'RDT_D2';
session_to_analyze_2 = 'RDT_D1'; %SHOCK TEST
Coor = final.(targetAnimal).(session_to_analyze_2).CNMFe_data.Coor;  

Cn_data = final.(targetAnimal).(session_to_analyze_2).CNMFe_data.Cn;
figure;
imagesc(Cn_data) %only one session will look good here typically



% these variables need to be customized based on what you want to plot
overlap_neurons = respClass_mouse.(targetAnimal).(session_to_analyze_1).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) == 1 & respClass_mouse.(targetAnimal).(session_to_analyze_2).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) == 1;
RDT_only_neurons = respClass_mouse.(targetAnimal).(session_to_analyze_1).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) == 1 & respClass_mouse.(targetAnimal).(session_to_analyze_2).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) ~= 1;
shock_test_only_neurons = respClass_mouse.(targetAnimal).(session_to_analyze_1).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) ~= 1 & respClass_mouse.(targetAnimal).(session_to_analyze_2).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) == 1;
neutral = respClass_mouse.(targetAnimal).(session_to_analyze_1).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) ~= 1 & respClass_mouse.(targetAnimal).(session_to_analyze_2).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) ~= 1;
% Create a combined array
Identity = zeros(size(Coor));

% Assign values based on conditions
Identity(overlap_neurons) = 1;
Identity(RDT_only_neurons) = 2;
Identity(shock_test_only_neurons) = 3;
Identity(neutral) = 4;
% 
% Coor = final.(targetAnimal).(session_to_analyze).CNMFe_data.Coor;  
% 
% % these variables need to be customized based on what you want to plot
% post_choice_reward = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) == 1;
% post_choice_reward_shk = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) == 1;
% both = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) == 1 & respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) == 1;
% neutral = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) ~= 1 & respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) ~= 1;
% 
% % Create a combined array
% Identity = zeros(size(Coor));
% 
% % Assign values based on conditions
% Identity(post_choice_reward) = 1;
% Identity(post_choice_reward_shk) = 2;
% Identity(both) = 3;
% Identity(neutral) = 4;
% 
% Cn_data = final.(targetAnimal).(session_to_analyze).CNMFe_data.Cn;
% figure;
% imagesc(Cn_data) %only one session will look good here typically



% Coor = final.(targetAnimal).(session_to_analyze).CNMFe_data.Coor;  
% 
% % these variables need to be customized based on what you want to plot
% post_choice_reward = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) == 1;
% post_choice_reward_shk = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) == 1;
% both = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) == 1 & respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) == 1;
% neutral = respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{1}).(identity_class_string_all{1}).(all_filter_args{1}) ~= 1 & respClass_mouse.(targetAnimal).(session_to_analyze).(epoc_to_align_all{2}).(identity_class_string_all{2}).(all_filter_args{2}) ~= 1;
% 
% % Create a combined array
% Identity = zeros(size(Coor));
% 
% % Assign values based on conditions
% Identity(post_choice_reward) = 1;
% Identity(post_choice_reward_shk) = 2;
% Identity(both) = 3;
% Identity(neutral) = 4;





%%
figure;
imagesc(Cn_data) %only one session will look good here typically
colormap gray;
hold on;
% Calculate the minimum and maximum values in the data
min_val = min(Cn_data(:));
max_val = max(Cn_data(:));

% Define the desired range for the color axis
new_min = .6;  % Set the minimum value you want to display
new_max = 1;  % Set the maximum value you want to display

% Update the color axis scaling
caxis([new_min, new_max]);



for i = 1:numel(Coor)
    % Get the coordinates of the current circle
    circleCoords = Coor{i};
    %
    % Extract x and y coordinates
    x = circleCoords(1, :);
    y = circleCoords(2, :);

    % Calculate centroid
    centroid = [mean(x), mean(y)];

    % % Store centroid in the array
    % centroids(i, :) = round(centroid);


    % if post_choice_reward(i) == 1
    %     plot_color = "red";
    % end
    % if post_choice_reward_shk(i) == 1
    %     plot_color = "blue";
    % end
    % if both(i) == 1
    %     plot_color = "green";
    % end
    % if neutral(i) == 1
    %     plot_color = "black";
    % end


    if Identity(i) == 1
        plot_color = "red";
    elseif Identity(i) == 2
        plot_color = "blue";
    elseif Identity(i) == 3
        plot_color = "green";
    elseif Identity(i) == 4
        plot_color = "black";
    end


    % Plot the circle
    fill(x, y, plot_color);


    % Add text next to circles #42, #46, and #54
    if i == 79 || i == 46 || i == 54
        text(centroid(1), centroid(2), num2str(i), 'Color', 'white', 'FontSize', 12, 'HorizontalAlignment', 'center');
    end


end

hold off;

% Add labels and title
xlabel('X');
ylabel('Y');
% title('Plot of Circles');

% Adjust the aspect ratio if needed
axis equal;


%%
figure;
hold on
width = 400; % Width of the figure
height = 550; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
ylim([-0.5 1.3]);
set(gca, 'XTick', [-8, 0, 8], 'YTick', [-0.5 0 0.5 1.0]);
shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 1}(shk_consistent==1, :)), nanmean(neuron_sem_array{1, 1}(shk_consistent==1, :)), 'lineProps', {'color','k'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_array{1, 2}(shk_consistent==1, :)), nanmean(neuron_sem_array{1, 2}(shk_consistent==1, :)), 'lineProps', {'color', 'r'});
% xtickformat('%.2f');
ytickformat('%.1f');
xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice/shock (s)');
% legend({'pre-choice active', 'post-choice reward active', 'consumption', 'not active'}, 'Location','northwest')

%% for checking whether unnormalized traces look similar - just in case the comparison as plotted above w/ zscores doesn't fly with reviewers? can use the data below, zscore all in one array, and plot separately (still need to do this)
plot_me = 3;

figure;
hold on
width = 400; % Width of the figure
height = 550; % Height of the figure (width is half of height)
set(gcf, 'Position', [100, 100, width, height]); % Set position and size [left, bottom, width, height]
xlim([-8 8]);
ylim([-0.5 1.2]);
set(gca, 'XTick', [-8, 0, 8], 'YTick', [-0.5 0 0.5 1.0]);
shadedErrorBar(ts1, nanmean(neuron_mean_mouse_unnormalized{plot_me , 1}(respClass_all_array_mouse{plot_me , 1}==1, :)), nanmean(neuron_sem_mouse_unnormalized{plot_me , 1} (respClass_all_array_mouse{plot_me , 1}==1, :)), 'lineProps', {'color','k'});
hold on;shadedErrorBar(ts1, nanmean(neuron_mean_mouse_unnormalized{plot_me , 2}(respClass_all_array_mouse{plot_me , 2}==1, :)), nanmean(neuron_sem_mouse_unnormalized{plot_me , 2} (respClass_all_array_mouse{plot_me , 2}==1, :)), 'lineProps', {'color', 'r'});
% xtickformat('%.2f');
ytickformat('%.2f');
xline(0);
% xline(median_start_time_from_choice, 'g', {'Median', 'start', 'time'})
% xline(median_collect_time_from_choice, 'r', {'Median', 'collect', 'latency'})
xlabel('Time from choice/shock (s)');
% legend({'pre-choice active', 'post-choice reward active', 'consumption', 'not active'}, 'Location','northwest')