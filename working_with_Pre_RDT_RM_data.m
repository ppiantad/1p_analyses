


excited_to_excited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1;
excited_to_excited_sum = sum(excited_to_excited);

exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} ~= 1;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1);
exclusive_activated_session_2 = respClass_all_array{1,1} ~= 1 & respClass_all_array{1,2} == 1;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2);

neutral = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3;



not_active = neuron_num - (exclusive_activated_session_1_sum + exclusive_activated_session_2_sum);

%%
% Example 2: Nested pie chart with custom colors for each wedge

% Initialize data points
inner_pie = [exclusive_activated_session_1_sum/neuron_num,...
            
            exclusive_activated_session_2_sum/neuron_num,...
           
            excited_to_excited_sum/neuron_num,...
            
            
            not_active/neuron_num];

figure; donutchart(inner_pie, 'InnerRadius', 0.5)
figure; pie(inner_pie)


%%
exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 3;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1);
exclusive_activated_session_2 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} == 3;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2);
exclusive_activated_session_3 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 1;
exclusive_activated_session_3_sum = sum(exclusive_activated_session_3);


not_active = neuron_num - (exclusive_activated_session_1_sum + exclusive_activated_session_2_sum + exclusive_activated_session_3_sum);



%% enter this section having saved respClass_all and respClass_mouse
%run data_loop with collectionTime large rew [-2 to 8]
large_consumption_ensemble_zall = zall_mean_all(respClass_all_array{1, 2}==1, :);
large_consumption_ensemble_sem = sem_all(respClass_all_array{1, 2}==1, :);


%run data_loop with collectionTime small rew [-2 to 8]
small_consumption_ensemble_zall = zall_mean_all(respClass_all_array{1, 2}==1, :);
small_consumption_ensemble_sem = sem_all(respClass_all_array{1, 2}==1, :);

figure; plot(ts1, mean(large_consumption_ensemble_zall));
hold on; plot(ts1, mean(small_consumption_ensemble_zall));
%%
figure;
shadedErrorBar(ts1, mean(small_consumption_ensemble_zall), mean(small_consumption_ensemble_sem));
hold on; shadedErrorBar(ts1, mean(large_consumption_ensemble_zall), mean(large_consumption_ensemble_sem));

%%
%run data_loop with choiceTime large rew [-10 to 5]
large_action_ensemble_zall = zall_mean_all(respClass_all_array{1, 1}==1, :);
large_action_ensemble_sem = sem_all(respClass_all_array{1, 1}==1, :);


%run data_loop with choiceTime small rew [-10 to 5]
small_action_ensemble_zall = zall_mean_all(respClass_all_array{1, 1}==1, :);
small_action_ensemble_sem = sem_all(respClass_all_array{1, 1}==1, :);


figure; plot(ts1, mean(large_action_ensemble_zall ));
hold on; plot(ts1, mean(small_action_ensemble_zall));

%%
figure;
shadedErrorBar(ts1, mean(small_action_ensemble_zall), mean(small_action_ensemble_sem));
hold on; shadedErrorBar(ts1, mean(large_action_ensemble_zall), mean(large_action_ensemble_sem));