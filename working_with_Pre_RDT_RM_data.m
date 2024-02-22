
not_active = neuron_num - (exclusive_activated_session_1_sum + exclusive_activated_session_2_sum);

excited_to_excited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1;
excited_to_excited_sum = sum(excited_to_excited);

exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1);
exclusive_activated_session_2 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2);



%%
% Example 2: Nested pie chart with custom colors for each wedge

% Initialize data points
inner_pie = [exclusive_activated_session_1_sum/neuron_num,...
            
            exclusive_activated_session_2_sum/neuron_num,...
           
            excited_to_excited_sum/neuron_num,...
            
            
            not_active/neuron_num];

figure; donutchart(inner_pie, 'InnerRadius', 0.8)
figure; pie(inner_pie)


%%
exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 3;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1);
exclusive_activated_session_2 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1 & respClass_all_array{1,3} == 3;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2);
exclusive_activated_session_3 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3 & respClass_all_array{1,3} == 1;
exclusive_activated_session_3_sum = sum(exclusive_activated_session_3);


not_active = neuron_num - (exclusive_activated_session_1_sum + exclusive_activated_session_2_sum + exclusive_activated_session_3_sum);




