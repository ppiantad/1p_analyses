%%
excited_to_excited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1;
inhibited_to_inhibited = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 2;
excited_to_excited_sum = sum(excited_to_excited);
indices_excited_to_excited = find(excited_to_excited == 1);
inhibited_to_inhibited_sum = sum(inhibited_to_inhibited);

excited_to_inhibited = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 2;
excited_to_inhibited_sum = sum(excited_to_inhibited);
excited_to_neutral = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3;
excited_to_neutral_sum = sum(excited_to_neutral);

inhibited_to_excited = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 1;
inhibited_to_excited_sum = sum(inhibited_to_excited);
inhibited_to_neutral = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 3;
inhibited_to_neutral_sum = sum(inhibited_to_neutral);

% above this line are good for single session stacked plots (e.g., how do
% neurons on session 1 that are activated (or inhibited) change to session
% 2?) 

neutral_to_excited = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1;
neutral_to_excited_sum = sum(neutral_to_excited);

neutral_to_inhibited = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 2;
neutral_to_inhibited_sum = sum(neutral_to_inhibited);

% need to add excited_to_excited twice because they are present in BOTH
% sessions
total_activated_possible = excited_to_excited_sum+ excited_to_excited_sum + excited_to_inhibited_sum+ excited_to_neutral_sum+ neutral_to_excited_sum + inhibited_to_excited_sum;
% need to add inhibited_to_inhibited twice because they are present in BOTH
% sessions
total_inhibited_possible = inhibited_to_inhibited_sum + inhibited_to_inhibited_sum + inhibited_to_excited_sum + inhibited_to_neutral_sum + neutral_to_inhibited_sum + excited_to_inhibited_sum;






neutral_to_neutral = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 3;
neutral_to_neutral_sum = sum(neutral_to_neutral);

exclusive_activated_session_1 = respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 3;
exclusive_activated_session_1_sum = sum(exclusive_activated_session_1);
exclusive_activated_session_2 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 1;
exclusive_activated_session_2_sum = sum(exclusive_activated_session_2);



exclusive_activated_session_1_indices = find(exclusive_activated_session_1(1,:) == 1);


exclusive_inhibited_session_1 = respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 3;
exclusive_inhibited_session_1_sum = sum(exclusive_inhibited_session_1);
exclusive_inhibited_session_2 = respClass_all_array{1,1} == 3 & respClass_all_array{1,2} == 2;
exclusive_inhibited_session_2_sum = sum(exclusive_inhibited_session_2);
% exclusive_activated_sum = sum(exclusive_activated);
% exclusive_inhibited_sum = sum(exclusive_inhibited);

%%
Block_1_activated = respClass_all_array{1,1} == 1;
Block_1_activated_sum = sum(Block_1_activated);
Block_2_activated = respClass_all_array{1,2} == 1;
Block_2_activated_sum = sum(Block_2_activated);
Block_2_consistent_activated = Block_1_activated & respClass_all_array{1,2} == 1;
disp(sum(Block_2_consistent_activated));

Block_1_responsive = respClass_all_array{1,1} == 1 | respClass_all_array{1,2} == 2;
disp(sum(Block_1_responsive));
Block_2_consistent_responsive = (respClass_all_array{1,1} == 1 & respClass_all_array{1,2} == 1) | (respClass_all_array{1,1} == 2 & respClass_all_array{1,2} == 2);
disp(sum(Block_2_consistent_responsive));




% this is the start of checking if neurons are MORE active than during
% other events, i.e. if you wanted to check if REW activated neurons are
% significantly differentially activated by SHK. preliminary poking around
% seems to suggest that few large reward active neurons have their activity
% increase in response to SHK
co_activated_indices = find(excited_to_excited(1,:) == 1);
co_activated_indices_sum = numel(co_activated_indices);
for qq = 1:size(co_activated_indices, 2)
    [h(qq),p(qq),ci{qq},stats{qq}] = ttest(neuron_mean_array{1, 1}(co_activated_indices(qq),sub_window_idx),neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx));
    mean_diff(qq) = mean(neuron_mean_array{1, 1}(co_activated_indices(qq),sub_window_idx) - mean(neuron_mean_array{1, 2}(co_activated_indices(qq),sub_window_idx)));
end

sig_increase_shk_from_large = co_activated_indices(p < 0.05 & mean_diff < 0);
sig_increase_shk_from_large_sum = numel(sig_increase_shk_from_large);
sig_increase_shk_from_large_ind = zeros(1, size(respClass_all_array{1,2}, 2));

sig_increase_shk_from_large_ind(:, sig_increase_shk_from_large) = 1;

shk_activated = respClass_all_array{1,2} == exclusive_activated_session_2 |  respClass_all_array{1,2} == sig_increase_shk_from_large_ind;
shk_activated_sum = sum(shk_activated)

figure; plot(ts1, mean(neuron_mean_array{1, 2}(shk_activated,:))); hold on; plot(ts1,  mean(neuron_mean_array{1, 1}(Block_1_activated,:)));



figure;
hold on; 
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(shk_activated,:)), mean(neuron_sem_array{1, 2}(shk_activated,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 1}(Block_1_activated,:)), mean(neuron_sem_array{1, 1}(Block_1_activated,:)), 'lineProps', {'color', batlowW(iter,:)});




% total_modulated = [sum_activated + sum_inhibited];
% total_co_modulated = [co_activated_sum + co_inhibited_sum];


total_modulated = [(Block_2_activated_sum/neuron_num)*100 (shk_activated_sum/neuron_num)*100];

A = total_modulated;
I = (sig_increase_shk_from_large_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end




total_modulated = [(Block_1_activated_sum/neuron_num)*100 (Block_2_activated_sum/neuron_num)*100];

A = total_modulated;
I = (co_activated_indices_sum/neuron_num)*100;
K = [A I];
figure; 
[H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
for i = 1:size(K, 2)
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),  [num2str(K(1,i))])
end


figure;
hold on; 
shadedErrorBar(ts1, mean(neuron_mean_array{1, 2}(exclusive_activated_session_2,:)), mean(neuron_sem_array{1, 2}(exclusive_activated_session_2,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 1}(exclusive_activated_session_1,:)), mean(neuron_sem_array{1, 1}(exclusive_activated_session_1,:)), 'lineProps', {'color', batlowW(iter,:)});
shadedErrorBar(ts1, mean(neuron_mean_array{1, 1}(excited_to_excited,:)), mean(neuron_sem_array{1, 1}(excited_to_excited,:)), 'lineProps', {'color', batlowW(iter,:)});

