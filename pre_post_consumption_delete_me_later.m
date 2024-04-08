pre_choice_only = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 3 & respClass_all_array{1, 3} == 3
sum(pre_choice_only)
post_choice_only = respClass_all_array{1, 1} == 3 & respClass_all_array{1, 2} == 1 & respClass_all_array{1, 3} == 3
sum(post_choice_only)
consumption_only = respClass_all_array{1, 1} == 3 & respClass_all_array{1, 2} == 3 & respClass_all_array{1, 3} == 1
sum(consumption_only)



figure; plot(ts1, neuron_mean_array{1, 1}(pre_choice_only, :))
figure; plot(ts1, mean(neuron_mean_array{1, 1}(pre_choice_only, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, 1}(post_choice_only, :)))
hold on; plot(ts1, mean(neuron_mean_array{1, 1}(consumption_only, :)))

