for ff = 1:size(respClass_all_array_mouse, 1)
    for dd = 1:size(respClass_all_array_mouse, 2)
        pos_modulated(ff, dd) = sum(respClass_all_array_mouse{ff, dd} == 1);
        neg_modulated(ff, dd) = sum(respClass_all_array_mouse{ff, dd} == 2);
        non_modulated(ff, dd) = sum(respClass_all_array_mouse{ff, dd} == 3);
    end

    


end

% plot these while using yoked trials to get reasonable estimates of small
% (from 4-5 trials in block 1) and large (from yoked 4-5 trials randomly
% selected in block 1)

small_but_not_large = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 3
sum(small_but_not_large)
figure; plot(ts1, mean(neuron_mean_array{1,1}(small_but_not_large  ==1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1,2}(small_but_not_large  ==1, :)))
large_but_not_small = respClass_all_array{1, 1} == 3 & respClass_all_array{1, 2} == 1
sum(large_but_not_small)
figure; plot(ts1, mean(neuron_mean_array{1,1}(large_but_not_small  ==1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1,2}(large_but_not_small  ==1, :)))
small_and_large = respClass_all_array{1, 1} == 1 & respClass_all_array{1, 2} == 1
sum(small_and_large)

figure; plot(ts1, mean(neuron_mean_array{1,1}(small_and_large  ==1, :)))
hold on; plot(ts1, mean(neuron_mean_array{1,2}(small_and_large  ==1, :)))

figure; imagesc(ts1, [], neuron_mean_array{1,1}(small_but_not_large  ==1, :))
figure; imagesc(ts1, [], neuron_mean_array{1,2}(small_but_not_large  ==1, :))

figure; imagesc(ts1, [], neuron_mean_array{1,1}(large_but_not_small  ==1, :))
figure; imagesc(ts1, [], neuron_mean_array{1,2}(large_but_not_small  ==1, :))

figure; imagesc(ts1, [], neuron_mean_array{1,1}(small_and_large  ==1, :))
figure; imagesc(ts1, [], neuron_mean_array{1,2}(small_and_large  ==1, :))