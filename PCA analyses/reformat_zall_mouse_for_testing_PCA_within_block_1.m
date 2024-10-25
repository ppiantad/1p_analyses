% array_to_reformat = unnormalized_by_mouse;
array_to_reformat = zall_mouse;

neuron_num = 0; 
for hh = 1:size(array_to_reformat, 1)
    current_array = array_to_reformat{hh, 1};
    trial_num_to_randomize = size(current_array{1, 1}, 1);
    trial_indices_to_randomize = randperm(trial_num_to_randomize);
    half = ceil(trial_num_to_randomize/2);
    indices_first = trial_indices_to_randomize(1:half); 
    indices_second = trial_indices_to_randomize(half+1:end); 
    for qq = 1:size(current_array, 2)
        
       
        
        
        neuron_num = neuron_num+1; 
        reformat_zall_mean_array{1, 1}(neuron_num, :) = mean(current_array{1, qq}(indices_first, :));
        reformat_zall_mean_array{1, 2}(neuron_num, :) = mean(current_array{1, qq}(indices_second, :));
        

    end




end