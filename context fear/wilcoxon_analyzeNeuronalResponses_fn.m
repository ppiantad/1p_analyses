% Analyze if neurons are significantly responsive to an event
% Input:
%   zall_array: Cell array where each cell contains a double array of neuronal activity
%               (rows = trials, columns = time samples)
%   ts1: Time array (1x120 double) from -2 to 10 in 0.1s increments

function modulation = wilcoxon_analyzeNeuronalResponses_fn(zall_array, ts1, pre_period_window, post_period_window)
    % Number of neurons
    numNeurons = length(zall_array);
    
    % Initialize output array: 0=no modulation, 1=positive, 2=negative
    modulation = zeros(numNeurons, 1);
    
    % Find indices for pre-event and post-event periods
    % Pre-event: 2 seconds before time=0
    % Post-event: 2 seconds after time=0
    zero_idx = find(ts1 >= 0, 1, 'first');
    pre_start_idx = find(ts1 >= pre_period_window(1), 1, 'first');
    pre_end_idx = find(ts1 >= pre_period_window(2), 1, 'first') - 1;
    post_start_idx = find(ts1 >= post_period_window(1), 1, 'first');
    post_end_idx = find(ts1 <= post_period_window(2), 1, 'last');
    
    % For each neuron
    for n = 1:numNeurons
        % Get neuronal data
        neuron_data = zall_array{n};
        
        % Calculate AUC for each trial before and after the event
        num_trials = size(neuron_data, 1);
        pre_auc = zeros(num_trials, 1);
        post_auc = zeros(num_trials, 1);
        
        for trial = 1:num_trials
            % Calculate AUC (Area Under Curve) using trapezoidal integration
            % Pre-event period
            pre_auc(trial) = trapz(ts1(pre_start_idx:pre_end_idx), neuron_data(trial, pre_start_idx:pre_end_idx));
            
            % Post-event period
            post_auc(trial) = trapz(ts1(post_start_idx:post_end_idx), neuron_data(trial, post_start_idx:post_end_idx));
        end
        
        % Perform Wilcoxon signed-rank test
        [p, h, stats] = signrank(pre_auc, post_auc);
        
        % Determine modulation based on p-value and sign of effect
        if h == 1  % Significant difference (p < 0.05)
            % Calculate the median difference to determine direction
            median_diff = median(post_auc - pre_auc);
            
            if median_diff > 0  % post > pre (positive modulation)
                modulation(n) = 1;
            else  % post < pre (negative modulation)
                modulation(n) = 2;
            end
        end
        % If not significant, modulation remains 0
    end
end