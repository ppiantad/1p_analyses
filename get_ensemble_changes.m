% Initialize empty arrays to store the indices
indices_2_1 = [];
indices_3_1 = [];

% Iterate over the cells in the 2nd column ({1,2})
for i = 1:numel(respClass_all_array{1,2})
    % Check if the value is 1 in the 2nd column, but not in the 1st column
    if respClass_all_array{1,2}(i) == 1 && respClass_all_array{1,1}(i) ~= 1
        indices_2_1 = [indices_2_1, i]; % Store the index
    end
end

% Iterate over the cells in the 3rd column ({1,3})
for i = 1:numel(respClass_all_array{1,3})
    % Check if the value is 1 in the 3rd column, but not in the 1st column
    if respClass_all_array{1,3}(i) == 1 && respClass_all_array{1,1}(i) ~= 1
        indices_3_1 = [indices_3_1, i]; % Store the index
    end
end

% % Display the indices for the comparison between {1,1} and {1,2}
% disp("Indices for comparison between {1,1} and {1,2}:");
% disp(indices_2_1);
% 
% % Display the indices for the comparison between {1,1} and {1,3}
% disp("Indices for comparison between {1,1} and {1,3}:");
% disp(indices_3_1);


%%
%classic loss of ensemble strength
figure; 
plot(nanmean(neuron_mean_array{1,1}(respClass_all_array{1,1} == 1,:)))
hold on; plot(nanmean(neuron_mean_array{1,2}(respClass_all_array{1,1} == 1,:)))
hold on; plot(nanmean(neuron_mean_array{1,3}(respClass_all_array{1,1} == 1,:)))

%uut actually ensemble strength is maintained, just by different neurons
figure; 
plot(nanmean(neuron_mean_array{1,1}(respClass_all_array{1,1} == 1,:)))
hold on; plot(nanmean(neuron_mean_array{1,2}(respClass_all_array{1,2} == 1,:)))
hold on; plot(nanmean(neuron_mean_array{1,3}(respClass_all_array{1,3} == 1,:)))

%here are the different neurons
figure; 
plot(nanmean(neuron_mean_array{1,1}(respClass_all_array{1,1} == 1,:)))
hold on; plot(nanmean(neuron_mean_array{1,2}(indices_2_1,:)))
hold on; plot(nanmean(neuron_mean_array{1,3}(indices_3_1,:)))



