
% this doesn't work particularly well because Ca is so sparse, selecting a
% random period within a block ends up with very low correlations

block_1_times = [BehavData.stTime(BehavData.Block == 1) BehavData.collectionTime(BehavData.Block == 1)]
block_2_times = [BehavData.stTime(BehavData.Block == 2) BehavData.collectionTime(BehavData.Block == 2)]
block_3_times = [BehavData.stTime(BehavData.Block == 3) BehavData.collectionTime(BehavData.Block == 3)]

block_1_session_ind = time_array > min(block_1_times(1, 1)) & time_array < max(block_1_times(end, 2));
block_2_session_ind = time_array > min(block_2_times(1, 1)) & time_array < max(block_2_times(end, 2));
block_3_session_ind = time_array > min(block_3_times(1, 1)) & time_array < max(block_3_times(end, 2));

% Ensure there are enough true values in the original indices
assert(sum(block_1_session_ind) >= 3000, 'Not enough true values in block_1_session_ind');
assert(sum(block_2_session_ind) >= 3000, 'Not enough true values in block_2_session_ind');
assert(sum(block_3_session_ind) >= 3000, 'Not enough true values in block_3_session_ind');

% Find the indices of true values for each block
block_1_true_ind = find(block_1_session_ind);
block_2_true_ind = find(block_2_session_ind);
block_3_true_ind = find(block_3_session_ind);

% Randomly select a starting index, ensuring there are at least 3000 consecutive true values
start_1 = randi([1, length(block_1_true_ind) - 3000 + 1]);
start_2 = randi([1, length(block_2_true_ind) - 3000 + 1]);
start_3 = randi([1, length(block_3_true_ind) - 3000 + 1]);

% Select 3000 consecutive true indices starting from the random index
block_1_sub_ind = block_1_true_ind(start_1:(start_1 + 2999));
block_2_sub_ind = block_2_true_ind(start_2:(start_2 + 2999));
block_3_sub_ind = block_3_true_ind(start_3:(start_3 + 2999));


% Create new logical indices for the sub-indices
new_block_1_session_ind = false(size(block_1_session_ind));
new_block_2_session_ind = false(size(block_2_session_ind));
new_block_3_session_ind = false(size(block_3_session_ind));

new_block_1_session_ind(block_1_sub_ind) = true;
new_block_2_session_ind(block_2_sub_ind) = true;
new_block_3_session_ind(block_3_sub_ind) = true;

% Display the new sub-indices (optional)
disp('New Block 1 Sub Index:');
disp(new_block_1_session_ind);
disp('New Block 2 Sub Index:');
disp(new_block_2_session_ind);
disp('New Block 3 Sub Index:');
disp(new_block_3_session_ind);

%%

ca_block_1_ind = ca(:, new_block_1_session_ind);
ca_block_2_ind = ca(:, new_block_2_session_ind);
ca_block_3_ind = ca(:, new_block_3_session_ind);



% Get the number of neurons (rows)
num_neurons = size(ca_block_1_ind, 1);

% Preallocate array to store correlation coefficients
corr_block_1_vs_block_2 = zeros(num_neurons, 1);

% Loop through each neuron (row) and calculate the correlation between block 1 and block 2
for neuron = 1:num_neurons
    % Get the calcium activity for the current neuron in both blocks
    ca_block_1_neuron = ca_block_1_ind(neuron, :);
    ca_block_2_neuron = ca_block_2_ind(neuron, :);
    
    % Compute the correlation coefficient between the neuron activity in block 1 and block 2
    R = corrcoef(ca_block_1_neuron, ca_block_2_neuron);
    
    % Store the correlation value (R(1,2) is the correlation between the two arrays)
    corr_block_1_vs_block_2(neuron) = R(1, 2);
end

% Display the resulting correlation coefficients
disp('Correlation coefficients between ca_block_1 and ca_block_2:');
disp(corr_block_1_vs_block_2);

% Preallocate array to store correlation coefficients between block 2 and block 3
corr_block_2_vs_block_3 = zeros(num_neurons, 1);

% Loop through each neuron and calculate the correlation between block 2 and block 3
for neuron = 1:num_neurons
    ca_block_2_neuron = ca_block_2_ind(neuron, :);
    ca_block_3_neuron = ca_block_3_ind(neuron, :);
    
    R = corrcoef(ca_block_2_neuron, ca_block_3_neuron);
    
    corr_block_2_vs_block_3(neuron) = R(1, 2);
end

% Display the resulting correlation coefficients
disp('Correlation coefficients between ca_block_2 and ca_block_3:');
disp(corr_block_2_vs_block_3);
