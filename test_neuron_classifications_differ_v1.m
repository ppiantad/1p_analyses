% Define time range
start_time = 1; % seconds
end_time = 3;   % seconds

% Convert time range to indices in ts1
start_index = find(ts1(1, :) >= start_time, 1);
end_index = find(ts1(1, :) <= end_time, 1, 'last');

% Initialize arrays to store comparison results
num_cells = size(zall_array, 2);
num_columns = size(zall_array{1}, 2);
comparison_results = zeros(num_cells, num_columns);

% Perform comparison for each column in each cell of zall_array
for cell_idx = 1:num_cells
    for col_idx = 1:num_columns
        % Extract data within the specified time range
        data_cell = zall_array{cell_idx, col_idx};
        data_time_range = data_cell(start_index:end_index);
        
        % Perform comparison (example: mean value)
        comparison_results(cell_idx, col_idx) = mean(data_time_range);
    end
end

% Now comparison_results contains the comparison result for each column
% in each cell within the specified time range (1 to 3 seconds)


differing_indexs = respClass_all_array{1, 1}  == 1 & respClass_all_array{1, 2} ~=1;

%%

plot_neuron_num = 909; 

first_array = caTraceTrials_unnormalized_array{1, plot_neuron_num};
second_array = caTraceTrials_unnormalized_array{2, plot_neuron_num};

array = [first_array; second_array]; 

events = [(ones(size(first_array, 1), 1)); (ones(size(second_array, 1), 1)*2)]; 

array = zscore(array, 0, 2); 


TrialWin_first_array = array(events == 1,evtWinIdx);
TrialWin_second_array = array(events == 2,evtWinIdx);

figure; plot(mean(TrialWin_first_array))
hold on; plot(mean(TrialWin_second_array))

figure; plot(mean(array(events == 1, :)))
hold on; plot(mean(array(events == 2, :)))

figure; plot(mean(zall_array{1, plot_neuron_num}))
hold on; plot(mean(zall_array{2, plot_neuron_num}))

% Perform two-sample t-test
[h, p, ci, stats] = ttest2(mean(TrialWin_first_array), mean(TrialWin_second_array));

% Display results
if h
    disp('The means of the two arrays are statistically different.');
else
    disp('The means of the two arrays are not statistically different.');
end
disp(['p-value: ', num2str(p)]);
disp(['Confidence interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['t-statistic: ', num2str(stats.tstat)]);
disp(['Degrees of freedom: ', num2str(stats.df)]);

%%

first_array = caTraceTrials_unnormalized_array{1, 1};
second_array = caTraceTrials_unnormalized_array{2, 1};

array = [first_array; second_array]; 

events = [(ones(size(first_array, 1), 1)); (ones(size(second_array, 1), 1)*2)]; 

array = zscore(array, 0, 2); 


TrialWin_first_array = array(events == 1,:);
TrialWin_second_array = array(events == 2,:);

% Perform two-sample t-test
[h, p, ci, stats] = ttest2(mean(TrialWin_first_array), mean(TrialWin_second_array));

% Display results
if h
    disp('The means of the two arrays are statistically different.');
else
    disp('The means of the two arrays are not statistically different.');
end
disp(['p-value: ', num2str(p)]);
disp(['Confidence interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['t-statistic: ', num2str(stats.tstat)]);
disp(['Degrees of freedom: ', num2str(stats.df)]);
