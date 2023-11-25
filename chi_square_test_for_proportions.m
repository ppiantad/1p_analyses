% must run relevant portions of
% eventRelatedActivityAndClassification_PTP_v2 first!!




% Define the observed data

total_neurons = neuron_num;

% comparison_data = sum_activated;

comparison_data = sum_inhibited;

observed_data = [comparison_data; total_neurons-sum_activated]; % First row: activated, Second row: not activated

% Calculate the expected data manually
row_totals = sum(observed_data, 2);
column_totals = sum(observed_data, 1);
total_sum = sum(observed_data(:));

% Calculate the expected data manually
expected_data = (row_totals * column_totals) / total_sum;


% Calculate the chi-squared statistic manually
chi2 = sum((observed_data(:) - expected_data(:)).^2 ./ expected_data(:));

% Calculate the degrees of freedom
degrees_of_freedom = (size(observed_data, 1) - 1) * (size(observed_data, 2) - 1);

% Calculate the p-value
p = 1 - chi2cdf(chi2, degrees_of_freedom);

% Display the results
disp(['Chi-Squared Statistic: ', num2str(chi2)]);
disp(['P-value: ', num2str(p)]);

% Check if the result is statistically significant at a chosen significance level (e.g., 0.05)
alpha = 0.05;
if p < alpha
    disp('The change in proportions is statistically significant.');
else
    disp('There is no significant change in proportions.');
end