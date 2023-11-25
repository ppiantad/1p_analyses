sum_activated_percent_sum = sum(sum_activated_percent);
data = [sum_activated_percent 100-sum_activated_percent_sum];

donutchart(data)