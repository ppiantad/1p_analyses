risk_percent = [60.60606	57.57576	65.15152	72.72727	75.75758	56.06061	62.12121	72.72727	78.78788	63.63636	59.09091		59.09091	71.21212];
shock_sens = [0.38	0.46	0.38	0.46	0.6	0.52	0.52	0.44	0.5	0.5	0.56		0.42	0.4];
free_feed_weight = [27.9 27.6 26.8 27.2 26.4 29 27.4 27 28.5 26.1 33.3 26.1 29];
session_weight_raw_grams = [22.6 22.5 21.9 22.3 22.5 24 21.1 21.6 23.1 22 25.9 21.5 23.7];
session_weight_percent_free_feed = [81.00358423 81.52173913 81.71641791 81.98529412 85.22727273 82.75862069 77.00729927 80 81.05263158 84.29118774 77.77777778 82.37547893 81.72413793];

%%
x = risk_percent;
y = shock_sens;

% Calculate the coefficient of determination (r-squared value)
p = polyfit(x, y, 1);
y_fit = polyval(p, x);
y_mean = mean(y);
SS_tot = sum((y - y_mean).^2);
SS_res = sum((y - y_fit).^2);
r_squared = 1 - (SS_res / SS_tot);
figure;
% Plot the scatter plot
scatter(x, y);
hold on;

% Plot the linear fit line
plot(x, y_fit, 'r');

% Add title and labels
title(['Scatter Plot with R-Squared = ', num2str(r_squared)]);
xlabel('RDT D1 AVG Risk %');
ylabel('Shock Sensitivity (mA)');

% Show the legend
legend('Data', 'Linear Fit', 'Location', 'best');

% Release hold
hold off;

%%
x = risk_percent;
y = session_weight_percent_free_feed;

% Calculate the coefficient of determination (r-squared value)
p = polyfit(x, y, 1);
y_fit = polyval(p, x);
y_mean = mean(y);
SS_tot = sum((y - y_mean).^2);
SS_res = sum((y - y_fit).^2);
r_squared = 1 - (SS_res / SS_tot);
figure;
% Plot the scatter plot
scatter(x, y);
hold on;

% Plot the linear fit line
plot(x, y_fit, 'r');

% Add title and labels
title(['Scatter Plot with R-Squared = ', num2str(r_squared)]);
xlabel('RDT D1 AVG Risk %');
ylabel('Session Weight (% free feed');

% Show the legend
legend('Data', 'Linear Fit', 'Location', 'best');

% Release hold
hold off;

%%
x = risk_percent;
y = session_weight_raw_grams;

% Calculate the coefficient of determination (r-squared value)
p = polyfit(x, y, 1);
y_fit = polyval(p, x);
y_mean = mean(y);
SS_tot = sum((y - y_mean).^2);
SS_res = sum((y - y_fit).^2);
r_squared = 1 - (SS_res / SS_tot);
figure;
% Plot the scatter plot
scatter(x, y);
hold on;

% Plot the linear fit line
plot(x, y_fit, 'r');

% Add title and labels
title(['Scatter Plot with R-Squared = ', num2str(r_squared)]);
xlabel('RDT D1 AVG Risk %');
ylabel('Session Weight (g)');

% Show the legend
legend('Data', 'Linear Fit', 'Location', 'best');

% Release hold
hold off;

%%
x = risk_percent;
y = free_feed_weight;

% Calculate the coefficient of determination (r-squared value)
p = polyfit(x, y, 1);
y_fit = polyval(p, x);
y_mean = mean(y);
SS_tot = sum((y - y_mean).^2);
SS_res = sum((y - y_fit).^2);
r_squared = 1 - (SS_res / SS_tot);
figure;
% Plot the scatter plot
scatter(x, y);
hold on;

% Plot the linear fit line
plot(x, y_fit, 'r');

% Add title and labels
title(['Scatter Plot with R-Squared = ', num2str(r_squared)]);
xlabel('RDT D1 AVG Risk %');
ylabel('Free Feed Weight (g)');

% Show the legend
legend('Data', 'Linear Fit', 'Location', 'best');

% Release hold
hold off;

