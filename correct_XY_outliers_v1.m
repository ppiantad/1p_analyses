function [X_data_corrected, Y_data_corrected] = correct_XY_outliers_v1(X_data, Y_data)


% % Define threshold for outlier detection
threshold = 2.5;
% 
% figure;
% plot(X_data); hold on; plot(Y_data); title('X_data and Y_data before')
% 
% % Detect outliers using z-score method
% z_scores = abs(zscore(X_data));
% outlier_indices = find(z_scores > threshold);
% 
% % Interpolate outliers using linear interpolation
% non_outlier_indices = setdiff(1:numel(X_data), outlier_indices);
% X_data(outlier_indices) = interp1(non_outlier_indices, X_data(non_outlier_indices), outlier_indices, 'linear', 'extrap');
% 
% 
% 
% 
% %%
% % Interpolate corresponding outliers in Y_data using linear interpolation
% Y_data(outlier_indices) = interp1(non_outlier_indices, Y_data(non_outlier_indices), outlier_indices, 'linear', 'extrap');
% 
% figure;
% plot(X_data); hold on; plot(Y_data); title('X_data and Y_data after correction')
% 
% X_data_corrected = X_data;
% Y_data_corrected = Y_data; 

%%


figure;
plot(X_data); hold on; plot(Y_data); title('X_data and Y_data before')

% Detect outliers in X_data using z-score method
z_scores_X = abs(zscore(X_data));
outlier_indices_X = find(z_scores_X > threshold);

% Detect outliers in Y_data using z-score method
z_scores_Y = abs(zscore(Y_data));
outlier_indices_Y = find(z_scores_Y > threshold);

% Combine outlier indices from X_data and Y_data
combined_outlier_indices = unique([outlier_indices_X(:); outlier_indices_Y(:)]);

% Interpolate outliers in X_data using linear interpolation
non_outlier_indices_X = setdiff(1:numel(X_data), outlier_indices_X);
X_data(outlier_indices_X) = interp1(non_outlier_indices_X, X_data(non_outlier_indices_X), outlier_indices_X, 'linear', 'extrap');

% Interpolate outliers in Y_data using linear interpolation
non_outlier_indices_Y = setdiff(1:numel(Y_data), outlier_indices_Y);
Y_data(outlier_indices_Y) = interp1(non_outlier_indices_Y, Y_data(non_outlier_indices_Y), outlier_indices_Y, 'linear', 'extrap');

% Interpolate both X_data and Y_data based on combined outlier array
non_outlier_indices_combined = setdiff(1:numel(X_data), combined_outlier_indices);
X_data(combined_outlier_indices) = interp1(non_outlier_indices_combined, X_data(non_outlier_indices_combined), combined_outlier_indices, 'linear', 'extrap');
Y_data(combined_outlier_indices) = interp1(non_outlier_indices_combined, Y_data(non_outlier_indices_combined), combined_outlier_indices, 'linear', 'extrap');


figure;
plot(X_data); hold on; plot(Y_data); title('X_data and Y_data after correction')

X_data_corrected = X_data;
Y_data_corrected = Y_data; 
